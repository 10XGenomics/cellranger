#![deny(missing_docs)]

use super::errors::DetectChemistryErrors;
use crate::detect_chemistry::chemistry_filter::{
    ChemistryFilter, DetectChemistryUnit, named_candidates,
};
use crate::detect_chemistry::length_filter::LengthFilter;
use anyhow::Result;
use barcode::BcSegSeq;
use cr_types::LibraryType;
use cr_types::chemistry::{BarcodeExtraction, ChemistryDef, ChemistryName};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use cr_types::rna_read::get_probe_barcode_range;
use fastq_set::read_pair::{ReadPair, ReadPart};
use itertools::Itertools;
use metric::{Histogram, PercentMetric, SimpleHistogram, TxHashSet, set};

pub(crate) const MIN_VALID_PROBE_BCS: usize = 1_000;

pub(crate) fn validate_no_probe_bc_mixture_in_sfrp(
    detected_chemistry: ChemistryName,
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
    feature_ref: Option<&FeatureReferenceFile>,
    feature_config: Option<&FeatureConfig>,
    min_major_probe_bc_frac: f64,
) -> Result<()> {
    #[allow(clippy::enum_glob_use)]
    use ChemistryName::*;

    for (unit, read_pairs) in units {
        // Subset allowed chems down to just the ones relevant for the library type.
        let allowed_chems = match (detected_chemistry, unit.library_type) {
            (SFRP, LibraryType::Gex) => set![SFRP, MFRP_RNA, MFRP_RNA_R1],
            (SFRP, LibraryType::Antibody) => set![SFRP, MFRP_Ab, MFRP_Ab_R1],
            (SFRP, LibraryType::Crispr) => set![SFRP, MFRP_CRISPR],
            (Flex_v2_singleplex, LibraryType::Gex) => {
                set![Flex_v2_singleplex, Flex_v2_R1, Flex_v2_RNA_R2]
            }
            (Flex_v2_singleplex, LibraryType::Antibody) => {
                set![Flex_v2_singleplex, Flex_v2_R1, Flex_v2_Ab_R2_64]
            }
            (Flex_v2_singleplex, LibraryType::Crispr) => {
                set![Flex_v2_singleplex, Flex_v2_R1, Flex_v2_RNA_R2]
            }
            _ => unreachable!(
                "chem={detected_chemistry} library_type={}",
                unit.library_type
            ),
        };
        // Check that read length is compatibe with MFRP i.e. includes probe barcodes
        let chemistries: Vec<ChemistryName> = LengthFilter::new(feature_ref, feature_config)?
            .process_unit(&named_candidates(&allowed_chems), unit, read_pairs)?
            .into_keys()
            .filter(|name| !name.is_sfrp())
            .collect();

        let chemistry = match chemistries.as_slice() {
            [] => {
                println!(
                    "Read length is incompatible with MFRP, skipping probe BC check in {unit}."
                );
                continue;
            }
            &[chemistry] => chemistry,
            // For Flex v1, use the R2 chemistry when both the R1 and R2 chemistry are compatible.
            [..] if chemistries.contains(&MFRP_RNA) => MFRP_RNA,
            [..] if chemistries.contains(&MFRP_Ab) => MFRP_Ab,
            // For Flex v2, use the R1 chemistry when both the R1 and R2 chemistry are compatible.
            [..] if chemistries.contains(&Flex_v2_R1) => Flex_v2_R1,
            [..] => unreachable!("{}", chemistries.iter().format(", ")),
        };

        let chem_def = ChemistryDef::named(chemistry);
        let barcode_components = chem_def.barcode_construct();
        let whitelist_source = chem_def.barcode_whitelist_source()?.probe();
        let id_map = whitelist_source.as_raw_seq_to_id()?;
        let whitelist = whitelist_source.as_whitelist()?;
        let (min_probe_bc_offset, max_probe_bc_offset) = match chem_def.barcode_extraction() {
            None => (0, 0),
            Some(BarcodeExtraction::VariableMultiplexingBarcode {
                min_offset,
                max_offset,
            }) => (*min_offset, *max_offset),
            Some(BarcodeExtraction::JointBc1Bc2 { .. }) => unreachable!(),
        };

        println!("\nChecking probe barcode mixtures in {unit}");
        let bc_counts: SimpleHistogram<_> = read_pairs
            .iter()
            .filter_map(|read_pair| {
                let bc_range = get_probe_barcode_range(
                    read_pair,
                    barcode_components,
                    &whitelist,
                    min_probe_bc_offset,
                    max_probe_bc_offset,
                );
                read_pair
                    .get_range(bc_range, ReadPart::Seq)
                    .and_then(|seq| {
                        whitelist.match_to_whitelist_allow_one_n(BcSegSeq::from_bytes(seq), false)
                    })
                    .map(|bc_in_wl| &id_map[&bc_in_wl])
            })
            .collect();

        let num_valid_bcs = bc_counts.raw_counts().sum();
        if usize::try_from(num_valid_bcs)? < MIN_VALID_PROBE_BCS {
            print!(
                "Skipping probe BC check in {unit}, because there were not enough reads with valid \
                 probe barcodes to confirm singleplex Flex chemistry.\n\
                 - Minimum number of required reads with valid probe barcodes = {MIN_VALID_PROBE_BCS}\n\
                 - Number of reads with valid probe barcodes available = {num_valid_bcs}\n",
            );
            continue;
        }

        let top_frac =
            PercentMetric::from_parts(*bc_counts.raw_counts().max().unwrap(), num_valid_bcs)
                .fraction()
                .unwrap_or(0.0f64);
        println!(
            "Top probe barcode was observed in {:.2}% of reads with valid barcodes.",
            top_frac * 100.0
        );

        if top_frac < min_major_probe_bc_frac {
            let mixture = bc_counts.top_n(3).into_iter().map(|(k, _)| k).collect();
            return Err(DetectChemistryErrors::ProbeBarcodeMixture {
                unit: Box::new(unit.clone()),
                mixture,
            }
            .into());
        }
    }
    Ok(())
}
