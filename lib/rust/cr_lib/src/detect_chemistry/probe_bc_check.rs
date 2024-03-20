use super::errors::DetectChemistryErrors;
use crate::detect_chemistry::chemistry_filter::{ChemistryFilter, DetectChemistryUnit};
use crate::detect_chemistry::length_filter::LengthFilter;
use anyhow::Result;
use barcode::BcSegSeq;
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use cr_types::LibraryType;
use fastq_set::read_pair::{ReadPair, ReadPart};
use metric::{set, PercentMetric, SimpleHistogram, TxHashSet};
use parameters_toml::min_major_probe_bc_frac;

pub(crate) const MIN_VALID_PROBE_BCS: usize = 1_000;

pub(crate) fn validate_no_probe_bc_mixture_in_sfrp(
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
    feature_ref: Option<&FeatureReferenceFile>,
    feature_config: Option<&FeatureConfig>,
) -> Result<()> {
    use ChemistryName::{MFRP_Ab, MFRP_Ab_R1, MFRP_CRISPR, MFRP_RNA, MFRP_RNA_R1, SFRP};

    for (unit, read_pairs) in units {
        // Subset allowed chems down to just the ones relevant for the library type.
        let allowed_chems = match unit.library_type {
            LibraryType::Gex => set![SFRP, MFRP_RNA, MFRP_RNA_R1],
            LibraryType::Antibody => set![SFRP, MFRP_Ab, MFRP_Ab_R1],
            LibraryType::Crispr => set![SFRP, MFRP_CRISPR],
            _ => unreachable!(),
        };
        // Check that read length is compatibe with MFRP i.e. includes probe barcodes
        let chemistries: Vec<ChemistryName> =
            LengthFilter::new(&allowed_chems, &allowed_chems, feature_ref, feature_config)?
                .process_unit(unit, read_pairs)?
                .into_iter()
                .filter(|&x| x != SFRP)
                .collect();

        let chemistry = match chemistries.as_slice() {
            [] => {
                println!(
                    "Read length is incompatible with MFRP, skipping probe BC check in {unit}."
                );
                continue;
            }
            [chemistry] => *chemistry,
            // Both MFRP and MFRP-R1 are compatible, so use MFRP.
            [..] if chemistries.contains(&MFRP_RNA) => MFRP_RNA,
            // Both MFRP-Ab and MFRP-Ab-R1 are compatible, so use MFRP.
            [..] if chemistries.contains(&MFRP_Ab) => MFRP_Ab,
            [..] => unreachable!(),
        };

        let chem_def = ChemistryDef::named(chemistry);
        let bc_range = chem_def.barcode_range().probe();
        let whitelist_source = chem_def.barcode_whitelist().probe().as_source(true)?;
        let id_map = whitelist_source.as_raw_seq_to_id()?;
        let whitelist = whitelist_source.as_whitelist()?;

        println!("\nChecking probe barcode mixtures in {unit}");
        let bc_counts: SimpleHistogram<_> = read_pairs
            .iter()
            .filter_map(|read_pair| {
                read_pair
                    .get_range(bc_range, ReadPart::Seq)
                    .and_then(|seq| whitelist.match_to_whitelist(BcSegSeq::from_bytes(seq)))
                    .map(|bc_in_wl| &id_map[&bc_in_wl])
            })
            .collect();

        let num_valid_bcs = bc_counts.raw_counts().sum();
        if usize::try_from(num_valid_bcs)? < MIN_VALID_PROBE_BCS {
            print!(
                "Skipping probe BC check in {unit}, because there were not enough reads with valid probe barcodes to confirm singleplex Fixed RNA Profiling chemistry.\n\
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

        let min_major_probe_bc_frac = *min_major_probe_bc_frac()?;
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
