use super::errors::DetectChemistryErrors;
use crate::detect_chemistry::chemistry_filter::{ChemistryFilter, DetectChemistryUnit};
use crate::detect_chemistry::length_filter::LengthFilter;
use anyhow::Result;
use barcode::{BcSegSeq, WhitelistSource};
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use fastq_set::read_pair::{ReadPair, ReadPart};
use metric::{set, PercentMetric, SimpleHistogram, TxHashSet};
use parameters_toml::min_major_probe_bc_frac;
use std::iter::zip;

pub(crate) const MIN_VALID_PROBE_BCS: usize = 1_000;

pub(crate) fn check_probe_bc(
    units: &[DetectChemistryUnit],
    read_pairs_all_units: &[Vec<ReadPair>],
    feature_ref: Option<&FeatureReferenceFile>,
    feature_config: Option<&FeatureConfig>,
) -> Result<()> {
    use ChemistryName::{MFRP, MFRP_R1, SFRP};
    let allowed_chems = set![SFRP, MFRP, MFRP_R1];

    for (unit, read_pairs) in zip(units, read_pairs_all_units) {
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
            [..] if chemistries.contains(&MFRP) => MFRP,
            [..] => unreachable!(),
        };

        let chem_def = ChemistryDef::named(chemistry);
        let whitelist_source =
            WhitelistSource::from_spec(chem_def.barcode_whitelist().probe(), true, None)?;
        let whitelist = whitelist_source.as_whitelist()?;
        let id_map = whitelist_source.as_raw_seq_to_id()?;
        let bc_range = chem_def.barcode_range().probe();

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
