//! Module for filtering raw reads based on adapter configurations to determine
//! compatible chemistries and trimming parameters.
//!
//! This filter operates before the standard read loading process. It
//! determines which adapter layout
//! (if any) matches the raw reads. This result is used to provide the necessary
//! configuration to split/trim the raw reads into a format compatible with
//! downstream processing, and to winnow down the possible chemistries.
use crate::detect_chemistry::chemistry_filter::LoadedUnit;
use anyhow::{Result, bail};
use cr_types::LibraryType;
use cr_types::chemistry::adapter::{trimmer_5p, trimmer_flex};
use cr_types::rna_read::AdapterTrimmingParameters;
use fastq_set::WhichRead;
use fastq_set::long_read_adapter_trimmer::{AdapterLayout, TrimmerConfig};
use fastq_set::read_pair::{ReadPair, ReadPart};
use itertools::Itertools;

/// Threshold of the fraction of reads that need to be successfully trimmed.
///
/// An adapter configuration that produces at least this fraction of valid
/// trimmed reads will be included in further chemistry detection.
const MIN_FRAC_TRIMMABLE: f64 = 0.1;

/// Process all units for one library type using adapter trimming.
///
/// Determine which adapter trimmer to use, then re-load the reads using the
/// trimmer to filter the reads. Returns the filtered units and the selected
/// adapter layout.
pub(crate) fn preprocess_units(
    loaded_units: Vec<LoadedUnit>,
    params: &AdapterTrimmingParameters,
    library_type: LibraryType,
    sample_reads: usize,
    total_reads: usize,
) -> Result<(Vec<LoadedUnit>, AdapterLayout)> {
    let selected = select_trimmer(params, library_type, &loaded_units)?;
    let selected_layout = selected.layout.clone();
    let trimmer = selected.into_trimmer();
    let units = loaded_units
        .into_iter()
        .map(|(unit, _)| {
            let reads = unit.sampled_read_pairs(sample_reads, total_reads, Some(&trimmer))?;
            anyhow::Ok((unit, reads))
        })
        .try_collect()?;
    Ok((units, selected_layout))
}

/// Score the trimmers on the units in order. Return the first trimmer that
/// exceeds the success threshold.
fn select_trimmer(
    params: &AdapterTrimmingParameters,
    library_type: LibraryType,
    units: &[LoadedUnit],
) -> Result<TrimmerConfig> {
    let trimmers = get_trimmers_for_detect_chem(params, library_type);
    if trimmers.is_empty() {
        bail!("unable to identify any candidate read trimmers for {library_type}");
    }
    let total_reads = units.iter().map(|(_unit, r)| r.len()).sum::<usize>();
    let mut bad_trims = vec![];
    for trimmer in trimmers {
        let count =
            count_trimmable_reads(trimmer.clone(), units.iter().flat_map(|(_, r)| r.iter()));
        let frac = count as f64 / total_reads as f64;
        println!("{frac} trimmable with {trimmer:?}");
        if frac >= MIN_FRAC_TRIMMABLE {
            return Ok(trimmer);
        }
        bad_trims.push((trimmer, frac));
    }

    bail!(
        "no trimmers handled more than {MIN_FRAC_TRIMMABLE} of reads:\n{}",
        bad_trims
            .into_iter()
            .map(|(trimmer, frac)| format!("{frac} trimmable with {trimmer:?}"))
            .format("\n")
    );
}

/// Get trimmers in priority order. Trimmers earlier in the list should take
/// precedence over trimmers later in the list (most to least specific).
fn get_trimmers_for_detect_chem(
    params: &AdapterTrimmingParameters,
    library_type: LibraryType,
) -> Vec<TrimmerConfig> {
    let mut trimmers = trimmer_flex(params, library_type);
    if let Some(trimmer) = trimmer_5p(params, library_type) {
        trimmers.push(trimmer);
    }
    trimmers
}

/// Count the number of reads that can be successfully trimmed using this config.
fn count_trimmable_reads<'a>(
    search_config: TrimmerConfig,
    reads: impl Iterator<Item = &'a ReadPair>,
) -> usize {
    let trimmer = search_config.into_trimmer();
    reads
        .filter(|read| {
            let Some(seq) = read.get(WhichRead::R1, ReadPart::Seq) else {
                return false;
            };
            let (inserts, _) = trimmer.boundary_detection(seq);
            !inserts.is_empty()
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect_chemistry::chemistry_filter::{DetectChemistryUnit, detect_chemistry_units};
    use cr_types::AlignerParam;
    use cr_types::chemistry::adapter::layouts;
    use cr_types::sample_def::{FastqMode, SampleDefBuilder};
    use fastq_set::long_read_adapter_trimmer::AdapterLayout;

    fn load_units(
        units: Vec<DetectChemistryUnit>,
        sample_reads: usize,
        total_reads: usize,
    ) -> Vec<LoadedUnit> {
        units
            .into_iter()
            .map(|u| {
                let reads = u
                    .sampled_read_pairs(sample_reads, total_reads, None)
                    .unwrap();
                (u, reads)
            })
            .collect()
    }

    /// Verify the preprocessing flow per library type:
    /// - We selected the right trimmer
    /// - Every returned read set has size that matches our regression expectation
    fn run_preprocess_units_test(
        path: &str,
        lib_type: LibraryType,
        expected_adapter_layout: AdapterLayout,
        expected_passing_reads: usize,
    ) {
        const NUM_READS: usize = 16_000;
        const SAMPLE_READS: usize = 8_000;
        let sdef = SampleDefBuilder::new(FastqMode::SBX, path.into())
            .library_type(lib_type)
            .build();
        let units = detect_chemistry_units(&[sdef], None, None)
            .unwrap()
            .into_values()
            .exactly_one()
            .unwrap();
        let params = units[0]
            .fastq_mode
            .adapter_trimming_parameters(Some(AlignerParam::Star))
            .unwrap();
        let loaded_units = load_units(units, SAMPLE_READS, NUM_READS);

        let selected_trimmer = select_trimmer(&params, lib_type, &loaded_units).unwrap();
        assert_eq!(expected_adapter_layout, selected_trimmer.layout);
        let (reloaded_units, _) =
            preprocess_units(loaded_units, &params, lib_type, NUM_READS, NUM_READS).unwrap();
        let (_, filtered_reads) = reloaded_units.into_iter().exactly_one().unwrap();

        assert_eq!(expected_passing_reads, filtered_reads.len());
    }

    // 3' and 5' Gex for the purposes of chemistry detection can use the
    // same trimming strategy.

    #[test]
    fn test_preprocess_units_sbx_threeprime_fulllen() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_threeprime_fulllen",
            LibraryType::Gex,
            layouts::FIVE_PRIME_GEX,
            14446,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_threeprime_frag() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_threeprime_frag",
            LibraryType::Gex,
            layouts::FIVE_PRIME_GEX,
            14272,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_fiveprime_fulllen() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_fiveprime_fulllen",
            LibraryType::Gex,
            layouts::FIVE_PRIME_GEX,
            14786,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_fiveprime_frag() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_fiveprime_frag",
            LibraryType::Gex,
            layouts::FIVE_PRIME_GEX,
            14538,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_flex_gex() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_flex_gex",
            LibraryType::Gex,
            layouts::FLEX_GEX,
            12612,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_flex_ab() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_flex_ab",
            LibraryType::Antibody,
            layouts::FLEX_ANTIBODY,
            12385,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_flex_crispr() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_flex_crispr",
            LibraryType::Crispr,
            layouts::FLEX_CRISPR,
            11880,
        );
    }

    #[test]
    fn test_preprocess_units_sbx_flex_hybrid_crispr() {
        run_preprocess_units_test(
            "test/fastq/sbx/sbx_flex_hybrid_crispr",
            LibraryType::Crispr,
            layouts::FLEX_HYBRID_CRISPR,
            12413,
        );
    }
}
