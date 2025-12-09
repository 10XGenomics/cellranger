//! Martian stage CHECK_BARCODES_COMPATIBILITY
//!
//! This stage checks that barcodes from different libraries are compatible.
//!
//! ## What is the meaning of compatibility?
//!
//! A set of cells processed through a single 10x gem well could be used to create multiple sequencing
//! libraries (GEX, Ab etc.). In order to analyze a sample with multiple library types, say GEX + Ab,
//! the user needs to specify both sets of fastqs to cellranger. This could lead to occassional
//! mix-ups where the GEX library from one gem well could be paired with Ab library which was
//! processed through a different gem well. This stage checks that the fastqs from different
//! libraries are compatible.
//!
//! <img src="v3-gel-bead.png" alt="3' v3 gel bead" width="400"/>
//!
//! The barcode associated with the Truseq handle (often called the Truseq barcode or the GEX barcode
//! informally) is not the same as the one associated with the Nextera handle (often called the
//! Nextera barcode). Hence we have two barcodes: a truseq barcode and a nextera barcode that refers
//! to the same partition. The gel beads are designed so that there is a known relationship between
//! each pair of barcodes - we call this a translation whitelist. Barcode translation refers to the
//! process of translating from the Nextera barcode to the Truseq barcode (by convention).
//!
//! However, there are also applications where the Antibody library is captured by polyA
//! (TotalSeq A, for example). Hence **for 3' v3** if we have GEX + FB, in some cases the barcodes
//! in the FB library need to be translated whereas in other cases they need not be. Chemistry
//! detection takes care of determining the correct whitelist for each library.
//!
//! ## How do we check for compatibility?
//!
//! If all the libraries came from the same set of cells, we expect to see a similar distribution of
//! the number of reads associated with each barcode. Conversely, if they came from different set of
//! cells, we expect a negligible overlap in this distribution (Given our barcode whitelist size,
//! the number of common cell associated barcodes that you get by chance between two gem wells is
//! very small for the cell loads we operate in). We compute the similarity using:
//! 1. Sample upto 1 million reads for each library and create a histogram of the number of reads
//!    associated with each valid barcode.
//! 2. Compute a modified [cosine similarity](https://en.wikipedia.org/wiki/Cosine_similarity) of
//!    the histogram of each library with the histogram of the gene expression library. Before applying
//!    the similarity measure, we first cap each count to a threshold value to reduces the impact of
//!    very high count outliers.
//!     - If translation is available (`3' v3`), we compute similarity using the translated barcode.
//!       If the similarity is improved, we flag the library for translation.
//! 3. The libraries are compatible if their modified cosine similarity is at least `0.1` (empirical
//!    threshold).
//!
//! ## NOTES
//! - No compatibility check is performed if we have only a single library. With the current
//!   whitelist setup, we can't know if a 3' v3 Antibody only sample will need to be translated or
//!   not. In order to do this correctly, we need to define a canonical whitelist instead of a mixed
//!   one containing both Truseq and Nextera barcodes. We made a deliberate choice to enable
//!   translation by default for 3' v3 antibody only samples. This is not the "correct" thing to do
//!   TotalSeqA. A consequence of this is that the cellranger output for TotalSeqA Antibody only 3' v3
//!   will **not be barcode compatible** with a GEX only run from the same gem well.
//! - There needs to be a GEX library if there are >1 library types.
#![deny(missing_docs)]

use anyhow::{Result, anyhow, ensure};
use barcode::{BarcodeConstruct, BcSegSeq, Whitelist};
use cr_types::chemistry::ChemistryDefs;
use cr_types::sample_def::SampleDef;
use cr_types::{ERROR_CODE_INFO, LibraryType};
use fastq_set::filenames::FindFastqs;
use fastq_set::read_pair::{ReadPart, RpRange};
use fastq_set::read_pair_iter::ReadPairIter;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use metric::{Histogram, Metric, SimpleHistogram, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};

const MAX_READS_BARCODE_COMPATIBILITY: usize = 1_000_000;
const ROBUST_FRACTION_THRESHOLD: f64 = 0.925;
const MIN_BARCODE_SIMILARITY: f64 = 0.1;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CheckBarcodesCompatibilityStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub sample_def: Vec<SampleDef>,
    pub check_library_compatibility: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct, PartialEq, Eq)]
pub struct CheckBarcodesCompatibilityStageOutputs {}

/// Martian stage CHECK_BARCODES_COMPATIBILITY
pub struct CheckBarcodesCompatibility;

pub fn sample_valid_barcodes(
    sample_def: &SampleDef,
    barcode_range: RpRange,
    wl: &Whitelist,
) -> Result<SimpleHistogram<BcSegSeq>> {
    let barcode_iter = sample_def
        .get_fastq_def()?
        .find_fastqs()?
        .into_iter()
        .map(|fastq| ReadPairIter::from_fastq_files(&fastq))
        .flatten_ok()
        .map(|x| x?)
        .filter_map_ok(|read| {
            read.get_range(barcode_range, ReadPart::Seq)
                .map(BcSegSeq::from_bytes)
        })
        .take(MAX_READS_BARCODE_COMPATIBILITY)
        .filter_map_ok(|seq| wl.match_to_whitelist_allow_one_n(seq, true));
    barcode_iter.process_results(|iter| SimpleHistogram::from_iter_owned(iter))
}

/// Calculate the robust cosine similarity of two barcode lists.
/// Implements the same basic cosine distance metric, but first
/// caps each count value to a threshold value such that robust_fracion_threshold
/// of the total counts are greater than or equal to the threshold.
/// This reduces the impact of very high count outliers.
///
/// Panics if the histograms is empty
fn robust_cosine_similarity(c1: &SimpleHistogram<BcSegSeq>, c2: &SimpleHistogram<BcSegSeq>) -> f64 {
    let Some(thresh1) = stats::compute_nx(c1.raw_counts(), ROBUST_FRACTION_THRESHOLD) else {
        // Empty counts => 0 similarity
        return 0.0;
    };
    let Some(thresh2) = stats::compute_nx(c2.raw_counts(), ROBUST_FRACTION_THRESHOLD) else {
        // Empty counts => 0 similarity
        return 0.0;
    };

    let mag1 = c1
        .raw_counts()
        .map(|c| (c.min(&thresh1) * c.min(&thresh1)) as f64)
        .sum::<f64>()
        .sqrt();

    let mag2 = c2
        .raw_counts()
        .map(|&c| (c.min(thresh2) * c.min(thresh2)) as f64)
        .sum::<f64>()
        .sqrt();

    let dot_prod: f64 = c1
        .iter()
        .map(|(bc, c)| (c.count().min(thresh1) * c2.get(bc).min(thresh2)) as f64)
        .sum();

    dot_prod / (mag1 * mag2)
}

#[make_mro(mem_gb = 2)]
impl MartianMain for CheckBarcodesCompatibility {
    type StageInputs = CheckBarcodesCompatibilityStageInputs;
    type StageOutputs = CheckBarcodesCompatibilityStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let unique_lib_types: TxHashSet<_> = args.chemistry_defs.keys().collect();

        // -----------------------------------------------------------------------------------------
        if args
            .chemistry_defs
            .values()
            .any(|x| matches!(x.barcode_whitelist_spec(), BarcodeConstruct::Segmented(_)))
        {
            return Ok(CheckBarcodesCompatibilityStageOutputs {});
        }

        // -----------------------------------------------------------------------------------------
        // Trivial case when we don't have more than 1 library type
        if unique_lib_types.len() < 2 {
            return Ok(CheckBarcodesCompatibilityStageOutputs {});
        }

        // -----------------------------------------------------------------------------------------
        // There needs to be a Gene expression library if there are >1 library type. This will
        // possibly be relaxed in the future
        ensure!(unique_lib_types.contains(&LibraryType::Gex), MISSING_GEX);

        // -----------------------------------------------------------------------------------------
        // Compute gel bead barcode histogram per library type
        let mut per_lib_bc_histogram = TxHashMap::default();
        for sdef in &args.sample_def {
            let chem = &args.chemistry_defs[&sdef.library_type.unwrap()];
            let whitelist = chem.barcode_whitelist_source()?.gel_bead().as_whitelist()?;
            let this_hist =
                sample_valid_barcodes(sdef, chem.barcode_range().gel_bead(), &whitelist)?;
            per_lib_bc_histogram
                .entry(sdef.library_type.unwrap())
                .or_insert_with(SimpleHistogram::default)
                .merge(this_hist);
        }

        let gex_hist = per_lib_bc_histogram.remove(&LibraryType::Gex).unwrap(); // This is guaranteed to succeed due to the check above

        // -----------------------------------------------------------------------------------------
        // Check similarity with the GEX library.
        for (lib_type, this_hist) in per_lib_bc_histogram {
            let similarity = robust_cosine_similarity(&gex_hist, &this_hist);
            println!("Similarity to GEX: {lib_type} - {similarity:?}");
            if args.check_library_compatibility {
                ensure!(
                    similarity >= MIN_BARCODE_SIMILARITY,
                    incompatible_message(LibraryType::Gex, lib_type, similarity),
                );
            }
        }

        Ok(CheckBarcodesCompatibilityStageOutputs {})
    }
}

const MISSING_GEX: &str = "Gene expression data is required if there are multiple library types.";

/// Return the error message for libraries with insufficient barcode overlap.
fn incompatible_message(
    lib_0_type: LibraryType,
    lib_1_type: LibraryType,
    similarity: f64,
) -> Error {
    anyhow!(
        "TXRNGR10006: Barcodes from the [{lib_0_type}] library and the [{lib_1_type}] library have \
         insufficient overlap: {similarity:.4} < {MIN_BARCODE_SIMILARITY}. \
         This usually indicates the libraries originated from different cells or samples. \
         This error can usually be fixed by providing correct FASTQ files from the same \
         sample. If you are certain the input libraries are matched, you can bypass \
         this check by adding `check-library-compatibility,false` to the \
         [gene-expression] section of your multi config CSV if using `cellranger multi` \
         or passing the --check-library-compatibility=false argument if using \
         `cellranger count`. If you have questions regarding this error or your results, \
         please contact support@10xgenomics.com. {ERROR_CODE_INFO}",
    )
}

#[cfg(test)]
mod barcode_compatibility_tests {
    use super::*;
    use cr_types::chemistry::{ChemistryDef, ChemistryName};
    use cr_types::sample_def::FastqMode;

    fn inputs(
        chemistry: ChemistryName,
        samples: Vec<SampleDef>,
    ) -> CheckBarcodesCompatibilityStageInputs {
        let chem = ChemistryDef::named(chemistry);
        CheckBarcodesCompatibilityStageInputs {
            chemistry_defs: samples
                .iter()
                .map(|s| (s.library_type.unwrap(), chem.clone()))
                .collect(),
            sample_def: samples,
            check_library_compatibility: true,
        }
    }

    #[test]
    fn test_single_lib_type() {
        let args = inputs(
            ChemistryName::ThreePrimeV3PolyA,
            vec![SampleDef {
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
        );
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_ab_and_gex_3pv3() {
        let chems: ChemistryDefs = [
            (
                LibraryType::Gex,
                ChemistryDef::named(ChemistryName::ThreePrimeV3PolyA),
            ),
            (
                LibraryType::Antibody,
                ChemistryDef::named(ChemistryName::ThreePrimeV3CS1),
            ),
        ]
        .into_iter()
        .collect();
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_defs: chems,
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_crispr_and_gex_3pv3() {
        let chems: ChemistryDefs = [
            (
                LibraryType::Gex,
                ChemistryDef::named(ChemistryName::ThreePrimeV3PolyA),
            ),
            (
                LibraryType::Crispr,
                ChemistryDef::named(ChemistryName::ThreePrimeV3CS1),
            ),
        ]
        .into_iter()
        .collect();
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_defs: chems,
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_crispr"
                            .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_crispr".into()]),
                    library_type: Some(LibraryType::Crispr),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex"
                        .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_gex".into()]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_incompatible() {
        // Antibody from pbmc_1k_protein_v3
        // GEX from CR-300-01
        let chems: ChemistryDefs = [
            (
                LibraryType::Gex,
                ChemistryDef::named(ChemistryName::ThreePrimeV3PolyA),
            ),
            (
                LibraryType::Antibody,
                ChemistryDef::named(ChemistryName::ThreePrimeV3CS1),
            ),
        ]
        .into_iter()
        .collect();
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_defs: chems,
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        assert_eq!(
            incompatible_message(
                LibraryType::Gex,
                LibraryType::Antibody,
                0.00013808299955230403,
            )
            .to_string(),
            CheckBarcodesCompatibility
                .test_run_tmpdir(args)
                .unwrap_err()
                .to_string(),
        );
    }

    #[test]
    fn test_incompatible_override() {
        // Same as test_incompatible but ensure it works with override
        let mut args = inputs(
            ChemistryName::ThreePrimeV3PolyA,
            vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
        );
        args.check_library_compatibility = false;
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_missing_gex() {
        let args = inputs(
            ChemistryName::ThreePrimeV3PolyA,
            vec![
                SampleDef {
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    library_type: Some(LibraryType::Crispr),
                    ..Default::default()
                },
            ],
        );
        assert_eq!(
            MISSING_GEX,
            &CheckBarcodesCompatibility
                .test_run_tmpdir(args)
                .unwrap_err()
                .to_string(),
        );
    }

    #[test]
    fn test_similarity() {
        let mut c1 = SimpleHistogram::default();
        c1.observe_by_owned(BcSegSeq::from_bytes(b"AA"), 10);
        c1.observe_by_owned(BcSegSeq::from_bytes(b"AC"), 20);

        let mut c2 = SimpleHistogram::default();
        c2.observe_by_owned(BcSegSeq::from_bytes(b"AA"), 5);
        c2.observe_by_owned(BcSegSeq::from_bytes(b"AG"), 10);

        assert!((robust_cosine_similarity(&c1, &c2) - 0.5f64).abs() < 10f64 * f64::EPSILON);

        let mut c3 = SimpleHistogram::default();
        c3.observe_by_owned(BcSegSeq::from_bytes(b"AC"), 5);
        c3.observe_by_owned(BcSegSeq::from_bytes(b"AG"), 10);

        assert!((robust_cosine_similarity(&c1, &c3) - 0.5f64).abs() < 10f64 * f64::EPSILON);
    }

    #[test]
    fn test_3pv3_total_seq_a() {
        let args = inputs(
            ChemistryName::ThreePrimeV3PolyA,
            vec![
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LibraryType::Antibody),
                    sample_indices: Some(vec!["ATTGGCAT".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-85_3pv3_totalseqA_antibody".into(),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LibraryType::Gex),
                    sample_indices: Some(vec!["SI-P2-A8".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-85_3pv3_totalseqA_gex".into(),
                    ..Default::default()
                },
            ],
        );
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_3pv2_total_seq_a() {
        let args = inputs(
            ChemistryName::ThreePrimeV2,
            vec![
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LibraryType::Antibody),
                    sample_indices: Some(vec!["GATCTGAT".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-82_3pv2_totalseqA_antibody".into(),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LibraryType::Gex),
                    sample_indices: Some(vec!["TTTGTACA".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-82_3pv2_totalseqA_gex".into(),
                    ..Default::default()
                },
            ],
        );
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_5p_ab_plus_gex() {
        let args = inputs(
            ChemistryName::FivePrimeR2,
            vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_antibody"
                            .into(),
                    sample_names: Some(vec!["vdj_v1_hs_pbmc2_antibody".into()]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_gex"
                        .into(),
                    sample_names: Some(vec!["vdj_v1_hs_pbmc2_5gex".into()]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
        );
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }

    #[test]
    fn test_cycle_failure() {
        // In this test, the first base of R1 is N in the antibody fastq
        let chems: ChemistryDefs = [
            (
                LibraryType::Gex,
                ChemistryDef::named(ChemistryName::ThreePrimeV3PolyA),
            ),
            (
                LibraryType::Antibody,
                ChemistryDef::named(ChemistryName::ThreePrimeV3CS1),
            ),
        ]
        .into_iter()
        .collect();
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_defs: chems,
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LibraryType::Antibody),
                    sample_indices: Some(vec!["SI-P2-C7".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-FZ-172_antibody"
                        .into(),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LibraryType::Gex),
                    sample_indices: Some(vec!["SI-P2-C5".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-FZ-172_gex".into(),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
    }
}
