//! Martian stage CHECK_BARCODES_COMPATIBILITY

use anyhow::{bail, Result};
use barcode::whitelist::find_whitelist;
use barcode::{BarcodeConstruct, BcSegSeq, Whitelist, WhitelistSource, WhitelistSpec};
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::rna_read::LegacyLibraryType;
use cr_types::sample_def::SampleDef;
use fastq_set::filenames::FindFastqs;
use fastq_set::read_pair::{ReadPart, RpRange};
use fastq_set::read_pair_iter::ReadPairIter;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use metric::{set, Metric, SimpleHistogram, TxHashMap, TxHashSet};
use parameters_toml::min_barcode_similarity;
use serde::{Deserialize, Serialize};

const MAX_READS_BARCODE_COMPATIBILITY: usize = 1_000_000;
const ROBUST_FRACTION_THRESHOLD: f64 = 0.925;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CheckBarcodesCompatibilityStageInputs {
    pub chemistry_def: ChemistryDef,
    pub sample_def: Vec<SampleDef>,
    pub check_library_compatibility: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct, PartialEq, Eq)]
#[cfg_attr(test, derive(Debug))]
pub struct CheckBarcodesCompatibilityStageOutputs {
    pub libraries_to_translate: TxHashSet<LegacyLibraryType>,
}

// This is our stage struct
pub struct CheckBarcodesCompatibility;

pub(crate) fn sample_valid_barcodes(
    sample_def: &SampleDef,
    barcode_range: RpRange,
    wl: &Whitelist,
) -> Result<SimpleHistogram<BcSegSeq>> {
    let mut histogram = SimpleHistogram::new();
    let mut num_reads = 0;
    for fastq in sample_def.get_fastq_def()?.find_fastqs()? {
        for read_pair in ReadPairIter::from_fastq_files(&fastq)? {
            if let Some(seq) = read_pair?.get_range(barcode_range, ReadPart::Seq) {
                // NOTE: This is robust to a single N cycle
                if let Some(bc_in_wl) = wl.match_to_whitelist(BcSegSeq::from_bytes(seq)) {
                    histogram.observe_owned(bc_in_wl);
                }
                num_reads += 1;
            }
            if num_reads >= MAX_READS_BARCODE_COMPATIBILITY {
                return Ok(histogram);
            }
        }
    }
    Ok(histogram)
}

/// Calculate the robust cosine similarity of two barcode lists.
/// Implements the same basic cosine distance metric, but first
/// caps each count value to a threshold value such that robust_fracion_threshold
/// of the total counts are greater than or equal to the threshold.
/// This reduces the impact of very high count outliers.
///
/// Panics if the histograms is empty
fn robust_cosine_similarity(c1: &SimpleHistogram<BcSegSeq>, c2: &SimpleHistogram<BcSegSeq>) -> f64 {
    let thresh1 = match stats::nx::nx(c1.raw_counts(), ROBUST_FRACTION_THRESHOLD) {
        Some(v) => v,
        None => return 0.0, // Empty counts => 0 similarity
    };
    let thresh2 = match stats::nx::nx(c2.raw_counts(), ROBUST_FRACTION_THRESHOLD) {
        Some(v) => v,
        None => return 0.0, // Empty counts => 0 similarity
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
        .distribution()
        .iter()
        .map(|(bc, c)| (c.count().min(thresh1) * c2.get(bc).min(thresh2)) as f64)
        .sum();

    dot_prod / (mag1 * mag2)
}

#[make_mro]
impl MartianMain for CheckBarcodesCompatibility {
    type StageInputs = CheckBarcodesCompatibilityStageInputs;
    type StageOutputs = CheckBarcodesCompatibilityStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let unique_lib_types: TxHashSet<_> = args
            .sample_def
            .iter()
            .map(|sdef| sdef.library_type.unwrap_or_default())
            .collect();

        // -----------------------------------------------------------------------------------------
        if let BarcodeConstruct::Segmented(_) = args.chemistry_def.barcode_whitelist() {
            let libraries_to_translate = set![];
            return Ok(CheckBarcodesCompatibilityStageOutputs {
                libraries_to_translate,
            });
        }

        // -----------------------------------------------------------------------------------------
        // Trivial case when we don't have more than 1 library type
        if unique_lib_types.len() < 2 {
            // We are making a choice to "translate" 3' v3 antibody only libraries by default. This
            // is not be the correct thing to do if the input is TotalSeqA, but we do it anyway,
            // because we do not have a notion of a "canonical" whitelist and we can only get one of
            // TotalSeqA or TotalSeqB correct until we fix that.
            let libraries_to_translate = if (args.chemistry_def.name == ChemistryName::ThreePrimeV3
                || args.chemistry_def.name == ChemistryName::ThreePrimeV3LT)
                && unique_lib_types.contains(&LegacyLibraryType::AntibodyCapture)
            {
                set![LegacyLibraryType::AntibodyCapture]
            } else {
                set![]
            };
            return Ok(CheckBarcodesCompatibilityStageOutputs {
                libraries_to_translate,
            });
        }

        // -----------------------------------------------------------------------------------------
        // There needs to be a Gene expression library if there are >1 library type. This will
        // possibly be relaxed in the future
        if !unique_lib_types.contains(&LegacyLibraryType::GeneExpression) {
            bail!("Gene expression data is required if there are multiple library types.");
        }

        // -----------------------------------------------------------------------------------------
        // Compute gel bead barcode histogram per library type
        let gb_barcode_whitelist = match args.chemistry_def.barcode_whitelist().gel_bead() {
            WhitelistSpec::TxtFile { name } => name,
            _ => unreachable!(),
        };
        let wl = Whitelist::new(&find_whitelist(gb_barcode_whitelist, false, None)?)?;
        let mut per_lib_bc_histogram = TxHashMap::default();
        for sdef in &args.sample_def {
            let this_hist =
                sample_valid_barcodes(sdef, args.chemistry_def.barcode_range().gel_bead(), &wl)?;
            per_lib_bc_histogram
                .entry(sdef.library_type.unwrap())
                .or_insert_with(SimpleHistogram::new)
                .merge(this_hist);
        }

        let gex_hist = per_lib_bc_histogram
            .remove(&LegacyLibraryType::GeneExpression)
            .unwrap(); // This is guaranteed to succeed due to the check above

        // -----------------------------------------------------------------------------------------
        // Check similarity with the GEX library and infer if we need translation.
        let translation_map = match find_whitelist(gb_barcode_whitelist, true, None) {
            Ok(path) => Some(WhitelistSource::txt_file(&path).as_translation()?),
            Err(_) => None,
        };
        let mut libraries_to_translate = set![];
        let min_barcode_similarity = *min_barcode_similarity()?;
        for (lib_type, this_hist) in per_lib_bc_histogram {
            let mut similarity = robust_cosine_similarity(&gex_hist, &this_hist);
            println!("Without translation: {lib_type} - {similarity:?}");
            if let Some(ref translate) = translation_map {
                let trans_similarity =
                    robust_cosine_similarity(&gex_hist, &this_hist.map_key(|key| translate[&key]));
                println!("With translation   : {lib_type} - {trans_similarity:?}");
                if trans_similarity > similarity {
                    libraries_to_translate.insert(lib_type);
                    similarity = trans_similarity;
                }
            }
            if args.check_library_compatibility && similarity < min_barcode_similarity {
                bail!(
                    "Barcodes from the [{}] library and the [{}] library have insufficient overlap. \
                     This usually indicates the libraries originated from different cells or samples. \
                     This error can usually be fixed by providing correct FASTQ files from the same \
                     sample. If you are certain the input libraries are matched, you can bypass \
                     this check by adding `check-library-compatibility,false` to the \
                     [gene-expression] section of your multi config CSV if using `cellranger multi` \
                     or passing the --check-library-compatibility=false argument if using \
                     `cellranger count`. If you have questions regarding this error or your results, \
                     please contact support@10xgenomics.com.",
                    LegacyLibraryType::GeneExpression,
                    lib_type
                );
            }
        }

        Ok(CheckBarcodesCompatibilityStageOutputs {
            libraries_to_translate,
        })
    }
}

#[cfg(test)]
mod barcode_compatibility_tests {
    use super::*;
    use cr_types::chemistry::ChemistryName;
    use cr_types::sample_def::FastqMode;
    use dui_tests::stage_test::StageFailTest;
    use dui_tests::{stage_fail_dui_test, DuiTest};

    #[test]
    fn test_single_lib_type() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![SampleDef {
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.libraries_to_translate, set![]);
    }

    #[test]
    fn test_ab_only_3pv3() {
        // By default we will translate the Antibody library in 3' v3. This is not the correct
        // thing to do for TotalSeqA, but we prefer TotalSeqB.
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                        .into(),
                sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                lanes: Some(vec![4]),
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.libraries_to_translate,
            set![LegacyLibraryType::AntibodyCapture]
        );
    }

    #[test]
    fn test_ab_only_3pv3lt() {
        // By default we will translate the Antibody library in 3' v3 LT. This is not the correct
        // thing to do for TotalSeqA, but we prefer TotalSeqB.
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3LT),
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry_ab".into(),
                sample_names: Some(vec!["1077130".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.libraries_to_translate,
            set![LegacyLibraryType::AntibodyCapture]
        );
    }

    #[test]
    fn test_5p_ab_only() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::FivePrimeR2),
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_antibody"
                    .into(),
                sample_names: Some(vec!["vdj_v1_hs_pbmc2_antibody".into()]),
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.libraries_to_translate, set![]);
    }

    #[test]
    fn test_ab_and_gex_3pv3() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.libraries_to_translate,
            set![LegacyLibraryType::AntibodyCapture]
        );
    }

    #[test]
    fn test_ab_and_gex_3pv3lt() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3LT),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry_gex_ab"
                        .into(),
                    sample_names: Some(vec!["1077130_ab".into()]),
                    lanes: Some(vec![1]),
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry_gex_ab"
                        .into(),
                    sample_names: Some(vec!["1077130_gex".into()]),
                    lanes: Some(vec![1]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.libraries_to_translate,
            set![LegacyLibraryType::AntibodyCapture]
        );
    }

    #[test]
    fn test_crispr_and_gex_3pv3() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_crispr"
                            .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_crispr".into()]),
                    library_type: Some(LegacyLibraryType::CrisprGuideCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex"
                        .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_gex".into()]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.libraries_to_translate,
            set![LegacyLibraryType::CrisprGuideCapture]
        );
    }

    #[test]
    fn test_incompatible() -> Result<()> {
        // Antibody from pbmc_1k_protein_v3
        // GEX from CR-300-01
        stage_fail_dui_test!(
            IncompatibleLibraries,
            description: "If the input data contains multiple library types, the pipeline produces \
            an error if the barcodes from the two libraries are not compatible (meaning those \
            two libraries very likely came from different set of cells).",
            stage: CheckBarcodesCompatibility,
            args: CheckBarcodesCompatibilityStageInputs {
                chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
                sample_def: vec![
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                                .into(),
                        sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                        lanes: Some(vec![4]),
                        library_type: Some(LegacyLibraryType::AntibodyCapture),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                            .into(),
                        sample_names: Some(vec!["CR-300-01".into()]),
                        library_type: Some(LegacyLibraryType::GeneExpression),
                        ..Default::default()
                    },
                ],
                check_library_compatibility: true,
            },
        );
        Ok(())
    }

    #[test]
    fn test_incompatible_override() -> Result<()> {
        // Same as test_incompatible but ensure it works with overid
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            check_library_compatibility: false,
        };
        CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        Ok(())
    }

    #[test]
    fn test_missing_gex() -> Result<()> {
        stage_fail_dui_test!(
            MissingGex,
            description: "If the input data contains multiple library types, at least one of them \
            needs to be Gene Expression. The test case contains two library types, Antibody \
            Capture and Crispr Guide Capture",
            stage: CheckBarcodesCompatibility,
            args: CheckBarcodesCompatibilityStageInputs {
                chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
                sample_def: vec![
                    SampleDef {
                        library_type: Some(LegacyLibraryType::AntibodyCapture),
                        ..Default::default()
                    },
                    SampleDef {
                        library_type: Some(LegacyLibraryType::CrisprGuideCapture),
                        ..Default::default()
                    },
                ],
                check_library_compatibility: true,
            },
        );
        Ok(())
    }

    #[test]
    fn test_similarity() {
        let mut c1 = SimpleHistogram::new();
        c1.observe_by_owned(BcSegSeq::from_bytes(b"AA"), 10);
        c1.observe_by_owned(BcSegSeq::from_bytes(b"AC"), 20);

        let mut c2 = SimpleHistogram::new();
        c2.observe_by_owned(BcSegSeq::from_bytes(b"AA"), 5);
        c2.observe_by_owned(BcSegSeq::from_bytes(b"AG"), 10);

        assert!((robust_cosine_similarity(&c1, &c2) - 0.5f64).abs() < 10f64 * std::f64::EPSILON);

        let mut c3 = SimpleHistogram::new();
        c3.observe_by_owned(BcSegSeq::from_bytes(b"AC"), 5);
        c3.observe_by_owned(BcSegSeq::from_bytes(b"AG"), 10);

        assert!((robust_cosine_similarity(&c1, &c3) - 0.5f64).abs() < 10f64 * std::f64::EPSILON);
    }

    #[test]
    fn test_3pv3_total_seq_a() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    sample_indices: Some(vec!["ATTGGCAT".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-85_3pv3_totalseqA_antibody".into(),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    sample_indices: Some(vec!["SI-P2-A8".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-85_3pv3_totalseqA_gex".into(),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.libraries_to_translate, set![]);
    }

    #[test]
    fn test_3pv2_total_seq_a() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV2),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    sample_indices: Some(vec!["GATCTGAT".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-82_3pv2_totalseqA_antibody".into(),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    sample_indices: Some(vec!["TTTGTACA".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-AB-82_3pv2_totalseqA_gex".into(),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.libraries_to_translate, set![]);
    }

    #[test]
    fn test_5p_ab_plus_gex() {
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::FivePrimeR2),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_antibody"
                            .into(),
                    sample_names: Some(vec!["vdj_v1_hs_pbmc2_antibody".into()]),
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_gex"
                        .into(),
                    sample_names: Some(vec!["vdj_v1_hs_pbmc2_5gex".into()]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.libraries_to_translate, set![]);
    }

    #[test]
    fn test_cycle_failure() {
        // In this test, the first base of R1 is N in the antibody fastq
        let args = CheckBarcodesCompatibilityStageInputs {
            chemistry_def: ChemistryDef::named(ChemistryName::ThreePrimeV3),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    sample_indices: Some(vec!["SI-P2-C7".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-FZ-172_antibody"
                        .into(),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    sample_indices: Some(vec!["SI-P2-C5".into()]),
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-FZ-172_gex".into(),
                    ..Default::default()
                },
            ],
            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibility.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.libraries_to_translate,
            set![LegacyLibraryType::AntibodyCapture]
        );
    }
}
