//! Martian stage CHECK_BARCODES_COMPATIBILITY_VDJ

use super::check_barcodes_compatibility::sample_valid_barcodes;
use anyhow::{bail, ensure, Result};
use barcode::{BcSegSeq, Whitelist};
use cr_types::chemistry::{ChemistryDef, ChemistryDefs, ChemistryName};
use cr_types::sample_def::SampleDef;
use cr_types::LibraryType;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use metric::{Metric, SimpleHistogram, TxHashSet};
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;

const MIN_SIMILARITY: f64 = 0.2;
const CELL_CALLING_THRESH: f64 = 0.9;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct CheckBarcodesCompatibilityVdjStageInputs {
    pub vdj_chemistry_def: ChemistryDef,
    pub vdj_sample_def: Vec<SampleDef>,
    pub count_chemistry_defs: ChemistryDefs,
    pub gex_sample_def: Vec<SampleDef>,
    pub check_library_compatibility: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct CheckBarcodesCompatibilityVdjStageOutputs {
    /// Computed only if there is a GEX library.
    pub similarity_score: Option<f64>,
}

// This is our stage struct
pub struct CheckBarcodesCompatibilityVdj;

/// Sort barcodes by read count and take the barcodes containing the top CELL_CALLING_THRESH percent
/// of reads as cells
fn approx_call_cells(c: &SimpleHistogram<BcSegSeq>) -> TxHashSet<BcSegSeq> {
    let thresh = c.raw_counts().map(|i| *(i) as f64).sum::<f64>() * CELL_CALLING_THRESH;
    let mut cells = TxHashSet::default();
    let mut acc = 0;
    for (bc, count) in c.distribution().iter().sorted_by_key(|(_, v)| Reverse(*v)) {
        if (acc as f64) < thresh {
            cells.insert(*bc);
        }
        acc += count.count();
    }
    cells
}

#[make_mro]
impl MartianMain for CheckBarcodesCompatibilityVdj {
    type StageInputs = CheckBarcodesCompatibilityVdjStageInputs;
    type StageOutputs = CheckBarcodesCompatibilityVdjStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let Some(gex_chemistry_def) = args.count_chemistry_defs.get(&LibraryType::Gex) else {
            // Antibody only case. We do not check for compatibility as of now. Something
            // we want to add in the future after looking at data.
            return Ok(CheckBarcodesCompatibilityVdjStageOutputs {
                similarity_score: None,
            });
        };

        let allowed_gex_chem = [
            ChemistryName::FivePrimePE,
            ChemistryName::FivePrimePEV3,
            ChemistryName::FivePrimeR1,
            ChemistryName::FivePrimeR2,
            ChemistryName::FivePrimeHT,
            ChemistryName::FivePrimeR2OH,
            ChemistryName::FivePrimeR2V3,
            ChemistryName::FivePrimeR2OHV3,
            ChemistryName::Custom,
        ];

        ensure!(
            allowed_gex_chem.contains(&gex_chemistry_def.name),
            "Only 5' gene expression library types are supported with VDJ. Got {:?}",
            gex_chemistry_def.name
        );

        let gex_sample_def: Vec<_> = args
            .gex_sample_def
            .into_iter()
            .filter(|sdef| sdef.library_type.unwrap_or_default() == LibraryType::Gex)
            .collect();

        assert!(!gex_sample_def.is_empty());

        // -----------------------------------------------------------------------------------------
        // Compute barcode histogram
        assert_eq!(
            gex_chemistry_def.barcode_whitelist(),
            args.vdj_chemistry_def.barcode_whitelist()
        );

        let wl = Whitelist::construct(gex_chemistry_def.barcode_whitelist(), false)?.gel_bead();

        let mut gex_bc_hist = SimpleHistogram::default();
        let mut vdj_bc_hist = SimpleHistogram::default();

        for sdef in &gex_sample_def {
            let this_hist =
                sample_valid_barcodes(sdef, gex_chemistry_def.barcode_range().gel_bead(), &wl)?;
            gex_bc_hist.merge(this_hist);
        }

        for sdef in &args.vdj_sample_def {
            let this_hist = sample_valid_barcodes(
                sdef,
                args.vdj_chemistry_def.barcode_range().gel_bead(),
                &wl,
            )?;
            vdj_bc_hist.merge(this_hist);
        }

        let vdj_cells = approx_call_cells(&vdj_bc_hist);
        let gex_cells = approx_call_cells(&gex_bc_hist);

        let similarity: f64 =
            (gex_cells.intersection(&vdj_cells).count() as f64) / (vdj_cells.len() as f64);

        if args.check_library_compatibility && similarity < MIN_SIMILARITY {
            bail!(
                "Barcodes from the [{}] library and the [VDJ] library have insufficient overlap. \
                 This usually indicates the libraries originated from different cells or samples. \
                 This error can usually be fixed by providing correct FASTQ files from the same \
                 sample. If you are certain the input libraries are matched, you can bypass \
                 this check by adding `check-library-compatibility,false` in the \
                 [gene-expression] section of your multi config CSV. \
                 If you have questions regarding this error or your results, please contact \
                 support@10xgenomics.com.",
                LibraryType::Gex
            );
        }

        Ok(CheckBarcodesCompatibilityVdjStageOutputs {
            similarity_score: Some(similarity),
        })
    }
}

#[cfg(test)]
mod barcode_compatibility_tests {
    use super::*;
    use cr_types::sample_def::FastqMode;

    fn gex_chem_map(name: ChemistryName) -> ChemistryDefs {
        [(LibraryType::Gex, ChemistryDef::named(name))]
            .into_iter()
            .collect()
    }

    #[test]
    fn test_invalid_gex_chemistry() {
        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::VdjAuto),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::ThreePrimeV1),
            gex_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::VdjAuto),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::ThreePrimeV2),
            gex_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::VdjAuto),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::ThreePrimeV3),
            gex_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());
    }

    #[test]
    fn test_fb_only() {
        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::VdjAuto),
                ..Default::default()
            }],
            count_chemistry_defs: [(
                LibraryType::Antibody,
                ChemistryDef::named(ChemistryName::FeatureBarcodingOnly),
            )]
            .into_iter()
            .collect(),
            gex_sample_def: vec![SampleDef {
                library_type: Some(LibraryType::Antibody),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args).unwrap();
        assert!(outs.similarity_score.is_none());
    }

    #[test]
    fn test_valid_pair() {
        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018272"
                        .into(),
                sample_names: Some(vec!["1018272_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["AGAATGGTTT".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017703"
                        .into(),
                sample_names: Some(vec!["1017703_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TATCAGCCTA".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        insta::assert_debug_snapshot!(outs.unwrap());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018288"
                        .into(),
                sample_names: Some(vec!["1018288_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["ATGGGTGAAA".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017703"
                        .into(),
                sample_names: Some(vec!["1017703_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TATCAGCCTA".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        insta::assert_debug_snapshot!(outs.unwrap());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018274"
                        .into(),
                sample_names: Some(vec!["1018274_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["CACAATCCCA".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017705"
                        .into(),
                sample_names: Some(vec!["1017705_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TGTCCCAACG".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        insta::assert_debug_snapshot!(outs.unwrap());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018290"
                        .into(),
                sample_names: Some(vec!["1018290_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["TCCGGGACAA".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017705"
                        .into(),
                sample_names: Some(vec!["1017705_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TGTCCCAACG".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        insta::assert_debug_snapshot!(outs.unwrap());
    }

    #[test]
    fn test_invalid_pair() {
        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018290"
                        .into(),
                sample_names: Some(vec!["1018290_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["TCCGGGACAA".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017703"
                        .into(),
                sample_names: Some(vec!["1017703_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TATCAGCCTA".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018274"
                        .into(),
                sample_names: Some(vec!["1018274_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["CACAATCCCA".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017703"
                        .into(),
                sample_names: Some(vec!["1017703_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TATCAGCCTA".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018272"
                        .into(),
                sample_names: Some(vec!["1018272_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["AGAATGGTTT".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017705"
                        .into(),
                sample_names: Some(vec!["1017705_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TGTCCCAACG".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());

        let args = CheckBarcodesCompatibilityVdjStageInputs {
            vdj_chemistry_def: ChemistryDef::named(ChemistryName::VdjPE),
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1018288"
                        .into(),
                sample_names: Some(vec!["1018288_PBMC_1k_v2_UserB_tcr".into()]),
                sample_indices: Some(vec!["ATGGGTGAAA".into()]),
                library_type: Some(LibraryType::VdjAuto),
                lanes: Some(vec![2]),
                ..Default::default()
            }],
            count_chemistry_defs: gex_chem_map(ChemistryName::FivePrimeR2),
            gex_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/barcode_matching_vdj_gex/1017705"
                        .into(),
                sample_names: Some(vec!["1017705_PBMC_1k_v2_UserB_gex".into()]),
                sample_indices: Some(vec!["TGTCCCAACG".into()]),
                library_type: Some(LibraryType::Gex),
                lanes: Some(vec![2]),
                ..Default::default()
            }],

            check_library_compatibility: true,
        };
        let outs = CheckBarcodesCompatibilityVdj.test_run_tmpdir(args);
        assert!(outs.is_err());
    }
}
