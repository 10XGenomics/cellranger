//! Martian stage DETECT_VDJ_RECEPTOR

use crate::detect_chemistry::chemistry_filter::detect_chemistry_units;
use anyhow::{bail, Result};
use cr_types::reference::feature_reference::{BeamMode, FeatureConfig, FeatureReferenceFile};
use cr_types::rna_read::LegacyLibraryType;
use cr_types::sample_def::SampleDef;
use fastq_set::read_pair::{ReadPart, WhichRead};
use fastq_set::read_pair_iter::ReadPairIter;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use metric::{Metric, TxHashSet};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::iter::zip;
use std::path::PathBuf;
use std::str::FromStr;
use vdj_reference::{KmerClassify, VdjReceptor, VdjReference};

/// How many reads to use to decide if the library is TCR or Ig
const MAX_READS_RECEPTOR_CLASSIFICATION: usize = 1_000_000;
const MIN_READS_RECEPTOR_CLASSIFICATION: usize = 10_000;
const MIN_FRAC_MAPPED_RECEPTOR_CLASSIFICATION: f64 = 0.05;
const MIN_MARGIN_RECEPTOR_CLASSIFICATION: f64 = 3.0;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct DetectVdjReceptorStageInputs {
    pub force_receptor: Option<String>,
    pub vdj_reference_path: Option<PathBuf>,
    pub feature_reference: Option<FeatureReferenceFile>,
    pub gex_sample_def: Option<Vec<SampleDef>>,
    pub vdj_sample_def: Vec<SampleDef>,
    pub is_multi: bool,
    pub feature_config: Option<FeatureConfig>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct DetectVdjReceptorStageOutputs {
    /// The output is an option because we could be in denovo mode
    /// without a VDJ reference
    pub receptor: Option<VdjReceptor>,
    pub beam_mode: Option<BeamMode>,
}

// This is our stage struct
pub struct DetectVdjReceptor;

#[derive(Debug, Metric, Serialize, Deserialize)]
struct ClassificationStats {
    tcr_reads: i64,
    ig_reads: i64,
    total_reads: i64,
}

impl fmt::Display for ClassificationStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:20} = {}", "Total Reads", self.total_reads)?;
        writeln!(f, "{:20} = {}", "Reads mapped to TR", self.tcr_reads)?;
        writeln!(f, "{:20} = {}", "Reads mapped to IG", self.ig_reads)
    }
}

impl ClassificationStats {
    fn compatible_receptor(&self) -> Option<VdjReceptor> {
        // Not enough reads
        if self.total_reads < (MIN_READS_RECEPTOR_CLASSIFICATION as i64) {
            return None;
        }
        let tcr_mapping_frac = (self.tcr_reads as f64) / (self.total_reads as f64);
        let ig_mapping_frac = (self.ig_reads as f64) / (self.total_reads as f64);
        // Not enough mapped reads
        if tcr_mapping_frac < MIN_FRAC_MAPPED_RECEPTOR_CLASSIFICATION
            && ig_mapping_frac < MIN_FRAC_MAPPED_RECEPTOR_CLASSIFICATION
        {
            return None;
        }
        if tcr_mapping_frac > MIN_MARGIN_RECEPTOR_CLASSIFICATION * ig_mapping_frac {
            Some(VdjReceptor::TR)
        } else if ig_mapping_frac > MIN_MARGIN_RECEPTOR_CLASSIFICATION * tcr_mapping_frac {
            Some(VdjReceptor::IG)
        } else {
            None
        }
    }

    fn help_text() -> String {
        format!(
            "In order to distinguish between the TR and the IG chain the following conditions \
            need to be satisfied:\n\
            - A minimum of {} total reads\n\
            - A minimum of {:.1}% of the total reads needs to map to TR or IG\n\
            - The number of reads mapped to TR should be at least {:.1}x compared to the number of \
            reads mapped to IG or vice versa",
            MIN_READS_RECEPTOR_CLASSIFICATION,
            100.0 * MIN_FRAC_MAPPED_RECEPTOR_CLASSIFICATION,
            MIN_MARGIN_RECEPTOR_CLASSIFICATION
        )
    }
}

fn check_feature_ref_and_config(
    feature_ref: Option<&FeatureReferenceFile>,
    vdj_receptor: Option<VdjReceptor>,
    gex_sample_def: Option<Vec<SampleDef>>,
    feature_config: Option<&FeatureConfig>,
) -> Result<Option<BeamMode>> {
    let has_antigen = match gex_sample_def {
        None => false,
        Some(sdef) => sdef
            .into_iter()
            .any(|s| s.library_type.unwrap_or_default() == LegacyLibraryType::AntigenCapture),
    };
    if has_antigen && (feature_ref.is_none() || vdj_receptor.is_none()) {
        panic!(
            "Expecting a valid feature reference and a specific VDJ receptor with Antigen Capture libraries!"
        )
    }
    let beam_mode = match (has_antigen, vdj_receptor.unwrap()) {
        (true, VdjReceptor::TR | VdjReceptor::TRGD) => Some(BeamMode::BeamT),
        (true, VdjReceptor::IG) => Some(BeamMode::BeamAB),
        (false, _) => None,
    };
    if let Some(beam) = beam_mode {
        if let Some(feature_config) = feature_config {
            if let Some(specificity_controls) = &feature_config.specificity_controls {
                if beam_mode == Some(BeamMode::BeamAB) && specificity_controls.has_mhc_allele_column
                {
                    bail!(
                        "[antigen-specificity] section of multi config CSV contains `mhc_allele` which is \
                    not supported by BCR Antigen Capture library (automatically detected based on VDJ chain detection)."
                    );
                }
                if beam_mode == Some(BeamMode::BeamT) && !specificity_controls.has_mhc_allele_column
                {
                    bail!(
                        "[antigen-specificity] section of multi config CSV does not contain `mhc_allele` which is \
                    required by TCR Antigen Capture library (automatically detected based on VDJ chain detection)."
                    );
                }
            }
        }

        feature_ref
            .unwrap()
            .read(feature_config)?
            .validate_beam_feature_ref(beam)?;
    }
    Ok(beam_mode)
}

// Note: We do not support auto-detection of G/D mode.
// User needs to specifically choose to run a G/D pipeline by selecting the correct chain_type
#[make_mro(volatile = strict)]
impl MartianMain for DetectVdjReceptor {
    type StageInputs = DetectVdjReceptorStageInputs;
    type StageOutputs = DetectVdjReceptorStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        // -----------------------------------------------------------------------------------------
        // A trivial case when we have no reference in denovo mode
        if args.vdj_reference_path.is_none() {
            return Ok(DetectVdjReceptorStageOutputs {
                receptor: None,
                beam_mode: None,
            });
        }

        // -----------------------------------------------------------------------------------------
        // Another trivial case when the receptor is explicitly specified. This is usually used as
        // a fallback when auto detection fails for some samples with low quality data.
        if let Some(force_receptor) = args.force_receptor {
            // "auto" and "all" will fail conversion here and fall-through
            if let Ok(force_receptor) = VdjReceptor::from_str(&force_receptor) {
                let beam_mode = check_feature_ref_and_config(
                    args.feature_reference.as_ref(),
                    Some(force_receptor),
                    args.gex_sample_def,
                    args.feature_config.as_ref(),
                )?;
                return Ok(DetectVdjReceptorStageOutputs {
                    receptor: Some(force_receptor),
                    beam_mode,
                });
            }
        }

        let vdj_ref = VdjReference::from_reference_folder(&args.vdj_reference_path.unwrap())?;
        let receptor_classifier = KmerClassify::<VdjReceptor>::new(&vdj_ref);

        let units = detect_chemistry_units(&args.vdj_sample_def, None, None)?;
        let resolution_text = match args.is_multi {
            true => {
                "Please specify the feature_types more specifically as either VDJ-T or VDJ-B. \
            If you are using our unsupported workflow for gamma/delta, then please use feature_type VDJ-T-GD."
            }
            false => {
                "Please check the input data and/or specify the chain via the --chain argument."
            }
        };
        let mut per_unit_receptors = Vec::with_capacity(units.len());
        for unit in &units {
            let mut stats = ClassificationStats::new();
            for fastq in &unit.fastqs {
                for read_pair in
                    ReadPairIter::from_fastq_files(fastq)?.take(MAX_READS_RECEPTOR_CLASSIFICATION)
                {
                    stats.total_reads += 1;
                    match receptor_classifier
                        .classify_rc(read_pair?.get(WhichRead::R2, ReadPart::Seq).unwrap())
                    {
                        Some(VdjReceptor::TR) => stats.tcr_reads += 1,
                        Some(VdjReceptor::TRGD) => stats.tcr_reads += 1,
                        Some(VdjReceptor::IG) => stats.ig_reads += 1,
                        None => {}
                    }
                }
                if stats.total_reads >= MAX_READS_RECEPTOR_CLASSIFICATION as i64 {
                    break;
                }
            }
            match stats.compatible_receptor() {
                Some(receptor) => per_unit_receptors.push(receptor),
                None => {
                    bail!(
                        "V(D)J Chain detection failed for {}.\n\n{}\n{}\n{}\n",
                        unit,
                        stats,
                        ClassificationStats::help_text(),
                        resolution_text,
                    );
                }
            }
            println!("{stats:?}");
        }

        let receptor_choices: TxHashSet<_> = per_unit_receptors.iter().collect();
        if receptor_choices.len() == 1 {
            let receptor = Some(*receptor_choices.into_iter().next().unwrap());
            let beam_mode = check_feature_ref_and_config(
                args.feature_reference.as_ref(),
                receptor,
                args.gex_sample_def,
                args.feature_config.as_ref(),
            )?;
            return Ok(DetectVdjReceptorStageOutputs {
                receptor,
                beam_mode,
            });
        }

        bail!(
            "Conflicting V(D)J chains detected.\n{}\nPlease check the input data.",
            zip(&units, &per_unit_receptors)
                .map(|(u, r)| format!("- {r} for {u}"))
                .join("\n")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::reference::feature_reference::SpecificityControls;
    use cr_types::sample_def::FastqMode;
    use dui_tests::stage_test::StageFailTest;
    use dui_tests::{stage_fail_dui_test, DuiTest};
    use std::collections::HashMap;
    use strum::IntoEnumIterator;

    #[test]
    fn test_denovo_no_ref() {
        let outs = DetectVdjReceptor
            .test_run_tmpdir(DetectVdjReceptorStageInputs {
                force_receptor: None,
                vdj_reference_path: None,
                feature_reference: None,
                gex_sample_def: None,
                vdj_sample_def: vec![],
                is_multi: false,
                feature_config: None,
            })
            .unwrap();
        assert_eq!(outs.receptor, None);
    }

    #[test]
    fn test_force_receptor() {
        for receptor in VdjReceptor::iter() {
            let outs = DetectVdjReceptor
                .test_run_tmpdir(DetectVdjReceptorStageInputs {
                    force_receptor: Some(receptor.to_string()),
                    vdj_reference_path: Some("/foo/bar".into()),
                    feature_reference: None,
                    gex_sample_def: None,
                    vdj_sample_def: vec![],
                    is_multi: false,
                    feature_config: None,
                })
                .unwrap();
            assert_eq!(outs.receptor, Some(receptor));
        }
    }

    #[test]
    fn test_vdj_ig() {
        let args = DetectVdjReceptorStageInputs {
            force_receptor: None,
            vdj_reference_path: Some(
                "../dui_tests/test_resources/reference/vdj_GRCm38_alts_ensembl-4.0.0/".into(),
            ),
            feature_reference: None,
            gex_sample_def: None,
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_mm_pbmc4_b_26x91"
                    .into(),
                sample_names: Some(vec!["vdj_v1_mm_pbmc4_b_26x91".into()]),
                ..Default::default()
            }],
            is_multi: false,
            feature_config: None,
        };
        let outs = DetectVdjReceptor.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.receptor, Some(VdjReceptor::IG));
    }

    #[test]
    fn test_vdj_tcr() {
        let args = DetectVdjReceptorStageInputs {
            force_receptor: None,
            vdj_reference_path: Some(
                "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0".into(),
            ),
            feature_reference: None,
            gex_sample_def: None,
            vdj_sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_t".into(),
                sample_names: Some(vec!["vdj_v1_hs_pbmc3_t".into()]),
                ..Default::default()
            }],
            is_multi: false,
            feature_config: None,
        };
        let outs = DetectVdjReceptor.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.receptor, Some(VdjReceptor::TR));
    }

    #[test]
    fn test_vdj_poor_mapping() -> Result<()> {
        // Uses 3'v2 data
        stage_fail_dui_test!(
            PoorMapping,
            description: "The input data is 3pv2, hence very few reads map to VDJ and the receptor \
            detection fails.",
            stage: DetectVdjReceptor,
            args: DetectVdjReceptorStageInputs {
                force_receptor: None,
                vdj_reference_path: Some(
                    "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0".into(),
                ),
                feature_reference: None,
                gex_sample_def: None,
                vdj_sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k"
                        .into(),
                    sample_names: Some(vec!["CR-120-01".into()]),
                    ..Default::default()
                }],
                is_multi: false,
                feature_config: None,
            },
        );
        Ok(())
    }

    #[test]
    fn test_vdj_poor_mapping_multi() -> Result<()> {
        // Multi has a different error message for poor mapping
        // Uses 3'v2 data
        stage_fail_dui_test!(
            PoorMappingMulti,
            description: "The input data is 3pv2, hence very few reads map to VDJ and the receptor \
            detection fails.",
            stage: DetectVdjReceptor,
            args: DetectVdjReceptorStageInputs {
                force_receptor: None,
                vdj_reference_path: Some(
                    "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0".into(),
                ),
                feature_reference: None,
                gex_sample_def: None,
                vdj_sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k"
                        .into(),
                    sample_names: Some(vec!["CR-120-01".into()]),
                    ..Default::default()
                }],
                is_multi: true,
                feature_config: None,
            },
        );
        Ok(())
    }

    #[test]
    fn test_vdj_mixed() -> Result<()> {
        stage_fail_dui_test!(
            MixedVdjReceptor,
            description: "The input data contains two sets of fastqs, one from TCR and the other \
            from IG.",
            stage: DetectVdjReceptor,
            args: DetectVdjReceptorStageInputs {
                force_receptor: None,
                vdj_reference_path: Some(
                    "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0".into(),
                ),
                feature_reference: None,
                gex_sample_def: None,
                vdj_sample_def: vec![
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_t"
                            .into(),
                        sample_names: Some(vec!["vdj_v1_hs_pbmc3_t".into()]),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_b"
                            .into(),
                        sample_names: Some(vec!["vdj_v1_hs_pbmc3_b".into()]),
                        ..Default::default()
                    },
                ],
                is_multi: false,
                feature_config: None,
            },
        );
        Ok(())
    }

    #[test]
    fn test_alleles_invalid_for_beam_ab() -> Result<()> {
        stage_fail_dui_test!(
            BeamAbMHCAllele,
            description: "The `mhc_allele` column is invalid for BCR Antigen Capture.",
            stage:DetectVdjReceptor,
            args: DetectVdjReceptorStageInputs {
                force_receptor: None,
                vdj_reference_path: Some(
                    "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0".into(),
                ),
                feature_reference: Some("test/multi/beam/feature_beam_t.csv".into()),
                gex_sample_def: Some(vec![SampleDef {
                    library_type: Some(LegacyLibraryType::AntigenCapture),
                    ..Default::default()
                }]),
                vdj_sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_b".into(),
                    sample_names: Some(vec!["vdj_v1_hs_pbmc3_b".into()]),
                    ..Default::default()
                }],
                is_multi: true,
                feature_config: Some(FeatureConfig {
                    specificity_controls: Some(SpecificityControls{
                        control_for_allele: HashMap::from([
                            ("a1".to_string(), "ag1".to_string()),
                        ]),
                        has_mhc_allele_column: true,}),
                    beam_mode: None,
                    functional_map: None,
                }),
            },
        );
        Ok(())
    }
}
