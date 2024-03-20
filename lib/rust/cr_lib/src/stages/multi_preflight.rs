//! Martian stage MULTI_PREFLIGHT
//! Run pre-flight checks.

use crate::preflight::{
    check_crispr_target_genes, check_resource_limits, check_target_panel,
    check_vdj_inner_enrichment_primers, check_vdj_known_enrichment_primers, hostname,
};
use anyhow::{bail, ensure, Result};
use cr_types::chemistry::ChemistryName;
use cr_types::reference::feature_reference::FeatureType;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::types::FileOrBytes;
use cr_types::FeatureBarcodeType;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use multi::config::preflight::{
    build_feature_reference_with_cmos, check_gex_reference, check_libraries, check_samples,
    check_vdj_reference, SectionCtx,
};
use multi::config::{multiconst, ChemistryParam, MultiConfigCsv, MultiConfigCsvFile};
use parameters_toml::max_multiplexing_tags;
use serde::{Deserialize, Serialize};
use std::io::Cursor;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MultiPreflightInputs {
    pub config: FileOrBytes,
    pub is_pd: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct MultiPreflightOutputs {}

// This is our stage struct
pub struct MultiPreflight;

#[make_mro(volatile = strict)]
impl MartianMain for MultiPreflight {
    type StageInputs = MultiPreflightInputs;
    type StageOutputs = MultiPreflightOutputs;

    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let hostname = hostname();
        check_resource_limits()?;
        // make a new chunk per sample_def to do chemistry detection
        let cfg = match args.config {
            FileOrBytes {
                bytes: Some(ref bytes),
                file: None,
            } => {
                let bytes = base64::decode(bytes)?;
                MultiConfigCsv::from_reader(Cursor::new(bytes), "bytes")
            }
            FileOrBytes {
                bytes: None,
                file: Some(ref file),
            } => MultiConfigCsvFile::from(file).read(),
            FileOrBytes {
                bytes: None,
                file: None,
            }
            | FileOrBytes {
                bytes: Some(_),
                file: Some(_),
            } => bail!("exactly one of either config file or config bytes must be provided"),
        }?;
        let gene_expression = cfg.gene_expression.as_ref();
        let (transcriptome, target_genes) = if let Some(gex) = gene_expression {
            ensure!(
                gex.chemistry != Some(ChemistryParam::Custom) || args.is_pd,
                "Unknown chemistry {} in [{}] section.",
                ChemistryName::Custom,
                multiconst::GENE_EXPRESSION,
            );

            let transcriptome = check_gex_reference(
                &SectionCtx {
                    section: "gene-expression",
                    field: "reference",
                },
                &gex.reference_path,
                &hostname,
            )?;
            let ref_info = ReferenceInfo::from_reference_path(&gex.reference_path)?;
            ref_info.validate()?;
            let genomes = ref_info.genomes;
            let target_genes = if let (Some(probe_set), Some(targeting_method)) =
                (gex.probe_set(), gex.targeting_method())
            {
                Some(check_target_panel(
                    &transcriptome,
                    genomes,
                    probe_set,
                    targeting_method,
                    args.is_pd,
                )?)
            } else {
                None
            };
            (Some(transcriptome), target_genes)
        } else {
            (None, None)
        };
        let max_multiplexing_tags = *max_multiplexing_tags()?;

        let (feature_reference, _tenx_cmos) =
            build_feature_reference_with_cmos(&cfg, args.is_pd, &hostname, max_multiplexing_tags)?;

        if let (Some(transcriptome), Some(feature_reference)) =
            (transcriptome.as_ref(), feature_reference.as_ref())
        {
            check_crispr_target_genes(transcriptome, feature_reference, target_genes.as_ref())?;
        }

        if let Some(ref vdj) = cfg.vdj {
            let vdj_ref = check_vdj_reference(
                &SectionCtx {
                    section: "vdj",
                    field: "reference",
                },
                &vdj.reference_path,
                &hostname,
            )?;
            if let Some(ref primers) = vdj.inner_enrichment_primers {
                check_vdj_inner_enrichment_primers(&vdj.reference_path, &vdj_ref, primers)?;
            } else {
                check_vdj_known_enrichment_primers(&vdj.reference_path, &vdj_ref)?;
            }
        }

        check_libraries(
            &cfg,
            feature_reference.clone(),
            args.is_pd,
            &hostname,
            max_multiplexing_tags,
        )?;
        check_samples(
            &cfg,
            feature_reference.clone(),
            args.is_pd,
            max_multiplexing_tags,
        )?;
        // TODO(CELLRANGER-7837) remove this check once PD support is removed
        if !args.is_pd
            && cfg
                .chemistry_specs()?
                .values()
                .any(|chem| chem.refined() == Some(ChemistryName::ThreePrimeV3LT))
        {
            bail!("The chemistry SC3Pv3LT (Single Cell 3'v3 LT) is no longer supported. To analyze this data, use Cell Ranger 7.2 or earlier.");
        }

        // Check for CRISPR guide capture libraries and v4 of 3pGEX, which are incompatible.
        if let Some(feature_reference) = feature_reference.as_ref() {
            for feat in &feature_reference.feature_defs {
                if feat.feature_type == FeatureType::Barcode(FeatureBarcodeType::Crispr)
                    && cfg.chemistry_specs()?.values().any(|chem| {
                        chem.refined() == Some(ChemistryName::ThreePrimeV4)
                            || chem.refined() == Some(ChemistryName::ThreePrimeV4OH)
                    })
                {
                    bail!("The chemistry SC3Pv4 (Single Cell 3'v4) is not supported with CRISPR Guide Capture libraries.")
                }
            }
        }

        Ok(MultiPreflightOutputs {})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glob::glob;
    use itertools::Itertools;
    use multi::config::preflight::check_feature_reference;
    use std::path::Path;
    use test_refdata::refdata_available;

    fn test_run_stage_is_pd(
        csv_file: impl AsRef<Path>,
        is_pd: bool,
    ) -> Result<MultiPreflightOutputs> {
        MultiPreflight.test_run_tmpdir(MultiPreflightInputs {
            config: FileOrBytes {
                bytes: None,
                file: Some(csv_file.as_ref().into()),
            },
            is_pd,
        })
    }

    fn test_run_stage(csv_file: impl AsRef<Path>) -> Result<MultiPreflightOutputs> {
        test_run_stage_is_pd(csv_file, false)
    }

    #[test]
    fn test_gex_missing_gex_section() {
        let outs = test_run_stage("test/multi/invalid_csvs/gex_missing_gex_section.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_3pgexv4_crispr_incompatible() {
        if !refdata_available() {
            return;
        }

        let outs = test_run_stage("test/multi/3pgexv4_crispr_incompatible.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_gex_multi_dup1() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_dup1.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_gex_multi_dup2() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_dup2.csv");
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_feature_references() -> Result<()> {
        let mut feat_csvs: Vec<_> = glob("test/feature/invalid_csvs/*.csv")?.try_collect()?;
        feat_csvs.sort();
        let outs: Vec<_> = feat_csvs
            .iter()
            .map(|feat_csv| {
                (
                    feat_csv,
                    check_feature_reference(&"feature reference", feat_csv, "Fake Host", None)
                        .expect_err(&format!(
                            "Expected error on invalid feature reference csv {feat_csv:?}"
                        ))
                        .to_string(),
                )
            })
            .collect();
        insta::assert_debug_snapshot!(outs);
        Ok(())
    }

    #[test]
    fn test_gex_multi_no_feature() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_no_feature.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_gex_multi_only() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_only.csv");
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_gex_multi_unsupp() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_unsupp.csv");
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_antigen_specificity_section() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/beamab_vdj_gex_invalid_antigen_spec1.csv");
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_antigen_specificity_section_allele_not_found() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/beamt_vdj_gex_invalid_antigen_spec1.csv");
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_beamab_with_allele() {
        let outs = test_run_stage("test/multi/beamab_vdj_gex_with_allele.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_invalid_beamab_config_w_allele_featureref_no_allele() {
        let outs = test_run_stage("test/multi/beamab_config_w_allele_featureref_no_allele.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamab_config_no_allele_featureref_w_allele() {
        let outs = test_run_stage("test/multi/beamab_config_no_allele_featureref_w_allele.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamt_config_no_allele_featureref_w_allele() {
        let outs = test_run_stage("test/multi/beamt_config_no_allele_featureref_w_allele.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamt_config_no_allele_featureref_no_allele() {
        let outs = test_run_stage("test/multi/beamt_config_no_allele_featureref_no_allele.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamt_config_incompatible_allele() {
        let outs = test_run_stage("test/multi/beamt_config_incompatible_allele.csv");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_invalid_feature_functional_map_section() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/beamab_vdj_gex_invalid_functional_map.csv");
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_mfrp_ab_multi() {
        if !refdata_available() {
            return;
        }
        assert!(test_run_stage_is_pd("test/multi/mfrp_ab_multi.csv", true).is_ok());
        let outs = test_run_stage("test/multi/mfrp_ab_multi.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_mfrp_ab_multi_invalid_probe_barcode() {
        if !refdata_available() {
            return;
        }
        let outs =
            test_run_stage("test/multi/invalid_csvs/mfrp_ab_multi_invalid_probe_barcode.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_mfrp_ab_multi_unknown_probe_barcode() {
        if !refdata_available() {
            return;
        }
        let outs =
            test_run_stage("test/multi/invalid_csvs/mfrp_ab_multi_unknown_probe_barcode.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_mfrp_ab_multi_without_ab_library() {
        if !refdata_available() {
            return;
        }
        assert!(
            test_run_stage_is_pd("test/multi/mfrp_ab_multi_without_ab_library.csv", true).is_ok()
        );
        let outs = test_run_stage("test/multi/mfrp_ab_multi_without_ab_library.csv");
        insta::assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_overhang_multiplexing() {
        if !refdata_available() {
            return;
        }

        let outs = test_run_stage("test/multi/invalid_csvs/invalid_overhang_ids.csv");
        insta::assert_debug_snapshot!(&outs.unwrap_err().to_string());
    }
}
