//! Martian stage MULTI_PREFLIGHT
//! Run pre-flight checks.
#![deny(missing_docs)]

use crate::preflight::{
    check_crispr_target_genes, check_resource_limits, check_vdj_inner_enrichment_primers,
    check_vdj_known_enrichment_primers, hostname, validate_target_panel,
};
use anyhow::{Context, Result, bail, ensure};
use base64::Engine;
use base64::prelude::BASE64_STANDARD;
use cr_types::AlignerParam;
use cr_types::chemistry::ChemistryName;
use cr_types::probe_set::merge_probe_set_csvs;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::types::FileOrBytes;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use multi::config::preflight::{
    SectionCtx, build_feature_reference_with_cmos, check_gex_reference, check_libraries,
    check_samples, check_vdj_reference,
};
use multi::config::{ChemistryParam, MultiConfigCsv, MultiConfigCsvFile, multiconst};
use multi::oscheck;
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

/// Martian stage MULTI_PREFLIGHT
pub struct MultiPreflight;

#[make_mro(volatile = strict, mem_gb = 2)]
impl MartianMain for MultiPreflight {
    type StageInputs = MultiPreflightInputs;
    type StageOutputs = MultiPreflightOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let hostname = hostname();
        check_resource_limits()?;
        if std::env::var("TENX_IGNORE_DEPRECATED_OS")
            .map(|s| s != "1")
            .unwrap_or(true)
        {
            oscheck(|warning| {
                rover.alarm(warning).unwrap();
            })?;
        }
        // make a new chunk per sample_def to do chemistry detection
        let cfg = match args.config {
            FileOrBytes {
                bytes: Some(ref bytes),
                file: None,
            } => {
                let bytes = BASE64_STANDARD.decode(bytes)?;
                MultiConfigCsv::from_reader(Cursor::new(bytes), "bytes")
            }
            FileOrBytes {
                bytes: None,
                file: Some(ref file),
            } => {
                if file.extension().and_then(|ext| ext.to_str()) != Some("csv") {
                    bail!("multi config file must have .csv extention")
                }
                MultiConfigCsvFile::from(file).read()
            }
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

            let transcriptome = gex
                .reference_path
                .as_deref()
                .map(|ref_path| {
                    check_gex_reference(
                        &SectionCtx {
                            section: "gene-expression",
                            field: "reference",
                        },
                        ref_path,
                        &hostname,
                    )
                })
                .transpose()?;

            let ref_info = gex
                .reference_path
                .as_deref()
                .map(ReferenceInfo::from_reference_path)
                .transpose()?;
            if let Some(ref_info) = &ref_info {
                ref_info.validate()?;
            }

            let probe_set_csv = match gex.probe_set() {
                [] => None,
                [probe_set] => Some(probe_set.clone()),
                probe_sets => {
                    // Ensure that the input probe CSVs themselves are valid before
                    // merging them.
                    for probe_set in probe_sets {
                        validate_target_panel(
                            &rover,
                            transcriptome.as_ref(),
                            ref_info.as_ref(),
                            probe_set,
                            gex.targeting_method().unwrap(),
                            false, // allow merging very small probe sets
                        )
                        .with_context(|| probe_set.to_string_lossy().to_string())?;
                    }

                    let probe_set_csv: TargetSetFile = rover.make_path("combined_probe_set");
                    merge_probe_set_csvs(
                        probe_sets,
                        probe_set_csv.buf_writer()?,
                        gex.reference_path.as_deref(),
                    )?;
                    Some(probe_set_csv)
                }
            };

            let target_genes = probe_set_csv
                .as_ref()
                .map(|probe_set_csv| {
                    validate_target_panel(
                        &rover,
                        transcriptome.as_ref(),
                        ref_info.as_ref(),
                        probe_set_csv,
                        gex.targeting_method().unwrap(),
                        true,
                    )
                })
                .transpose()?;
            (transcriptome, target_genes)
        } else {
            (None, None)
        };

        let aligner = gene_expression.and_then(|gex| gex.aligner.as_ref());
        if !args.is_pd && aligner == Some(&AlignerParam::Minimap2) {
            bail!("minimap2 is an unsupported aligner");
        }

        let (feature_reference, _tenx_cmos) =
            build_feature_reference_with_cmos(&cfg, args.is_pd, &hostname)?;

        if let (Some(transcriptome), Some(feature_reference)) = (&transcriptome, &feature_reference)
        {
            check_crispr_target_genes(transcriptome, feature_reference, target_genes.as_ref())?;
        }

        if let Some(ref vdj) = cfg.vdj {
            if let Some(ref reference_path) = vdj.reference_path {
                let vdj_ref = check_vdj_reference(
                    &SectionCtx {
                        section: "vdj",
                        field: "reference",
                    },
                    reference_path,
                    &hostname,
                )?;
                if let Some(ref primers) = vdj.inner_enrichment_primers {
                    check_vdj_inner_enrichment_primers(
                        Some(reference_path),
                        Some(&vdj_ref),
                        primers,
                    )?;
                } else {
                    check_vdj_known_enrichment_primers(reference_path, &vdj_ref)?;
                }
            }
            // denovo mode
            if let Some(ref primers) = vdj.inner_enrichment_primers {
                check_vdj_inner_enrichment_primers(None, None, primers)?;
            }
            if vdj.denovo == Some(true) {
                ensure!(
                    cfg.samples.is_none() || args.is_pd,
                    "Multiplexing with denovo mode is unsupported.",
                );
            }
        }

        check_libraries(&cfg, feature_reference.clone(), args.is_pd, &hostname)?;
        check_samples(&cfg, feature_reference.clone(), args.is_pd)?;
        // TODO(CELLRANGER-7837) remove this check once PD support is removed
        if !args.is_pd
            && cfg
                .chemistry_specs()?
                .values()
                .any(|chem| chem.refined() == Some(ChemistryName::ThreePrimeV3LT))
        {
            bail!(
                "The chemistry SC3Pv3LT (Single Cell 3'v3 LT) is no longer supported. To analyze this data, use Cell Ranger 7.2 or earlier."
            );
        }

        Ok(MultiPreflightOutputs {})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glob::glob;
    use insta::{assert_debug_snapshot, assert_snapshot};
    use itertools::Itertools;
    use multi::config::preflight::check_feature_reference;
    use std::io::Write;
    use std::path::Path;
    use tempfile::NamedTempFile;
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
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_gex_multi_dup1() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_dup1.csv");
        assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_gex_multi_dup2() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_dup2.csv");
        assert_snapshot!(&outs.unwrap_err());
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
        assert_debug_snapshot!(outs);
        Ok(())
    }

    #[test]
    fn test_gex_multi_no_feature() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_no_feature.csv");
        assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_gex_multi_only() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_only.csv");
        assert_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_gex_multi_unsupp() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/gex_multi_unsupp.csv");
        assert_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_antigen_specificity_section() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/beamab_vdj_gex_invalid_antigen_spec1.csv");
        assert_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_antigen_specificity_section_allele_not_found() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/beamt_vdj_gex_invalid_antigen_spec1.csv");
        assert_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_invalid_beamab_with_allele() {
        let outs = test_run_stage("test/multi/beamab_vdj_gex_with_allele.csv");
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_invalid_beamab_config_w_allele_featureref_no_allele() {
        let outs = test_run_stage("test/multi/beamab_config_w_allele_featureref_no_allele.csv");
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamab_config_no_allele_featureref_w_allele() {
        let outs = test_run_stage("test/multi/beamab_config_no_allele_featureref_w_allele.csv");
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamt_config_no_allele_featureref_w_allele() {
        let outs = test_run_stage("test/multi/beamt_config_no_allele_featureref_w_allele.csv");
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamt_config_no_allele_featureref_no_allele() {
        let outs = test_run_stage("test/multi/beamt_config_no_allele_featureref_no_allele.csv");
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn beamt_config_incompatible_allele() {
        let outs = test_run_stage("test/multi/beamt_config_incompatible_allele.csv");
        assert_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_invalid_feature_functional_map_section() {
        if !refdata_available() {
            return;
        }
        let outs = test_run_stage("test/multi/beamab_vdj_gex_invalid_functional_map.csv");
        assert_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_mfrp_ab_multi() {
        if !refdata_available() {
            return;
        }
        assert!(test_run_stage_is_pd("test/multi/mfrp_ab_multi.csv", true).is_ok());
        let outs = test_run_stage("test/multi/mfrp_ab_multi.csv");
        assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_mfrp_ab_multi_invalid_probe_barcode() {
        if !refdata_available() {
            return;
        }
        let outs =
            test_run_stage("test/multi/invalid_csvs/mfrp_ab_multi_invalid_probe_barcode.csv");
        assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_mfrp_ab_multi_unknown_probe_barcode() {
        if !refdata_available() {
            return;
        }
        let outs =
            test_run_stage("test/multi/invalid_csvs/mfrp_ab_multi_unknown_probe_barcode.csv");
        assert_debug_snapshot!(&outs);
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
        assert_debug_snapshot!(&outs);
    }

    #[test]
    fn test_overhang_multiplexing() {
        if !refdata_available() {
            return;
        }

        let outs = test_run_stage("test/multi/invalid_csvs/invalid_overhang_ids.csv");
        assert_debug_snapshot!(&outs.unwrap_err().to_string());
    }

    #[test]
    fn test_misspelled_hash_tag() {
        let feature_ref = r"id,name,read,pattern,sequence,feature_type
TotalSeqB_Hashtag_1,TotalB_HashTag_1,R2,^NNNNNNNNNN(BC)NNNNNNNNN,GTCAACTCTTTAGCG,FEATURETEST
TotalSeqB_Hashtag_2,TotalB_HashTag_2,R2,^NNNNNNNNNN(BC)NNNNNNNNN,TGATGGCCTATTGGG,FEATURETEST
TotalSeqB_Hashtag_3,TotalB_HashTag_3,R2,^NNNNNNNNNN(BC)NNNNNNNNN,TTCCGCCTCTCTTTG,FEATURETEST
TotalSeqB_Hashtag_4,TotalB_HashTag_4,R2,^NNNNNNNNNN(BC)NNNNNNNNN,AGTAAGTTCAGCGTA,FEATURETEST";
        let mut temp_file = NamedTempFile::new().expect("Failed to create temp file");
        write!(temp_file, "{feature_ref}").expect("Failed to write to temp file");

        let csv = format!(
            r"[gene-expression]
ref,/blah/reporters
create-bam,false

[feature]
ref,{}

[libraries]
fastq_path,sample_indices,lanes,physical_library_id,feature_types,subsample_rate,chemistry
blah/fastq_path,CGTGACATGC,1|2|3|4|5|6|7|8,1384285,Gene Expression,,
blah/fastq_path,AGATGAGAAT,1|2|3|4|5|6|7|8,1384321,Multiplexing Capture,,

[samples]
sample_id,cmo_ids,description
Samp1,TotalB_HashTag_1,Samp1,,,,
Samp2,TotalB_HashTag_2,Samp2,,,,
Samp3,TotalB_HashTag_3,Samp3,,,,
Samp4,TotalB_HashTag_4,Samp4,,,,",
            temp_file.path().display()
        );
        let mut csv_file = NamedTempFile::new().expect("Failed to make temp file");
        write!(csv_file, "{csv}").expect("Failed write");

        let cfg = MultiConfigCsv::from_csv(csv_file.path()).unwrap();
        let (feature_reference, _tenx_cmos) =
            build_feature_reference_with_cmos(&cfg, true, &hostname()).unwrap();

        let error_msg = r"Unknown cmo_ids ('TotalB_HashTag_1') provided for sample 'Samp1', please ensure you are either using valid 10x CMO IDs or are providing the correct [gene-expression] cmo-set.

Valid IDs are currently:
TotalSeqB_Hashtag_1
TotalSeqB_Hashtag_2
TotalSeqB_Hashtag_3
TotalSeqB_Hashtag_4

Did you perhaps mean any of the following? (Input -> Intended Value)?
TotalB_HashTag_1 -> TotalSeqB_Hashtag_1
";

        let result = check_samples(&cfg, feature_reference, true);
        match result {
            Err(e) => assert_eq!(e.to_string(), error_msg),
            Ok(()) => panic!("Failed to flag invalid sequence."),
        };
    }
    #[test]
    fn test_invalid_hashtag_ids() {
        if !refdata_available() {
            return;
        }

        let outs = test_run_stage("test/multi/invalid_csvs/invalid_hashtag_ids.csv");
        assert_debug_snapshot!(&outs.unwrap_err().to_string());
    }
}
