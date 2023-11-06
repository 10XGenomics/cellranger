//! Martian stage PARSE_MULTI_CONFIG

use crate::preflight::hostname;
use anyhow::{anyhow, bail, Result};
use cr_types::chemistry::{ChemistryDef, IndexScheme};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use cr_types::rna_read::LegacyLibraryType;
use cr_types::sample_def::SampleDef;
use cr_types::types::FileOrBytes;
use cr_types::{AlignerParam, TargetingMethod};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct, MartianType};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use metric::TxHashMap;
use multi::barcode_sample_assignment::SampleAssignmentCsv;
use multi::config::preflight::build_feature_reference_with_cmos;
use multi::config::{create_feature_config, FeatureType, MultiConfigCsvFile};
use parameters_toml::max_multiplexing_tags;
use serde::{Deserialize, Serialize};
use serde_json::Value as JValue;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::PathBuf;

#[derive(Clone, Default, Serialize, Deserialize, MartianType)]
pub struct MultiParams {
    #[serde(default)]
    initial_reads: Option<usize>,
    #[serde(default)]
    subsample_rate: Option<f64>,
    #[serde(default)]
    primers: Vec<Primers>,
    #[serde(default)]
    index_scheme: Option<IndexScheme>,
    #[serde(default)]
    barcode_whitelist: Option<String>,
    #[serde(default)]
    special_genomic_regions: Option<Vec<String>>,
}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct ParseMultiConfigStageInputs {
    pub sample_id: String,
    pub sample_desc: String,
    pub config: FileOrBytes,
    pub config_hash: Option<String>,
    pub params: Option<MultiParams>,
    pub is_pd: bool,
}

/// Copy of `struct CellCallingParam`, defined in _basic_sc_rna_counter_stages.mro
/// The `per_gem_well` parameter is applied at the library-level (CMO-multiplexing, standard GEX).
/// The `per_sample` field is only valid in case of RTL multiplexed inputs.
/// The fields `per_gem_well` and `per_sample` are mutually exclusive.
/// Martian will unify this copy of the type and the copy
/// in _basic_sc_rna_counter_stages.mro
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct CellCallingParam {
    pub per_gem_well: Option<usize>,
    pub per_sample: Option<TxHashMap<String, Option<usize>>>,
}

impl CellCallingParam {
    // Return total values for parameters that can be input either at the per_gem_well or per_sample level
    pub fn sum(&self) -> Option<usize> {
        match (self.per_gem_well, self.per_sample.as_ref()) {
            (Some(x), Some(xs)) => {
                assert!(xs.values().all(Option::is_none));
                Some(x)
            }
            (Some(x), None) => Some(x),
            (None, Some(xs)) => xs.values().fold(Some(0), |acc, x| match (acc, x) {
                (Some(acc), Some(x)) => Some(acc + x),
                _ => None,
            }),
            (None, None) => None,
        }
    }
}

/// Copy of `struct CellCalling`, defined in _basic_sc_rna_counter_stages.mro
/// Carries options to customize the cell calling mode.
/// Martian will unify this copy of the type and the copy
/// in _basic_sc_rna_counter_stages.mro
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct CellCalling {
    pub recovered_cells: CellCallingParam,
    pub force_cells: CellCallingParam,
    pub emptydrops_minimum_umis: CellCallingParam,
    pub cell_barcodes: Option<JsonFile<()>>,
    pub override_mode: Option<String>,
    pub override_library_types: Option<Vec<String>>,
    pub disable_ab_aggregate_detection: bool,
    pub disable_high_occupancy_gem_detection: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct CommonInputs {
    pub sample_id: String,
    pub sample_desc: String,
    pub multi_config_sha: String,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct BarcodeAssignments {
    pub sample_barcodes: Option<JsonFile<()>>,
    pub non_singlet_barcodes: Option<JsonFile<()>>,
    pub cells_per_tag: Option<JsonFile<()>>,
}

pub type Primers = JValue;
pub type GeneticDemuxParams = JValue;

/// CountInputs struct
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct CountInputs {
    pub sample_def: Vec<SampleDef>,
    pub chemistry: Option<String>,
    pub custom_chemistry_def: Option<ChemistryDef>,
    pub reference_path: PathBuf,
    pub gene_index: JsonFile<()>,
    #[mro_type = "map[]"]
    pub primers: Vec<Primers>,
    pub cell_calling_config: CellCalling,
    pub subsample_rate: Option<f64>,
    pub initial_reads: Option<usize>,
    pub primer_initial_reads: Option<usize>,
    pub special_genomic_regions: Option<Vec<String>>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub trim_polya_min_score: Option<i64>,
    pub trim_tso_min_score: Option<i64>,
    pub no_secondary_analysis: bool,
    pub no_target_umi_filter: bool,
    pub filter_probes: Option<bool>,
    pub feature_reference: Option<FeatureReferenceFile>,
    pub include_exons: bool,
    pub include_introns: bool,
    pub targeting_method: Option<TargetingMethod>,
    pub aligner: Option<AlignerParam>,
    #[mro_type = "map"]
    pub genetic_demux_params: Option<GeneticDemuxParams>,
    pub throughput: Option<String>,
    pub check_library_compatibility: bool,
    pub no_bam: bool,
    pub force_sample_barcodes: BarcodeAssignments,
    pub tenx_cmos: Option<bool>,
    pub min_assignment_confidence: Option<f64>,
    pub annotations: Option<Vec<CsvFile<()>>>,
}

/// General VdjInputs which are not chain specific
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct VdjGenInputs {
    pub reference_path: Option<PathBuf>,
    pub vdj_reference_path: Option<PathBuf>,
    pub filter_flags: VdjFilterFlags,
}

#[derive(Default, Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct VdjFilterFlags {
    pub multiplet_filter: Option<bool>,
    pub shared_contig_filter: Option<bool>,
    pub umi_baseline_filter: Option<bool>,
}

/// VDJInputs struct
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct VdjInputs {
    pub sample_def: Vec<SampleDef>,
    pub chemistry: Option<String>,
    pub custom_chemistry_def: Option<ChemistryDef>,
    #[mro_type = "map[]"]
    pub primers: Vec<Primers>,
    pub subsample_rate: Option<f64>,
    pub initial_reads: Option<usize>,
    pub primer_initial_reads: Option<usize>,
    pub special_genomic_regions: Option<Vec<String>>,
    pub denovo: bool,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub ground_truth_clonotype_path: Option<PathBuf>,
    pub inner_enrichment_primers: Option<PathBuf>,
    pub chain_type: Option<String>,
    pub physical_library_id: Option<String>,
    pub r2_revcomp: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct BasicPipelineConfig {
    pub disable_count: bool,
    pub disable_vdj: bool,
    /// boolean to disable stages that are only needed in the multi pipeline
    pub disable_multi: bool,
    /// boolean to disable stages that are only needed when count libraries are
    /// present in the multi pipeline
    pub disable_multi_count: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct ParseMultiConfigStageOutputs {
    pub common_input: CommonInputs,
    pub count_input: Option<CountInputs>,
    pub vdj_inputs: Vec<VdjInputs>, // or just JSON w/ outer array?
    pub vdj_gen_inputs: Option<VdjGenInputs>,
    pub basic_config: BasicPipelineConfig,
    pub config_file: MultiConfigCsvFile,
    pub feature_config: Option<FeatureConfig>,
    #[mro_retain]
    pub feature_ref: Option<FeatureReferenceFile>,
    #[mro_retain]
    pub cell_barcodes: Option<JsonFile<()>>,
    #[mro_retain]
    pub sample_barcodes: Option<JsonFile<()>>,
    #[mro_retain]
    pub non_singlet_barcodes: Option<JsonFile<()>>,
    #[mro_retain]
    pub cells_per_tag: Option<JsonFile<()>>,
    #[mro_retain]
    pub barcode_sample_assignments: Option<CsvFile<()>>,
}

// This is our stage struct
pub struct ParseMultiConfig;

#[make_mro(mem_gb = 6, volatile = strict)]
impl MartianMain for ParseMultiConfig {
    type StageInputs = ParseMultiConfigStageInputs;
    type StageOutputs = ParseMultiConfigStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // make a new chunk per sample_def to do chemistry detection
        let config_bytes = match args.config {
            FileOrBytes {
                file: None,
                bytes: Some(ref bytes),
            } => base64::decode(bytes)?,
            FileOrBytes {
                bytes: None,
                file: Some(ref file),
            } => {
                let mut bytes = vec![];
                let mut handle = File::open(file)?;
                let _ = handle.read_to_end(&mut bytes)?;
                bytes
            }
            _ => {
                bail!("exactly one of config file or config bytes must be provided");
            }
        };
        if !cfg!(test) {
            std::io::stdout().write_all(&config_bytes[..])?;
        }
        let config_file = {
            let file: MultiConfigCsvFile = rover.make_path("config");
            let mut writer = file.buf_writer()?;
            writer.write_all(&config_bytes[..])?;
            writer.flush()?;
            file
        };
        let cfg = config_file.read()?;

        let MultiParams {
            initial_reads,
            subsample_rate,
            primers,
            index_scheme,
            barcode_whitelist,
            special_genomic_regions,
        } = args.params.unwrap_or_default();

        let (
            count_input,
            feature_ref,
            cell_barcodes,
            sample_barcodes,
            non_singlet_barcodes,
            cells_per_tag,
            barcode_sample_assignments,
        ) = {
            // TODO: need to march over distinct GemWells --
            //   so I need to provide an API for that over on MultiConfigCsv,
            //   which is good, because I probably want it for validating [gem-wells]
            //   misc thought, do deep validation of libraries and gem-wells settings,
            //   e.g. vdj_force_cells vs (gex_)force_cells
            let mut sample_def = cfg
                .libraries
                .0
                .iter()
                .fold(Ok(vec![]), |acc: Result<_>, x| {
                    let mut acc = acc?;
                    if x.is_count()? {
                        acc.push(x.to_sample_def(&cfg)?);
                    }
                    Ok(acc)
                })?;
            if sample_def.is_empty() {
                (None, None, None, None, None, None, None)
            } else {
                let gex = cfg.gene_expression.as_ref().ok_or_else(
                    #[cold]
                    || {
                        anyhow!(
                            "[gene-expression] section with a path to the transcriptome reference \
                             is a required input (even if only antibody-data is present)."
                        )
                    },
                )?;
                let feature = cfg.feature.as_ref();
                for sample in &mut sample_def {
                    if let Some(LegacyLibraryType::GeneExpression) = sample.library_type {
                        sample.r1_length = gex.r1_length;
                        sample.r2_length = gex.r2_length;
                    } else {
                        sample.r1_length = feature.and_then(|f| f.r1_length);
                        sample.r2_length = feature.and_then(|f| f.r2_length);
                    }
                }
                // write the gene_index json if we need to
                let gene_index: JsonFile<()> = rover.make_path("gene_index");
                transcriptome::python_gene_index::write_gene_index(
                    &gex.reference_path,
                    gene_index.as_ref(),
                )?;
                let chemistry = gex.chemistry.name().unwrap_or("auto").to_string();
                let chemistry = Some(chemistry);
                /* TODO: when we divvy up CountInputs into a vec by gem-well,
                //   thread in gem-wells params for force_cells, etc
                let force_cells = cfg
                    .gem_wells
                    .as_ref()
                    .map(|gw| gw.0.get(&GemWell(1)).map(|p| p.force_cells))
                    .flatten()
                    .flatten();
                */

                let (
                    cell_barcodes,
                    sample_barcodes,
                    non_singlet_barcodes,
                    cells_per_tag,
                    barcode_sample_assignments,
                ) = if let Some(bsa) = gex.barcode_sample_assignment.as_ref() {
                    let sample_assignments = SampleAssignmentCsv::from_file(bsa, &cfg)?;
                    let cell_barcodes: JsonFile<()> = rover.make_path("cell_barcodes");
                    let sample_barcodes: JsonFile<()> = rover.make_path("sample_barcodes");
                    let non_singlet_barcodes: JsonFile<()> =
                        rover.make_path("non_singlet_barcodes");
                    let cells_per_tag: JsonFile<()> = rover.make_path("cells_per_tag");
                    let barcode_sample_assignment_csv: CsvFile<()> =
                        rover.make_path("barcode_sample_assignment");

                    std::fs::copy(bsa, &barcode_sample_assignment_csv)?;
                    sample_assignments.to_cell_barcodes_json(&cell_barcodes)?;
                    sample_assignments.to_sample_barcodes_json(&sample_barcodes)?;
                    sample_assignments.to_non_singlet_barcodes_json(&non_singlet_barcodes)?;
                    let cells_per_tag =
                        if sample_assignments.to_cells_per_tag_json(&cells_per_tag)? {
                            Some(cells_per_tag)
                        } else {
                            None
                        };
                    (
                        Some(cell_barcodes),
                        Some(sample_barcodes),
                        Some(non_singlet_barcodes),
                        cells_per_tag,
                        Some(barcode_sample_assignment_csv),
                    )
                } else {
                    (None, None, None, None, None)
                };

                let cell_calling_config = CellCalling {
                    force_cells: CellCallingParam {
                        per_gem_well: gex.force_cells,
                        per_sample: cfg
                            .samples
                            .as_ref()
                            .map(multi::config::SamplesCsv::get_force_cells),
                    },
                    recovered_cells: CellCallingParam {
                        per_gem_well: gex.expect_cells,
                        per_sample: cfg
                            .samples
                            .as_ref()
                            .map(multi::config::SamplesCsv::get_expect_cells),
                    },
                    emptydrops_minimum_umis: CellCallingParam {
                        per_gem_well: gex.emptydrops_minimum_umis,
                        per_sample: cfg
                            .samples
                            .as_ref()
                            .map(multi::config::SamplesCsv::get_emptydrops_minimum_umis),
                    },
                    cell_barcodes: cell_barcodes.clone(),
                    override_mode: None,
                    override_library_types: None,
                    disable_ab_aggregate_detection: false,
                    disable_high_occupancy_gem_detection: false,
                };

                let (feature_reference_file, tenx_cmos) = match build_feature_reference_with_cmos(
                    &cfg,
                    args.is_pd,
                    &hostname(),
                    *max_multiplexing_tags()?,
                )? {
                    (Some(_), None) => {
                        // Not using CMO multiplexing; use the original feature reference file.
                        let src = cfg
                            .feature
                            .as_ref()
                            .unwrap()
                            .reference_path
                            .as_ref()
                            .unwrap();
                        let dst: FeatureReferenceFile = rover.make_path("feature_ref");
                        if std::fs::hard_link(src, &dst).is_err() {
                            std::fs::copy(src, &dst)?;
                        }
                        (Some(dst), None)
                    }
                    (Some(feature_ref), Some(tenx_cmos)) => {
                        // Using CMO, need to write out the constructed feature reference file.
                        let dst: FeatureReferenceFile = rover.make_path("feature_ref");
                        let mut w = BufWriter::new(File::create(&dst)?);
                        feature_ref.to_csv(&mut w)?;
                        (Some(dst), Some(tenx_cmos))
                    }
                    (None, _) => (None, None),
                };

                let custom_chemistry_def = index_scheme.map(|is| {
                    is.to_chemistry_def(
                        barcode_whitelist
                            .expect("Barcode set is required for custom chemistry")
                            .as_ref(),
                    )
                });

                let count_input = Some(CountInputs {
                    sample_def,
                    chemistry,
                    custom_chemistry_def,
                    reference_path: PathBuf::from(&gex.reference_path),
                    gene_index,
                    primers: primers.clone(),
                    cell_calling_config,
                    subsample_rate,
                    initial_reads,
                    primer_initial_reads: Some(1000000),
                    special_genomic_regions: special_genomic_regions.clone(),
                    r1_length: None,
                    r2_length: None,
                    trim_polya_min_score: None,
                    trim_tso_min_score: None,
                    no_secondary_analysis: gex.no_secondary_analysis,
                    no_target_umi_filter: false,
                    filter_probes: gex.filter_probes,
                    feature_reference: feature_reference_file.clone(),
                    include_exons: true,
                    include_introns: gex.include_introns,
                    targeting_method: gex.targeting_method(),
                    aligner: gex.aligner,
                    genetic_demux_params: None,
                    throughput: None,
                    check_library_compatibility: gex.check_library_compatibility,
                    no_bam: gex.no_bam,
                    force_sample_barcodes: BarcodeAssignments {
                        sample_barcodes: sample_barcodes.clone(),
                        non_singlet_barcodes: non_singlet_barcodes.clone(),
                        cells_per_tag: cells_per_tag.clone(),
                    },
                    tenx_cmos,
                    min_assignment_confidence: gex.min_assignment_confidence,
                    annotations: None,
                });

                (
                    count_input,
                    feature_reference_file,
                    cell_barcodes,
                    sample_barcodes,
                    non_singlet_barcodes,
                    cells_per_tag,
                    barcode_sample_assignments,
                )
            }
        };
        let (vdj_inputs, vdj_gen_inputs) = {
            let mut per_lib_vdj_sample_def = BTreeMap::new();
            for lib in &cfg.libraries.0 {
                // TODO: Pass only a single VdjInput with all sample_defs and
                // vector-types for other per-chain configurables
                if !lib.is_count()? {
                    per_lib_vdj_sample_def
                        .entry((lib.physical_library_id(), lib.feature_types()))
                        .or_insert_with(Vec::new)
                        .push(lib.to_sample_def(&cfg)?);
                }
            }
            if per_lib_vdj_sample_def.len() > 3 {
                // TODO: This needs to be relaxed with multi gem well
                bail!(
                    "Found {} VDJ libraries, but we expect at most 3 libraries",
                    per_lib_vdj_sample_def.len()
                );
            }
            // TODO: we need to partition sample_def by gem-well and vdj_b or vdj_t,
            // but only know the first now...
            if per_lib_vdj_sample_def.is_empty() {
                (vec![], None)
            } else {
                let vdj = cfg.vdj.as_ref().ok_or_else(
                    #[cold]
                    || anyhow!("missing [vdj] table for VDJ libraries"),
                )?;
                let mut vdj_inputs = Vec::new();
                for ((physical_library_id, feature_types), sample_def) in per_lib_vdj_sample_def {
                    use FeatureType::{VDJ_B, VDJ_T, VDJ_T_GD};
                    let chain_type = match feature_types {
                        [VDJ_T] => Some("TR".to_string()),
                        [VDJ_T_GD] => {
                            // In gamma/delta mode we need inner-enrichment-primers to be present in multi config
                            if vdj.inner_enrichment_primers.is_none() {
                                bail!(
                                    "VDJ-T-GD library requires inner enrichment primers to be specified in the multi config file."
                                );
                            }
                            Some("TR_GD".to_string())
                        }
                        [VDJ_B] => Some("IG".to_string()),
                        _ => None,
                    };
                    vdj_inputs.push(VdjInputs {
                        sample_def,
                        chemistry: Some("SCVDJ_auto".to_string()),
                        custom_chemistry_def: None,

                        primers: primers.clone(),
                        subsample_rate,
                        initial_reads,
                        primer_initial_reads: Some(1000000),
                        special_genomic_regions: special_genomic_regions.clone(),
                        denovo: false,
                        r1_length: vdj.r1_length,
                        r2_length: vdj.r2_length,
                        ground_truth_clonotype_path: None,
                        inner_enrichment_primers: vdj.inner_enrichment_primers.clone(),
                        chain_type,
                        physical_library_id: Some(physical_library_id.to_string()),
                        r2_revcomp: vdj.r2_revcomp.unwrap_or_default(),
                    });
                }
                let vdj_gen_inputs = VdjGenInputs {
                    reference_path: cfg
                        .gene_expression
                        .as_ref()
                        .map(|gex| PathBuf::from(&gex.reference_path)),
                    vdj_reference_path: Some(vdj.reference_path.clone()),
                    filter_flags: VdjFilterFlags {
                        multiplet_filter: vdj.multiplet_filter,
                        shared_contig_filter: vdj.shared_contig_filter,
                        umi_baseline_filter: vdj.umi_baseline_filter,
                    },
                };
                (vdj_inputs, Some(vdj_gen_inputs))
            }
        };
        let basic_config = BasicPipelineConfig {
            disable_count: count_input.is_none(),
            disable_vdj: vdj_inputs.is_empty(),
            disable_multi: false,
            disable_multi_count: count_input.is_none(),
        };

        let common_input = CommonInputs {
            sample_id: args.sample_id,
            sample_desc: args.sample_desc,
            multi_config_sha: {
                let mut hasher = Sha256::new();
                hasher.update(&std::fs::read_to_string(&config_file)?);
                format!("{:x}", hasher.finalize())
            },
        };

        let feature_config = create_feature_config(
            cfg.antigen_specificity.as_ref(),
            cfg.functional_map.as_ref(),
            cfg.libraries.beam_mode(),
        );

        // Create feature reference object to make sure features are continous
        // And validate feature ref if beam mode can be deciphered (VDJ chain type provided)
        if let Some(ref feature_ref) = feature_ref {
            let feature_reference = feature_ref.read(feature_config.as_ref())?;
            if let Some(beam_mode) = feature_config.as_ref().and_then(|c| c.beam_mode) {
                feature_reference.validate_beam_feature_ref(beam_mode)?;
            }
        }

        Ok(ParseMultiConfigStageOutputs {
            common_input,
            count_input,
            vdj_inputs,
            vdj_gen_inputs,
            basic_config,
            config_file,
            feature_config,
            feature_ref,
            cell_barcodes,
            sample_barcodes,
            non_singlet_barcodes,
            cells_per_tag,
            barcode_sample_assignments,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multi::config::ProbeBarcodeIterationMode;
    use std::path::Path;

    fn test_run_stage_is_pd(
        csv_file: impl AsRef<Path>,
        sample_id: &str,
        is_pd: bool,
    ) -> Result<ParseMultiConfigStageOutputs> {
        ParseMultiConfig.test_run_tmpdir(ParseMultiConfigStageInputs {
            sample_id: sample_id.into(),
            sample_desc: String::new(),
            config: FileOrBytes {
                bytes: None,
                file: Some(csv_file.as_ref().into()),
            },
            config_hash: None,
            params: None,
            is_pd,
        })
    }

    fn test_run_stage(
        csv_file: impl AsRef<Path>,
        sample_id: &str,
    ) -> Result<ParseMultiConfigStageOutputs> {
        test_run_stage_is_pd(csv_file, sample_id, false)
    }

    fn insta_settings(test_csv: &str) -> insta::Settings {
        use insta::dynamic_redaction;
        let replace_if_not_none = |value, key| {
            if matches!(value, insta::internals::Content::None) {
                value
            } else {
                format!("[{key}:redacted]").into()
            }
        };
        let mut settings = insta::Settings::clone_current();
        settings.add_redaction(".config_file", test_csv);
        settings.add_redaction(
            ".count_input.gene_index",
            dynamic_redaction(move |value, _| replace_if_not_none(value, "gene_index")),
        );
        settings.add_redaction(
            ".count_input.feature_reference",
            dynamic_redaction(move |value, _| replace_if_not_none(value, "feature_reference")),
        );
        settings.add_redaction(
            ".feature_ref",
            dynamic_redaction(move |value, _| replace_if_not_none(value, "feature_ref")),
        );
        settings.add_redaction(".multi_graph", "[multi_graph:redacted]");
        settings
    }

    fn json_snapshot(sample_id: &str) {
        let test_csv = &format!("test/multi/{sample_id}.csv");
        insta_settings(test_csv).bind(|| {
            let outs = test_run_stage(test_csv, sample_id).unwrap();
            insta::assert_json_snapshot!(sample_id, &outs);
        });
    }

    #[test]
    fn test_parse_cmos() {
        let cfg = MultiConfigCsvFile::new("test/multi", "gex_multi_cmos.csv")
            .read()
            .unwrap();
        let cmos = cfg.sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::All);
        assert_eq!(cmos.len(), 3);
        assert!(cmos.contains("CMO1"));
        assert!(cmos.contains("CMO2"));
        assert!(cmos.contains("CMO3"));

        assert!(cfg
            .samples
            .unwrap()
            .get_translated_probe_barcodes()
            .is_empty())
    }

    #[test]
    fn test_parse_probe_barcode_ids() {
        let cfg = MultiConfigCsvFile::new("test/multi", "mfrp_multi.csv")
            .read()
            .unwrap();
        let probe_barcode_ids =
            cfg.sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::All);
        assert_eq!(probe_barcode_ids.len(), 2);
        assert!(probe_barcode_ids.contains("BC001"));
        assert!(probe_barcode_ids.contains("BC002"));
    }

    #[test]
    fn test_parse_mapped_probe_barcode_ids() {
        let cfg = MultiConfigCsvFile::from("test/multi/mfrp_ab_multi.csv")
            .read()
            .unwrap();
        let probe_barcode_ids =
            cfg.sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::All);
        assert_eq!(probe_barcode_ids.len(), 4);
        assert!(probe_barcode_ids.contains("BC001"));
        assert!(probe_barcode_ids.contains("BC002"));
        assert!(probe_barcode_ids.contains("AB002"));
        assert!(probe_barcode_ids.contains("BC004"));
        let probe_barcode_ids =
            cfg.sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::Mapped);
        assert_eq!(probe_barcode_ids.len(), 3);
        assert!(probe_barcode_ids.contains("BC001"));
        assert!(probe_barcode_ids.contains("BC002"));
        assert!(probe_barcode_ids.contains("BC004"));
    }

    #[test]
    fn test_vdj_internal() {
        json_snapshot("vdj_micro")
    }

    #[test]
    fn test_vdj_gd_missing_primers() {
        let outs = test_run_stage(
            "test/multi/invalid_csvs/vdj_micro_gd_no_primer.csv",
            "vdj_micro_gd_noprimer",
        );
        insta::assert_display_snapshot!(&outs.unwrap_err());
    }

    #[test]
    fn test_vdj_gex_internal() {
        json_snapshot("vdj_gex_micro")
    }

    #[test]
    fn test_gex_fbc_internal() {
        json_snapshot("gex_fbc_micro")
    }

    #[test]
    fn test_gex_fbc_dos_internal() {
        json_snapshot("gex_fbc_micro_dos")
    }

    #[test]
    fn test_gex_fbc_dos_utf8_internal() {
        json_snapshot("gex_fbc_micro_dos_utf8")
    }

    #[test]
    fn test_gex_fbc_mac_utf8_internal() {
        json_snapshot("gex_fbc_micro_mac_utf8")
    }

    #[test]
    fn test_fb_only_missing_gex_section() {
        let outs = test_run_stage(
            "test/multi/invalid_csvs/fb_only_missing_gex_section.csv",
            "fb",
        );
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_gex_missing_gex_section() {
        let outs = test_run_stage("test/multi/invalid_csvs/gex_missing_gex_section.csv", "gex");
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_gex_vdj_beamab_internal() {
        json_snapshot("beamab_vdj_gex")
    }

    #[test]
    fn test_gex_vdj_beamt_internal() {
        json_snapshot("beamt_vdj_gex")
    }

    #[test]
    fn test_gex_vdj_beamab_with_antigen_specificity_internal() {
        json_snapshot("beamab_vdj_gex_antigen_spec")
    }

    #[test]
    fn test_non_continous_feature_ref() {
        let outs = test_run_stage(
            "test/multi/non_continous_feature_ref.csv",
            "non_continous_feature_ref",
        );
        insta::assert_display_snapshot!(outs.unwrap_err());
    }

    #[test]
    fn test_vdj_filters() {
        json_snapshot("vdj_micro_filters")
    }
}
