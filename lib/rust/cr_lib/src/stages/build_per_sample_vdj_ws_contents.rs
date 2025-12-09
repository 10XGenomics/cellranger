//! Martian stage BUILD_PER_SAMPLE_VDJ_WS_CONTENTS
#![deny(missing_docs)]

use crate::SequencingMetricsFormat;
use crate::stages::parse_multi_config::VdjGenInputs;
use anyhow::Result;
use cr_types::{
    BarcodeMultiplexingType, CellLevel, CrMultiGraph, LibraryType, MetricsFile, ReadLevel,
    SampleAssignment,
};
use cr_websummary::multi::metrics::{
    ActiveConditions, CountAndPercentTransformer, MetricsProcessor, load_metrics_etl_special,
    load_metrics_etl_vdj,
};
use cr_websummary::multi::websummary::{JsonMetricSummary, MetricsTraitWrapper, Section};
#[allow(clippy::wildcard_imports)]
use cr_websummary::multi::websummary_vdj::*;
use cr_websummary::{AlertContext, ChartWithHelp, TitleWithHelp};
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use martian_filetypes::json_file::{JsonFile, JsonFormat};
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::iter::zip;
use vdj_reference::VdjReceptor;

const VDJ_BARCODE_RANK_PLOT_TITLE: &str = "V(D)J Barcode Rank Plot";
const VDJ_BARCODE_RANK_PLOT_HELP: &str = "The plot shows the count of filtered UMIs mapped to each barcode. A barcode must have a contig that aligns to a V segment to be identified as a targeted cell. (In the denovo case, the only requirement is a contig's presence.) There must also be at least three filtered UMIs with at least two read pairs each. It is possible that a barcode with at least as many filtered UMIs as another cell-associated barcode is not identified as a targeted cell. The color of the graph is based on the local density of cell-associated barcodes. Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts.";

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsLibraryContents {
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub metrics: Vec<JsonMetricSummary>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsSampleContents {
    pub clonotype_info: Option<ClonotypeInfo>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub metrics: Vec<JsonMetricSummary>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsContents {
    pub receptor: VdjReceptor,
    pub denovo: bool,
    pub chemistry: String,
    pub vdj_reference: Option<String>,
    pub vdj_reference_path: Option<String>,
    pub lib_contents: VdjWsLibraryContents,
    pub sample_contents: VdjWsSampleContents,
    pub filter_metrics: Option<Value>,
}

impl VdjWsContents {
    pub fn into_library_ws(self) -> LibraryVdjWebSummary {
        LibraryVdjWebSummary {
            parameters_table: VdjParametersTable {
                chemistry: self.chemistry,
                vdj_reference: self.vdj_reference,
                vdj_reference_path: self.vdj_reference_path,
                gamma_delta: self.receptor == VdjReceptor::TRGD,
                denovo: self.denovo,
            },
            barcode_rank_plot: self.lib_contents.barcode_rank_plot,
            metrics: MetricsTraitWrapper(self.lib_contents.metrics),
        }
    }

    pub fn diagnostics(&self) -> VdjDiagnostics {
        VdjDiagnostics {
            filter_metrics: self.filter_metrics.clone(),
        }
    }

    pub fn into_sample_ws(self) -> SampleVdjWebSummary {
        SampleVdjWebSummary {
            clonotype_info: self.sample_contents.clonotype_info,
            barcode_rank_plot: self.sample_contents.barcode_rank_plot,
            metrics: MetricsTraitWrapper(self.sample_contents.metrics),
        }
    }
}

martian_filetype!(_VdjWsContentsFile, "vwc");
pub type VdjWsContentsFormat = JsonFormat<_VdjWsContentsFile, VdjWsContents>;

#[derive(Deserialize, MartianStruct)]
pub struct BuildPerSampleVdjWsContentsStageInputs {
    pub receptor: VdjReceptor,
    pub denovo: bool,
    pub physical_library_id: String,
    pub lib_level_metrics: MetricsFile,
    pub per_sample_metrics: HashMap<SampleAssignment, MetricsFile>,
    pub vdj_gen_inputs: VdjGenInputs,
    pub sequencing_metrics: SequencingMetricsFormat,
    pub lib_level_vdj_ws_json: JsonFile<Value>, // Web summary JSON from cellranger vdj pipeline
    pub per_sample_vdj_ws_json: HashMap<SampleAssignment, JsonFile<Value>>, // Web summary JSON from cellranger vdj pipeline
    pub filter_metrics: HashMap<SampleAssignment, Option<JsonFile<Value>>>,
    pub multi_graph: JsonFile<CrMultiGraph>,
    pub vdj_cells_per_tag_json: Option<JsonFile<HashMap<String, usize>>>,
    pub clonotype_info_json: HashMap<SampleAssignment, Option<JsonFile<ClonotypeInfo>>>,
}

#[derive(Debug, Serialize, Deserialize, Clone, MartianStruct)]
pub struct BuildPerSampleVdjWsContentsChunkInputs {
    sample: SampleAssignment,
}

#[derive(Debug, Serialize, Deserialize, Clone, MartianStruct)]
pub struct BuildPerSampleVdjWsContentsChunkOutputs {
    vdj_ws: VdjWsContentsFormat,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct BuildPerSampleVdjWsContentsStageOutputs {
    vdj_ws_contents: HashMap<SampleAssignment, VdjWsContentsFormat>,
}

/// Martian stage BUILD_PER_SAMPLE_VDJ_WS_CONTENTS
pub struct BuildPerSampleVdjWsContents;

fn generate_barcore_rank_plot(value: &Value) -> Option<ChartWithHelp> {
    let cells = &value["summary"]["summary_tab"]["cells"];
    cells.get("barcode_knee_plot").map(|plot| ChartWithHelp {
        plot: {
            let mut plot = plot.clone();
            plot["layout"].as_object_mut().unwrap().remove("title");
            serde_json::from_value(plot).unwrap()
        },
        help: TitleWithHelp {
            help: VDJ_BARCODE_RANK_PLOT_HELP.to_string(),
            title: VDJ_BARCODE_RANK_PLOT_TITLE.to_string(),
        },
    })
}

/// Construct the metrics per tag table.
fn build_metrics_per_tag(
    metrics_proc: &MetricsProcessor,
    active_conditions: &ActiveConditions,
    alert_context: &AlertContext,
    section: Section,
    multi_graph: &CrMultiGraph,
    vdj_cells_per_tag: Option<&JsonFile<HashMap<String, usize>>>,
) -> Result<Vec<JsonMetricSummary>> {
    let special_metrics_group = match multi_graph.barcode_multiplexing_type() {
        Some(BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag)) => {
            "vdj_library_metrics_per_hashtag_id"
        }
        Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::OH)) => {
            "vdj_library_metrics_per_ocm_barcode"
        }
        _ => {
            return Ok(vec![]);
        }
    };

    let Some(vdj_cells_per_tag) = vdj_cells_per_tag else {
        panic!("vdj_cells_per_tag file is unexpectedly null");
    };

    let tag_id_to_sample_id_data = multi_graph.get_tag_name_to_sample_id_map()?;
    let vdj_cells_per_tag = vdj_cells_per_tag.read()?;

    let mut metrics_proc = metrics_proc.with_default_transformers();

    metrics_proc.add_transformer(
        "CellsFraction",
        CountAndPercentTransformer::new(vdj_cells_per_tag.values().sum()),
    );

    #[derive(JsonReport)]
    struct Metrics<'a> {
        tag_id: &'a str,
        sample_id: &'a str,
        vdj_cells_per_tag: usize,
    }

    let mut metrics = vec![];
    for (tag_id, (sample_id, _)) in tag_id_to_sample_id_data {
        let cells = vdj_cells_per_tag[tag_id];
        metrics.extend(metrics_proc.process_group(
            special_metrics_group,
            section,
            active_conditions,
            alert_context,
            &Metrics {
                tag_id,
                sample_id,
                vdj_cells_per_tag: cells,
            },
        )?);
    }
    Ok(metrics)
}

#[make_mro(volatile = strict)]
impl MartianStage for BuildPerSampleVdjWsContents {
    type StageInputs = BuildPerSampleVdjWsContentsStageInputs;
    type StageOutputs = BuildPerSampleVdjWsContentsStageOutputs;
    type ChunkInputs = BuildPerSampleVdjWsContentsChunkInputs;
    type ChunkOutputs = BuildPerSampleVdjWsContentsChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        Ok(args
            .per_sample_metrics
            .keys()
            .map(|sample| BuildPerSampleVdjWsContentsChunkInputs {
                sample: sample.clone(),
            })
            .collect())
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let sample = chunk_args.sample;
        // For now, we will extract the bacode knee plot and the clonotype table/hist
        // from the vdj web summary json.
        // TODO: This is not very pretty. Cleanup in the future.
        let lib_vdj_ws_json = args.lib_level_vdj_ws_json.read()?;
        let lib_bc_rank_plot = generate_barcore_rank_plot(&lib_vdj_ws_json);

        let sample_vdj_ws_json = args.per_sample_vdj_ws_json.get(&sample).unwrap().read()?;
        let sample_bc_rank_plot = generate_barcore_rank_plot(&sample_vdj_ws_json);

        let clonotype_info_json = args.clonotype_info_json.get(&sample).unwrap().clone();
        let clonotype_info_plot: Option<ClonotypeInfo> = match clonotype_info_json {
            Some(clonotype_info_json) => Some(clonotype_info_json.read()?),
            None => None,
        };

        // End of extraction from vdj web summary json

        let (lib_metrics_proc, sample_metrics_proc) = load_metrics_etl_vdj()?;

        let mut library_metrics: TxHashMap<String, Value> = args.lib_level_metrics.read()?;
        library_metrics.insert(
            "physical_library_id".to_string(),
            args.physical_library_id.into(),
        );

        let sample_metrics: TxHashMap<String, Value> =
            args.per_sample_metrics.get(&sample).unwrap().read()?;

        let vdj_reference = if let Some(vdj_reference) =
            library_metrics.get("vdj_reference_genomes")
        {
            let vdj_ref_version = library_metrics
                .get("vdj_reference_version")
                .and_then(|v| v.as_str())
                .map(|x| format!("-{x}"))
                .unwrap_or_default();
            let vdj_ref_string = format!("{}{}", vdj_reference.as_str().unwrap(), vdj_ref_version);
            Some(vdj_ref_string)
        } else {
            None
        };

        let multi_graph = args.multi_graph.read()?;

        let section = match args.receptor {
            VdjReceptor::TR => Section::VdjT,
            VdjReceptor::TRGD => Section::VdjTGd,
            VdjReceptor::IG => Section::VdjB,
        };

        let active_conditions = ActiveConditions {
            section,
            vdj_receptor: Some(args.receptor),
            is_multiplexed: multi_graph.barcode_multiplexing_type().is_some(),
            is_cell_multiplexed: multi_graph.is_cell_level_multiplexed(),
            is_read_multiplexed: multi_graph.is_read_level_multiplexed(),
            is_rtl: false,          //FIXME CELLRANGER-8444
            has_gdna: false,        //FIXME CELLRANGER-8444
            include_introns: false, //FIXME CELLRANGER-8444
            has_vdj_reference: vdj_reference.is_some(),
        };

        let alert_context = AlertContext {
            is_rtl: active_conditions.is_rtl,
            is_arc_chemistry: false,
            library_types: multi_graph.library_types().collect(),
            multiplexing_method: multi_graph.barcode_multiplexing_type(),
            is_fiveprime: false, // FIXME CELLRANGER-8444
            include_introns: active_conditions.include_introns,
            no_preflight: false, // FIXME CELLRANGER-8444
        };

        let mut lib_metrics_out = lib_metrics_proc.process(
            section,
            &active_conditions,
            &alert_context,
            &library_metrics,
        )?;

        let special_metrics_proc = load_metrics_etl_special()?;

        lib_metrics_out.extend(build_metrics_per_tag(
            &special_metrics_proc,
            &active_conditions,
            &alert_context,
            section,
            &multi_graph,
            args.vdj_cells_per_tag_json.as_ref(),
        )?);

        for seq_metrics in args
            .sequencing_metrics
            .read()?
            // CELLRANGER-7889 this may need to be updated eventually
            .remove(&LibraryType::VdjAuto)
            .unwrap()
        {
            lib_metrics_out.extend(special_metrics_proc.process_group(
                "sequencing_metrics",
                section,
                &active_conditions,
                &alert_context,
                &seq_metrics,
            )?);
        }

        let contents = VdjWsContents {
            receptor: args.receptor,
            denovo: args.denovo,
            chemistry: library_metrics["chemistry_description"]
                .as_str()
                .unwrap()
                .to_string(),
            vdj_reference,
            vdj_reference_path: args
                .vdj_gen_inputs
                .vdj_reference_path
                .map(|p| p.display().to_string()),
            lib_contents: VdjWsLibraryContents {
                metrics: lib_metrics_out,
                barcode_rank_plot: lib_bc_rank_plot,
            },
            sample_contents: VdjWsSampleContents {
                metrics: sample_metrics_proc.process(
                    section,
                    &active_conditions,
                    &alert_context,
                    &sample_metrics,
                )?,
                clonotype_info: clonotype_info_plot,
                barcode_rank_plot: sample_bc_rank_plot,
            },
            filter_metrics: args
                .filter_metrics
                .get(&sample)
                .unwrap()
                .as_ref()
                .map(FileTypeRead::read)
                .transpose()?,
        };

        let vdj_ws_contents: VdjWsContentsFormat = rover.make_path("vdj_ws_contents");
        vdj_ws_contents.write(&contents)?;

        Ok(BuildPerSampleVdjWsContentsChunkOutputs {
            vdj_ws: vdj_ws_contents,
        })
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let vdj_ws_contents = zip(chunk_defs, chunk_outs)
            .map(|(defs, outs)| (defs.sample, outs.vdj_ws))
            .collect();
        Ok(BuildPerSampleVdjWsContentsStageOutputs { vdj_ws_contents })
    }
}
