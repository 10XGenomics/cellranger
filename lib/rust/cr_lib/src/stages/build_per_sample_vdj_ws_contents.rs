//! Martian stage BUILD_PER_SAMPLE_VDJ_WS_CONTENTS
use crate::stages::parse_multi_config::VdjGenInputs;
use crate::SequencingMetricsFormat;
use anyhow::Result;
use cr_types::{
    BarcodeMultiplexingType, CellLevel, CrMultiGraph, LibraryType, MetricsFile, ReadLevel,
    SampleAssignment,
};
use cr_websummary::multi::tables::{
    SequencingMetricsTable, VdjBEnrichmentMetricsTable, VdjBLibraryCellMetricsTable,
    VdjBSampleAnnotationMetricsTable, VdjBSampleHeroMetricsTable, VdjLibraryMetricsPerHashtagIdRow,
    VdjLibraryMetricsPerHashtagIdTable, VdjLibraryMetricsPerOcmBarcodeRow,
    VdjLibraryMetricsPerOcmBarcodeTable, VdjPhysicalLibraryMetricsTable,
    VdjTEnrichmentMetricsTable, VdjTLibraryCellMetricsTable, VdjTSampleAnnotationMetricsTable,
    VdjTSampleHeroMetricsTable, VdjTgdEnrichmentMetricsTable, VdjTgdLibraryCellMetricsTable,
    VdjTgdSampleAnnotationMetricsTable, VdjTgdSampleHeroMetricsTable,
};
use cr_websummary::multi::websummary::{
    ClonotypeInfo, LibraryVdjWebSummary, SampleVdjWebSummary, VdjChainTypeSpecific, VdjDiagnostics,
    VdjEnrichmentMetricsTable, VdjLibraryCellMetricsTable, VdjParametersTable,
    VdjSampleAnnotationMetricsTable, VdjSampleHeroMetricsTable,
};
use cr_websummary::{ChartWithHelp, CountAndPercent, TitleWithHelp};
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::{JsonFile, JsonFormat};
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{PercentMetric, TxHashMap};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::iter::zip;
use vdj_reference::VdjReceptor;

const VDJ_BARCODE_RANK_PLOT_TITLE: &str = "V(D)J Barcode Rank Plot";
const VDJ_BARCODE_RANK_PLOT_HELP: &str = "The plot shows the count of filtered UMIs mapped to each barcode. A barcode must have a contig that aligns to a V segment to be identified as a targeted cell. (In the denovo case, the only requirement is a contig's presence.) There must also be at least three filtered UMIs with at least two read pairs each. It is possible that a barcode with at least as many filtered UMIs as another cell-associated barcode is not identified as a targeted cell. The color of the graph is based on the local density of cell-associated barcodes. Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts.";
const VDJ_CLONOTYPE_INFO_PLOT_TITLE: &str = "Top 10 Clonotypes";
const VDJ_CLONOTYPE_INFO_PLOT_HELP: &str = r#"The histogram displays the fraction of cells (percentage of cells) occupied by the 10 most abundant clonotypes in this sample. The clonotype IDs on the X axis correspond to the clonotype IDs listed in the table. The table lists the CDR3 sequence of the first exact subclonotype of the 10 most abundant clonotypes in this sample. For each of the top 10 clonotypes, the constant region, number of cells (frequency), and what percentage of the dataset those cells occupy (proportion) are also displayed. For the full table and more details, please refer to the "clonotypes.csv" and "consensus_annotations.csv" files produced by the pipeline."#;

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsLibraryContents {
    pub sequencing_metrics_table: SequencingMetricsTable,
    pub physical_library_metrics_table: VdjPhysicalLibraryMetricsTable,
    pub cell_metrics_table: VdjLibraryCellMetricsTable,
    pub enrichment_metrics_table: VdjEnrichmentMetricsTable,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub metrics_per_ocm_barcode_table: Option<VdjLibraryMetricsPerOcmBarcodeTable>,
    pub metrics_per_hashtag_id_table: Option<VdjLibraryMetricsPerHashtagIdTable>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsSampleContents {
    pub hero_metrics: VdjSampleHeroMetricsTable,
    pub annotation_metrics_table: VdjSampleAnnotationMetricsTable,
    pub enrichment_metrics_table: Option<VdjEnrichmentMetricsTable>,
    pub clonotype_info: Option<ClonotypeInfo>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsContents {
    pub receptor: VdjReceptor,
    pub chemistry: String,
    pub vdj_reference: String,
    pub vdj_reference_path: String,
    pub lib_contents: VdjWsLibraryContents,
    pub sample_contents: VdjWsSampleContents,
    pub filter_metrics: Option<Value>,
}

impl VdjWsContents {
    pub fn to_library_ws(self) -> LibraryVdjWebSummary {
        LibraryVdjWebSummary {
            parameters_table: VdjParametersTable {
                chemistry: self.chemistry,
                vdj_reference: self.vdj_reference,
                vdj_reference_path: self.vdj_reference_path,
                gamma_delta: self.receptor == VdjReceptor::TRGD,
            },
            sequencing_metrics_table: self.lib_contents.sequencing_metrics_table.into(),
            cell_metrics_table: self.lib_contents.cell_metrics_table.into(),
            enrichment_metrics_table: self.lib_contents.enrichment_metrics_table.into(),
            physical_library_metrics_table: self.lib_contents.physical_library_metrics_table.into(),
            barcode_rank_plot: self.lib_contents.barcode_rank_plot,
            metrics_per_ocm_barcode_table: self
                .lib_contents
                .metrics_per_ocm_barcode_table
                .map(Into::into),
            metrics_per_hashtag_id_table: self
                .lib_contents
                .metrics_per_hashtag_id_table
                .map(Into::into),
        }
    }

    pub fn diagnostics(&self) -> VdjDiagnostics {
        VdjDiagnostics {
            filter_metrics: self.filter_metrics.clone(),
        }
    }

    pub fn to_sample_ws(self) -> SampleVdjWebSummary {
        SampleVdjWebSummary {
            parameters_table: VdjParametersTable {
                chemistry: self.chemistry,
                vdj_reference: self.vdj_reference,
                vdj_reference_path: self.vdj_reference_path,
                gamma_delta: self.receptor == VdjReceptor::TRGD,
            },
            hero_metrics: self.sample_contents.hero_metrics.into(),
            annotation_metrics_table: self.sample_contents.annotation_metrics_table.into(),
            enrichment_metrics_table: self
                .sample_contents
                .enrichment_metrics_table
                .map(Into::into),
            clonotype_info: self.sample_contents.clonotype_info,
            barcode_rank_plot: self.sample_contents.barcode_rank_plot,
        }
    }
}

martian_filetype!(_VdjWsContentsFile, "vwc");
pub type VdjWsContentsFormat = JsonFormat<_VdjWsContentsFile, VdjWsContents>;

#[derive(Deserialize, MartianStruct)]
pub struct BuildPerSampleVdjWsContentsStageInputs {
    pub receptor: VdjReceptor,
    pub physical_library_id: String,
    pub multiplexing_method: Option<BarcodeMultiplexingType>,
    pub lib_level_metrics: MetricsFile,
    pub per_sample_metrics: HashMap<SampleAssignment, MetricsFile>,
    pub vdj_gen_inputs: VdjGenInputs,
    pub sequencing_metrics: SequencingMetricsFormat,
    pub lib_level_vdj_ws_json: JsonFile<Value>, // Web summary JSON from cellranger vdj pipeline
    pub per_sample_vdj_ws_json: HashMap<SampleAssignment, JsonFile<Value>>, // Web summary JSON from cellranger vdj pipeline
    pub filter_metrics: HashMap<SampleAssignment, Option<JsonFile<Value>>>,
    pub multi_graph: JsonFile<CrMultiGraph>,
    pub vdj_cells_per_tag_json: Option<JsonFile<HashMap<String, usize>>>,
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

fn generate_clonotype_info_plot(value: &Value) -> Option<ClonotypeInfo> {
    if let Some(analysis_tab) = value["summary"].get("vdj_analysis_tab") {
        let clonotype_info = ClonotypeInfo {
            table: serde_json::from_value(analysis_tab["vdj_clonotype"]["table"].clone()).unwrap(),
            plot: serde_json::from_value(analysis_tab["vdj_clonotype_hist"]["plot"].clone())
                .unwrap(),
            help: TitleWithHelp {
                help: VDJ_CLONOTYPE_INFO_PLOT_HELP.to_string(),
                title: VDJ_CLONOTYPE_INFO_PLOT_TITLE.to_string(),
            },
        };
        Some(clonotype_info)
    } else {
        None
    }
}

fn build_metrics_per_tag(
    multi_graph: &CrMultiGraph,
    vdj_cells_per_tag: Option<JsonFile<HashMap<String, usize>>>,
) -> (TxHashMap<String, String>, TxHashMap<String, i64>, i64) {
    let tag_id_to_sample_id = multi_graph
        .get_tag_name_to_sample_id_map()
        .unwrap()
        .into_iter()
        .map(|(u, v)| (u.to_string(), v.0.to_string()))
        .collect();

    let vdj_cells_per_tag: TxHashMap<String, i64> = vdj_cells_per_tag
        .unwrap()
        .read()
        .unwrap()
        .into_iter()
        .map(|(id, n)| (id.clone(), n as i64))
        .collect();

    let vdj_cells_in_library: i64 = vdj_cells_per_tag.values().sum();
    (tag_id_to_sample_id, vdj_cells_per_tag, vdj_cells_in_library)
}
fn build_metrics_per_hashtag_id_table(
    multi_graph: &CrMultiGraph,
    vdj_cells_per_tag: Option<JsonFile<HashMap<String, usize>>>,
) -> Option<VdjLibraryMetricsPerHashtagIdTable> {
    if multi_graph.barcode_multiplexing_type()
        != Some(BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag))
    {
        return None;
    }
    let (tag_id_to_sample_id, vdj_cells_per_tag, vdj_cells_in_library) =
        build_metrics_per_tag(multi_graph, vdj_cells_per_tag);

    let rows: Vec<VdjLibraryMetricsPerHashtagIdRow> = vdj_cells_per_tag
        .into_iter()
        .map(|(tag_id, cells)| {
            let sample_id = tag_id_to_sample_id.get(&tag_id).cloned();
            let fraction_cells_per_tag = Some(CountAndPercent(PercentMetric::from_parts(
                cells,
                vdj_cells_in_library,
            )));
            match sample_id {
                Some(x) => VdjLibraryMetricsPerHashtagIdRow {
                    hashtag_id: Some(tag_id.to_string()),
                    sample_id: Some(x.to_string()),
                    vdj_cells_per_tag: fraction_cells_per_tag,
                },
                None => VdjLibraryMetricsPerHashtagIdRow {
                    hashtag_id: Some(tag_id.to_string()),
                    sample_id: None,
                    vdj_cells_per_tag: fraction_cells_per_tag,
                },
            }
        })
        .collect();
    Some(VdjLibraryMetricsPerHashtagIdTable(rows))
}

fn build_metrics_per_ocm_barcode_table(
    multi_graph: &CrMultiGraph,
    vdj_cells_per_tag: Option<JsonFile<HashMap<String, usize>>>,
) -> Option<VdjLibraryMetricsPerOcmBarcodeTable> {
    if multi_graph.barcode_multiplexing_type()
        != Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::OH))
    {
        return None;
    }
    let (tag_id_to_sample_id, vdj_cells_per_tag, vdj_cells_in_library) =
        build_metrics_per_tag(multi_graph, vdj_cells_per_tag);

    let rows: Vec<VdjLibraryMetricsPerOcmBarcodeRow> = vdj_cells_per_tag
        .into_iter()
        .map(|(tag_id, cells)| {
            let sample_id = tag_id_to_sample_id.get(&tag_id).cloned();
            let fraction_cells_per_tag = Some(CountAndPercent(PercentMetric::from_parts(
                cells,
                vdj_cells_in_library,
            )));
            match sample_id {
                Some(x) => VdjLibraryMetricsPerOcmBarcodeRow {
                    ocm_barcode_id: Some(tag_id.to_string()),
                    sample_id: Some(x.to_string()),
                    vdj_cells_per_tag: fraction_cells_per_tag,
                },
                None => VdjLibraryMetricsPerOcmBarcodeRow {
                    ocm_barcode_id: Some(tag_id.to_string()),
                    sample_id: None,
                    vdj_cells_per_tag: fraction_cells_per_tag,
                },
            }
        })
        .collect();
    Some(VdjLibraryMetricsPerOcmBarcodeTable(rows))
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
        let clonotype_info_plot = generate_clonotype_info_plot(&sample_vdj_ws_json);

        // End of extraction from vdj web summary json

        let mut library_metrics: HashMap<String, Value> = args.lib_level_metrics.read()?;
        library_metrics.insert(
            "physical_library_id".to_string(),
            args.physical_library_id.into(),
        );

        let sample_metrics: HashMap<String, Value> =
            args.per_sample_metrics.get(&sample).unwrap().read()?;

        let vdj_ref_version = match library_metrics["vdj_reference_version"].as_str() {
            Some(x) => format!("-{x}"),
            None => String::default(),
        };

        let is_read_level_multiplexed = matches!(
            args.multiplexing_method,
            Some(BarcodeMultiplexingType::ReadLevel(_))
        );

        let multi_graph = args.multi_graph.read()?;

        let contents = VdjWsContents {
            receptor: args.receptor,
            chemistry: library_metrics["chemistry_description"]
                .as_str()
                .unwrap()
                .to_string(),
            vdj_reference: format!(
                "{}{}",
                library_metrics["vdj_reference_genomes"].as_str().unwrap(),
                vdj_ref_version
            ),
            vdj_reference_path: args
                .vdj_gen_inputs
                .vdj_reference_path
                .unwrap()
                .display()
                .to_string(),
            lib_contents: VdjWsLibraryContents {
                sequencing_metrics_table: SequencingMetricsTable(
                    args.sequencing_metrics
                        .read()?
                        // CELLRANGER-7889 this may need to be updated eventually
                        .remove(&LibraryType::VdjAuto)
                        .unwrap()
                        .into_iter()
                        .map(Into::into)
                        .collect(),
                ),
                physical_library_metrics_table: VdjPhysicalLibraryMetricsTable::from_metrics(
                    &library_metrics,
                )?,
                cell_metrics_table: match args.receptor {
                    VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                        VdjTLibraryCellMetricsTable::from_metrics(&library_metrics)?,
                    ),
                    VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                        VdjTgdLibraryCellMetricsTable::from_metrics(&library_metrics)?,
                    ),
                    VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                        VdjBLibraryCellMetricsTable::from_metrics(&library_metrics)?,
                    ),
                },
                enrichment_metrics_table: match args.receptor {
                    VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                        VdjTEnrichmentMetricsTable::from_metrics(&library_metrics)?,
                    ),
                    VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                        VdjTgdEnrichmentMetricsTable::from_metrics(&library_metrics)?,
                    ),
                    VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                        VdjBEnrichmentMetricsTable::from_metrics(&library_metrics)?,
                    ),
                },
                barcode_rank_plot: lib_bc_rank_plot,
                metrics_per_ocm_barcode_table: build_metrics_per_ocm_barcode_table(
                    &multi_graph,
                    args.vdj_cells_per_tag_json.clone(),
                ),
                metrics_per_hashtag_id_table: build_metrics_per_hashtag_id_table(
                    &multi_graph,
                    args.vdj_cells_per_tag_json.clone(),
                ),
            },
            sample_contents: VdjWsSampleContents {
                hero_metrics: match args.receptor {
                    VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                        VdjTSampleHeroMetricsTable::from_metrics(&sample_metrics)?,
                    ),
                    VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                        VdjTgdSampleHeroMetricsTable::from_metrics(&sample_metrics)?,
                    ),
                    VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                        VdjBSampleHeroMetricsTable::from_metrics(&sample_metrics)?,
                    ),
                },
                annotation_metrics_table: match args.receptor {
                    VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                        VdjTSampleAnnotationMetricsTable::from_metrics(&sample_metrics)?,
                    ),
                    VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                        VdjTgdSampleAnnotationMetricsTable::from_metrics(&sample_metrics)?,
                    ),
                    VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                        VdjBSampleAnnotationMetricsTable::from_metrics(&sample_metrics)?,
                    ),
                },
                enrichment_metrics_table: is_read_level_multiplexed.then_some(
                    match args.receptor {
                        VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                            VdjTEnrichmentMetricsTable::from_metrics(&sample_metrics).unwrap(),
                        ),
                        VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                            VdjTgdEnrichmentMetricsTable::from_metrics(&sample_metrics).unwrap(),
                        ),
                        VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                            VdjBEnrichmentMetricsTable::from_metrics(&sample_metrics).unwrap(),
                        ),
                    },
                ),
                clonotype_info: clonotype_info_plot,
                barcode_rank_plot: sample_bc_rank_plot,
            },
            filter_metrics: args
                .filter_metrics
                .get(&sample)
                .unwrap()
                .as_ref()
                .map(martian_filetypes::FileTypeRead::read)
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
