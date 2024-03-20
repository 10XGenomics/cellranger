//! Martian stage BUILD_PER_SAMPLE_VDJ_WS_CONTENTS

use crate::stages::parse_multi_config::VdjGenInputs;
use crate::SequencingMetricsFormat;
use anyhow::Result;
use cr_types::{LibraryType, MetricsFile, SampleAssignment};
use cr_websummary::multi::tables::{
    SequencingMetricsTable, VdjBEnrichmentMetricsTable, VdjBSampleAnnotationMetricsTable,
    VdjBSampleHeroMetricsTable, VdjLibraryCellMetricsTable, VdjPhysicalLibraryMetricsTable,
    VdjTEnrichmentMetricsTable, VdjTSampleAnnotationMetricsTable, VdjTSampleHeroMetricsTable,
    VdjTgdEnrichmentMetricsTable, VdjTgdSampleAnnotationMetricsTable, VdjTgdSampleHeroMetricsTable,
};
use cr_websummary::multi::websummary::{
    ClonotypeInfo, LibraryVdjWebSummary, SampleVdjWebSummary, VdjChainTypeSpecific, VdjDiagnostics,
    VdjEnrichmentMetricsTable, VdjParametersTable, VdjSampleAnnotationMetricsTable,
    VdjSampleHeroMetricsTable,
};
use cr_websummary::{ChartWithHelp, TitleWithHelp};
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::{JsonFile, JsonFormat};
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::iter::zip;
use vdj_reference::VdjReceptor;

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsContents {
    pub chemistry: String,
    pub receptor: VdjReceptor,
    pub vdj_reference: String,
    pub vdj_reference_path: String,
    pub sequencing_metrics_table: SequencingMetricsTable,
    pub cell_metrics_table: VdjLibraryCellMetricsTable,
    pub enrichment_metrics_table: VdjEnrichmentMetricsTable,
    pub physical_library_metrics_table: VdjPhysicalLibraryMetricsTable,
    pub hero_metrics: VdjSampleHeroMetricsTable,
    pub annotation_metrics_table: VdjSampleAnnotationMetricsTable,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub clonotype_info: Option<ClonotypeInfo>,
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
            sequencing_metrics_table: self.sequencing_metrics_table.into(),
            cell_metrics_table: self.cell_metrics_table.into(),
            enrichment_metrics_table: self.enrichment_metrics_table.into(),
            physical_library_metrics_table: self.physical_library_metrics_table.into(),
            barcode_rank_plot: self.barcode_rank_plot,
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
            hero_metrics: self.hero_metrics.into(),
            annotation_metrics_table: self.annotation_metrics_table.into(),
            clonotype_info: self.clonotype_info,
            barcode_rank_plot: None,
        }
    }
}

martian_filetype!(_VdjWsContentsFile, "vwc");
pub type VdjWsContentsFormat = JsonFormat<_VdjWsContentsFile, VdjWsContents>;

#[derive(Deserialize, MartianStruct)]
pub struct BuildPerSampleVdjWsContentsStageInputs {
    pub lib_level_metrics: MetricsFile,
    pub per_sample_metrics: HashMap<SampleAssignment, MetricsFile>,
    pub receptor: VdjReceptor,
    pub physical_library_id: String,
    pub vdj_gen_inputs: VdjGenInputs,
    pub sequencing_metrics: SequencingMetricsFormat,
    pub vdj_ws_json: HashMap<SampleAssignment, JsonFile<Value>>, // Web summary JSON for the vdj pipeline
    pub filter_metrics: HashMap<SampleAssignment, Option<JsonFile<Value>>>,
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
        let vdj_ws_json = args.vdj_ws_json.get(&sample).unwrap().read()?;
        let cells_value = &vdj_ws_json["summary"]["summary_tab"]["cells"];
        let barcode_rank_plot = cells_value.get("barcode_knee_plot").map(|plot| ChartWithHelp {
            plot: {
                let mut plot = plot.clone();
                plot["layout"].as_object_mut().unwrap().remove("title");
                serde_json::from_value(plot).unwrap()
            },
            help: TitleWithHelp {
                help: "The plot shows the count of filtered UMIs mapped to each barcode. A barcode must have a contig that aligns to a V segment to be identified as a targeted cell. (In the denovo case, the only requirement is a contig's presence.) There must also be at least three filtered UMIs with at least two read pairs each. It is possible that a barcode with at least as many filtered UMIs as another cell-associated barcode is not identified as a targeted cell. The color of the graph is based on the local density of cell-associated barcodes. Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts.".into(),
                title: "V(D)J Barcode Rank Plot".into(),
            }
        });
        let clonotype_info = if let Some(analysis_tab) =
            &vdj_ws_json["summary"].get("vdj_analysis_tab")
        {
            let clonotype_info = ClonotypeInfo {
                plot: serde_json::from_value(analysis_tab["vdj_clonotype_hist"]["plot"].clone())?,
                table: serde_json::from_value(analysis_tab["vdj_clonotype"]["table"].clone())?,
                help: TitleWithHelp {
                    help: r#"The histogram displays the fraction of cells (percentage of cells) occupied by the 10 most abundant clonotypes in this sample. The clonotype IDs on the X axis correspond to the clonotype IDs listed in the table. The table lists the CDR3 sequence of the first exact subclonotype of the 10 most abundant clonotypes in this sample. For each of the top 10 clonotypes, the constant region, number of cells (frequency), and what percentage of the dataset those cells occupy (proportion) are also displayed. For the full table and more details, please refer to the "clonotypes.csv" and "consensus_annotations.csv" files produced by the pipeline."#.into(),
                    title: "Top 10 Clonotypes".into(),
                }
            };
            Some(clonotype_info)
        } else {
            None
        };

        // End of extraction from vdj web summary json

        let mut library_metrics: Value = args.lib_level_metrics.read()?;
        library_metrics.as_object_mut().unwrap().insert(
            "physical_library_id".to_string(),
            args.physical_library_id.into(),
        );

        let sample_metrics: Value = args.per_sample_metrics.get(&sample).unwrap().read()?;

        let vdj_ref_version = match library_metrics["vdj_reference_version"].as_str() {
            Some(x) => format!("-{x}"),
            None => String::default(),
        };

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
            cell_metrics_table: VdjLibraryCellMetricsTable::from_json_value(&library_metrics),
            enrichment_metrics_table: match args.receptor {
                VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                    VdjTEnrichmentMetricsTable::from_json_value(&library_metrics),
                ),
                VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                    VdjTgdEnrichmentMetricsTable::from_json_value(&library_metrics),
                ),
                VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                    VdjBEnrichmentMetricsTable::from_json_value(&library_metrics),
                ),
            },
            physical_library_metrics_table: VdjPhysicalLibraryMetricsTable::from_json_value(
                &library_metrics,
            ),
            hero_metrics: match args.receptor {
                VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                    VdjTSampleHeroMetricsTable::from_json_value(&sample_metrics),
                ),
                VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                    VdjTgdSampleHeroMetricsTable::from_json_value(&sample_metrics),
                ),
                VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                    VdjBSampleHeroMetricsTable::from_json_value(&sample_metrics),
                ),
            },
            annotation_metrics_table: match args.receptor {
                VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                    VdjTSampleAnnotationMetricsTable::from_json_value(&sample_metrics),
                ),
                VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                    VdjTgdSampleAnnotationMetricsTable::from_json_value(&sample_metrics),
                ),
                VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                    VdjBSampleAnnotationMetricsTable::from_json_value(&sample_metrics),
                ),
            },
            sequencing_metrics_table: args
                .sequencing_metrics
                .read()?
                // CELLRANGER-7889 this may need to be updated eventually
                .remove(&LibraryType::VdjAuto)
                .unwrap(),
            barcode_rank_plot,
            clonotype_info,
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
