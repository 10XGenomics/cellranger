//! Martian stage BUILD_VDJ_WS_CONTENTS

use crate::stages::parse_multi_config::{VdjGenInputs, VdjInputs};
use crate::SequencingMetricsFormat;
use anyhow::Result;
use cr_types::rna_read::LegacyLibraryType;
use cr_types::MetricsFile;
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
use vdj_reference::VdjReceptor;

#[derive(Clone, Serialize, Deserialize)]
pub struct VdjWsContents {
    chemistry: String,
    receptor: VdjReceptor,
    vdj_reference: String,
    vdj_reference_path: String,
    sequencing_metrics_table: SequencingMetricsTable,
    cell_metrics_table: VdjLibraryCellMetricsTable,
    enrichment_metrics_table: VdjEnrichmentMetricsTable,
    physical_library_metrics_table: VdjPhysicalLibraryMetricsTable,
    hero_metrics: VdjSampleHeroMetricsTable,
    annotation_metrics_table: VdjSampleAnnotationMetricsTable,
    barcode_rank_plot: Option<ChartWithHelp>,
    clonotype_info: ClonotypeInfo,
    filter_metrics: Option<Value>,
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
pub struct BuildVdjWsContentsStageInputs {
    pub metrics_summary: MetricsFile,
    pub receptor: VdjReceptor,
    pub vdj_inputs: VdjInputs,
    pub vdj_gen_inputs: VdjGenInputs,
    pub sequencing_metrics: SequencingMetricsFormat,
    pub vdj_ws_json: JsonFile<Value>, // Web summary JSON for the vdj pipeline
    pub filter_metrics: Option<JsonFile<Value>>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct BuildVdjWsContentsStageOutputs {
    vdj_ws_contents: VdjWsContentsFormat,
}

pub struct BuildVdjWsContents;

#[make_mro(volatile = strict)]
impl MartianMain for BuildVdjWsContents {
    type StageInputs = BuildVdjWsContentsStageInputs;
    type StageOutputs = BuildVdjWsContentsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // For now, we will extract the bacode knee plot and the clonotype table/hist
        // from the vdj web summary json.
        // TODO: This is not very pretty. Cleanup in the future.
        let vdj_ws_json = args.vdj_ws_json.read()?;
        let cells_value = &vdj_ws_json["summary"]["summary_tab"]["cells"];
        let barcode_rank_plot = cells_value.get("barcode_knee_plot").map(|plot| ChartWithHelp {
            plot: {
                let mut plot = plot.clone();
                plot["layout"].as_object_mut().unwrap().remove("title");
                serde_json::from_value(plot).unwrap()
            },
            help: TitleWithHelp {
                help: "The plot shows the count of filtered UMIs mapped to each barcode. A barcode must have a contig that aligns to a V segment to be identified as a targeted cell. (In the denovo case, the only requirement is a contig's presence.) There must also be at least three filtered UMIs with at least two read pairs each. It is possible that a barcode with at least as many filtered UMIs as another cell-associated barcode is not identified as a targeted cell. The color of the graph is based on the local density of cell-associated barcodes.".into(),
                title: "V(D)J Barcode Rank Plot".into(),
            }
        });
        let analysis_tab = &vdj_ws_json["summary"]["vdj_analysis_tab"];
        let clonotype_info = ClonotypeInfo {
            plot: serde_json::from_value(analysis_tab["vdj_clonotype_hist"]["plot"].clone())?,
            table: serde_json::from_value(analysis_tab["vdj_clonotype"]["table"].clone())?,
            help: TitleWithHelp {
                help: r#"The histogram displays the fraction of cells (percentage of cells) occupied by the 10 most abundant clonotypes in this sample. The clonotype IDs on the X axis correspond to the clonotype IDs listed in the table. The table lists the CDR3 sequence of the first exact subclonotype of the 10 most abundant clonotypes in this sample. For each of the top 10 clonotypes, the constant region, number of cells (frequency), and what percentage of the dataset those cells occupy (proportion) are also displayed. For the full table and more details, please refer to the "clonotypes.csv" and "consensus_annotations.csv" files produced by the pipeline."#.into(),
                title: "Top 10 Clonotypes".into(),
            }
        };

        // End of extraction from vdj web summary json

        let mut metric_json_val: Value = args.metrics_summary.read()?;
        metric_json_val.as_object_mut().unwrap().insert(
            "physical_library_id".to_string(),
            args.vdj_inputs.physical_library_id.clone().unwrap().into(),
        );

        let vdj_ref_version = match metric_json_val["vdj_reference_version"].as_str() {
            Some(x) => format!("-{x}"),
            None => String::default(),
        };

        let contents = VdjWsContents {
            receptor: args.receptor,
            chemistry: metric_json_val["chemistry_description"]
                .as_str()
                .unwrap()
                .to_string(),
            vdj_reference: format!(
                "{}{}",
                metric_json_val["vdj_reference_genomes"].as_str().unwrap(),
                vdj_ref_version
            ),
            vdj_reference_path: args
                .vdj_gen_inputs
                .vdj_reference_path
                .unwrap()
                .display()
                .to_string(),
            cell_metrics_table: VdjLibraryCellMetricsTable::from_json_value(&metric_json_val),
            enrichment_metrics_table: match args.receptor {
                VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                    VdjTEnrichmentMetricsTable::from_json_value(&metric_json_val),
                ),
                VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                    VdjTgdEnrichmentMetricsTable::from_json_value(&metric_json_val),
                ),
                VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                    VdjBEnrichmentMetricsTable::from_json_value(&metric_json_val),
                ),
            },
            physical_library_metrics_table: VdjPhysicalLibraryMetricsTable::from_json_value(
                &metric_json_val,
            ),
            hero_metrics: match args.receptor {
                VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                    VdjTSampleHeroMetricsTable::from_json_value(&metric_json_val),
                ),
                VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                    VdjTgdSampleHeroMetricsTable::from_json_value(&metric_json_val),
                ),
                VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                    VdjBSampleHeroMetricsTable::from_json_value(&metric_json_val),
                ),
            },
            annotation_metrics_table: match args.receptor {
                VdjReceptor::TR => VdjChainTypeSpecific::VdjT(
                    VdjTSampleAnnotationMetricsTable::from_json_value(&metric_json_val),
                ),
                VdjReceptor::TRGD => VdjChainTypeSpecific::VdjTgd(
                    VdjTgdSampleAnnotationMetricsTable::from_json_value(&metric_json_val),
                ),
                VdjReceptor::IG => VdjChainTypeSpecific::VdjB(
                    VdjBSampleAnnotationMetricsTable::from_json_value(&metric_json_val),
                ),
            },
            sequencing_metrics_table: args
                .sequencing_metrics
                .read()?
                .remove(&LegacyLibraryType::Vdj)
                .unwrap(),
            barcode_rank_plot,
            clonotype_info,
            filter_metrics: args.filter_metrics.map(|f| f.read()).transpose()?,
        };

        let vdj_ws_contents: VdjWsContentsFormat = rover.make_path("vdj_ws_contents");
        vdj_ws_contents.write(&contents)?;

        Ok(BuildVdjWsContentsStageOutputs { vdj_ws_contents })
    }
}
