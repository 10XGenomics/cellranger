//! Martian stage GENERATE_CAS_METRICS
//! Generate cell annotation metrics json
#![deny(missing_docs)]
use crate::cell_annotation_ws_parameters::{CellAnnotationMetrics, CellTypeInfo, MetadataStruct};
use anyhow::{Ok, Result};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct GenerateCasMetricsStageInputs {
    cell_type_info: CellTypeInfo,
    metadata: Option<JsonFile<MetadataStruct>>,
    aggregate_and_annotate: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct GenerateCasMetricsStageOutputs {
    cell_annotation_metrics: Option<JsonFile<CellAnnotationMetrics>>,
    cell_annotation_cloupe_name: Option<String>,
}

/// Martian stage GENERATE_CAS_METRICS
pub struct GenerateCasMetrics;

#[make_mro(volatile = strict)]
impl MartianMain for GenerateCasMetrics {
    type StageInputs = GenerateCasMetricsStageInputs;
    type StageOutputs = GenerateCasMetricsStageOutputs;
    fn main(
        &self,
        args: Self::StageInputs,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        if args.cell_type_info.disable_cas_ws.unwrap_or(false) {
            return Ok(GenerateCasMetricsStageOutputs {
                cell_annotation_metrics: None,
                cell_annotation_cloupe_name: None,
            });
        }

        let metadata_metrics = args
            .metadata
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?
            .unwrap_or_default();

        let cell_annotation_metrics = CellAnnotationMetrics {
            cell_annotation_beta_model: metadata_metrics.is_beta_model != Some(false),
            cell_annotation_model_developer: metadata_metrics
                .developer
                .clone()
                .unwrap_or_else(|| "Unknown".to_string()),
            cell_annotation_model: args
                .cell_type_info
                .cell_annotation_model
                .clone()
                .unwrap_or_else(|| "NOT FOUND".to_string()),
            cell_annotation_tree_version_used: metadata_metrics
                .tree_version_used
                .clone()
                .unwrap_or_else(|| "NOT FOUND".to_string()),
            cell_annotation_display_map_version_used: metadata_metrics
                .display_map_version_used
                .unwrap_or_else(|| "NOT FOUND".to_string()),
            cell_annotation_frac_returned_bcs: args.cell_type_info.cas_frac_returned_bcs,
            cell_annotation_success: args.cell_type_info.cas_success,
            cell_annotation_differential_expression: args
                .cell_type_info
                .disable_differential_expression,
            pipeline_version: rover.pipelines_version(),
            aggregate_and_annotate: args.aggregate_and_annotate,
        };
        let cell_annotation_cloupe_name = Some(cell_annotation_metrics.get_cloupe_track_name());

        let json: JsonFile<CellAnnotationMetrics> = rover.make_path("cell_annotation_metrics");
        json.write(&cell_annotation_metrics)?;

        Ok(GenerateCasMetricsStageOutputs {
            cell_annotation_metrics: Some(json),
            cell_annotation_cloupe_name,
        })
    }
}
