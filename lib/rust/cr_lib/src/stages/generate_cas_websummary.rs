//! Martian stage GENERATE_CAS_WEBSUMMARY
//! Dump vega plots into cell annotation websummary
#![deny(missing_docs)]
use crate::HtmlFile;
use crate::cell_annotation_ws_parameters::{
    CellAnnotationMetrics, CellTypeWebSummaryBundle, generate_cas_bc_mismatch_alert,
    generate_cas_de_warn_alert, generate_cas_failure_alert, generate_cell_type_barcharts_from_json,
    generate_cell_type_diffexp_from_json, generate_cell_type_metrics,
    generate_cell_type_umap_plot_from_json, generate_cell_type_violin_plot_from_json,
};
use anyhow::{Ok, Result};
use cr_types::constants::{COMMAND_LINE_ENV_DEFAULT_VALUE, COMMAND_LINE_ENV_VARIABLE_NAME};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::env;
use tenx_websummary::components::{
    Card, CommandLine, DifferentialExpressionTable, InlineHelp, TableMetric, Title, VegaLitePlot,
    WithTitle, WsNavBar,
};
use tenx_websummary::{HtmlTemplate, SinglePageHtml};
use websummary_build::build_files;

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct MetadataStruct {
    tree_version_used: Option<String>,
    display_map_version_used: Option<String>,
    fraction_non_informative_annotations: Option<f64>,
    is_beta_model: Option<bool>,
    developer: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct GenerateCasWebsummaryStageInputs {
    sample_id: Option<String>,
    sample_desc: Option<String>,
    cas_frac_returned_bcs: Option<f64>,
    metadata: Option<JsonFile<MetadataStruct>>,
    cell_annotation_model: Option<String>,
    cell_type_bar_chart: Option<JsonFile<Value>>,
    spatial_cell_types_chart: Option<JsonFile<Value>>,
    cell_type_websummary_bundle: CellTypeWebSummaryBundle,
    cas_success: Option<bool>,
    disable_differential_expression: Option<bool>,
    alert_string: Option<String>,
    disable_cas_ws: Option<bool>,
    pipestance_type: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct GenerateCasWebsummaryStageOutputs {
    summary: Option<HtmlFile>,
    cell_annotation_metrics: Option<JsonFile<CellAnnotationMetrics>>,
    cell_annotation_cloupe_name: Option<String>,
}

#[derive(Serialize, HtmlTemplate)]
pub struct WebSummaryContent {
    disclaimer_banner: Option<InlineHelp>,
    metrics_table: Card<WithTitle<TableMetric>>,
    cell_type_interactive_bar_chart: Option<Card<WithTitle<VegaLitePlot>>>,
    cell_types_box_plot: Option<Card<WithTitle<VegaLitePlot>>>,
    cell_types_umap_plot: Option<Card<WithTitle<VegaLitePlot>>>,
    diffexp_table: Option<Card<WithTitle<DifferentialExpressionTable>>>,
    cell_type_bar_chart: Option<WithTitle<VegaLitePlot>>,
    spatial_cell_types_chart: Option<WithTitle<VegaLitePlot>>,
    command_line_card: Option<Card<CommandLine>>,
}

pub const NO_ANNOTATE_WS_PIPESTANCE_TYPES: [&str; 1] = ["SC_MULTI_CORE_SAMPLE"];

impl WebSummaryContent {
    fn from_stage_inputs(
        args: &GenerateCasWebsummaryStageInputs,
        cell_annotation_metrics: CellAnnotationMetrics,
        cmdline: &str,
    ) -> Result<Self> {
        Ok(WebSummaryContent {
            disclaimer_banner: cell_annotation_metrics.generate_disclaiming_banner(),
            metrics_table: Card::full_width(generate_cell_type_metrics(cell_annotation_metrics)?),
            cell_type_interactive_bar_chart: args.cell_type_websummary_bundle.cell_type_interactive_bar_chart
                .as_ref()
                .map(|x| Ok(Card::full_width(generate_cell_type_barcharts_from_json(x)?)))
                .transpose()?,
            cell_types_box_plot: args.cell_type_websummary_bundle.cell_types_box_plot
                .as_ref()
                .map(|x| {
                    Ok(Card::full_width(generate_cell_type_violin_plot_from_json(
                        x,
                    )?))
                })
                .transpose()?,
            cell_types_umap_plot: args.cell_type_websummary_bundle.cell_types_umap_plot
                .as_ref()
                .map(|x| Ok(Card::full_width(generate_cell_type_umap_plot_from_json(x)?)))
                .transpose()?,
            diffexp_table: args.cell_type_websummary_bundle.diffexp
                .as_ref()
                .map(|x| Ok(Card::full_width(generate_cell_type_diffexp_from_json(x)?)))
                .transpose()?,
            cell_type_bar_chart: args
                .cell_type_bar_chart
                .as_ref()
                .map(|x| {
                    Ok(WithTitle {
                        title: Title::new(
                            "cell annotation cell types and clusters discovered by unbiased clustering"
                                .to_string(),
                        ),
                        inner: VegaLitePlot {
                            spec: x.read()?,
                            actions: None,
                            renderer: None,
                        },
                    })
                })
                .transpose()?,
            spatial_cell_types_chart: args
                .spatial_cell_types_chart
                .as_ref()
                .map(|x| {
                    Ok(WithTitle {
                        title: Title::new("Spatial distribution of annotated spots".to_string()),
                        inner: VegaLitePlot {
                            spec: x.read()?,
                            actions: None,
                            renderer: None,
                        },
                    })
                })
                .transpose()?,
            command_line_card: Some(Card::full_width(CommandLine::new(cmdline)?)),
        })
    }
}

/// Martian stage GENERATE_CAS_WEBSUMMARY
pub struct GenerateCasWebsummary;

#[make_mro]
impl MartianMain for GenerateCasWebsummary {
    type StageInputs = GenerateCasWebsummaryStageInputs;
    type StageOutputs = GenerateCasWebsummaryStageOutputs;
    fn main(
        &self,
        args: Self::StageInputs,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        if args.disable_cas_ws.unwrap_or(false) {
            return Ok(GenerateCasWebsummaryStageOutputs {
                summary: None,
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
            cell_annotation_frac_returned_bcs: args.cas_frac_returned_bcs,
            cell_annotation_success: args.cas_success,
            cell_annotation_differential_expression: args.disable_differential_expression,
            pipeline_version: rover.pipelines_version(),
        };

        let cell_annotation_cloupe_name = Some(cell_annotation_metrics.get_cloupe_track_name());

        let cmdline = env::var(COMMAND_LINE_ENV_VARIABLE_NAME)
            .unwrap_or_else(|_| COMMAND_LINE_ENV_DEFAULT_VALUE.to_string());

        let json: JsonFile<CellAnnotationMetrics> = rover.make_path("cell_annotation_metrics");
        json.write(&cell_annotation_metrics)?;

        let web_summary_content =
            WebSummaryContent::from_stage_inputs(&args, cell_annotation_metrics, &cmdline)?;

        let nav_bar = WsNavBar {
            pipeline: "Cell Annotation".to_string(),
            id: args.sample_id.unwrap_or_default(),
            description: args.sample_desc.unwrap_or_default(),
        };
        let summary_html: HtmlFile = rover.make_path("minimal_websummary");
        let initial_alerts = match (args.cas_success, args.disable_differential_expression) {
            (Some(false), _) => Some(vec![generate_cas_failure_alert()]),
            (Some(true), Some(true)) => Some(vec![generate_cas_de_warn_alert()]),
            _ => None,
        };

        let final_alerts = match (initial_alerts, args.alert_string) {
            (Some(alert_vec), None) => Some(alert_vec),
            (Some(mut alert_vec), Some(alert_str)) => {
                alert_vec.push(generate_cas_bc_mismatch_alert(alert_str));
                Some(alert_vec)
            }
            (None, Some(alert_str)) => Some(vec![generate_cas_bc_mismatch_alert(alert_str)]),
            _ => None,
        };

        let pipestance_type = args.pipestance_type.unwrap_or_default();

        let summary_output = if !(NO_ANNOTATE_WS_PIPESTANCE_TYPES
            .iter()
            .any(|&ps_type| pipestance_type.contains(ps_type)))
        {
            let html = SinglePageHtml::new(nav_bar, web_summary_content, final_alerts);
            html.generate_html_file_with_build_files(&summary_html, build_files()?)?;
            Some(summary_html)
        } else {
            None
        };

        Ok(GenerateCasWebsummaryStageOutputs {
            summary: summary_output,
            cell_annotation_metrics: Some(json),
            cell_annotation_cloupe_name,
        })
    }
}
