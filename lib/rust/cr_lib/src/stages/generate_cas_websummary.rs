//! Martian stage GENERATE_CAS_WEBSUMMARY
//! Dump vega plots into cell annotation websummary
#![deny(missing_docs)]
use crate::HtmlFile;
use crate::cell_annotation_ws_parameters::{
    CellAnnotationMetrics, CellTypeWebSummaryBundle, CellTypeWebSummaryBundleValue,
    WebSummaryCellTypeInputs, generate_cas_bc_mismatch_alert, generate_cas_de_warn_alert,
    generate_cas_failure_alert, generate_cell_type_barcharts_from_json,
    generate_cell_type_diffexp_from_json, generate_cell_type_metrics,
    generate_cell_type_umap_plot_from_json, generate_cell_type_violin_plot_from_json,
};
use anyhow::{Ok, Result, bail};
use cr_types::constants::{COMMAND_LINE_ENV_DEFAULT_VALUE, COMMAND_LINE_ENV_VARIABLE_NAME};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};
use std::env;
use tenx_websummary::components::{
    Card, CommandLine, DifferentialExpressionTable, InlineHelp, MultiLayerImages, TableMetric,
    Tabs, TermDesc, Title, TitleWithTermDesc, VegaLitePlot, WithTitle, WsNavBar,
};
use tenx_websummary::{Alert, HtmlTemplate, SinglePageHtml};
use websummary_build::build_files;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct GenerateCasWebsummaryStageInputs {
    sample_id: Option<String>,
    sample_desc: Option<String>,
    alert_string: Option<String>,
    cell_type_websummary_bundle: Option<CellTypeWebSummaryBundle>,
    cell_type_websummary_bundles: Option<WebSummaryCellTypeInputs>,
    pipestance_type: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct GenerateCasWebsummaryStageOutputs {
    summary: Option<HtmlFile>,
    /// Flattened collection of all websummary data.
    ///
    /// Only produced when a single bundle was provided in the inputs.
    cell_type_websummary_bundle_value: Option<JsonFile<CellTypeWebSummaryBundleValue>>,
}

#[derive(Serialize, HtmlTemplate)]
pub struct WebSummaryContent {
    disclaimer_banner: Option<InlineHelp>,
    metrics_table: Card<WithTitle<TableMetric>>,
    cell_type_spatial_plot: Option<Card<WithTitle<MultiLayerImages>>>,
    cell_type_interactive_bar_chart: Option<Card<WithTitle<VegaLitePlot>>>,
    cell_types_box_plot: Option<Card<WithTitle<VegaLitePlot>>>,
    cell_types_umap_plot: Option<Card<WithTitle<VegaLitePlot>>>,
    diffexp_table: Option<Card<WithTitle<DifferentialExpressionTable>>>,
    cell_type_bar_chart: Option<WithTitle<VegaLitePlot>>,
    command_line_card: Option<Card<CommandLine>>,
}

pub const NO_ANNOTATE_WS_PIPESTANCE_TYPES: [&str; 4] = [
    "SC_MULTI_CORE_SAMPLE",
    "SPATIAL_RNA_COUNTER_HD_SEGMENTATION",
    "CELLRANGER_ANNOTATE_CS",
    "SPACERANGER_ANNOTATE_CS",
];

impl WebSummaryContent {
    fn from_stage_inputs(
        cell_type_websummary_bundle: &CellTypeWebSummaryBundle,
        cell_annotation_metrics: CellAnnotationMetrics,
        cmdline: &str,
    ) -> Result<Self> {
        // Use cell_type_spatial_plot from bundle (MultiLayerImages)
        let spatial_chart = cell_type_websummary_bundle
            .cell_type_spatial_plot
            .as_ref()
            .map(|chart| {
                Ok(Card::full_width(WithTitle {
                    title: TitleWithTermDesc {
                        title: "Cell Segmentation".into(),
                        data: vec![TermDesc::with_one_desc(
                            "Cell Segmentation With Cell Type Annotations",
                            "Cell boundaries are colored and filled by cell type assignment. \
                                    The tissue image is cropped to the capture area + 30 pixels.",
                        )],
                    }
                    .into(),
                    inner: chart.read()?,
                }))
            })
            .transpose()?;

        Ok(WebSummaryContent {
            disclaimer_banner: cell_annotation_metrics.generate_disclaiming_banner(),
            metrics_table: Card::full_width(generate_cell_type_metrics(cell_annotation_metrics)?),
            cell_type_interactive_bar_chart: cell_type_websummary_bundle.cell_type_interactive_bar_chart
                .as_ref()
                .map(|x| Ok(Card::full_width(generate_cell_type_barcharts_from_json(x)?)))
                .transpose()?,
            cell_types_box_plot: cell_type_websummary_bundle.cell_types_box_plot
                .as_ref()
                .map(|x| {
                    Ok(Card::full_width(generate_cell_type_violin_plot_from_json(
                        x,
                    )?))
                })
                .transpose()?,
            cell_types_umap_plot: cell_type_websummary_bundle.cell_types_umap_plot
                .as_ref()
                .map(|x| Ok(Card::full_width(generate_cell_type_umap_plot_from_json(x)?)))
                .transpose()?,
            diffexp_table: cell_type_websummary_bundle.diffexp
                .as_ref()
                .map(|x| Ok(Card::full_width(generate_cell_type_diffexp_from_json(x)?)))
                .transpose()?,
            cell_type_bar_chart: cell_type_websummary_bundle.cell_type_pd_plots.cell_type_bar_chart
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
            cell_type_spatial_plot: spatial_chart,
            command_line_card: Some(Card::full_width(CommandLine::new(cmdline)?)),
        })
    }
}

fn process_websummary_bundle(
    cell_type_websummary_bundle: &CellTypeWebSummaryBundle,
    cmdline: &str,
    pipeline_version: &str,
) -> Result<(impl HtmlTemplate + Serialize, Option<Vec<Alert>>)> {
    let metadata_metrics = cell_type_websummary_bundle
        .cell_annotation_metadata
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
        cell_annotation_model: cell_type_websummary_bundle
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
        cell_annotation_frac_returned_bcs: cell_type_websummary_bundle
            .cell_type_info
            .cas_frac_returned_bcs,
        cell_annotation_success: cell_type_websummary_bundle.cell_type_info.cas_success,
        cell_annotation_differential_expression: cell_type_websummary_bundle
            .cell_type_info
            .disable_differential_expression,
        pipeline_version: pipeline_version.to_string(),
        aggregate_and_annotate: cell_type_websummary_bundle.aggregate_and_annotate,
    };

    let web_summary_content = WebSummaryContent::from_stage_inputs(
        cell_type_websummary_bundle,
        cell_annotation_metrics,
        cmdline,
    )?;

    let initial_alerts = match (
        cell_type_websummary_bundle.cell_type_info.cas_success,
        cell_type_websummary_bundle
            .cell_type_info
            .disable_differential_expression,
    ) {
        (Some(false), _) => Some(vec![generate_cas_failure_alert()]),
        (Some(true), Some(true)) => Some(vec![generate_cas_de_warn_alert()]),
        _ => None,
    };

    Ok((web_summary_content, initial_alerts))
}

fn process_websummary_bundles(
    cell_type_websummary_bundle: &WebSummaryCellTypeInputs,
    cmdline: &str,
    pipeline_version: &str,
) -> Result<(impl HtmlTemplate + Serialize, Option<Vec<Alert>>)> {
    let mut tabs = Tabs::new();
    let mut cumulative_alerts = vec![];
    for (title, websummary_bundle) in &cell_type_websummary_bundle.cell_type_websummary_bundles {
        if let Some(websummary_bundle) = websummary_bundle {
            if websummary_bundle.cell_type_info.cas_success == Some(false) {
                continue;
            }
            let (websummary_content, alerts) =
                process_websummary_bundle(websummary_bundle, cmdline, pipeline_version)?;
            tabs.push(title, websummary_content);
            if let Some(mut alert_vec) = alerts {
                cumulative_alerts.append(&mut alert_vec);
            }
        }
    }
    Ok((
        tabs,
        if cumulative_alerts.is_empty() {
            None
        } else {
            Some(cumulative_alerts)
        },
    ))
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
        // Produce the flattened version of the input data.
        let output_bundle = args
            .cell_type_websummary_bundle
            .as_ref()
            .map(|bundle| {
                let val = bundle.read()?;
                rover
                    .make_path::<JsonFile<_>>("cell_type_websummary_bundle")
                    .with_content(&val)
            })
            .transpose()?;

        if let Some(cell_type_websummary_bundle) = args.cell_type_websummary_bundle.as_ref()
            && cell_type_websummary_bundle
                .cell_type_info
                .disable_cas_ws
                .unwrap_or(false)
        {
            return Ok(GenerateCasWebsummaryStageOutputs {
                summary: None,
                cell_type_websummary_bundle_value: output_bundle,
            });
        }
        if let Some(cell_type_websummary_bundles) = args.cell_type_websummary_bundles.as_ref()
            && cell_type_websummary_bundles
                .cell_type_websummary_bundles
                .iter()
                .all(|(_, bundle)| {
                    bundle
                        .as_ref()
                        .is_none_or(|b| b.cell_type_info.disable_cas_ws.unwrap_or(false))
                })
        {
            return Ok(GenerateCasWebsummaryStageOutputs {
                summary: None,
                cell_type_websummary_bundle_value: output_bundle,
            });
        }

        let pipeline_version = rover.pipelines_version();
        let cmdline = env::var(COMMAND_LINE_ENV_VARIABLE_NAME)
            .unwrap_or_else(|_| COMMAND_LINE_ENV_DEFAULT_VALUE.to_string());
        let nav_bar = WsNavBar {
            pipeline: "Cell Annotation".to_string(),
            id: args.sample_id.unwrap_or_default(),
            description: args.sample_desc.unwrap_or_default(),
        };
        let summary_html: HtmlFile = rover.make_path("minimal_websummary");
        let pipestance_type = args.pipestance_type.unwrap_or_default();

        let summary_output = match (
            args.cell_type_websummary_bundle,
            args.cell_type_websummary_bundles,
        ) {
            (Some(bundle), None) => {
                let (websummary_content, initial_alerts) =
                    process_websummary_bundle(&bundle, &cmdline, &pipeline_version)?;
                let final_alerts = match (initial_alerts, args.alert_string) {
                    (Some(alert_vec), None) => Some(alert_vec),
                    (Some(mut alert_vec), Some(alert_str)) => {
                        alert_vec.push(generate_cas_bc_mismatch_alert(alert_str));
                        Some(alert_vec)
                    }
                    (None, Some(alert_str)) => {
                        Some(vec![generate_cas_bc_mismatch_alert(alert_str)])
                    }
                    _ => None,
                };
                if !(NO_ANNOTATE_WS_PIPESTANCE_TYPES
                    .iter()
                    .any(|&ps_type| pipestance_type.contains(ps_type)))
                {
                    let html = SinglePageHtml::new(nav_bar, websummary_content, final_alerts);
                    html.generate_html_file_with_build_files(&summary_html, build_files()?)?;
                    Some(summary_html)
                } else {
                    None
                }
            }
            (None, Some(bundles)) => {
                let (websummary_content, initial_alerts) =
                    process_websummary_bundles(&bundles, &cmdline, &pipeline_version)?;
                let final_alerts = match (initial_alerts, args.alert_string) {
                    (Some(alert_vec), None) => Some(alert_vec),
                    (Some(mut alert_vec), Some(alert_str)) => {
                        alert_vec.push(generate_cas_bc_mismatch_alert(alert_str));
                        Some(alert_vec)
                    }
                    (None, Some(alert_str)) => {
                        Some(vec![generate_cas_bc_mismatch_alert(alert_str)])
                    }
                    _ => None,
                };

                let html = SinglePageHtml::new(nav_bar, websummary_content, final_alerts);
                html.generate_html_file_with_build_files(&summary_html, build_files()?)?;
                Some(summary_html)
            }

            _ => bail!(
                "Either cell_type_websummary_bundle or cell_type_websummary_bundles with at least one element must be provided."
            ),
        };

        Ok(GenerateCasWebsummaryStageOutputs {
            summary: summary_output,
            cell_type_websummary_bundle_value: output_bundle,
        })
    }
}
