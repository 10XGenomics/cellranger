// Cell Annotation Web Summary Parameters used across multi, count, and spatial
#![deny(missing_docs)]
use anyhow::{Ok, Result};
use martian_derive::MartianStruct;
use martian_filetypes::FileTypeRead;
use martian_filetypes::json_file::JsonFile;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use tenx_websummary::components::{
    DifferentialExpressionTable, InlineHelp, MultiLayerImages, TableMetric, TermDesc,
    TitleWithHelp, TitleWithTermDesc, VegaLitePlot, WithTitle,
};
use tenx_websummary::{Alert, AlertLevel};

const CELL_ANNOTATION_DISCLAIMER_STRING: &str = r"<p><b>Important Note</b><br>
Cell type annotation is currently in beta and relies on the Chan Zuckerberg CELL by GENE reference.
As this reference is community-driven it may not cover all tissue types, which could impact your results. For further details, please visit the
<a href='https://10xgen.com/cell-annotation' target='_blank' title='10x Support: Cell Annotation' rel='noopener noreferrer'>10x support site</a>.
</p>";

const CELL_ANNOTATION_BROAD_MODEL: &str = "Annotation was performed using a model that was co-developed by 10x Genomics and the Cellarium AI Lab at the Data Sciences Platform of the Broad Institute.";

const CELL_ANNOTATION_DISCLAIMER_STRING_AZIMUTH: &str = r"<p>
<b>Important Note</b>
<br>
Cell type annotation using Pan-human Azimuth is currently in beta and relies on <a href='https://github.com/satijalab/panhumanpy' target='_blank' title='panhumanpy github' rel='noopener noreferrer'>panhumanpy</a>.
This model may not cover all tissue types, which could impact your results. Please review the documentation and <a href='https://10xgen.com/1We' target='_blank' title='Pan Human Azimuth Publication' rel='noopener noreferrer'>publication</a>
for supported data types and examine results from other tissues carefully. For further details, please visit the
<a href='https://10xgen.com/cell-annotation' target='_blank' title='10x Support: Cell Annotation' rel='noopener noreferrer'>10x support site</a>.
</p>";

const CELL_ANNOTATION_FAILURE_TITLE: &str = "No cell annotations produced!";

const CELL_ANNOTATION_FAILURE_MESSAGE: &str = r"<p>Please check your cellranger annotate or spaceranger annotate logs.
If you wish to attempt cell type annotation again please use
<a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-cell-annotation-pipeline'
target='_blank' title='Space Ranger Annotate' rel='noopener noreferrer'>cellranger annotate</a> or <a href='https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/space-ranger-annotate'
target='_blank' title='Cell Ranger Annotate' rel='noopener noreferrer'>spaceranger annotate</a> .</p>";

const CELL_ANNOTATION_DE_WARN_TITLE: &str = "Cell type differential expression not run";
const CELL_ANNOTATION_DE_WARN_MESSAGE: &str = "Too few cell types to run differential expression.";
const CELL_ANNOTATION_BC_MISMATCH_WARN_TITLE: &str =
    "Barcodes in cloupe file and feature barcode matrix do not match";

const CELL_ANNOTATION_METRICS_TABLE_TITLE: &str = "Cell Annotation Summary";
const CELL_ANNOTATION_METRICS_TABLE_HELP: &str = "<b>Pipeline version</b><br>
The version of Cell Ranger or Space Ranger used to generate this annotation.<br>
<b>Cell annotation model</b><br>
The model used to generate cell type annotations.<br>
<b>Annotate meta-cells</b><br>
For annotation, whether to group cells with at least 30 UMIs and at most 500 UMIs with similar gene expression into meta-cells of at least 500 UMIs.";

const CELL_ANNOTATION_INTERACTIVE_BARCHART_TITLE: &str = "Cell Type Composition";
const CELL_ANNOTATION_INTERACTIVE_BARCHART_HELP: &str = "This plot shows the major cell types in your sample. To view sub-type annotations, click on a bar.
Only major classifications with at least 10 cells are shown.<br>
<b>Left-hand bar chart</b>: This bar chart shows the annotated cell types on the y-axis and the number of barcodes annotated with them on the x-axis.<br>
<b>Right-hand table</b>: This table shows the sub-types within a selected major cell type. Fractions of barcodes shown are the percentage of cells relative to
this selected cell type—not the entire sample.";

const CELL_ANNOTATION_VIOLIN_PLOT_TITLE: &str = "UMI distribution by cell type";
const CELL_ANNOTATION_VIOLIN_PLOT_HELP: &str =
    "Box plot showing the distribution of (UMI+1) by annotated cell types.<br>
    Only cell types with >= 10 barcodes included. If a cell type has more<br>
    than 1000 barcodes associated with it the group is subsampled to 1000 barcodes.<br>
    Minimum and maximum values are kept. Y-axis is log scale";

const CELL_ANNOTATION_UMAP_PLOT_TITLE: &str = "UMAP projection of cell types";
const CELL_ANNOTATION_UMAP_PLOT_HELP: &str =
    "UMAP projection of annotated cells. Click legend to highlight<br>
     a single cell type. Shift+Click to highlight multiple cell types.<br>
     Double-click to reset plot.";

const CELL_ANNOTATION_DIFFEXP_TITLE: &str = "Top Features by Cell Type";
const CELL_ANNOTATION_DIFFEXP_HELP: &str =
    "Differential expression analysis identifies, for each cell type,  genes that are more highly expressed in that cell type relative to the rest of the sample.<br>
    <b>Log2 Fold-Change (L2FC)</b>: An estimate of the log2 ratio of expression of a gene in a cell type to the mean of all other cell types. A value of 1 indicates
    2-fold greater expression in the cell type of interest.<br>
    <b>P-Value</b>: A measure of the statistical significance of the expression difference, after correcting for multiple tests to avoid false positives.<br>
    <b>How to use this table</b>: Click on any column to sort data by that feature. Genes not meeting our criteria (L2FC < 0 or FDR-adjusted p-value >= 0.10)
    are shown in gray. For a complete look at the data, check the 'differential_expression.csv' files from our analysis pipeline.";

const CELL_TYPING_CLOUPE_NAME: &str = "Cell Types";
const CELL_TYPING_BETA_SUFFIX: &str = "(beta)";

/// Disclaimer banner content.
#[derive(Serialize, Deserialize)]
pub struct DisclaimerBanner {
    /// Raw HTML content of the disclaimer.
    pub html: String,
}

/// Interactive cell types plot
#[derive(Serialize, Deserialize)]
pub struct CellAnnotationPlot {
    /// Interactive cell types plot help text
    pub help: TitleWithHelp,
    /// Vega-Lite specification.
    pub spec: Value,
}

/// Cell Annotation metrics table
#[derive(Serialize, Deserialize)]
pub struct CellAnnotationMetricsTableData {
    /// Metrics help text
    pub help: TitleWithTermDesc,
    /// Two-column key/value table data.
    pub table: TableMetric,
}

/// Differential expression table with a help header, used as a prop in the CellTypesTab React component.
///
/// The `table_data` fields (`table`, `default_page_size`, `allow_fixed_columns`) are flattened
/// alongside `help` so JavaScript can destructure them directly.
#[derive(Serialize, Deserialize)]
pub struct CellAnnotationDiffexpTableData {
    /// Title and HTML help text rendered by `HeaderWithHelp.js`.
    pub help: TitleWithHelp,
    /// DifferentialExpressionTable props, flattened alongside `help`.
    #[serde(flatten)]
    pub table_data: DifferentialExpressionTable,
}

/// Pre-formatted props for the `CellTypesTab` React component in the count websummary.
///
/// Fields correspond directly to the props consumed by `CellTypesTab.js`.
/// Python passes them through without any reformatting.
#[derive(Serialize, Deserialize)]
pub struct CellTypesTabData {
    /// Tab title shown in the tab bar.
    pub title: String,
    /// Optional disclaimer banner rendered above the metrics table for beta models.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_annotation_disclaimer_banner: Option<DisclaimerBanner>,
    /// Cell annotation summary metrics table.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_annotation_metrics_table: Option<CellAnnotationMetricsTableData>,
    /// Interactive bar chart of cell type composition.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_type_bar_chart: Option<CellAnnotationPlot>,
    /// UMAP projection of annotated cells.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_type_umap_plot: Option<CellAnnotationPlot>,
    /// UMI distribution box plot by cell type.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_cells_bar_chart: Option<CellAnnotationPlot>,
    /// Differential expression table.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diffexp_table: Option<CellAnnotationDiffexpTableData>,
}

/// Builds the pre-formatted tab props for the `CellTypesTab` React component.
fn build_cell_types_tab_data(
    model_family: &str,
    metrics: Option<&CellAnnotationMetrics>,
    bar_chart: Option<&Value>,
    umap_plot: Option<&Value>,
    box_plot: Option<&Value>,
    diffexp: Option<&DifferentialExpressionTable>,
) -> Option<CellTypesTabData> {
    if metrics.is_none()
        && bar_chart.is_none()
        && umap_plot.is_none()
        && box_plot.is_none()
        && diffexp.is_none()
    {
        return None;
    }
    let model_family_display = if model_family.is_empty() {
        "Cell Types"
    } else {
        model_family
    };
    let title = format!("Cell Annotation: {model_family_display}");
    let cell_annotation_disclaimer_banner = metrics.and_then(|m| {
        m.generate_disclaimer_html_fragment()
            .map(|html| DisclaimerBanner { html })
    });
    let cell_annotation_metrics_table = metrics.map(|m| CellAnnotationMetricsTableData {
        help: TitleWithTermDesc {
            title: CELL_ANNOTATION_METRICS_TABLE_TITLE.to_string(),
            data: vec![
                TermDesc(
                    "Pipeline version".to_string(),
                    vec![
                        "The version of Cell Ranger used to generate this annotation.".to_string(),
                    ],
                ),
                TermDesc(
                    "Cell annotation model".to_string(),
                    vec!["The model used to generate cell type annotations.".to_string()],
                ),
            ],
        },
        table: TableMetric {
            rows: m.get_table(),
        },
    });
    let cell_type_bar_chart = bar_chart.map(|spec| CellAnnotationPlot {
        help: TitleWithHelp {
            title: CELL_ANNOTATION_INTERACTIVE_BARCHART_TITLE.to_string(),
            help: CELL_ANNOTATION_INTERACTIVE_BARCHART_HELP.to_string(),
        },
        spec: spec.clone(),
    });
    let cell_type_umap_plot = umap_plot.map(|spec| CellAnnotationPlot {
        help: TitleWithHelp {
            title: CELL_ANNOTATION_UMAP_PLOT_TITLE.to_string(),
            help: CELL_ANNOTATION_UMAP_PLOT_HELP.to_string(),
        },
        spec: spec.clone(),
    });
    let total_cells_bar_chart = box_plot.map(|spec| CellAnnotationPlot {
        help: TitleWithHelp {
            title: CELL_ANNOTATION_VIOLIN_PLOT_TITLE.to_string(),
            help: CELL_ANNOTATION_VIOLIN_PLOT_HELP.to_string(),
        },
        spec: spec.clone(),
    });
    let diffexp_table = diffexp.map(|d| CellAnnotationDiffexpTableData {
        help: TitleWithHelp {
            title: CELL_ANNOTATION_DIFFEXP_TITLE.to_string(),
            help: CELL_ANNOTATION_DIFFEXP_HELP.to_string(),
        },
        table_data: d.clone(),
    });
    Some(CellTypesTabData {
        title,
        cell_annotation_disclaimer_banner,
        cell_annotation_metrics_table,
        cell_type_bar_chart,
        cell_type_umap_plot,
        total_cells_bar_chart,
        diffexp_table,
    })
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
/// Plots for cell type annotation results.
pub struct CellTypePdPlots {
    /// Interactive bar chart of cell type composition.
    pub cell_type_bar_chart: Option<JsonFile<Value>>,
    /// Spatial cell type chart.
    pub spatial_cell_types_chart: Option<JsonFile<Value>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
/// Information about cell type annotation run used in websummary or downstream analysis.
pub struct CellTypeInfo {
    /// Fraction of barcodes that were successfully annotated.
    pub cas_frac_returned_bcs: Option<f64>,
    /// Cell annotation model used.
    pub cell_annotation_model: Option<String>,
    /// If cell annotation run was successful.
    pub cas_success: Option<bool>,
    /// If differential should be disabled.
    pub disable_differential_expression: Option<bool>,
    /// If cell annotation websummary generation should be disabled.
    pub disable_cas_ws: Option<bool>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
/// Metadata related to cell type annotation.
pub struct MetadataStruct {
    /// Cell annotation tree used.
    pub tree_version_used: Option<String>,
    /// Cell annotation display map used.
    pub display_map_version_used: Option<String>,
    /// Fraction of non-informative annotations.
    pub fraction_non_informative_annotations: Option<f64>,
    /// Whether beta model was used.
    pub is_beta_model: Option<bool>,
    /// Model developer 10x, Broad, etc.
    pub developer: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
/// Inputs for generating cell type websummary across multiple samples.
pub struct WebSummaryCellTypeInputs {
    /// Cell type websummary bundles for each sample.
    pub cell_type_websummary_bundles: TxHashMap<String, Option<CellTypeWebSummaryBundle>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
/// Bundle of files for plots and tables used to generate cell type websummary.
pub struct CellTypeWebSummaryBundle {
    /// Cell annotation model family used.
    pub cell_annotation_model_family: String,
    /// Interactive bar chart of cell type composition.
    pub cell_type_interactive_bar_chart: Option<JsonFile<Value>>,
    /// Box plot of cell type UMI counts.
    pub cell_types_box_plot: Option<JsonFile<Value>>,
    /// UMAP plot of cell types.
    pub cell_types_umap_plot: Option<JsonFile<Value>>,
    /// Differential expression table.
    pub diffexp: Option<JsonFile<DifferentialExpressionTable>>,
    /// Cell annotation metrics.
    pub cell_annotation_metrics: Option<JsonFile<CellAnnotationMetrics>>,
    /// Cell annotation PD plots.
    pub cell_type_pd_plots: CellTypePdPlots,
    /// Cell annotation metadata.
    pub cell_annotation_metadata: Option<JsonFile<MetadataStruct>>,
    /// Cell type information.
    pub cell_type_info: CellTypeInfo,
    /// Spatial plot of cell types.
    pub cell_type_spatial_plot: Option<JsonFile<MultiLayerImages>>,
    /// Whether metacell aggregation was used before annotation.
    pub aggregate_and_annotate: Option<bool>,
}

impl CellTypeWebSummaryBundle {
    /// Reads the referenced files in the bundle and loads their content into memory.
    /// Constructs a `CellTypeWebSummaryBundleValue` containing the deserialized data.
    pub fn read(&self) -> Result<CellTypeWebSummaryBundleValue> {
        let cell_type_interactive_bar_chart = self
            .cell_type_interactive_bar_chart
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?;
        let cell_types_box_plot = self
            .cell_types_box_plot
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?;
        let cell_types_umap_plot = self
            .cell_types_umap_plot
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?;
        let diffexp = self.diffexp.as_ref().map(FileTypeRead::read).transpose()?;
        let cell_annotation_metrics = self
            .cell_annotation_metrics
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?;
        let cell_type_spatial_plot = self
            .cell_type_spatial_plot
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?;
        let cell_types_tab_data = if self.cell_type_info.disable_cas_ws.unwrap_or(false) {
            None
        } else {
            build_cell_types_tab_data(
                &self.cell_annotation_model_family,
                cell_annotation_metrics.as_ref(),
                cell_type_interactive_bar_chart.as_ref(),
                cell_types_umap_plot.as_ref(),
                cell_types_box_plot.as_ref(),
                diffexp.as_ref(),
            )
        };
        Ok(CellTypeWebSummaryBundleValue {
            cell_type_interactive_bar_chart,
            cell_types_box_plot,
            cell_types_umap_plot,
            diffexp,
            cell_annotation_metrics,
            cell_annotation_model_family: self.cell_annotation_model_family.clone(),
            cell_type_spatial_plot,
            cell_types_tab_data,
        })
    }
}

#[derive(Serialize, Deserialize)]
/// Values read from CellTypeWebSummaryBundle files.
/// A collection of web summary visualization data and metrics for cell type annotation.
pub struct CellTypeWebSummaryBundleValue {
    /// The name or identifier of the cell annotation model family used
    pub cell_annotation_model_family: String,
    /// JSON value containing data for an interactive bar chart visualization of cell types
    pub cell_type_interactive_bar_chart: Option<Value>,
    /// JSON value containing data for a box plot visualization of cell types
    pub cell_types_box_plot: Option<Value>,
    /// JSON value containing data for a UMAP plot visualization with cell type annotations
    pub cell_types_umap_plot: Option<Value>,
    /// Differential expression analysis results table comparing cell types
    pub diffexp: Option<DifferentialExpressionTable>,
    /// Metrics summarizing the performance and quality of cell type annotations
    pub cell_annotation_metrics: Option<CellAnnotationMetrics>,
    /// Spatial segmentation plot with cell type annotations
    pub cell_type_spatial_plot: Option<MultiLayerImages>,
    /// Pre-formatted props for the CellTypesTab React component; `None` when websummary is disabled or nothing to display
    pub cell_types_tab_data: Option<CellTypesTabData>,
}

impl CellTypeWebSummaryBundleValue {
    /// Checks if any of the visualization or metric components in the bundle are present.
    ///
    /// Returns `true` if at least one of the major components (bar chart, box plot,
    /// UMAP, differential expression, or metrics) is `Some`.
    pub fn has_any_value(&self) -> bool {
        self.cell_type_interactive_bar_chart.is_some()
            || self.cell_types_box_plot.is_some()
            || self.cell_types_umap_plot.is_some()
            || self.diffexp.is_some()
            || self.cell_annotation_metrics.is_some()
            || self.cell_type_spatial_plot.is_some()
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
/// Metrics and models related to cell type annotation used in websummary and downstream analysis.
pub struct CellAnnotationMetrics {
    /// The specific name of the cell annotation model used.
    pub cell_annotation_model: String,
    /// The version of the cell ontology tree used for classification.
    pub cell_annotation_tree_version_used: String,
    /// The version of the display map configuration used.
    pub cell_annotation_display_map_version_used: String,
    /// The fraction of input barcodes that were successfully assigned a cell type.
    pub cell_annotation_frac_returned_bcs: Option<f64>,
    /// Whether the cell annotation process completed successfully.
    pub cell_annotation_success: Option<bool>,
    /// Whether differential expression analysis was performed on the annotated groups.
    pub cell_annotation_differential_expression: Option<bool>,
    /// Indicates if the model used is considered a beta/experimental version.
    pub cell_annotation_beta_model: bool,
    /// The organization or entity that developed the model (e.g., "10x", "Broad").
    pub cell_annotation_model_developer: String,
    /// The version of the analysis pipeline used to generate these metrics.
    pub pipeline_version: String,
    /// Whether metacell aggregation was used before annotation.
    pub aggregate_and_annotate: Option<bool>,
}

impl CellAnnotationMetrics {
    /// Formats the key metrics into a vector of label-value pairs for display in a summary table.
    fn get_table(&self) -> Vec<(String, String)> {
        let mut rows = vec![
            (
                "Pipeline version".to_string(),
                self.pipeline_version.clone(),
            ),
            (
                "Cell annotation model".to_string(),
                self.cell_annotation_model.clone(),
            ),
        ];
        if let Some(agg) = self.aggregate_and_annotate {
            rows.push((
                "Annotate meta-cells".to_string(),
                if agg {
                    "Yes".to_string()
                } else {
                    "No".to_string()
                },
            ));
        }
        rows
    }
    /// Generates an inline help banner containing a disclaimer if a beta model was used.
    ///
    /// Returns `Some(InlineHelp)` if `is_beta_model` is true, otherwise returns `None`.
    /// Customized based on the model developer (e.g., specific text for Broad Institute models).
    pub fn generate_disclaiming_banner(&self) -> Option<InlineHelp> {
        if self.cell_annotation_beta_model {
            let disclaimer = if self.cell_annotation_model_developer.contains("Broad") {
                format!("{CELL_ANNOTATION_DISCLAIMER_STRING}\n{CELL_ANNOTATION_BROAD_MODEL}")
            } else {
                CELL_ANNOTATION_DISCLAIMER_STRING_AZIMUTH.to_string()
            };

            Some(InlineHelp::with_content(disclaimer))
        } else {
            None
        }
    }
    /// Generates a raw HTML fragment string for the disclaimer if a beta model was used.
    ///
    /// This is similar to `generate_disclaiming_banner` but returns a raw String
    /// suitable for embedding directly in HTML templates.
    pub fn generate_disclaimer_html_fragment(&self) -> Option<String> {
        if self.cell_annotation_beta_model {
            let disclaimer = if self.cell_annotation_model_developer.contains("Broad") {
                format!("{CELL_ANNOTATION_DISCLAIMER_STRING}\n{CELL_ANNOTATION_BROAD_MODEL}")
            } else {
                CELL_ANNOTATION_DISCLAIMER_STRING_AZIMUTH.to_string()
            };

            Some(disclaimer)
        } else {
            None
        }
    }
    /// Returns the appropriate track name for Cloupe visualization.
    ///
    /// Appends a "Beta" suffix to the track name if `is_beta_model` is true.
    pub fn get_cloupe_track_name(&self) -> String {
        if self.cell_annotation_beta_model {
            format!("{CELL_TYPING_CLOUPE_NAME} {CELL_TYPING_BETA_SUFFIX}")
        } else {
            CELL_TYPING_CLOUPE_NAME.to_string()
        }
    }
}
/// Boilerplate build vega plot
fn build_vega_plot(title: &str, help: &str, spec: Value) -> Result<WithTitle<VegaLitePlot>> {
    Ok(WithTitle {
        title: TitleWithHelp {
            title: title.to_string(),
            help: help.to_string(),
        }
        .into(),
        inner: VegaLitePlot {
            spec,
            actions: Some(Value::Bool(false)),
            renderer: None,
        },
    })
}

/// Reads a JSON file containing bar chart data and converts it into a VegaLite plot component with a title.
pub fn generate_cell_type_barcharts_from_json(
    json_file: &JsonFile<Value>,
) -> Result<WithTitle<VegaLitePlot>> {
    generate_cell_type_barchart_from_value(json_file.read()?)
}
/// Wraps a raw JSON Value for a bar chart into a VegaLite plot component with a standard title and help text.
pub fn generate_cell_type_barchart_from_value(value: Value) -> Result<WithTitle<VegaLitePlot>> {
    build_vega_plot(
        CELL_ANNOTATION_INTERACTIVE_BARCHART_TITLE,
        CELL_ANNOTATION_INTERACTIVE_BARCHART_HELP,
        value,
    )
}
/// Reads a JSON file containing violin plot data and converts it into a VegaLite plot component with a title.
pub fn generate_cell_type_violin_plot_from_json(
    json_file: &JsonFile<Value>,
) -> Result<WithTitle<VegaLitePlot>> {
    generate_cell_type_violin_plot_from_value(json_file.read()?)
}
/// Wraps a raw JSON Value for a violin plot into a VegaLite plot component with a standard title and help text.
pub fn generate_cell_type_violin_plot_from_value(value: Value) -> Result<WithTitle<VegaLitePlot>> {
    build_vega_plot(
        CELL_ANNOTATION_VIOLIN_PLOT_TITLE,
        CELL_ANNOTATION_VIOLIN_PLOT_HELP,
        value,
    )
}
/// Reads a JSON file containing UMAP plot data and converts it into a VegaLite plot component with a title.
pub fn generate_cell_type_umap_plot_from_json(
    json_file: &JsonFile<Value>,
) -> Result<WithTitle<VegaLitePlot>> {
    generate_cell_type_umap_plot_from_value(json_file.read()?)
}
/// Wraps a raw JSON Value for a UMAP plot into a VegaLite plot component with a standard title and help text.
pub fn generate_cell_type_umap_plot_from_value(value: Value) -> Result<WithTitle<VegaLitePlot>> {
    build_vega_plot(
        CELL_ANNOTATION_UMAP_PLOT_TITLE,
        CELL_ANNOTATION_UMAP_PLOT_HELP,
        value,
    )
}
/// Reads a JSON file containing differential expression data and converts it into a table component with a title.
pub fn generate_cell_type_diffexp_from_json(
    json_file: &JsonFile<DifferentialExpressionTable>,
) -> Result<WithTitle<DifferentialExpressionTable>> {
    generate_cell_type_diffexp_from_value(json_file.read()?)
}
/// Wraps a raw DifferentialExpressionTable into a component with a standard title and help text.
pub fn generate_cell_type_diffexp_from_value(
    value: DifferentialExpressionTable,
) -> Result<WithTitle<DifferentialExpressionTable>> {
    Ok(WithTitle {
        title: TitleWithHelp {
            title: CELL_ANNOTATION_DIFFEXP_TITLE.to_string(),
            help: CELL_ANNOTATION_DIFFEXP_HELP.to_string(),
        }
        .into(),
        inner: value,
    })
}
/// Converts cell annotation metrics into a table component with a standard title and help text.
pub fn generate_cell_type_metrics(
    cell_annotation_metrics: CellAnnotationMetrics,
) -> Result<WithTitle<TableMetric>> {
    Ok(WithTitle {
        title: TitleWithHelp {
            title: CELL_ANNOTATION_METRICS_TABLE_TITLE.to_string(),
            help: CELL_ANNOTATION_METRICS_TABLE_HELP.to_string(),
        }
        .into(),
        inner: TableMetric {
            rows: cell_annotation_metrics.get_table(),
        },
    })
}
/// Converts cell annotation metrics into a raw TableMetric struct without a title wrapper.
pub fn generate_cell_type_parameter_table(
    cell_annotation_metrics: CellAnnotationMetrics,
) -> Result<TableMetric> {
    Ok(TableMetric {
        rows: cell_annotation_metrics.get_table(),
    })
}
/// Generates a standard error alert for cell annotation failures.
pub fn generate_cas_failure_alert() -> Alert {
    Alert {
        level: AlertLevel::Error,
        title: CELL_ANNOTATION_FAILURE_TITLE.to_string(),
        formatted_value: None,
        message: CELL_ANNOTATION_FAILURE_MESSAGE.to_string(),
    }
}
/// Generates a warning alert regarding differential expression issues in cell annotation.
pub fn generate_cas_de_warn_alert() -> Alert {
    Alert {
        level: AlertLevel::Error,
        title: CELL_ANNOTATION_DE_WARN_TITLE.to_string(),
        formatted_value: None,
        message: CELL_ANNOTATION_DE_WARN_MESSAGE.to_string(),
    }
}
/// Generates a warning alert for barcode mismatches during cell annotation.
///
/// Takes a custom `alert_string` message to display specific mismatch details.
pub fn generate_cas_bc_mismatch_alert(alert_string: String) -> Alert {
    Alert {
        level: AlertLevel::Warn,
        title: CELL_ANNOTATION_BC_MISMATCH_WARN_TITLE.to_string(),
        formatted_value: None,
        message: alert_string,
    }
}
