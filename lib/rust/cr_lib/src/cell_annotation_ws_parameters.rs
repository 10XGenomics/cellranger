use anyhow::{Ok, Result};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use tenx_websummary::components::{
    DifferentialExpressionTable, InlineHelp, TableMetric, TitleWithHelp, VegaLitePlot, WithTitle,
};
use tenx_websummary::{Alert, AlertLevel};

const CELL_ANNOTATION_DISCLAIMER_STRING: &str = r#"<p><b>Important Note</b><br>
Cell type annotation is currently in beta and relies on the Chan Zuckerberg CELL by GENE reference. 
As this reference is community-driven it may not cover all tissue types, which could impact your results. For further details, please visit the
<a href='https://10xgen.com/cell-annotation' target='_blank' title='10x Support: Cell Annotation' rel='noopener noreferrer'>10x support site</a>.
</p>"#;

const CELL_ANNOTATION_BROAD_MODEL: &str =
    "Annotation was performed using a model that was co-developed by 10x Genomics and the Cellarium AI Lab at the Data Sciences Platform of the Broad Institute.";

const CELL_ANNOTATION_FAILURE_TITLE: &str = "No cell annotations produced!";

const CELL_ANNOTATION_FAILURE_MESSAGE: &str = r#"<p>Please check your cellranger annotate logs. 
If you wish to attempt cell type annotation again please use 
<a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-cell-ranger' 
target='_blank' title='What is Cell Ranger' rel='noopener noreferrer'>cellranger annotate</a>.</p>"#;

const CELL_ANNOTATION_DE_WARN_TITLE: &str = "Cell type differential expression not run";
const CELL_ANNOTATION_DE_WARN_MESSAGE: &str = "Too few cell types to run differential expression.";

const CELL_ANNOTATION_BC_MISMATCH_WARN_TITLE: &str =
    "Barcodes in cloupe file and feature barcode matrix do not match";

const CELL_ANNOTATION_METRICS_TABLE_TITLE: &str = "Cell Annotation Summary";
const CELL_ANNOTATION_METRICS_TABLE_HELP: &str = "<b>Pipeline version</b><br>
The version of Cell Ranger used to generate this annotation.<br>
<b>Cell annotation model</b><br>
The model used to generate cell type annotations.";

const CELL_ANNOTATION_INTERACTIVE_BARCHART_TITLE: &str = "Cell Type Composition";
const CELL_ANNOTATION_INTERACTIVE_BARCHART_HELP: &str = "This plot shows the major cell types in your sample. To view sub-type annotations, click on a bar. 
Only major classifications with at least 10 cells are shown.<br>
<b>Left-hand bar chart</b>: This bar chart shows the annotated cell types on the y-axis and the number of barcodes annotated with them on the x-axis.<br>
<b>Right-hand table</b>: This table shows the sub-types within a selected major cell type. Fractions of barcodes shown are the percentage of cells relative to 
this selected cell type—not the entire sample.";

const CELL_ANNOTATION_VIOLIN_PLOT_TITLE: &str = "UMI distribution by cell type";
const CELL_ANNOTATION_VIOLIN_PLOT_HELP: &str =
    "Box plot showing the distribution of (UMI+1) by annotated cell types.<br>
    Only cell types with >= 10 barcodes included. Y-axis is log scale";

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
#[derive(Serialize, Deserialize, Clone)]
pub struct CellAnnotationMetrics {
    pub cell_annotation_model: String,
    pub cell_annotation_tree_version_used: String,
    pub cell_annotation_display_map_version_used: String,
    pub cell_annotation_frac_returned_bcs: Option<f64>,
    pub cell_annotation_success: Option<bool>,
    pub cell_annotation_differential_expression: Option<bool>,
    pub cell_annotation_beta_model: bool,
    pub cell_annotation_model_developer: String,
    pub pipeline_version: String,
}

impl CellAnnotationMetrics {
    pub fn get_table(&self) -> Vec<(String, String)> {
        vec![
            (
                "Pipeline version".to_string(),
                self.pipeline_version.clone(),
            ),
            (
                "Cell annotation model".to_string(),
                self.cell_annotation_model.clone(),
            ),
        ]
    }

    pub(crate) fn generate_disclaiming_banner(&self) -> Option<InlineHelp> {
        if self.cell_annotation_beta_model {
            let disclaimer = if self.cell_annotation_model_developer.contains("Broad") {
                format!("{CELL_ANNOTATION_DISCLAIMER_STRING}\n{CELL_ANNOTATION_BROAD_MODEL}")
            } else {
                CELL_ANNOTATION_DISCLAIMER_STRING.to_string()
            };

            Some(InlineHelp::with_content(disclaimer))
        } else {
            None
        }
    }

    pub(crate) fn generate_disclaimer_html_fragment(&self) -> Option<String> {
        if self.cell_annotation_beta_model {
            let disclaimer = if self.cell_annotation_model_developer.contains("Broad") {
                format!("{CELL_ANNOTATION_DISCLAIMER_STRING}\n{CELL_ANNOTATION_BROAD_MODEL}")
            } else {
                CELL_ANNOTATION_DISCLAIMER_STRING.to_string()
            };

            Some(disclaimer)
        } else {
            None
        }
    }

    pub(crate) fn get_cloupe_track_name(&self) -> String {
        if self.cell_annotation_beta_model {
            format!("{CELL_TYPING_CLOUPE_NAME} {CELL_TYPING_BETA_SUFFIX}")
        } else {
            CELL_TYPING_CLOUPE_NAME.to_string()
        }
    }
}

pub(crate) fn generate_cell_type_barcharts_from_json(
    json_file: &JsonFile<Value>,
) -> Result<WithTitle<VegaLitePlot>> {
    generate_cell_type_barchart_from_value(json_file.read()?)
}
pub(crate) fn generate_cell_type_barchart_from_value(
    value: Value,
) -> Result<WithTitle<VegaLitePlot>> {
    Ok(WithTitle {
        title: TitleWithHelp {
            title: CELL_ANNOTATION_INTERACTIVE_BARCHART_TITLE.to_string(),
            help: CELL_ANNOTATION_INTERACTIVE_BARCHART_HELP.to_string(),
        }
        .into(),
        inner: VegaLitePlot {
            spec: value,
            actions: Some(Value::Bool(false)),
            renderer: None,
        },
    })
}

pub(crate) fn generate_cell_type_violin_plot_from_json(
    json_file: &JsonFile<Value>,
) -> Result<WithTitle<VegaLitePlot>> {
    generate_cell_type_violin_plot_from_value(json_file.read()?)
}
pub(crate) fn generate_cell_type_violin_plot_from_value(
    value: Value,
) -> Result<WithTitle<VegaLitePlot>> {
    Ok(WithTitle {
        title: TitleWithHelp {
            title: CELL_ANNOTATION_VIOLIN_PLOT_TITLE.to_string(),
            help: CELL_ANNOTATION_VIOLIN_PLOT_HELP.to_string(),
        }
        .into(),
        inner: VegaLitePlot {
            spec: value,
            actions: Some(Value::Bool(false)),
            renderer: None,
        },
    })
}

pub(crate) fn generate_cell_type_umap_plot_from_json(
    json_file: &JsonFile<Value>,
) -> Result<WithTitle<VegaLitePlot>> {
    generate_cell_type_umap_plot_from_value(json_file.read()?)
}
pub(crate) fn generate_cell_type_umap_plot_from_value(
    value: Value,
) -> Result<WithTitle<VegaLitePlot>> {
    Ok(WithTitle {
        title: TitleWithHelp {
            title: CELL_ANNOTATION_UMAP_PLOT_TITLE.to_string(),
            help: CELL_ANNOTATION_UMAP_PLOT_HELP.to_string(),
        }
        .into(),
        inner: VegaLitePlot {
            spec: value,
            actions: Some(Value::Bool(false)),
            renderer: None,
        },
    })
}

pub(crate) fn generate_cell_type_diffexp_from_json(
    json_file: &JsonFile<DifferentialExpressionTable>,
) -> Result<WithTitle<DifferentialExpressionTable>> {
    generate_cell_type_diffexp_from_value(json_file.read()?)
}
pub(crate) fn generate_cell_type_diffexp_from_value(
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

pub(crate) fn generate_cell_type_metrics(
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

pub(crate) fn generate_cell_type_parameter_table(
    cell_annotation_metrics: CellAnnotationMetrics,
) -> Result<TableMetric> {
    Ok(TableMetric {
        rows: cell_annotation_metrics.get_table(),
    })
}

pub(crate) fn generate_cas_failure_alert() -> Alert {
    Alert {
        level: AlertLevel::Error,
        title: CELL_ANNOTATION_FAILURE_TITLE.to_string(),
        formatted_value: None,
        message: CELL_ANNOTATION_FAILURE_MESSAGE.to_string(),
    }
}

pub(crate) fn generate_cas_de_warn_alert() -> Alert {
    Alert {
        level: AlertLevel::Error,
        title: CELL_ANNOTATION_DE_WARN_TITLE.to_string(),
        formatted_value: None,
        message: CELL_ANNOTATION_DE_WARN_MESSAGE.to_string(),
    }
}

pub(crate) fn generate_cas_bc_mismatch_alert(alert_string: String) -> Alert {
    Alert {
        level: AlertLevel::Warn,
        title: CELL_ANNOTATION_BC_MISMATCH_WARN_TITLE.to_string(),
        formatted_value: None,
        message: alert_string,
    }
}
