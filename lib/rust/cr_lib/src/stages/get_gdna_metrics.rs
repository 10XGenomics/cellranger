//! Martian stage GET_GDNA_METRICS

// Other imports for stage
use crate::gdna_utils::compute_gdna_metrics;
// Bring the procedural macros in scope:
// #[derive(MartianStruct)], #[derive(MartianType)], #[make_mro], martian_filetype!
use crate::H5File;
use anyhow::Result;
use json_report_derive::JsonReport;
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeWrite;
use metric::{JsonReport, Metric};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Metric, JsonReport, Serialize, Deserialize)]
struct GdnaMetric {
    estimated_gdna_content: f64,
    estimated_gdna_unspliced_threshold: f64,
}

#[derive(Serialize, Deserialize)]
pub struct GdnaPlottingSummary {
    unspliced_counts: Vec<f64>,
    spliced_counts: Vec<f64>,
    model_crit_point: f64,
    model_constant: f64,
    model_slope: f64,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct GetGdnaMetricsStageInputs {
    pub molecule_info: H5File,
    pub reference_path: PathBuf,
    pub probe_set: CsvFile<()>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct GetGdnaMetricsStageOutputs {
    pub summary: JsonFile<()>,
    pub gdna_plot_sufficient_stats: JsonFile<GdnaPlottingSummary>,
}

// This is our stage struct
pub struct GetGdnaMetrics;

#[make_mro(mem_gb = 4)]
impl MartianMain for GetGdnaMetrics {
    type StageInputs = GetGdnaMetricsStageInputs;
    type StageOutputs = GetGdnaMetricsStageOutputs; // Use `MartianVoid` if empty
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let gdna_metrics =
            compute_gdna_metrics(&args.molecule_info, &args.probe_set, &args.reference_path);

        let metrics = GdnaMetric {
            estimated_gdna_content: gdna_metrics
                .estimated_percentage_of_gdna_umi
                .fraction()
                .unwrap_or(f64::NAN),
            estimated_gdna_unspliced_threshold: gdna_metrics.estimated_gdna_per_probe.ceil(),
        };
        let summary: JsonFile<()> = rover.make_path("summary");
        let reporter = metrics.to_json_reporter();
        reporter.report(&summary)?;

        let plotly_summary = GdnaPlottingSummary {
            unspliced_counts: gdna_metrics.unspliced_counts.clone(),
            spliced_counts: gdna_metrics.spliced_counts.clone(),
            model_constant: gdna_metrics.estimated_model.model.constant,
            model_crit_point: gdna_metrics.estimated_model.model.critical_point,
            model_slope: gdna_metrics.estimated_model.model.slope,
        };
        let gdna_plot_sufficient_stats: JsonFile<_> = rover.make_path("gdna_plot_sufficient_stats");
        gdna_plot_sufficient_stats.write(&plotly_summary).unwrap();

        Ok(GetGdnaMetricsStageOutputs {
            summary,
            gdna_plot_sufficient_stats,
        })
    }
}
