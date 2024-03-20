//! Martian stage MergeMetrics
//! Combine the summary.json files of the MatrixComputer substages.

use anyhow::Result;
use cr_types::MetricsFile;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::FileTypeRead;
use metric::{JsonReporter, Metric};
use serde::{Deserialize, Serialize};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct MergeMetricsStageInputs {
    pub summaries: Vec<Option<MetricsFile>>,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MergeMetricsStageOutputs {
    pub summary: MetricsFile,
}

/// Output the summary metrics CSV file `metrics_summary.csv`.
pub struct MergeMetrics;

/// Combine the summary.json files of the MatrixComputer substages.
#[make_mro(volatile = strict)]
impl MartianMain for MergeMetrics {
    type StageInputs = MergeMetricsStageInputs;
    type StageOutputs = MergeMetricsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let metrics = args.summaries.iter().flatten().try_fold(
            JsonReporter::default(),
            |mut metrics, json_file| -> Result<_> {
                metrics.merge(json_file.read()?);
                Ok(metrics)
            },
        )?;
        Ok(MergeMetricsStageOutputs {
            summary: MetricsFile::from_reporter(&rover, "summary", &metrics)?,
        })
    }
}
