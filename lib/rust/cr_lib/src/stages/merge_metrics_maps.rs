//! Martian stage MERGE_METRICS_MAPS
//! Merge metrics JSON files.
#![deny(missing_docs)]

use anyhow::Result;
use cr_types::MetricsFile;
use itertools::Itertools;
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use metric::JsonReporter;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

/// The Martian stage inputs.
#[derive(Deserialize, MartianStruct)]
pub struct MergeMetricsMapsStageInputs {
    top_level_name: String,
    summaries: HashMap<String, Option<MetricsFile>>,
}

/// The Martian stage outputs.
#[derive(Serialize, Deserialize, MartianStruct)]
pub struct MergeMetricsMapsStageOutputs {
    summary: Option<MetricsFile>,
}

/// Martian stage MERGE_METRICS_MAPS
pub struct MergeMetricsMaps;

/// Combine the summary.json files of map calls.
#[make_mro(volatile = strict)]
impl MartianMain for MergeMetricsMaps {
    type StageInputs = MergeMetricsMapsStageInputs;
    type StageOutputs = MergeMetricsMapsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let metrics: serde_json::Map<String, Value> = args
            .summaries
            .into_iter()
            .filter_map(|(key, filename)| Some((key, filename?)))
            .map(|(key, filename)| anyhow::Ok((key, filename.read()?)))
            .try_collect()?;
        let summary = if metrics.is_empty() {
            None
        } else {
            Some(MetricsFile::from_reporter(
                &rover,
                "summary",
                &JsonReporter::from((args.top_level_name, metrics)),
            )?)
        };
        Ok(Self::StageOutputs { summary })
    }
}
