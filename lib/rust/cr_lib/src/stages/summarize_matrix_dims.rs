//! Martian stage SUMMARIZE_MATRIX_DIMS
//! Summarize the dimensions of an H5 matrix.
#![deny(missing_docs)]

use anyhow::Result;
use cr_h5::count_matrix::CountMatrixFile;
use cr_types::MetricsFile;
use cr_types::probe_set::ProbeType;
use cr_types::reference::feature_reference::FeatureType;
use cr_types::reference::probe_set_reference::TargetSetFile;
use itertools::chain;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use metric::{JsonReporter, join_metric_name};
use serde::{Deserialize, Serialize};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct SummarizeMatrixDimsStageInputs {
    /// Metric name prefix
    pub prefix: String,
    /// Feature-barcode matrix
    pub matrix: CountMatrixFile,
    /// Probe set reference CSV for counting gapfill probes
    pub probe_set: Option<TargetSetFile>,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct SummarizeMatrixDimsStageOutputs {
    /// JSON summary of matrix dimensions
    pub summary: MetricsFile,
}

/// Martian stage SUMMARIZE_MATRIX_DIMS
pub struct SummarizeMatrixDims;

#[make_mro(volatile = strict)]
impl MartianMain for SummarizeMatrixDims {
    type StageInputs = SummarizeMatrixDimsStageInputs;
    type StageOutputs = SummarizeMatrixDimsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let dims = args.matrix.load_dimensions()?;
        let feature_type_counts = args.matrix.load_feature_type_counts()?;

        let prefix = &args.prefix;
        let gapfill_probes_count = if let Some(probe_set) = args.probe_set {
            assert_eq!(prefix, "raw_probe_matrix");
            let count = probe_set
                .read(None)?
                .into_iter()
                .filter(|x| {
                    matches!(
                        x.probe_type,
                        ProbeType::PairedGapAlign | ProbeType::UnpairedGapAlign
                    )
                })
                .count();
            Some((format!("{prefix}_gapfill_probes_count"), count))
        } else {
            None
        };
        let metrics: JsonReporter = chain!(
            [
                (format!("{prefix}_features_count"), dims.num_features),
                (format!("{prefix}_barcodes_count"), dims.num_barcodes),
                (
                    format!("{prefix}_nonzero_entries_count"),
                    dims.num_non_zeros,
                ),
            ],
            feature_type_counts
                .into_iter()
                .filter(|(feature_type, _)| *feature_type != FeatureType::Peaks)
                .map(|(feature_type, count)| {
                    (
                        join_metric_name(feature_type, &format!("{prefix}_feature_type_count")),
                        count,
                    )
                }),
            gapfill_probes_count,
        )
        .collect();

        Ok(SummarizeMatrixDimsStageOutputs {
            summary: MetricsFile::from_reporter(&rover, "summary", &metrics)?,
        })
    }
}
