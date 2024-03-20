//! Martian stage SET_TARGETED_UMI_FILTER

use crate::types::FeatureReferenceFormat;
use crate::BcUmiInfoShardFile;
use anyhow::Result;
use cr_types::types::BcUmiInfo;
use cr_types::MetricsFile;
use json_report_derive::JsonReport;
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::FileTypeRead;
use metric::{PercentMetric, SimpleHistogram};
use serde::{Deserialize, Serialize};
use shardio::ShardReader;

const QUANTILE: f64 = 0.9;
const NUM_ORDMAG: f64 = 2.0;

/// Set the parameter `umi_read_count_threshold` based on a random sample of barcodes.
/// Receive UMI counts from a set of barcodes and construct a histogram of the
/// reads/UMI distribution. Compute a threshold proportional to some particular
/// quantile of that distribution, such that UMIs with fewer read counts than the
/// threshold are predicted to be spurious and can be filtered out without
/// sacrificing many true positive detections.
pub struct SetTargetedUmiFilter;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct StageInputs {
    pub bc_umi_info: Vec<BcUmiInfoShardFile>,
    pub feature_reference: FeatureReferenceFormat,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    pub umi_read_count_threshold: u64,
    pub summary: MetricsFile,
}

#[derive(JsonReport)]
struct UmiFilteringMetric {
    filtered_target_umi_count_threshold: i64,
    initial_filtered_target_umis: PercentMetric,
    initial_filtered_target_umis_reads: PercentMetric,
}

#[make_mro(mem_gb = 8, volatile = strict)]
impl MartianMain for SetTargetedUmiFilter {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // load target set information
        let target_set = args
            .feature_reference
            .read()?
            .target_set
            .expect("No target set found in feature_reference input.");

        // open shard reader over barcode-UMI counts
        let reader = ShardReader::<BcUmiInfo, BcUmiInfo>::open_set(&args.bc_umi_info)?;

        // Track a histogram of reads per UMI over targeted features
        let mut rpu_histogram = SimpleHistogram::default();
        for bc_umi_info in reader.iter()? {
            let bc_umi_info = bc_umi_info?;
            for c in bc_umi_info.umi_counts {
                if target_set.is_on_target(c.feature_idx) {
                    rpu_histogram.observe(&c.read_count);
                }
            }
        }

        // Compute RPU threshold using quantile + order of magnitude adjustment
        let threshold = match rpu_histogram.quantiles([QUANTILE]) {
            Some(quantiles) => {
                let threshold_float = quantiles[0] / (10_f64).powf(NUM_ORDMAG);
                threshold_float.ceil() as u64
            }
            // No counts
            None => 1_u64,
        };

        // Compute metrics
        // - The denominators here are the number of on-target reads (UMIs) counted
        //   if there were no targeted UMI filtering.
        let mut frac_umis_filtered = PercentMetric::default();
        let mut frac_reads_filtered = PercentMetric::default();
        for (&read_count, umis) in rpu_histogram.distribution() {
            let read_count = read_count as i64;
            let is_filtered = read_count < threshold as i64;
            frac_umis_filtered.increment_by(umis.count(), is_filtered);
            frac_reads_filtered.increment_by(read_count * umis.count(), is_filtered);
        }

        // Write metrics to summary
        let metrics = UmiFilteringMetric {
            filtered_target_umi_count_threshold: threshold as i64,
            initial_filtered_target_umis: frac_umis_filtered,
            initial_filtered_target_umis_reads: frac_reads_filtered,
        };

        Ok(StageOutputs {
            umi_read_count_threshold: threshold,
            summary: MetricsFile::from_reporter(&rover, "summary", &metrics)?,
        })
    }
}
