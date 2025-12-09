//! Martian stage WRITE_BARCODE_INDEX
//! Assign a distinct integer to each barcode sequence.
#![deny(missing_docs)]

use anyhow::Result;
use cr_types::{
    BarcodeIndexFormat, BarcodeIndexOutput, MetricsFile, PerLibrarySortedBarcodeCounts,
};
use itertools::Itertools;
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use metric::JsonReport;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteBarcodeIndexStageInputs {
    pub barcode_counts: PerLibrarySortedBarcodeCounts,
    pub barcode_index_override: Option<BarcodeIndexOutput>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteBarcodeIndexStageOutputs {
    pub barcode_index_output: BarcodeIndexOutput,
    pub summary: MetricsFile,
}

#[derive(Serialize, Deserialize, JsonReport)]
struct BarcodeIndexSummary {
    total_barcodes_detected: i64,
}

/// Martian stage WRITE_BARCODE_INDEX
pub struct WriteBarcodeIndex;

#[make_mro(volatile = strict, mem_gb = 1)]
impl MartianMain for WriteBarcodeIndex {
    type StageInputs = WriteBarcodeIndexStageInputs;
    type StageOutputs = WriteBarcodeIndexStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let path: PathBuf = rover.make_path("barcode_index");

        let make_metrics_summary = |num_barcodes| {
            MetricsFile::from_reporter(
                &rover,
                "summary.json",
                &BarcodeIndexSummary {
                    total_barcodes_detected: num_barcodes as i64,
                }
                .to_json_reporter(),
            )
        };

        if let Some(index_override) = args.barcode_index_override {
            return Ok(WriteBarcodeIndexStageOutputs {
                summary: make_metrics_summary(index_override.num_barcodes)?,
                barcode_index_output: index_override.hardlink_or_copy(&path)?,
            });
        }
        let barcode_index_output = args
            .barcode_counts
            .iter_barcodes()?
            .process_results(|bc_iter| BarcodeIndexFormat::write(&path, bc_iter))??;

        Ok(WriteBarcodeIndexStageOutputs {
            summary: make_metrics_summary(barcode_index_output.num_barcodes)?,
            barcode_index_output,
        })
    }
}
