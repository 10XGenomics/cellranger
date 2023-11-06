//! Martian stage COMPUTE_ANTIGEN_VDJ_METRICS

use crate::align_metrics::{BarcodeKind, BarcodeMetrics, LibFeatThenBarcodeOrder, VisitorMetrics};
use crate::BarcodeMetricsShardFile;
use anyhow::Result;
use barcode::Barcode;
use cr_types::FeatureType;
use itertools::Itertools;
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{JsonReport, MeanMetric, Metric, SimpleHistogram, TxHashSet};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};
use shardio::ShardReader;

martian_filetype! { AntigenVdjMetricsFile, "ag.vdj" }
pub type AntigenVdjMetricsFormat = BinaryFormat<AntigenVdjMetricsFile, AntigenVdjMetrics>;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ComputeAntigenVdjMetricsStageInputs {
    vdj_cell_barcodes: JsonFile<Vec<String>>,
    per_barcode_count_metrics: Vec<BarcodeMetricsShardFile>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ComputeAntigenVdjMetricsStageOutputs {
    metrics_json: JsonFile<()>,
    metrics_bin: AntigenVdjMetricsFormat,
}

pub struct ComputeAntigenVdjMetrics;

#[derive(Serialize, Deserialize, Metric, JsonReport)]
pub struct AntigenVdjMetrics {
    pub num_cells: i64,
    pub mean_usable_reads_per_cell: f64,
    pub median_umis_per_cell: f64,
}

#[derive(Serialize, Deserialize, Metric)]
struct RawMetrics {
    usable_reads_per_cell: MeanMetric,
    umis_per_cell: SimpleHistogram<i64>,
}

impl RawMetrics {
    fn median_umis_per_cell(&self) -> f64 {
        self.umis_per_cell.quantiles([0.5]).map_or(0.0, |x| x[0])
    }
}

impl RawMetrics {
    fn observe(&mut self, metrics: VisitorMetrics) {
        self.usable_reads_per_cell
            .record(metrics.usable_reads as f64);
        self.umis_per_cell.observe(&metrics.umi_counts);
    }
}

#[make_mro(mem_gb = 4, volatile = strict)]
impl MartianMain for ComputeAntigenVdjMetrics {
    type StageInputs = ComputeAntigenVdjMetricsStageInputs;
    type StageOutputs = ComputeAntigenVdjMetricsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let vdj_cells: TxHashSet<Barcode> = {
            let cells = args.vdj_cell_barcodes.read()?;
            cells.into_iter().map(|bc| bc.parse()).try_collect()
        }?;

        let metrics_reader: ShardReader<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardReader::open_set(&args.per_barcode_count_metrics)?;

        let mut raw_metrics = RawMetrics::new();
        for barcode_metrics in metrics_reader.iter()? {
            let BarcodeMetrics {
                barcode,
                library_feats,
                metrics,
            } = barcode_metrics?;

            match barcode {
                BarcodeKind::Invalid => {}
                BarcodeKind::Valid(bc) => {
                    if library_feats.has_feature(FeatureType::Antigen) && vdj_cells.contains(&bc) {
                        raw_metrics.observe(metrics);
                    }
                }
            }
        }

        let summary_json: JsonFile<()> = rover.make_path("metrics_json");
        let antigen_vdj_metrics = AntigenVdjMetrics {
            num_cells: vdj_cells.len() as i64,
            mean_usable_reads_per_cell: raw_metrics.usable_reads_per_cell.mean(),
            median_umis_per_cell: raw_metrics.median_umis_per_cell(),
        };
        let summary_bin: AntigenVdjMetricsFormat = rover.make_path("metrics_bin");
        summary_bin.write(&antigen_vdj_metrics)?;

        let mut reporter = antigen_vdj_metrics.to_json_reporter();
        reporter.add_prefix("ANTIGEN_vdj");
        reporter.report(&summary_json)?;

        Ok(ComputeAntigenVdjMetricsStageOutputs {
            metrics_json: summary_json,
            metrics_bin: summary_bin,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::rna_read::LegacyLibraryType;
    use shardio::ShardWriter;

    #[test]
    fn test_run_compute_antigen_vdj_stage() -> Result<()> {
        let tmp_dir = tempfile::tempdir()?;

        let barcodes = vec![
            "AAACCTGAGAAACCAT-1".to_string(),
            "AAACCTGAGAAACCGC-1".to_string(),
            "AAACCTGAGAAACCTA-1".to_string(),
            "AAACCTGAGAAACGAG-1".to_string(),
            "AAACCTGAGAAACGCC-1".to_string(),
            "AAACCTGAGAAAGTGG-1".to_string(),
        ];
        let vdj_cells: Vec<_> = barcodes.iter().take(3).cloned().collect();
        let vdj_cell_barcodes = JsonFile::new(&tmp_dir, "vdj_cell_barcodes");
        vdj_cell_barcodes.write(&vdj_cells)?;

        let metrics_shard = BarcodeMetricsShardFile::new(&tmp_dir, "metrics_shard");
        let mut metrics_writer: ShardWriter<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardWriter::new(&metrics_shard, 100, 100, 100)?;
        let mut sender = metrics_writer.get_sender();

        for (i, bc) in barcodes.iter().enumerate() {
            let barcode: Barcode = bc.parse()?;
            for (j, library_type) in [
                LegacyLibraryType::GeneExpression,
                LegacyLibraryType::AntigenCapture,
            ]
            .into_iter()
            .enumerate()
            {
                let barcode_kinds = if i == 0 {
                    vec![
                        Some(BarcodeKind::Valid(barcode)),
                        Some(BarcodeKind::Invalid),
                    ]
                } else {
                    vec![Some(BarcodeKind::Valid(barcode)), None]
                };
                for (k, kind) in barcode_kinds.into_iter().flatten().enumerate() {
                    sender.send(BarcodeMetrics {
                        barcode: kind,
                        library_feats: library_type.into(),
                        metrics: {
                            let mut m = VisitorMetrics::new();
                            m.usable_reads = 10 * (i + j + k + 1) as i64;
                            m.umi_counts = (i + j + k + 1) as i64;
                            m
                        },
                    })?;
                }
            }
        }
        sender.finished()?;
        metrics_writer.finish()?;

        let outs = ComputeAntigenVdjMetrics.test_run(
            &tmp_dir,
            ComputeAntigenVdjMetricsStageInputs {
                vdj_cell_barcodes,
                per_barcode_count_metrics: vec![metrics_shard],
            },
        )?;

        let actual_metrics = outs.metrics_bin.read()?;
        assert_eq!(actual_metrics.num_cells, 3);
        assert_eq!(actual_metrics.mean_usable_reads_per_cell, 30.0);
        assert_eq!(actual_metrics.median_umis_per_cell, 3.0);

        drop(tmp_dir);
        Ok(())
    }
}
