//! Martian stage COMPUTE_ANTIGEN_VDJ_METRICS

use crate::align_metrics::{BarcodeKind, BarcodeMetrics, LibFeatThenBarcodeOrder, VisitorMetrics};
use crate::BarcodeMetricsShardFile;
use anyhow::Result;
use barcode::Barcode;
use cr_types::FeatureBarcodeType;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{MeanMetric, SimpleHistogram, TxHashSet};
use serde::{Deserialize, Serialize};
use shardio::ShardReader;

martian_filetype! { AntigenVdjMetricsFile, "ag.vdj" }
pub type AntigenVdjMetricsFormat = BinaryFormat<AntigenVdjMetricsFile, AntigenVdjMetrics>;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct ComputeAntigenVdjMetricsStageInputs {
    vdj_cell_barcodes: JsonFile<Vec<String>>,
    per_barcode_count_metrics: Vec<BarcodeMetricsShardFile>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct ComputeAntigenVdjMetricsStageOutputs {
    metrics_json: JsonFile<AntigenVdjMetrics>,
    metrics_bin: AntigenVdjMetricsFormat,
}

#[derive(Serialize, Deserialize)]
pub struct AntigenVdjMetrics {
    #[serde(rename = "ANTIGEN_vdj_num_cells")]
    pub num_cells: i64,
    #[serde(rename = "ANTIGEN_vdj_mean_usable_reads_per_cell")]
    pub mean_usable_reads_per_cell: f64,
    #[serde(rename = "ANTIGEN_vdj_median_umis_per_cell")]
    pub median_umis_per_cell: f64,
}

#[derive(Default)]
struct RawMetrics {
    usable_reads_per_cell: MeanMetric,
    umis_per_cell: SimpleHistogram<i64>,
}

impl RawMetrics {
    fn median_umis_per_cell(&self) -> f64 {
        self.umis_per_cell.quantiles([0.5]).map_or(0.0, |x| x[0])
    }

    fn observe(&mut self, metrics: VisitorMetrics) {
        self.usable_reads_per_cell
            .record(metrics.usable_reads as f64);
        self.umis_per_cell.observe(&metrics.umi_counts);
    }
}

pub struct ComputeAntigenVdjMetrics;

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

        let mut raw_metrics = RawMetrics::default();
        for barcode_metrics in metrics_reader.iter()? {
            let BarcodeMetrics {
                barcode,
                library_type,
                metrics,
            } = barcode_metrics?;

            match barcode {
                BarcodeKind::Invalid => {}
                BarcodeKind::Valid(bc) => {
                    if library_type.is_fb_type(FeatureBarcodeType::Antigen)
                        && vdj_cells.contains(&bc)
                    {
                        raw_metrics.observe(metrics);
                    }
                }
            }
        }

        let antigen_vdj_metrics = AntigenVdjMetrics {
            num_cells: vdj_cells.len() as i64,
            mean_usable_reads_per_cell: raw_metrics.usable_reads_per_cell.mean(),
            median_umis_per_cell: raw_metrics.median_umis_per_cell(),
        };
        let metrics_bin: AntigenVdjMetricsFormat = rover.make_path("metrics_bin");
        metrics_bin.write(&antigen_vdj_metrics)?;

        let metrics_json: JsonFile<_> = rover.make_path("metrics_json");
        metrics_json.write(&antigen_vdj_metrics)?;

        Ok(ComputeAntigenVdjMetricsStageOutputs {
            metrics_json,
            metrics_bin,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::LibraryType;
    use shardio::ShardWriter;

    #[test]
    fn test_run_compute_antigen_vdj_stage() -> Result<()> {
        let tmp_dir = tempfile::tempdir()?;

        let barcodes = [
            "AAACCTGAGAAACCAT-1",
            "AAACCTGAGAAACCGC-1",
            "AAACCTGAGAAACCTA-1",
            "AAACCTGAGAAACGAG-1",
            "AAACCTGAGAAACGCC-1",
            "AAACCTGAGAAAGTGG-1",
        ];
        let vdj_cells: Vec<_> = barcodes
            .iter()
            .take(3)
            .copied()
            .map(str::to_string)
            .collect();
        let vdj_cell_barcodes = JsonFile::new(&tmp_dir, "vdj_cell_barcodes");
        vdj_cell_barcodes.write(&vdj_cells)?;

        let metrics_shard = BarcodeMetricsShardFile::new(&tmp_dir, "metrics_shard");
        let mut metrics_writer: ShardWriter<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardWriter::new(&metrics_shard, 100, 100, 100)?;
        let mut sender = metrics_writer.get_sender();

        for (i, bc) in barcodes.into_iter().enumerate() {
            let barcode: Barcode = bc.parse()?;
            for (j, library_type) in [LibraryType::Gex, LibraryType::Antigen]
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
                        library_type,
                        metrics: {
                            let mut m = VisitorMetrics::default();
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
