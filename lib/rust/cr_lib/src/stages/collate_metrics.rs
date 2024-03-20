//! Martian stage COLLATE_METRICS
//! This stage receives per barcode metrics from ALIGN_AND_COUNT and creates the summary JSON.

use crate::align_metrics::{BarcodeKind, BarcodeMetrics, LibFeatThenBarcodeOrder, VisitorMetrics};
use crate::types::{BarcodeMetricsShardFile, FeatureReferenceFormat};
use crate::AggregateBarcode;
use anyhow::Result;
use barcode::Barcode;
use cr_types::filtered_barcodes::{read_filtered_barcodes_set, FilteredBarcodesCsv};
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::{
    GenomeName, LibraryType, MetricsFile, SampleAssignment, SampleBarcodes, SampleBarcodesFile,
};
use itertools::{process_results, Itertools};
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeRead;
use metric::{JsonReport, JsonReporter, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter::zip;
use std::ops::Add;
use std::path::PathBuf;

const FILTERED_BARCODES_SUFFIX: &str = "_in_filtered_barcodes";

/// Aggregate barcode metrics for one sample.
#[allow(non_snake_case)]
#[derive(JsonReport, derive_more::Add)]
struct AggregateBarcodesMetrics {
    /// Number of GEMs that contain protein aggregates.
    number_aggregate_GEMs: i64,

    /// Fraction antibody reads in aggregate barcodes.
    reads_lost_to_aggregate_GEMs: f64,
}

impl From<&AggregateBarcode> for AggregateBarcodesMetrics {
    /// Convert one AggregateBarcode to an AggregateBarcodesMetrics
    fn from(aggregate_barcode: &AggregateBarcode) -> Self {
        Self {
            number_aggregate_GEMs: 1,
            reads_lost_to_aggregate_GEMs: aggregate_barcode.frac_sample_reads,
        }
    }
}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct CollateMetricsStageInputs {
    pub per_barcode_metrics: Vec<BarcodeMetricsShardFile>,
    pub reference_path: PathBuf,
    pub feature_reference: FeatureReferenceFormat,
    pub filtered_barcodes: Option<FilteredBarcodesCsv>,
    pub aggregate_barcodes: Option<CsvFile<AggregateBarcode>>,
    pub sample_barcodes: Option<SampleBarcodesFile>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct SampleMetrics {
    pub sample: SampleAssignment,
    pub summary: MetricsFile,
    pub per_barcode_metrics: Option<CsvFile<()>>,
}

#[derive(Serialize, MartianStruct)]
pub struct CollateMetricsStageOutputs {
    pub summary: Option<MetricsFile>,
    pub per_barcode_metrics: Option<CsvFile<()>>,
    pub multi_metrics: Option<Vec<SampleMetrics>>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CollateMetricsChunkInputs {
    pub sample: SampleAssignment,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct CollateMetricsChunkOutputs {
    pub summary: MetricsFile,
    pub per_barcode_metrics: Option<CsvFile<()>>,
}

/// VisitorMetrics for both all barcodes and filtered barcodes.
#[derive(Default)]
struct FilteredBarcodesMetrics {
    /// Metrics for all barcodes.
    all_barcodes_metrics: VisitorMetrics,

    /// Metrics for filtered barcodes. None when filtered_barcodes.csv is not available.
    filtered_barcodes_metrics: Option<VisitorMetrics>,
}

impl std::ops::Add for FilteredBarcodesMetrics {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        FilteredBarcodesMetrics {
            all_barcodes_metrics: self.all_barcodes_metrics + rhs.all_barcodes_metrics,
            filtered_barcodes_metrics: match (
                self.filtered_barcodes_metrics,
                rhs.filtered_barcodes_metrics,
            ) {
                (None, None) => None,
                (Some(lhs), None) => Some(lhs),
                (None, Some(rhs)) => Some(rhs),
                (Some(lhs), Some(rhs)) => Some(lhs + rhs),
            },
        }
    }
}

impl std::iter::Sum for FilteredBarcodesMetrics {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), Self::add)
    }
}

impl FilteredBarcodesMetrics {
    /// Construct FilteredBarcodesMetrics from BarcodeMetrics.
    /// Set filtered_barcodes_metrics to Some when the barcode is found in filtered_barcodes.
    fn new(
        barcode_metrics: BarcodeMetrics,
        filtered_barcodes: Option<&TxHashSet<Barcode>>,
    ) -> Self {
        let in_filtered_barcodes = if let (Some(filtered_barcodes), BarcodeKind::Valid(barcode)) =
            (filtered_barcodes, &barcode_metrics.barcode)
        {
            filtered_barcodes.contains(barcode)
        } else {
            false
        };

        FilteredBarcodesMetrics {
            filtered_barcodes_metrics: in_filtered_barcodes
                .then(|| barcode_metrics.metrics.clone()),
            all_barcodes_metrics: barcode_metrics.metrics,
        }
    }

    /// Construct a JsonReport of metrics in all barcodes and in filtered barcodes.
    fn make_report(
        self,
        library_type: LibraryType,
        genomes: &[GenomeName],
        has_targeted_gex: bool,
    ) -> JsonReporter {
        let mut reporter =
            self.all_barcodes_metrics
                .make_report(library_type, genomes, has_targeted_gex);

        // Construct metrics for filtered barcodes.
        // Exclude *_barcoded_* barcodes, because all reads in filtered barcodes are barcoded.
        // Append to the metric key a suffix FILTERED_BARCODES_SUFFIX.
        if let Some(metrics) = self.filtered_barcodes_metrics {
            reporter.extend(
                metrics
                    .make_report(library_type, genomes, has_targeted_gex)
                    .into_iter()
                    .filter_map(|(key, value)| {
                        if key.contains("_barcoded_") {
                            None
                        } else {
                            Some((key + FILTERED_BARCODES_SUFFIX, value))
                        }
                    }),
            );
        }

        reporter
    }
}

/// Calculate metrics per library.
fn calculate_metrics_per_library(
    iter: impl Iterator<Item = BarcodeMetrics>,
    barcodes_for_sample: Option<&TxHashSet<Barcode>>,
    filtered_barcodes: Option<&TxHashSet<Barcode>>,
    genomes: &[GenomeName],
    csv_writer: &mut BufWriter<File>,
) -> Result<TxHashMap<LibraryType, FilteredBarcodesMetrics>> {
    iter.group_by(|x| x.library_type)
        .into_iter()
        .map(|(library_type, bc_metrics_iter)| {
            if library_type == LibraryType::Gex {
                BarcodeMetrics::to_csv_header(genomes, csv_writer)?;
            }

            let library_metrics = bc_metrics_iter
                .filter(|bcm| {
                    match (barcodes_for_sample, bcm.barcode) {
                        // Skip invalid barcodes of multiplexed samples.
                        (Some(_), BarcodeKind::Invalid) => false,
                        // Skip barcodes not in this sample.
                        (Some(barcodes), BarcodeKind::Valid(bc)) => barcodes.contains(&bc),
                        // Not multiplexed, keep all barcodes.
                        (None, _) => true,
                    }
                })
                .group_by(|x| x.barcode)
                .into_iter()
                .map(|(barcode, mut group)| {
                    let metrics = match barcode {
                        BarcodeKind::Invalid => group.map(|x| x.metrics).sum(),
                        BarcodeKind::Valid(barcode) => {
                            let metrics = group.next().unwrap().metrics;
                            assert!(group.next().is_none(), "Duplicate barcode {barcode}");
                            metrics
                        }
                    };
                    let bcm = BarcodeMetrics {
                        library_type,
                        barcode,
                        metrics,
                    };
                    if library_type == LibraryType::Gex {
                        bcm.to_csv_row(genomes, csv_writer)?;
                    }
                    Ok(FilteredBarcodesMetrics::new(bcm, filtered_barcodes))
                })
                .sum::<Result<_>>()?;

            Ok((library_type, library_metrics))
        })
        .collect()
}

/// Calculate aggregate barcodes metrics.
fn calculate_aggregate_barcodes_metrics(
    aggregate_metrics: &[AggregateBarcode],
    barcodes_for_sample: &TxHashSet<Barcode>,
) -> HashMap<LibraryType, AggregateBarcodesMetrics> {
    aggregate_metrics
        .iter()
        .filter(|x| barcodes_for_sample.contains(&x.barcode))
        .map(|x| (x.library_type, AggregateBarcodesMetrics::from(x)))
        .into_grouping_map()
        .sum()
}

/// Output the metrics per barcode CSV file `per_barcode_metrics.csv`.
pub struct CollateMetrics;

#[make_mro(volatile = strict)]
impl MartianStage for CollateMetrics {
    type StageInputs = CollateMetricsStageInputs;
    type StageOutputs = CollateMetricsStageOutputs;
    type ChunkInputs = CollateMetricsChunkInputs;
    type ChunkOutputs = CollateMetricsChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let mem_gb = if args.sample_barcodes.is_some() { 5 } else { 1 };
        Ok(SampleBarcodes::read_samples(args.sample_barcodes.as_ref())?
            .into_iter()
            .map(|sample| {
                (
                    CollateMetricsChunkInputs { sample },
                    Resource::with_mem_gb(mem_gb),
                )
            })
            .collect())
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let filtered_barcodes = if let Some(filtered_barcodes) = &args.filtered_barcodes {
            Some(read_filtered_barcodes_set(filtered_barcodes)?)
        } else {
            None
        };

        let barcodes_for_sample = SampleBarcodes::read_from_json(args.sample_barcodes.as_ref())?
            .into_barcodes(&chunk_args.sample)
            .map(TxHashSet::from_iter);

        let ref_info = ReferenceInfo::from_reference_path(&args.reference_path)?;

        let has_targeted_gex = args.feature_reference.read()?.target_set.is_some();

        let per_barcode_metrics: CsvFile<_> = rover.make_path("per_barcode_metrics");
        let mut csv_writer = per_barcode_metrics.buf_writer()?;

        // We will accumulate metrics per library type and report the metrics with
        // the appropriate prefix. The shard files are sorted by (library_type, barcode).
        let metrics_reader: ShardReader<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardReader::open_set(&args.per_barcode_metrics)?;
        let per_lib_type_metrics = process_results(metrics_reader.iter()?, |iter| {
            calculate_metrics_per_library(
                iter,
                barcodes_for_sample.as_ref(),
                filtered_barcodes.as_ref(),
                &ref_info.genomes,
                &mut csv_writer,
            )
        })??;
        csv_writer.flush()?;

        let per_barcode_metrics = if per_lib_type_metrics.contains_key(&LibraryType::Gex) {
            assert_ne!(per_barcode_metrics.metadata()?.len(), 0);
            Some(per_barcode_metrics)
        } else {
            assert_eq!(per_barcode_metrics.metadata()?.len(), 0);
            None
        };

        let per_lib_type_reports: TxHashMap<_, _> = per_lib_type_metrics
            .into_iter()
            .map(|(library_type, metrics)| {
                (
                    library_type,
                    metrics.make_report(library_type, &ref_info.genomes, has_targeted_gex),
                )
            })
            .collect();

        let aggregate_barcode_metrics =
            if let (Some(aggeregate_barcodes), Some(barcodes_for_sample)) =
                &(args.aggregate_barcodes, barcodes_for_sample)
            {
                calculate_aggregate_barcodes_metrics(
                    &aggeregate_barcodes.read()?,
                    barcodes_for_sample,
                )
                .to_json_reporter()
            } else {
                JsonReporter::default()
            };

        let metrics: TxHashMap<_, _> = chain!(
            per_lib_type_reports.to_json_reporter(),
            ref_info.into_json_report(),
            aggregate_barcode_metrics,
        )
        .collect();

        let summary: MetricsFile = rover.make_path("summary");
        summary.write(&metrics)?;

        Ok(CollateMetricsChunkOutputs {
            summary,
            per_barcode_metrics,
        })
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        if chunk_defs.len() == 1 && chunk_defs[0].sample == SampleAssignment::NonMultiplexed {
            let Self::ChunkOutputs {
                summary,
                per_barcode_metrics,
            } = chunk_outs.into_iter().next().unwrap();
            return Ok(CollateMetricsStageOutputs {
                summary: Some(summary),
                per_barcode_metrics,
                multi_metrics: None,
            });
        }

        let multi_metrics = zip(chunk_defs, chunk_outs)
            .map(|(chunk_def, chunk_out)| SampleMetrics {
                sample: chunk_def.sample,
                summary: chunk_out.summary,
                per_barcode_metrics: chunk_out.per_barcode_metrics,
            })
            .collect();

        Ok(CollateMetricsStageOutputs {
            summary: None,
            per_barcode_metrics: None,
            multi_metrics: Some(multi_metrics),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Seek};
    use std::iter;

    fn make_metrics_per_library(
        library_type: LibraryType,
        csv_file: File,
    ) -> Result<TxHashMap<LibraryType, FilteredBarcodesMetrics>> {
        let mut metrics = BarcodeMetrics {
            barcode: BarcodeKind::Invalid,
            library_type,
            metrics: VisitorMetrics::default(),
        };
        metrics.metrics.total_reads.increment();

        calculate_metrics_per_library(
            iter::repeat(metrics).take(3),
            None,
            None,
            &["GRCh38".into()],
            &mut BufWriter::new(csv_file),
        )
    }

    #[test]
    fn test_calculate_metrics_per_library_gex() -> Result<()> {
        let mut csv_file = tempfile::tempfile()?;
        let metrics_per_library =
            make_metrics_per_library(LibraryType::Gex, csv_file.try_clone()?)?;
        assert_eq!(metrics_per_library.len(), 1);

        let gex_metrics = &metrics_per_library[&LibraryType::Gex];
        assert_eq!(gex_metrics.all_barcodes_metrics.total_reads.count(), 3);
        assert!(gex_metrics.filtered_barcodes_metrics.is_none());

        let mut csv_string = String::new();
        csv_file.rewind()?;
        csv_file.read_to_string(&mut csv_string)?;
        assert_eq!(csv_string.matches('\n').count(), 2);
        assert_eq!(csv_string.matches("NO_BARCODE,3").count(), 1);

        Ok(())
    }

    #[test]
    fn test_calculate_metrics_per_library_antibody() -> Result<()> {
        let antibody = LibraryType::Antibody;
        let csv_file = tempfile::tempfile()?;
        let metrics_per_library = make_metrics_per_library(antibody, csv_file.try_clone()?)?;
        assert_eq!(metrics_per_library.len(), 1);

        let antibody_metrics = &metrics_per_library[&antibody];
        assert_eq!(antibody_metrics.all_barcodes_metrics.total_reads.count(), 3);
        assert!(antibody_metrics.filtered_barcodes_metrics.is_none());
        assert_eq!(csv_file.metadata()?.len(), 0);

        Ok(())
    }
}
