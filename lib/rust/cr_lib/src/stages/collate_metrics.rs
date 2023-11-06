//! Martian stage COLLATE_METRICS
//! This stage receives per barcode metrics from ALIGN_AND_COUNT and creates the summary JSON.

use crate::align_metrics::{BarcodeKind, BarcodeMetrics, LibFeatThenBarcodeOrder, VisitorMetrics};
use crate::types::{BarcodeMetricsShardFile, FeatureReferenceFormat};
use crate::AggregateBarcode;
use anyhow::Result;
use barcode::Barcode;
use cr_types::reference::genome_of_chrom::GenomeName;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::types::{LibraryFeatures, SampleAssignment, SampleBarcodes};
use cr_types::{MetricsFile, SampleBarcodesFile};
use itertools::Itertools;
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::tabular_file::{CsvFile, CsvFileNoHeader};
use martian_filetypes::FileTypeRead;
use metric::{JsonReport, JsonReporter, Metric, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::collections::HashMap;
use std::io::Write;
use std::iter::zip;
use std::path::PathBuf;
use std::str::FromStr;

const FILTERED_BARCODES_SUFFIX: &str = "_in_filtered_barcodes";

/// Aggregate barcode metrics for one sample.
#[allow(non_snake_case)]
#[derive(Clone, JsonReport, Serialize, derive_more::Add)]
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
    pub filtered_barcodes: Option<CsvFileNoHeader<(String, String)>>,
    pub aggregate_barcodes: Option<CsvFile<AggregateBarcode>>,
    pub sample_barcodes: Option<SampleBarcodesFile>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct SampleMetrics {
    pub sample: SampleAssignment,
    pub summary: MetricsFile,
    pub per_barcode_metrics: CsvFile<()>,
}

#[derive(Clone, Serialize, MartianStruct)]
pub struct CollateMetricsStageOutputs {
    pub summary: Option<MetricsFile>,
    pub per_barcode_metrics: Option<CsvFile<()>>,
    pub multi_metrics: Option<Vec<SampleMetrics>>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CollateMetricsChunkInputs {
    pub sample: SampleAssignment,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CollateMetricsChunkOutputs {
    pub summary: MetricsFile,
    pub per_barcode_metrics: CsvFile<()>,
}

/// VisitorMetrics for both all barcodes and filtered barcodes.
struct FilteredBarcodesMetrics {
    /// Metrics for all barcodes.
    all_barcodes_metrics: VisitorMetrics,

    /// Metrics for filtered barcodes. None when filtered_barcodes.csv is not available.
    filtered_barcodes_metrics: Option<VisitorMetrics>,
}

impl FilteredBarcodesMetrics {
    /// Construct a new FilteredBarcodesMetrics.
    /// Calculte metrics in filtered barcodes if have_filtered_barcodes is true.
    fn new(have_filtered_barcodes: bool) -> FilteredBarcodesMetrics {
        FilteredBarcodesMetrics {
            all_barcodes_metrics: VisitorMetrics::new(),
            filtered_barcodes_metrics: if have_filtered_barcodes {
                Some(VisitorMetrics::new())
            } else {
                None
            },
        }
    }

    /// Construct a JsonReport of metrics in all barcodes and in filtered barcodes.
    fn make_report(
        self,
        library_feats: LibraryFeatures,
        genomes: &[GenomeName],
        has_targeted_gex: bool,
    ) -> JsonReporter {
        let mut reporter =
            self.all_barcodes_metrics
                .make_report(library_feats, genomes, has_targeted_gex);

        // Construct metrics for filtered barcodes.
        // Exclude *_barcoded_* barcodes, because all reads in filtered barcodes are barcoded.
        // Append to the metric key a suffix FILTERED_BARCODES_SUFFIX.
        if let Some(metrics) = self.filtered_barcodes_metrics {
            reporter.extend(
                metrics
                    .make_report(library_feats, genomes, has_targeted_gex)
                    .into_iter()
                    .filter_map(|(key, value)| {
                        if key.contains("_barcoded_") {
                            None
                        } else {
                            Some((key + FILTERED_BARCODES_SUFFIX, value))
                        }
                    }),
            )
        }

        reporter
    }
}

/// Read the filtered_barcodes CSV file.
/// Its two columns are genome and barcode.
fn read_filtered_barcodes_set(
    filtered_barcodes_filename: &CsvFileNoHeader<(String, String)>,
) -> Result<TxHashSet<Barcode>> {
    filtered_barcodes_filename
        .read()?
        .into_iter()
        .map(|(_genome, barcode): (String, String)| Barcode::from_str(&barcode))
        .collect()
}

/// Calculate aggregate barcodes metrics.
fn calculate_aggregate_barcodes_metrics(
    aggregate_metrics: &[AggregateBarcode],
    barcodes_for_sample: &TxHashSet<Barcode>,
) -> HashMap<&'static str, AggregateBarcodesMetrics> {
    aggregate_metrics
        .iter()
        .filter(|x| barcodes_for_sample.contains(&x.barcode))
        .map(|x| {
            (
                x.library_type.metric_prefix().unwrap(),
                AggregateBarcodesMetrics::from(x),
            )
        })
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

        let metrics_reader: ShardReader<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardReader::open_set(&args.per_barcode_metrics)?;

        let ref_info = ReferenceInfo::from_reference_path(&args.reference_path)?;

        let has_targeted_gex = args.feature_reference.read()?.target_set.is_some();

        let per_barcode_metrics: CsvFile<_> = rover.make_path("per_barcode_metrics");
        let mut csv_writer = per_barcode_metrics.buf_writer()?;
        // We will accumulate metrics per library type and report the metrics with
        // the appropriate prefix. The shard files are sorted by (library_feats, barcode).
        let mut per_lib_type_metrics = TxHashMap::default();
        for (lib_feats, bc_metrics_iter) in &metrics_reader
            .iter()?
            .group_by(|m| m.as_ref().ok().map(|x| x.library_feats))
        {
            let mut lib_type_metrics = FilteredBarcodesMetrics::new(filtered_barcodes.is_some());
            if lib_feats == Some(LibraryFeatures::gex()) {
                BarcodeMetrics::to_csv_header(&ref_info.genomes, &mut csv_writer);
            }
            for bc_metric in bc_metrics_iter {
                let bcm = bc_metric?;

                // skip this barcode if the sample doesn't contain it
                match (&barcodes_for_sample, bcm.barcode) {
                    (&Some(_), BarcodeKind::Invalid) => {
                        // sample is multiplexed so invalid barcode is no good
                        continue;
                    }
                    (Some(bcs), BarcodeKind::Valid(bc)) => {
                        if !bcs.contains(&bc) {
                            continue;
                        }
                    }
                    _ => {}
                }
                if lib_feats == Some(LibraryFeatures::gex()) {
                    bcm.to_csv_row(&ref_info.genomes, &mut csv_writer);
                }

                if let Some(metrics) = lib_type_metrics.filtered_barcodes_metrics.as_mut() {
                    let in_filtered_barcodes = match &bcm.barcode {
                        BarcodeKind::Invalid => false,
                        BarcodeKind::Valid(bc) => filtered_barcodes.as_ref().unwrap().contains(bc),
                    };
                    if in_filtered_barcodes {
                        *metrics += bcm.metrics.clone();
                    }
                }

                lib_type_metrics.all_barcodes_metrics += bcm.metrics;
            }
            per_lib_type_metrics.insert(lib_feats.unwrap(), lib_type_metrics);
        }
        csv_writer.flush()?;

        let mut reporter = per_lib_type_metrics
            .into_iter()
            .map(|(library_features, metrics)| {
                (
                    library_features
                        .legacy_library_type()
                        .metric_prefix()
                        .unwrap_or("_"), // "_" means no prefix, applies to GEX library
                    metrics.make_report(library_features, &ref_info.genomes, has_targeted_gex),
                )
            })
            .collect::<TxHashMap<_, _>>()
            .to_json_reporter();

        if let (Some(aggeregate_barcodes), Some(barcodes_for_sample)) =
            &(args.aggregate_barcodes, barcodes_for_sample)
        {
            reporter.merge(
                calculate_aggregate_barcodes_metrics(
                    &aggeregate_barcodes.read()?,
                    barcodes_for_sample,
                )
                .to_json_reporter(),
            );
        }

        reporter.merge(ref_info.to_json_reporter());

        Ok(CollateMetricsChunkOutputs {
            summary: MetricsFile::from_reporter(&rover, "summary", &reporter)?,
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
                per_barcode_metrics: Some(per_barcode_metrics),
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
