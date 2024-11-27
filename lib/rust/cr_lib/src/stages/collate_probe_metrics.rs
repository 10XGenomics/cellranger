//! Martian stage COLLATE_PROBE_METRICS

// Bring the procedural macros in scope:
// Other imports for stage
use crate::gdna_utils::{
    compute_gdna_corrected_median_genes_per_spot, get_filtered_per_probe_metrics,
    GdnaCorrectedMetrics,
};
use crate::probe_barcode_matrix::{write_probe_bc_matrix, ProbeCounts};
use anyhow::Result;
use cr_types::filtered_barcodes::{read_filtered_barcodes_set, FilteredBarcodesCsv};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::{BarcodeIndexFormat, CountShardFile, H5File};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::TxHashSet;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Deserialize, MartianStruct)]
pub struct CollateProbeMetricsStageInputs {
    pub probe_barcode_counts: Vec<CountShardFile>,
    pub reference_path: Option<PathBuf>,
    pub probe_set: TargetSetFile,
    pub filtered_barcodes: FilteredBarcodesCsv,
    pub probe_set_name: String,
    pub barcode_index_path: BarcodeIndexFormat,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct CollateProbeMetricsStageOutputs {
    pub per_probe_metrics: CsvFile<ProbeCounts>,
    pub raw_probe_bc_matrix: H5File,
    pub estimated_gdna_metrics: JsonFile<GdnaCorrectedMetrics>,
}

// This is our stage struct
pub struct CollateProbeMetrics;

#[make_mro(mem_gb = 8, volatile = strict)]
impl MartianMain for CollateProbeMetrics {
    type StageInputs = CollateProbeMetricsStageInputs;
    type StageOutputs = CollateProbeMetricsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Construct a set of filtered barcode to check membership in
        let filtered_barcodes = read_filtered_barcodes_set(&args.filtered_barcodes)?;

        // Decide if gDNA analysis should be run. If so run it and get if a probe is filtered. Everything
        // written out in the form of a vec<ProbeCounts>
        let filtered_per_probe_metrics: Vec<ProbeCounts> = get_filtered_per_probe_metrics(
            &args.probe_barcode_counts,
            &filtered_barcodes,
            &args.probe_set,
            args.reference_path.as_deref(),
        )?;

        // Compute gDNA corrected metrics and write them out out to JSON.
        let estimated_gdna_metrics: JsonFile<_> = rover.make_path("estimated_gdna_metrics");
        estimated_gdna_metrics.write(&compute_gdna_corrected_median_genes_per_spot(
            &args.probe_barcode_counts,
            &filtered_per_probe_metrics,
            &filtered_barcodes,
        )?)?;

        // Write out per probe metrics
        let per_probe_metrics: CsvFile<_> = rover.make_path("filtered_probes");
        let mut filtered_probes_csv_writer = per_probe_metrics.lazy_writer()?;
        for row in &filtered_per_probe_metrics {
            filtered_probes_csv_writer.write_item(row)?;
        }

        // Get filtered probeset
        let filtered_probe_set: TxHashSet<_> = filtered_per_probe_metrics
            .into_iter()
            .filter(|x| x.pass_filter)
            .map(|x| x.probe_id)
            .collect();

        // Read barcode index
        let barcode_index = args.barcode_index_path.read()?;

        // Write out raw-probe-filtered-barcode matrix. Filtered probes are annotated.
        let raw_probe_bc_matrix: H5File = rover.make_path("raw_probe_bc_matrix");
        write_probe_bc_matrix(
            &args.probe_set,
            args.reference_path.as_deref(),
            &args.probe_barcode_counts,
            &raw_probe_bc_matrix,
            &barcode_index,
            &filtered_probe_set,
            &barcode_index.into_indicator_vec(&filtered_barcodes),
            &args.probe_set_name,
        )?;

        Ok(Self::StageOutputs {
            per_probe_metrics,
            raw_probe_bc_matrix,
            estimated_gdna_metrics,
        })
    }
}
