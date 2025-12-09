//! Martian stage COLLATE_PROBE_METRICS
#![deny(missing_docs)]

// Bring the procedural macros in scope:
// Other imports for stage
use crate::gdna_utils::{
    GdnaCorrectedMetrics, compute_gdna_corrected_median_genes_per_spot,
    get_filtered_per_probe_metrics,
};
use crate::probe_barcode_matrix::{ProbeCounts, write_probe_bc_matrix};
use crate::utils::estimate_mem::barcode_mem_gib;
use anyhow::Result;
use cr_h5::count_matrix::CountMatrixFile;
use cr_types::filtered_barcodes::{FilteredBarcodesCsv, read_filtered_barcodes_set};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::{BarcodeIndexFormat, CountShardFile, H5File};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::TxHashSet;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, MartianStruct)]
pub struct CollateProbeMetricsStageInputs {
    pub probe_barcode_counts: Vec<CountShardFile>,
    pub reference_info: Option<ReferenceInfo>,
    pub probe_set: TargetSetFile,
    pub filtered_barcodes: FilteredBarcodesCsv,
    pub probe_set_name: String,
    pub barcode_index_path: BarcodeIndexFormat,
    pub raw_feature_bc_matrix: CountMatrixFile,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct CollateProbeMetricsStageOutputs {
    pub per_probe_metrics: CsvFile<ProbeCounts>,
    pub raw_probe_bc_matrix: H5File,
    pub estimated_gdna_metrics: JsonFile<GdnaCorrectedMetrics>,
}

/// Martian stage COLLATE_PROBE_METRICS
pub struct CollateProbeMetrics;

#[make_mro(volatile = strict)]
impl MartianStage for CollateProbeMetrics {
    type StageInputs = CollateProbeMetricsStageInputs;
    type StageOutputs = CollateProbeMetricsStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;
    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let total_barcodes_detected = args.raw_feature_bc_matrix.load_dimensions()?.num_barcodes;
        // barcode_mem_gib params are empirically determined.
        let barcode_index_mem_gib = barcode_mem_gib(total_barcodes_detected, 99, 6);
        println!("total_barcodes_detected={total_barcodes_detected}");
        println!("barcode_index_mem_gib={barcode_index_mem_gib}");

        Ok(StageDef::with_join_resource(Resource::with_mem_gb(
            barcode_index_mem_gib,
        )))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: MartianVoid,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<MartianVoid>,
        _chunk_outs: Vec<MartianVoid>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        // Construct a set of filtered barcode to check membership in
        let filtered_barcodes = read_filtered_barcodes_set(&args.filtered_barcodes)?;

        let reference_path = args
            .reference_info
            .as_ref()
            .and_then(|x| x.get_reference_path());
        // Decide if gDNA analysis should be run. If so run it and get if a probe is filtered. Everything
        // written out in the form of a vec<ProbeCounts>
        let filtered_per_probe_metrics: Vec<ProbeCounts> = get_filtered_per_probe_metrics(
            &args.probe_barcode_counts,
            &filtered_barcodes,
            &args.probe_set,
            reference_path,
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
            reference_path,
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
