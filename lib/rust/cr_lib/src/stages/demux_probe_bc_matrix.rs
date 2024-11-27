//! Martian stage DEMUX_PROBE_BC_MATRIX

use crate::gdna_utils::get_filtered_per_probe_metrics;
use crate::probe_barcode_matrix::{read_bc_json, write_probe_bc_matrix, ProbeCounts};
use barcode::Barcode;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::{BarcodeIndex, CountShardFile, H5File};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::TxHashSet;
use serde::de::IgnoredAny;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::iter::zip;
use std::path::PathBuf;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixStageInputs {
    pub probe_barcode_counts: Vec<CountShardFile>,
    pub reference_path: Option<PathBuf>,
    pub probe_set: TargetSetFile,
    pub probe_set_name: String,
    pub sample_barcodes: JsonFile<HashMap<String, Vec<String>>>,
    pub sample_cell_barcodes: JsonFile<HashMap<String, Vec<String>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixStageOutputs {
    pub sample_raw_probe_bc_matrices: HashMap<String, H5File>,
    pub samples_per_probe_metrics: HashMap<String, CsvFile<ProbeCounts>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixChunkInputs {
    pub sample_name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixChunkOutputs {
    pub sample_per_probe_metrics: CsvFile<ProbeCounts>,
    pub sample_raw_probe_bc_matrix: H5File,
}

// This is our stage struct
pub struct DemuxProbeBcMatrix;

#[make_mro(volatile = strict)]
impl MartianStage for DemuxProbeBcMatrix {
    type StageInputs = DemuxProbeBcMatrixStageInputs;
    type StageOutputs = DemuxProbeBcMatrixStageOutputs;
    type ChunkInputs = DemuxProbeBcMatrixChunkInputs;
    type ChunkOutputs = DemuxProbeBcMatrixChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>, Error> {
        // Reading the sample barcodes
        let sample_cell_barcodes: HashMap<String, IgnoredAny> =
            JsonFile::from_path(&args.sample_cell_barcodes).read()?;
        // Split into one sample per chunk.
        Ok(sample_cell_barcodes
            .into_keys()
            .map(|sample_name| {
                (
                    DemuxProbeBcMatrixChunkInputs { sample_name },
                    Resource::with_mem_gb(5),
                )
            })
            .collect())
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs, Error> {
        // re-reading the sample barcodes. Because can not use barcodes in split
        let sample_bcs = read_bc_json(&args.sample_barcodes)?
            .remove(&chunk_args.sample_name)
            .unwrap();
        let sample_cell_bcs: TxHashSet<Barcode> = read_bc_json(&args.sample_cell_barcodes)?
            .remove(&chunk_args.sample_name)
            .unwrap()
            .into_iter()
            .collect();

        // Decide if gDNA analysis should be run. If so run it and get if a probe is filtered. Everything
        // written out in the form of a vec<ProbeCounts>
        let filtered_per_probe_metrics = get_filtered_per_probe_metrics(
            &args.probe_barcode_counts,
            &sample_cell_bcs,
            &args.probe_set,
            args.reference_path.as_deref(),
        )?;

        // Write out sample per probe metrics
        let sample_per_probe_metrics: CsvFile<_> = rover.make_path("sample_per_probe_metrics");
        sample_per_probe_metrics.write(&filtered_per_probe_metrics)?;

        // Get filtered probeset
        let filtered_probe_set: TxHashSet<_> = filtered_per_probe_metrics
            .into_iter()
            .filter(|x| x.pass_filter)
            .map(|x| x.probe_id)
            .collect();

        // Write out raw-probe-raw-barcode matrix. Filtered probes are annotated. Filtered barcodes are not
        let sample_raw_probe_bc_matrix: H5File = rover.make_path("sample_raw_probe_bc_matrix");
        let bc_index = BarcodeIndex::from_iter(sample_bcs);
        let cell_bc_indicator_vec = bc_index.into_indicator_vec(&sample_cell_bcs);
        write_probe_bc_matrix(
            &args.probe_set,
            args.reference_path.as_deref(),
            &args.probe_barcode_counts,
            &sample_raw_probe_bc_matrix,
            &bc_index,
            &filtered_probe_set,
            &cell_bc_indicator_vec,
            &args.probe_set_name,
        )?;

        Ok(DemuxProbeBcMatrixChunkOutputs {
            sample_raw_probe_bc_matrix,
            sample_per_probe_metrics,
        })
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs, Error> {
        let sample_names: Vec<_> = chunk_defs.into_iter().map(|x| x.sample_name).collect();
        let samples_per_probe_metrics: HashMap<_, _> = zip(&sample_names, &chunk_outs)
            .map(|(sample_name, out)| (sample_name.clone(), out.sample_per_probe_metrics.clone()))
            .collect();
        let sample_raw_probe_bc_matrices: HashMap<_, _> = zip(sample_names, chunk_outs)
            .map(|(sample_name, out)| (sample_name, out.sample_raw_probe_bc_matrix))
            .collect();
        Ok(DemuxProbeBcMatrixStageOutputs {
            sample_raw_probe_bc_matrices,
            samples_per_probe_metrics,
        })
    }
}
