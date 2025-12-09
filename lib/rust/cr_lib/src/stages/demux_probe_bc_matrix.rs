//! Martian stage DEMUX_PROBE_BC_MATRIX
#![deny(missing_docs)]

use crate::gdna_utils::get_filtered_per_probe_metrics;
use crate::probe_barcode_matrix::{ProbeCounts, write_probe_bc_matrix};
use barcode::Barcode;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::{
    BarcodeIndex, CountShardFile, H5File, SampleAssignment, SampleBarcodes, SampleBarcodesFile,
    SampleId,
};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeWrite;
use martian_filetypes::tabular_file::CsvFile;
use metric::TxHashSet;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::iter::zip;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixStageInputs {
    pub probe_barcode_counts: Vec<CountShardFile>,
    pub reference_info: Option<ReferenceInfo>,
    pub probe_set: TargetSetFile,
    pub probe_set_name: String,
    pub sample_barcodes: SampleBarcodesFile,
    pub sample_cell_barcodes: SampleBarcodesFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixStageOutputs {
    pub sample_raw_probe_bc_matrices: HashMap<SampleId, H5File>,
    pub samples_per_probe_metrics: HashMap<SampleId, CsvFile<ProbeCounts>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixChunkInputs {
    pub sample_name: SampleId,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct DemuxProbeBcMatrixChunkOutputs {
    pub sample_per_probe_metrics: CsvFile<ProbeCounts>,
    pub sample_raw_probe_bc_matrix: H5File,
}

/// Martian stage DEMUX_PROBE_BC_MATRIX
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
        let sample_barcodes_lines = crate::io::count_lines(&args.sample_barcodes)?;
        let mem_gib = 1 + (120 * sample_barcodes_lines).div_ceil(1024 * 1024 * 1024);
        println!("sample_barcodes_lines={sample_barcodes_lines},mem_gib={mem_gib}");

        // Split into one sample per chunk.
        Ok(
            SampleBarcodes::read_samples(Some(&args.sample_cell_barcodes))?
                .into_iter()
                .map(|sample_name| {
                    let SampleAssignment::Assigned(sample_name) = sample_name else {
                        unreachable!("unexpected {sample_name}");
                    };
                    (
                        DemuxProbeBcMatrixChunkInputs { sample_name },
                        Resource::with_mem_gb(mem_gib as isize),
                    )
                })
                .collect(),
        )
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs, Error> {
        let sample_cell_bcs: TxHashSet<Barcode> =
            SampleBarcodes::read_from_json(Some(&args.sample_cell_barcodes))?
                .into_barcodes(&SampleAssignment::from(chunk_args.sample_name.clone()))
                .unwrap()
                .into_iter()
                .collect();

        let reference_path = args
            .reference_info
            .as_ref()
            .and_then(ReferenceInfo::get_reference_path);

        // Decide if gDNA analysis should be run. If so run it and get if a probe is filtered. Everything
        // written out in the form of a vec<ProbeCounts>
        let filtered_per_probe_metrics = get_filtered_per_probe_metrics(
            &args.probe_barcode_counts,
            &sample_cell_bcs,
            &args.probe_set,
            reference_path,
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
        // re-reading the sample barcodes. Because can not use barcodes in split
        let bc_index = {
            BarcodeIndex::from_sorted(
                SampleBarcodes::read_from_json(Some(&args.sample_barcodes))?
                    .into_barcodes(&SampleAssignment::from(chunk_args.sample_name))
                    .unwrap()
                    .into_iter()
                    .sorted(),
            )
        };
        let cell_bc_indicator_vec = bc_index.into_indicator_vec(&sample_cell_bcs);
        write_probe_bc_matrix(
            &args.probe_set,
            reference_path,
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
