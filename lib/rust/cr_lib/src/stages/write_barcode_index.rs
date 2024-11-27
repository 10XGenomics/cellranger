//! Martian stage WRITE_BARCODE_INDEX
//! Assign a distinct integer to each barcode sequence.

use crate::utils::estimate_mem::{barcode_mem_gib, get_total_barcodes_detected};
use anyhow::Result;
use cr_types::barcode_index::BarcodeIndex;
use cr_types::types::{BarcodeIndexFormat, BcCountFormat};
use cr_types::MetricsFile;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteBarcodeIndexStageInputs {
    pub barcode_counts: BcCountFormat,
    pub barcodes_under_tissue: Option<JsonFile<Vec<String>>>,
    pub barcode_correction_summary: MetricsFile,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteBarcodeIndexStageOutputs {
    pub barcode_index: BarcodeIndexFormat,
}

pub struct WriteBarcodeIndex;

#[make_mro(volatile = strict)]
impl MartianStage for WriteBarcodeIndex {
    type StageInputs = WriteBarcodeIndexStageInputs;
    type StageOutputs = WriteBarcodeIndexStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // Multi uses more memory for sample_barcodes and other data.
        let barcodes_count = get_total_barcodes_detected(&args.barcode_correction_summary.read()?);
        // bytes_per_barcode and offset_gib are empirically determined.
        let mem_gib = barcode_mem_gib(barcodes_count, 220, 2);
        println!("barcode_count={barcodes_count},mem_gib={mem_gib}");
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(mem_gib)))
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
        {
            let barcode_index: BarcodeIndexFormat = rover.make_path("barcode_index");
            barcode_index.write(&BarcodeIndex::new(
                &args.barcode_counts,
                args.barcodes_under_tissue.as_ref(),
            )?)?;
            Ok(WriteBarcodeIndexStageOutputs { barcode_index })
        }
    }
}
