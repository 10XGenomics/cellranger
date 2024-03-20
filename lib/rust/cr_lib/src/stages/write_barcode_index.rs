//! Martian stage WRITE_BARCODE_INDEX
//! Assign a distinct integer to each barcode sequence.

use anyhow::Result;
use cr_types::barcode_index::BarcodeIndex;
use cr_types::types::{BarcodeIndexFormat, BcCountFormat};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeWrite;
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteBarcodeIndexStageInputs {
    pub barcode_counts: BcCountFormat,
    pub barcodes_under_tissue: Option<JsonFile<Vec<String>>>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteBarcodeIndexStageOutputs {
    pub barcode_index: BarcodeIndexFormat,
}

pub struct WriteBarcodeIndex;

#[make_mro(mem_gb = 6, volatile = strict)]
impl MartianMain for WriteBarcodeIndex {
    type StageInputs = WriteBarcodeIndexStageInputs;
    type StageOutputs = WriteBarcodeIndexStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let barcode_index: BarcodeIndexFormat = rover.make_path("barcode_index");
        barcode_index.write(&BarcodeIndex::new(
            &args.barcode_counts,
            args.barcodes_under_tissue.as_ref(),
        )?)?;
        Ok(WriteBarcodeIndexStageOutputs { barcode_index })
    }
}
