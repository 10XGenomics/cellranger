#![deny(missing_docs)]
use anyhow::Result;
use barcode::Barcode;
use martian_derive::MartianStruct;
use martian_filetypes::json_file::JsonFile;
use metric::TxHashSet;
use serde::Deserialize;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct V1PatternFixParams {
    #[expect(dead_code)]
    affected_barcodes: JsonFile<Vec<String>>,
    #[expect(dead_code)]
    correction_factor: f64,
}

impl V1PatternFixParams {
    pub(super) fn barcode_subsampling(&self) -> Result<(TxHashSet<Barcode>, Option<f64>)> {
        Ok((TxHashSet::default(), None))
    }
}
