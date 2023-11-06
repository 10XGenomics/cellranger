use crate::stages::align_and_count::StageInputs;
use anyhow::Result;
use barcode::Barcode;
use metric::TxHashSet;

pub fn get_barcode_subsampling(_args: &StageInputs) -> Result<(TxHashSet<Barcode>, Option<f64>)> {
    Ok((TxHashSet::default(), None))
}
