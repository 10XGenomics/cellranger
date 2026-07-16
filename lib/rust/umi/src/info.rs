#![expect(missing_docs)]
use crate::{UmiQual, UmiSeq, is_valid};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct UmiInfo {
    pub seq: UmiSeq,
    pub qual: UmiQual,
    pub is_valid: bool,
}

impl UmiInfo {
    pub fn new(seq: UmiSeq, qual: UmiQual) -> UmiInfo {
        UmiInfo {
            is_valid: is_valid(&seq, &qual),
            seq,
            qual,
        }
    }
}
