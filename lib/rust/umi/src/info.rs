use crate::{UmiQual, UmiSeq};
use fastq_set::metric_utils::ILLUMINA_QUAL_OFFSET;
use serde::{Deserialize, Serialize};

const UMI_POLYT_SUFFIX_LENGTH: usize = 5;
const UMI_MIN_QV: u8 = 10;

#[derive(Serialize, Deserialize, Debug)]
pub struct UmiInfo {
    pub seq: UmiSeq,
    pub qual: UmiQual,
    pub is_valid: bool,
    pub low_min_qual: bool,
    pub has_n: bool,
    pub is_homopolymer: bool,
    pub has_polyt: bool,
}

impl UmiInfo {
    pub fn new(seq: UmiSeq, qual: UmiQual) -> UmiInfo {
        let has_n = has_n(seq.as_bytes());
        let has_polyt = has_polyt(seq.as_bytes());
        let is_homopolymer = is_homopolymer(seq.as_bytes());
        let low_min_qual = low_min_qual(qual.as_bytes(), UMI_MIN_QV);
        let is_valid = !(has_n || is_homopolymer || low_min_qual);

        UmiInfo {
            seq,
            qual,
            is_valid,
            low_min_qual,
            has_n,
            is_homopolymer,
            has_polyt,
        }
    }
}

fn has_n(seq: &[u8]) -> bool {
    seq.contains(&b'N')
}

fn has_polyt(seq: &[u8]) -> bool {
    let end_idx = if UMI_POLYT_SUFFIX_LENGTH <= seq.len() {
        seq.len() - UMI_POLYT_SUFFIX_LENGTH
    } else {
        seq.len()
    };
    for base in &seq[..end_idx] {
        if *base != b'T' {
            return false;
        }
    }
    true
}

fn is_homopolymer(seq: &[u8]) -> bool {
    // TODO this would call NNNNN as a homopolymer
    for i in 1..seq.len() {
        if seq[i - 1] != seq[i] {
            return false;
        }
    }
    true
}

fn low_min_qual(qual: &[u8], min_qv: u8) -> bool {
    for qv in qual {
        if qv - ILLUMINA_QUAL_OFFSET < min_qv {
            return true;
        }
    }
    false
}
