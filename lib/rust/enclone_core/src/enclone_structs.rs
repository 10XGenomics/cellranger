// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::barcode_fate::BarcodeFate;
use crate::defs::{CloneInfo, ExactClonotype};
use enclone_proto::types::DonorReferenceItem;
use qd::Double;
use std::collections::HashMap;

pub type BarcodeFates = HashMap<String, BarcodeFate>;

// FIXME: this being i32 is legacy and we should replace this with either u32,
// usize, or a newtype index.
pub type CloneInfoIndex = i32;

pub type CandidateClonotype = Vec<CloneInfoIndex>;

#[derive(Default, Clone)]
pub struct EncloneExacts {
    pub to_bc: HashMap<(usize, usize), Vec<String>>,
    pub exact_clonotypes: Vec<ExactClonotype>,
    pub raw_joins: Vec<Vec<usize>>,
    pub info: Vec<CloneInfo>,
    pub candidate_clonotypes: Vec<CandidateClonotype>,
    pub drefs: Vec<DonorReferenceItem>,
    pub sr: Vec<Vec<Double>>,
}
