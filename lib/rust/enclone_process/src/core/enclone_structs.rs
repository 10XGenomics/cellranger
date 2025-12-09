// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use super::barcode_fate::BarcodeFate;
use super::defs::{CloneInfo, ExactClonotype};
use enclone_proto::types::DonorReferenceItem;
use qd::Double;
use std::collections::HashMap;

pub(crate) type BarcodeFates = HashMap<String, BarcodeFate>;

// FIXME: this being i32 is legacy and we should replace this with either u32,
// usize, or a newtype index.
pub(crate) type CloneInfoIndex = i32;

pub(crate) type CandidateClonotype = Vec<CloneInfoIndex>;

#[derive(Default, Clone)]
pub(crate) struct EncloneExacts {
    pub(crate) to_bc: HashMap<(usize, usize), Vec<String>>,
    pub(crate) exact_clonotypes: Vec<ExactClonotype>,
    pub(crate) raw_joins: Vec<Vec<usize>>,
    pub(crate) info: Vec<CloneInfo>,
    pub(crate) candidate_clonotypes: Vec<CandidateClonotype>,
    pub(crate) drefs: Vec<DonorReferenceItem>,
    pub(crate) sr: Vec<Vec<Double>>,
}
