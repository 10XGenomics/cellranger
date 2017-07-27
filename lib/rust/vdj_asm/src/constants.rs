//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

pub const P: usize = 8;
pub const K: usize = 47;
pub const READ_LEN: usize = 150;

pub const QUAL_OFFSET: u8 = 33;

pub const MATCH_SCORE: f32 = 2.0;
pub const MISMATCH_SCORE: f32 = 2.0;
pub const GAP_OPEN: f32 = 5.0;
pub const GAP_EXTEND: f32 = 2.0;
pub const CLIP: f32 = 1.0;
pub const SEED_LEN: usize = 20;
pub const MAX_NUM_READPAIRS: usize = 2e5 as usize;
pub const MAX_NUM_KMERS: usize = 1e6 as usize;

pub const DUMMY_CONTIG_NAME: &'static str = "";

pub const RAW_UMI_TAG: &'static str = "UR";
pub const PROCESSED_UMI_TAG: &'static str = "UB";
pub const PROCESSED_BC_TAG: &'static str = "CB";

pub type UmiType = u32;
pub type ReadType = u32;

pub const EPSILON : f64 = 1.0e-10;

pub static NUCS : &'static [u8] = b"ACGT";
