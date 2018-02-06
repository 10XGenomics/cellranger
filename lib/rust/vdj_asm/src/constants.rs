//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

pub const QUAL_OFFSET: u8 = 33;

pub const MATCH_SCORE: i32 = 2;
pub const MISMATCH_SCORE: i32 = 4;
pub const GAP_OPEN: i32 = 10;
pub const GAP_EXTEND: i32 = 1;
pub const CLIP: i32 = 10;
pub const SEED_LEN: usize = 20;
pub const MAX_NUM_READPAIRS: usize = 1e5 as usize;
pub const MAX_NUM_KMERS: usize = 150000 as usize;

pub const DUMMY_CONTIG_NAME: &'static str = "";

pub const RAW_UMI_TAG: &'static str = "UR";
pub const PROCESSED_UMI_TAG: &'static str = "UB";
pub const PROCESSED_BC_TAG: &'static str = "CB";

pub type UmiType = u32;
pub type ReadType = u32;

pub const EPSILON : f64 = 1.0e-10;

pub static NUCS : &'static [u8] = b"ACGT";

pub const KMER_LEN_BANDED_ALIGN: usize = 11;
pub const WINDOW_SIZE_BANDED_ALIGN: usize = 10;

pub const LSE_COEFF: [[f64; 3]; 6] = [[0.00000000, 1.00000000, 0.00000000],[0.02572060, 0.97427940, 0.11953626],[0.08771177, 0.91228823, 0.29721147],[0.18094899, 0.81905101, 0.47282944],[0.29835882, 0.70164118, 0.60946732],[0.43104496, 0.56895504, 0.68360721],];
