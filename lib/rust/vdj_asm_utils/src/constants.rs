#![expect(missing_docs)]

pub const MATCH_SCORE: i32 = 2;
pub const MISMATCH_SCORE: i32 = 4;
pub const GAP_OPEN: i32 = 10;
pub const GAP_EXTEND: i32 = 1;
pub const CLIP: i32 = 10;

pub type UmiType = u32;
pub type ReadType = u32;

pub const KMER_LEN_BANDED_ALIGN: usize = 11;
pub const WINDOW_SIZE_BANDED_ALIGN: usize = 10;

#[cfg(test)]
pub const LSE_COEFF: [[f64; 3]; 6] = [
    [0.00000000, 1.00000000, 0.00000000],
    [0.02572060, 0.97427940, 0.11953626],
    [0.08771177, 0.91228823, 0.29721147],
    [0.18094899, 0.81905101, 0.47282944],
    [0.29835882, 0.70164118, 0.60946732],
    [0.43104496, 0.56895504, 0.68360721],
];

// Chain type constants.

pub const CHAIN_TYPES: [&str; 8] = ["TRB", "None", "IGH", "TRA", "TRG", "TRD", "IGK", "IGL"];
pub const CHAIN_NONE_POS: i8 = 1;

pub const CHAIN_TYPESX: [&str; 7] = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];

pub const OR_CHAIN_TYPES: [&str; 14] = [
    "fw.IGH", "fw.IGK", "fw.IGL", "fw.TRA", "fw.TRB", "fw.TRD", "fw.TRG", "rc.IGH", "rc.IGK",
    "rc.IGL", "rc.TRA", "rc.TRB", "rc.TRD", "rc.TRG",
];

pub const PRIMER_EXT_LEN: usize = 40;
