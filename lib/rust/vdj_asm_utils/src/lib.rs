//! vdj_asm_utils
#![deny(missing_docs)]
pub mod asm;
pub mod bam_utils;
pub mod barcode_data;
pub mod constants;
mod contigs;
pub mod exact_clonotyping;
mod fastq;
pub mod graph_read;
pub mod heuristics;
pub mod log_opts;
pub mod primers;
pub mod process;
mod ref_free;
pub mod sw;
pub mod umi_data;
mod utils;

pub(crate) fn reverse_complement(x: &mut [u8]) {
    x.reverse();
    for entry in x {
        *entry = match entry {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => *entry,
        }
    }
}
