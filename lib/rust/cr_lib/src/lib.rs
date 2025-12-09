//! cr_lib
#![deny(missing_docs)]

mod align_and_count_metrics;
mod align_metrics;
pub(crate) mod io;

/// Align homopolymer sequences.
pub mod align_homopolymer;

/// Align reads, annotate, deduplicate UMI and count genes for
/// one group of reads from a single barcode
mod aligner;

/// Barcode correction stage (and metrics)
mod barcode_correction_metrics;

mod barcode_overlap;

/// Barcode sorting workflow used by MAKE_SHARD
mod barcode_sort;

/// Struct containing cell annotation metrics for all figures
mod cell_annotation_ws_parameters;

mod detect_chemistry;

/// Read environment variables.
mod env;

/// Functions to do piecewise linear fitting for gDNA
mod fit_piecewise_linear_model;

/// gDNA utils
mod gdna_utils;

/// Miscellaneous macros
mod macros;

/// Metrics for MAKE_SHARD
mod make_shard_metrics;

/// Minimap2 aligner
mod minimap2;

/// Parquet file IO
pub mod parquet_file;

/// Preflight checks and utilities
mod preflight;

/// Probe barcode matrix I/O
mod probe_barcode_matrix;

/// Shared code for handling read-level multiplexing.
mod read_level_multiplexing;

/// Martian stages.
pub mod stages;

/// Testing utilities
pub mod testing;

/// Datatypes and martian file declarations
mod types;

/// Utility code: should be refactored to other places gradually
mod utils;

pub use align_metrics::MULTI_GENOME;
pub use aligner::BarcodeSummary;
pub use barcode_sort::BarcodeOrder;
pub use read_level_multiplexing::map_multiplexing_seq_to_id;
pub use stages::{correct_barcode_in_read, load_dist, make_bam_header, write_web_summary_html};
pub use types::*;

// initialize insta test harness
#[cfg(test)]
#[ctor::ctor]
fn init() {
    // when we need things like samtools or bgzip, put them on path
    dui_tests::bazel_utils::set_env_vars(&[("PATH", "lib/bin")], &[]).unwrap_or_default();
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    unsafe { std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root) }
}
