// Warning groups (as of rust 1.55)
#![deny(
    future_incompatible,
    nonstandard_style,
    rust_2018_compatibility,
    rust_2021_compatibility,
    rust_2018_idioms,
    unused
)]
// Other warnings (as of rust 1.55)
#![deny(
    asm_sub_register,
    bad_asm_style,
    bindings_with_variant_name,
    clashing_extern_declarations,
    confusable_idents,
    const_item_mutation,
    deprecated,
    deref_nullptr,
    drop_bounds,
    dyn_drop,
    elided_lifetimes_in_paths,
    exported_private_dependencies,
    function_item_references,
    improper_ctypes,
    improper_ctypes_definitions,
    incomplete_features,
    inline_no_sanitize,
    invalid_value,
    irrefutable_let_patterns,
    large_assignments,
    mixed_script_confusables,
    non_shorthand_field_patterns,
    no_mangle_generic_items,
    overlapping_range_endpoints,
    renamed_and_removed_lints,
    stable_features,
    temporary_cstring_as_ptr,
    trivial_bounds,
    type_alias_bounds,
    uncommon_codepoints,
    unconditional_recursion,
    unknown_lints,
    unnameable_test_items,
    unused_comparisons,
    while_true
)]

#[macro_use]
extern crate itertools;

use martian_derive::martian_filetype;
use serde::{Deserialize, Serialize};

martian_filetype! {HtmlFile, "html"}

mod align_and_count_metrics;
pub mod align_metrics;

/// Align homopolymer sequences.
pub mod align_homopolymer;

/// Align reads, annotate, deduplicate UMI and count genes for
/// one group of reads from a single barcode
pub mod aligner;

/// Barcode correction stage (and metrics)
mod barcode_correction_metrics;

mod barcode_overlap;

/// Barcode sorting workflow used by MAKE_SHARD
pub mod barcode_sort;

/// Struct containing cell annotation metrics for all figures
mod cell_annotation_ws_parameters;

pub mod detect_chemistry;

/// Read environment variables.
mod env;

/// Functions to do piecewise linear fitting for gDNA
mod fit_piecewise_linear_model;

/// gDNA utils
mod gdna_utils;

/// Miscellaneous macros
mod macros;

/// Metrics for MAKE_SHARD
pub mod make_shard_metrics;

/// Parquet file IO
pub mod parquet_file;

/// Preflight checks and utilities
pub mod preflight;

/// Probe barcode matrix I/O
mod probe_barcode_matrix;

/// Shared code for handling read-level multiplexing.
pub mod read_level_multiplexing;

/// Martian stages.
pub mod stages;

/// Testing utilities
pub mod testing;

/// Datatypes and martian file declarations
mod types;
pub use types::*;

/// Utility code: should be refactored to other places gradually
mod utils;

// initialize insta test harness
#[cfg(test)]
#[ctor::ctor]
fn init() {
    // when we need things like samtools or bgzip, put them on path
    dui_tests::bazel_utils::set_env_vars(&[("PATH", "lib/bin")], &[]).unwrap_or_default();
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
}
