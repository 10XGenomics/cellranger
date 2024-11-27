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
extern crate io_utils;

#[macro_use]
extern crate lazy_static;

// Modules introducing macros need to come before the modules that use them.

mod analyze_vs;
pub mod asm;
pub mod bam_utils;
pub mod barcode_data;
pub mod constants;
mod contigs;
pub mod exact_clonotyping;
mod fastq;
pub mod graph_read;
pub mod heuristics;
pub mod hops;
pub mod log_opts;
mod metrics;
pub mod primers;
pub mod process;
mod ref_free;
pub mod sw;
mod utils;
