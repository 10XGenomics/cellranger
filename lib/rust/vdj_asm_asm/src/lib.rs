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
pub mod airrfilter;
pub mod asm_call_cells;
mod assembly;
pub use assembly::{Assembly, BarcodeDataBriefFile};
mod adapter;
pub mod assembly_types;
mod contig_aligner;
pub mod handle_gex_cells;
pub mod make_filter_switch;
mod translator;
pub mod write_ann_csv;
pub mod write_contig_outs;

// initialize insta test harness
#[cfg(test)]
#[ctor::ctor]
fn init() {
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
}
