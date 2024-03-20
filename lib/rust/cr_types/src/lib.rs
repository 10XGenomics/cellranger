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

use anyhow::{ensure, Context, Result};
use csv_parser::CsvParser;
use std::path::{Path, PathBuf};

pub mod aggr;
pub mod barcode_index;
pub use barcode_index::*;
mod bit_encode;
pub mod chemistry;
pub mod clonotype;
pub mod constants;
pub mod csv_parser;
pub mod probe_set;
pub mod reference;
pub mod rna_read;
pub mod sample_def;
mod serde_helpers;
pub mod target_panel_summary;
pub mod types;
pub use types::*;
mod metrics_file;
pub mod spill_vec;
pub mod utils;
pub use metrics_file::*;
pub mod filtered_barcodes;
pub mod websummary;

#[derive(Debug)]
pub struct LegacyLibrariesEntry {
    pub fastqs: PathBuf,
    pub library_type: LibraryType,
    pub sample: String,
    pub project: Option<String>,
}

pub fn parse_legacy_libraries_csv(file: &Path) -> Result<Vec<LegacyLibrariesEntry>> {
    let required_fields = ["fastqs", "library_type", "sample"];

    let mut parser = CsvParser::new(file, required_fields, "Libraries CSV")
        .with_context(|| format!("Error opening --libraries file: {}", file.display()))?;

    let mut libraries: Vec<LegacyLibrariesEntry> = Vec::new();

    for i in 0..parser.len() {
        parser.set_line(i);

        let fastqs = parser.require_string("fastqs")?;
        let sample = parser.require_string("sample")?;
        let library_type: LibraryType = parser.parse_field("library_type", "valid library_type")?;
        // CELLRANGER-7889: this previously only parsed "auto VDJ", so maintain this invariant.
        if let Some(chain_type) = library_type.vdj_chain_type() {
            ensure!(
                chain_type == VdjChainType::Auto,
                "invalid library_type {library_type}"
            );
        }

        let project = parser.get_extra_data().get("project").cloned();

        let def = LegacyLibrariesEntry {
            fastqs: fastqs.into(),
            library_type,
            sample,
            project,
        };

        libraries.push(def);
    }

    Ok(libraries)
}

// initialize insta test harness
#[cfg(test)]
#[ctor::ctor]
fn init() {
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
}
