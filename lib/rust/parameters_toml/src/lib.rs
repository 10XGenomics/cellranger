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

use anyhow::{Context, Result};
use log::warn;
use serde::Deserialize;
use std::sync::OnceLock;

#[derive(Debug, Deserialize, Clone)]
struct Parameters {
    /// DETECT_CHEMISTRY samples detect_chemistry_sample_reads from detect_chemistry_total_reads.
    detect_chemistry_sample_reads: usize,
    /// DETECT_CHEMISTRY samples detect_chemistry_sample_reads from detect_chemistry_total_reads.
    detect_chemistry_total_reads: usize,
    /// DETECT_CHEMISTRY ensures that the reads agree with the barcode whitelist.
    min_fraction_whitelist_match: f64,
    /// CHECK_BARCODES_COMPATIBILITY ensures that two libraries have sufficient barcode overlap.
    min_barcode_similarity: f64,
    /// Maximum number of reads to use per barcode. In paired
    /// end mode, one read pair counts as two reads.
    vdj_max_reads_per_barcode: usize,
    /// Command line parameters of STAR.
    star_parameters: String,
    /// Alternative minimum UMI read length
    umi_min_read_length: Option<usize>,
    /// Maximum number of multiplexing capture tags allowed
    max_multiplexing_tags: usize,
    /// 5' chemistries with multiplexing
    fiveprime_multiplexing: bool,
    /// 3' LT chemistries with multiplexing
    threeprime_lt_multiplexing: bool,
    /// Minimum fraction of the single major probe barcode for
    /// singleplex FRP libraries
    min_major_probe_bc_frac: f64,
}

const DEFAULT_PARAMETERS: Parameters = Parameters {
    detect_chemistry_sample_reads: 100_000,
    detect_chemistry_total_reads: 2_000_000,
    min_fraction_whitelist_match: 0.1,
    min_barcode_similarity: 0.1,
    vdj_max_reads_per_barcode: 80_000,
    star_parameters: String::new(),
    umi_min_read_length: None,
    max_multiplexing_tags: 12,
    fiveprime_multiplexing: true,
    threeprime_lt_multiplexing: false,
    min_major_probe_bc_frac: 0.7,
};
static PARAMETERS: OnceLock<Result<Parameters>> = OnceLock::new();

/// Return a reference to the global parameters.
/// The parameters may need to be loaded; if loading fails, return Err.
fn parameters() -> &'static Result<Parameters> {
    // TODO: use get_or_try_init once [#109737](https://github.com/rust-lang/rust/issues/109737) is stabilized
    PARAMETERS.get_or_init(|| {
        let path = bazel_utils::current_exe()
            .context("Unable to locate the running executable")?
            .with_file_name("parameters.toml");
        if !path.exists() {
            warn!(
                "could not find parameters.toml at {}, falling back to defaults",
                path.display()
            );
            Ok(DEFAULT_PARAMETERS)
        } else {
            let s = std::fs::read_to_string(&path).with_context(|| path.display().to_string())?;
            Ok(toml::from_str(&s).with_context(|| path.display().to_string())?)
        }
    })
}

macro_rules! parameter_getter {
    ($a:ident, $t:ty) => {
        pub fn $a() -> Result<&'static $t> {
            let val = match parameters() {
                Err(e) => return Err(anyhow::anyhow!(e)),
                Ok(p) => &p.$a,
            };
            if DEFAULT_PARAMETERS.$a != *val {
                warn!("using non-default {} = {:?}", stringify!($a), val);
            }
            Ok(val)
        }
    };
}

parameter_getter!(detect_chemistry_sample_reads, usize);
parameter_getter!(detect_chemistry_total_reads, usize);
parameter_getter!(min_fraction_whitelist_match, f64);
parameter_getter!(min_barcode_similarity, f64);
parameter_getter!(vdj_max_reads_per_barcode, usize);
parameter_getter!(umi_min_read_length, Option<usize>);
parameter_getter!(max_multiplexing_tags, usize);
parameter_getter!(fiveprime_multiplexing, bool);
parameter_getter!(threeprime_lt_multiplexing, bool);
parameter_getter!(min_major_probe_bc_frac, f64);
parameter_getter!(star_parameters, str);
