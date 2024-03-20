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

use anyhow::Result;
use docopt::Docopt;
use martian::prelude::*;
use serde::Deserialize;

const HEADER: &str = "# Copyright 2023 10x Genomics, Inc. All rights reserved.";

const USAGE: &str = "
Test Rust stage
Usage:
  cr_vdj martian <adapter>...
  cr_vdj mro [--file=<filename>] [--rewrite]
  cr_vdj --help
Options:
     --help            Show this screen.
";

#[derive(Debug, Deserialize)]
struct Args {
    // Martian interface
    cmd_martian: bool,
    cmd_mro: bool,
    arg_adapter: Vec<String>,
    flag_file: Option<String>,
    flag_rewrite: bool,
}

fn main() -> Result<()> {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    let (stage_registry, mro_registry) = martian_stages![
        cr_vdj::copy_vdj_reference::CopyVdjReference,
        cr_vdj::summarize_vdj_filters::SummarizeVdjFilters,
        cr_vdj::create_barcode_csv::CreateBarcodeCsv,
        clonotype_assigner::Assigner,
        clonotype_assigner::write_clonotype_outs::WriteClonotypeOuts,
        clonotype_assigner::fill_clonotype_info::FillClonotypeInfo,
        clonotype_assigner::handle_no_clonotyping::HandleNoClonotyping,
        clonotype_assigner::write_concat_ref_outs::WriteConcatRefOuts,
        clonotype_assigner::write_consensus_bam::WriteConsensusBam,
        clonotype_assigner::write_consensus_txt::WriteConsensusTxt,
        vdj_asm_asm::Assembly,
        vdj_asm_asm::asm_call_cells::AsmCallCells,
        vdj_asm_asm::airrfilter::AirrFilter,
        vdj_asm_asm::write_contig_outs::WriteContigOuts,
        vdj_asm_asm::handle_gex_cells::HandleGexCells,
        vdj_asm_asm::make_filter_switch::MakeFilterSwitch,
        vdj_asm_asm::merge_per_sample_annotations::MergePerSampleAnnotations,
        vdj_asm_asm::write_ann_csv::WriteAnnCsv,
        vdj_asm_asm::subset_assembly_outs::SubsetAssemblyOuts,
    ];

    if args.cmd_martian {
        // Call the martian adapter
        let adapter = MartianAdapter::new(stage_registry);

        // Suppress any logging that would be emitted via crate log.
        let adapter = adapter.log_level(martian::LevelFilter::Warn);

        let retcode = adapter.run(args.arg_adapter);
        std::process::exit(retcode);
    } else if args.cmd_mro {
        // Create the mro for all the stages in this adapter
        martian_make_mro(HEADER, args.flag_file, args.flag_rewrite, mro_registry)?;
    } else {
        // If you need custom commands, implement them here
        unimplemented!()
    }
    Ok(())
}
