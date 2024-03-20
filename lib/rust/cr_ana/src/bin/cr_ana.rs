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
use clap::Parser;
use martian::prelude::*;

#[derive(Debug, Parser)]
enum Args {
    Martian {
        adapter_args: Vec<String>,
    },
    Mro {
        #[clap(long)]
        file: Option<String>,
        #[clap(long)]
        rewrite: bool,
    },
}

const HEADER: &str = "# Copyright 2023 10x Genomics, Inc. All rights reserved.";

fn main() -> Result<()> {
    let args: Args = Args::parse();

    let (stage_registry, mro_registry) = martian_stages![
        cr_ana::stages::diff_exp_stage::DiffExpStage,
        cr_ana::stages::graph_clustering::GraphClusteringStage,
        cr_ana::stages::hierarchical_clustering::HierarchicalClusteringStage,
        cr_ana::stages::pca::PcaStage,
        cr_ana::stages::pca2::Pca2Stage,
        cr_ana::stages::tsne::TsneStage,
        cr_ana::stages::umap::UmapStage,
    ];

    match args {
        Args::Martian { adapter_args } => {
            let adapter = MartianAdapter::new(stage_registry).log_level(martian::LevelFilter::Info);
            let retcode = adapter.run(adapter_args);
            std::process::exit(retcode);
        }
        Args::Mro { file, rewrite } => {
            martian_make_mro(HEADER, file, rewrite, mro_registry)?;
        }
    }

    Ok(())
}
