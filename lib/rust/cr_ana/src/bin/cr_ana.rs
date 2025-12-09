//! cr_ana
#![deny(missing_docs)]

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
        cr_ana::stages::pca2::RunPca2,
        cr_ana::stages::tsne::RunTsne,
        cr_ana::stages::umap::UmapStage,
    ];

    match args {
        Args::Martian { adapter_args } => {
            let adapter = MartianAdapter::new(stage_registry).log_level(LevelFilter::Info);
            let retcode = adapter.run(adapter_args);
            std::process::exit(retcode);
        }
        Args::Mro { file, rewrite } => {
            martian_make_mro(HEADER, file, rewrite, mro_registry)?;
        }
    }

    Ok(())
}
