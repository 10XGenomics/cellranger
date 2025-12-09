//! cr_aggr
#![deny(missing_docs)]

use anyhow::Result;
use docopt::Docopt;
use martian::prelude::*;
use serde::Deserialize;

const HEADER: &str = "# Copyright 2023 10x Genomics, Inc. All rights reserved.";

const USAGE: &str = "
Test Rust stage
Usage:
  cr_aggr martian <adapter>...
  cr_aggr mro [--file=<filename>] [--rewrite]
  cr_aggr --help
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
        cr_aggr::merge_molecules::MergeMolecules,
        cr_aggr::process_vdj_proto::ProcessVdjProto,
        cr_aggr::setup_vdj_aggr::SetupVdjAggr,
        cr_aggr::parse_aggr_csv::ParseAggrCsv,
        cr_aggr::write_contig_proto::WriteContigProto,
        cr_aggr::match_vdj_outs::MatchVdjOuts,
        cr_aggr::write_aggr_ann::WriteAggrAnn,
        cr_aggr::write_ws_json::WriteWsJson,
        cr_aggr::create_antigen_clonotype_clustermap::CreateClonotypeClustermap,
        cr_aggr::run_enclone_aggr::RunEncloneAggr,
    ];

    if args.cmd_martian {
        // Call the martian adapter
        let adapter = MartianAdapter::new(stage_registry);

        // Suppress any logging that would be emitted via crate log.
        let adapter = adapter.log_level(LevelFilter::Warn);

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
