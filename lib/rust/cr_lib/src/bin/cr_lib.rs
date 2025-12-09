//! cr_lib
#![deny(missing_docs)]

use anyhow::{Context, Result, bail};
use cr_lib::{load_dist, write_web_summary_html};
use docopt::Docopt;
use martian::prelude::*;
use serde::Deserialize;
use std::fs::File;

const HEADER: &str = "# Copyright 2023 10x Genomics, Inc. All rights reserved.";

const USAGE: &str = "
Test Rust stage
Usage:
  cr_lib martian <adapter>...
  cr_lib mro [--file=<filename>] [--rewrite]
  cr_lib rebuild-multi-ws <json-in> <html-out>
  cr_lib --help
Options:
     --help            Show this screen.
";

#[derive(Deserialize)]
struct Args {
    // Martian interface
    cmd_martian: bool,
    cmd_mro: bool,
    cmd_rebuild_multi_ws: bool,
    arg_adapter: Vec<String>,
    flag_file: Option<String>,
    flag_rewrite: bool,
    arg_json_in: Option<String>,
    arg_html_out: Option<String>,
}

fn main() -> Result<()> {
    #[allow(clippy::wildcard_imports)]
    use cr_lib::stages::*;

    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    let (stage_registry, mro_registry) = martian_stages![
        AlignAndCount,
        BarcodeCorrection,
        BuildPerSampleVdjWsContents,
        CallTagsGenetic,
        CallTagsOH,
        CallTagsRTL,
        CheckBarcodesCompatibility,
        CheckBarcodesCompatibilityVdj,
        CheckSingleBeamMode,
        CollateMetrics,
        CollateProbeMetrics,
        ComputeAntigenVdjMetrics,
        CopyChemistrySpec,
        CountAlleles,
        CreateMultiGraph,
        DemuxProbeBcMatrix,
        DetectChemistry,
        DetectVdjReceptor,
        ExpectSingleBarcodeWhitelist,
        ExtractSingleChemistry,
        GenerateCasWebsummary,
        GetChemistryDef,
        GetGdnaMetrics,
        LogicNot,
        MakeCorrectionMap,
        MakeShard,
        MergeGemWellFiles,
        MergeMetrics,
        MultiPreflight,
        MultiSetupChunks,
        ParseMultiConfig,
        PickBeamAnalyzer,
        RustBridge,
        SetupReferenceInfo,
        SetupVdjAnalysis,
        SetupVDJDemux,
        WriteBarcodeIndex,
        WriteBarcodeSummary,
        WriteGeneIndex,
        WriteH5Matrix,
        WriteMinimapIndex,
        WriteMatrixMarket,
        WriteMoleculeInfo,
        WriteMultiWebSummary,
        WritePosBam,
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
    } else if args.cmd_rebuild_multi_ws {
        let Some(json_path) = args.arg_json_in else {
            bail!("provide JSON path as first argument");
        };
        let Some(output_path) = args.arg_html_out else {
            bail!("provide HTML output path as second argument");
        };
        rebuild_multi_ws(&json_path, &output_path)?;
    } else {
        // If you need custom commands, implement them here
        unimplemented!()
    }
    Ok(())
}

/// Helper tool to inject a pre-built JSON data payload into the current multi websummary.
/// To extract the data payload from an existing websummary, use the extract-web-summary.sh
/// tool in the root of the websummary submodule.
///
/// Run this via the Bazel target to include the websummary build dependencies.
///
/// Note that websummaries will not be byte-for-byte identical, because the
/// extracted and re-written JSON payload will not be in the same key order,
/// and floating-point numbers may change in the lowest digit.
fn rebuild_multi_ws(json_path: &str, output_path: &str) -> Result<()> {
    let json_file = File::open(json_path).with_context(|| json_path.to_string())?;
    let parsed_data: serde_json::Value = serde_json::from_reader(&json_file)?;
    let mut output_file = File::create(output_path).with_context(|| output_path.to_string())?;
    write_web_summary_html(&parsed_data, &load_dist()?, &mut output_file)?;
    Ok(())
}
