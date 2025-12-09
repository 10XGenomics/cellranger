#![expect(missing_docs)]
use crate::IntoExitCode;
use crate::utils::AllArgs;
use anyhow::{Context, Result};
use std::process::{Command, ExitCode};

/// Currently we capture all the arguments and pass them through to
/// the old shell wrapper code for mkfastq.
/// (full Rust-ification is an early work in progress for now)
pub fn run_mkfastq(args: &AllArgs, legacy_wrapper: &str) -> Result<ExitCode> {
    let mut warning_msg = r"
The `cellranger mkfastq` pipeline is deprecated and will be removed in a future release.
Please use Illumina's BCL Convert to generate Cell Ranger-compatible FASTQ files.
For detailed guidance, refer to the Generating FASTQs support page:
https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-direct-demultiplexing
    ";

    if legacy_wrapper == "_cellranger_arc_internal" {
        // Swap the warning message to point to the correct support page for CellRanger ARC
        warning_msg = r"
The `cellranger-arc mkfastq` pipeline is deprecated and will be removed in a future release.
Please use Illumina's BCL Convert to generate Cell Ranger ARC-compatible FASTQ files.
For detailed guidance, refer to the Generating FASTQs support page:
https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/inputs/generating-fastqs-mkfastq
        ";
    }

    eprintln!("{warning_msg}");

    let mut cmd = Command::new(legacy_wrapper);
    cmd.arg("mkfastq");

    for r in &args.input {
        cmd.arg(r);
    }

    let result = cmd
        .status()
        .with_context(|| format!("Running {legacy_wrapper}"));

    eprintln!("{warning_msg}");

    Ok(result?.into_exit_code())
}
