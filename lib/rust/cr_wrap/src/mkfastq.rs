use crate::utils::AllArgs;
use crate::IntoExitCode;
use anyhow::{Context, Result};
use std::process::{Command, ExitCode};

/// Currently we capture all the arguments and pass them through to
/// the old shell wrapper code for mkfastq.
/// (full Rust-ification is an early work in progress for now)
pub fn run_mkfastq(args: &AllArgs, legacy_wrapper: &str) -> Result<ExitCode> {
    let mut cmd = Command::new(legacy_wrapper);
    cmd.arg("mkfastq");

    for r in &args.input {
        cmd.arg(r);
    }

    Ok(cmd
        .status()
        .with_context(|| format!("Running {legacy_wrapper}"))?
        .into_exit_code())
}
