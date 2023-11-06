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

// handling pipeline environments
pub mod arc;
pub mod chemistry_arg;
mod deprecated_os;
pub mod env;
pub mod fastqs;
pub mod mkfastq;
pub mod mkref;
pub mod mrp_args;
pub mod shared_cmd;
pub mod targeted_compare;
pub mod utils;

use anyhow::{ensure, Context, Result};
use itertools::Itertools;
use mrp_args::MrpArgs;
use serde::Serialize;
use std::fs::{create_dir, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, ExitStatus, Stdio};

/// Convert something to an ExitCode.
trait IntoExitCode {
    fn into_exit_code(self) -> ExitCode;
}

impl IntoExitCode for ExitStatus {
    /// Convert an ExitStatus to an ExitCode.
    fn into_exit_code(self) -> ExitCode {
        self.code()
            .map_or(ExitCode::FAILURE, |x| ExitCode::from(x as u8))
    }
}

/// Generate an MRO invocation string
/// Args:
///  - `call`: the pipeline to invoke.
///  - `args`: the pipeline argument - will be serialized to json and passed to mrg
///  - `pipeline_mro_file`: the MRO file declaring the pipeline. Must be on the MROPATH.
pub fn make_mro<T: Serialize>(call: &str, args: &T, pipeline_mro_file: &str) -> Result<String> {
    let args =
        serde_json::to_value(args).with_context(|| "error serializing pipeline args to json")?;

    let json = serde_json::json!({
        "call": call,
        "args": args,
        "mro_file": pipeline_mro_file,
    });

    let msg = serde_json::to_string_pretty(&json)?;

    let mro_string =
        call_mrg(&msg).with_context(|| format!("failure calling mrg on json:\n {msg}"))?;
    Ok(mro_string)
}

pub fn make_mro_with_comment<T: Serialize>(
    call: &str,
    args: &T,
    pipeline_mro_file: &str,
    comment: &str,
) -> Result<String> {
    let comment_with_hashes: String = comment
        .lines()
        .map(|line| "# ".to_string() + line + "\n")
        .collect();
    Ok(comment_with_hashes + &make_mro(call, args, pipeline_mro_file)?)
}

fn call_mrg(msg: &str) -> Result<String> {
    let mut child = Command::new("mrg")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| "Failed to run mrg")?;

    {
        let stdin = child.stdin.as_mut().expect("Failed to open stdin of mrg");
        stdin
            .write_all(msg.as_bytes())
            .with_context(|| "Failed to write to stdin of mrg")?;
    }

    let invocation = child
        .wait_with_output()
        .with_context(|| "Failed to read stdout of mrg")?;
    let result = String::from_utf8(invocation.stdout)?;

    ensure!(
        invocation.status.success(),
        "Creating MRO invocation: {result}"
    );

    Ok(result)
}

pub enum MroInvocation {
    MroString(String),
    File(PathBuf),
}

/// Execute a Martian pipeline with mrp and return an ExitCode.
/// Args:
///  - `job_id` the job name. Job results will be written to this directory.
///  - `invocation`: string of the MRO invocation to run.
///  - `mrg_args`: additional parameters to pass to mrp controlling the job execution
///  - `dry_run`: emit the MRO that would be run but don't invoke with mrp.
pub fn execute(
    job_id: &str,
    invocation: &str,
    mrp_args: &MrpArgs,
    dry_run: bool,
) -> Result<ExitCode> {
    Ok(execute_to_status(job_id, invocation, mrp_args, dry_run)?.into_exit_code())
}

/// Execute a Martian pipeline with mrp and return an ExitStatus.
pub fn execute_to_status(
    job_id: &str,
    invocation: &str,
    mrp_args: &MrpArgs,
    dry_run: bool,
) -> Result<ExitStatus> {
    let inv = MroInvocation::MroString(invocation.to_string());
    execute_any(job_id, inv, mrp_args, dry_run)
}

pub fn execute_any(
    job_id: &str,
    invocation: MroInvocation,
    mrp_args: &MrpArgs,
    dry_run: bool,
) -> Result<ExitStatus> {
    let (mro_file, tmp) = match invocation {
        MroInvocation::MroString(mro) => {
            let filename: PathBuf = format!("__{job_id}.mro").into();
            let mut f = File::create(&filename).with_context(|| "couldn't open MRO file")?;
            f.write_all(mro.as_bytes())?;
            (filename, true)
        }
        MroInvocation::File(f) => (f, false),
    };

    if dry_run {
        println!("Dry Run Mode");
        println!();
        println!(
            "mrp command: {:?} {} {}",
            mro_file,
            job_id,
            mrp_args.get_args().join(" ")
        );
        println!("mro file: {mro_file:?}");

        return Ok(Command::new("mro")
            .arg("check")
            .arg(mro_file)
            .output()?
            .status);
    }

    let exit_status = run_mrp(job_id, &mro_file, mrp_args)?;
    if tmp {
        std::fs::remove_file(mro_file)?;
    }
    let output_dir = mrp_args.output_dir.as_deref().unwrap_or(job_id);
    let _ = run_tarmri(output_dir, exit_status)?;
    Ok(exit_status)
}

/// Set CMDLINE to the command line arguments, needed to create the pipeline output file _cmdline.
fn set_env_cmdline() {
    std::env::set_var("CMDLINE", std::env::args().join(" "));
}

// Set COLUMNS to 80 when the terminal size is unknown.
// Wrap the output of --help to 80 columns when the terminal size is unknown.
// The default value of clap is 100.
fn set_env_columns() {
    if terminal_size::terminal_size().is_none() && std::env::var_os("COLUMNS").is_none() {
        std::env::set_var("COLUMNS", "80");
    }
}

// Set environment variables.
pub fn set_env_vars() {
    set_env_cmdline();
    set_env_columns();
}

fn run_mrp(job_id: &str, mro_path: &Path, mrp_args: &MrpArgs) -> Result<ExitStatus> {
    if let Some(output_dir) = &mrp_args.output_dir {
        // Create output_dir to ensure that it is created successfully.
        match create_dir(output_dir) {
            Ok(_) => (),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(error) => return Err(error).context(output_dir.clone()),
        }
    };

    let args = mrp_args.get_args();
    Command::new("mrp")
        .arg(mro_path)
        .arg(job_id)
        .args(&args)
        .status()
        .context(format!(
            "running mrp {} {job_id} {}",
            mro_path.display(),
            args.join(" ")
        ))
}

fn run_tarmri(output_dir: &str, exit_status: ExitStatus) -> Result<ExitCode> {
    let mut cmd = Command::new("tarmri");

    // mro and output path
    cmd.arg(output_dir);
    cmd.arg(exit_status.code().unwrap_or(1).to_string());

    let cmdline: Vec<String> = std::env::args().collect();
    let cmdline = cmdline.join(" ");
    cmd.arg(cmdline);

    Ok(cmd
        .status()
        .with_context(|| "Error reading from tarmri")?
        .into_exit_code())
}

pub fn check_deprecated_os() -> Result<()> {
    if std::env::var("TENX_IGNORE_DEPRECATED_OS")
        .map(|s| s != "0")
        .unwrap_or(false)
    {
        return Ok(());
    }
    // by running this check in a separate subprocess we can do so inside a
    // sanitized environment sans e.g. LD_PRELOAD or LD_LIBRARY_PATH,
    // which could render these checks unreliable
    // if we're running `cellranger oscheck`, don't recurse
    let output = Command::new(std::env::current_exe()?)
        .arg("oscheck")
        .env("TENX_IGNORE_DEPRECATED_OS", "1")
        .output()
        .expect("failed to run oscheck");

    let msg = format!(
        "{}\n{}\n",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    ensure!(output.status.success(), "{}", msg.trim_end_matches('\n'));
    eprint!("{msg}");
    Ok(())
}
