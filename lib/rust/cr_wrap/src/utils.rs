use crate::IntoExitCode;
use anyhow::{bail, ensure, Context, Result};
use clap::{self, Parser};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display, Formatter};
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::str::FromStr;

/// Convert an io::error to a string and strip "(os error 4)" from the end.
fn io_error_to_string(err: &std::io::Error) -> String {
    let s = err.to_string();
    s.strip_suffix(&format!(" (os error {})", err.raw_os_error().unwrap_or(0)))
        .unwrap_or(&s)
        .to_string()
}

/// Print an error chain.
pub fn print_error_chain(err: &anyhow::Error) {
    let error_chain = err.chain().join("\n\tCaused by: ");
    if let Some(io_err) = err.downcast_ref::<std::io::Error>() {
        let io_err_str = io_error_to_string(io_err);
        match err.chain().len() {
            1 => println!("ERROR: {io_err_str}"),
            2 => println!("ERROR: {io_err_str}: {err}"),
            _ => println!("ERROR: {error_chain}"),
        };
    } else {
        println!("ERROR: {error_chain}");
    };
}

/// Use this type for input paths that should
/// be canonicalized to a fully qualified path
/// before being passed to the MRO.
#[derive(Serialize, Deserialize, Clone)]
#[serde(transparent)]
pub struct CliPath {
    path: PathBuf,
}

impl From<PathBuf> for CliPath {
    fn from(path: PathBuf) -> Self {
        CliPath { path }
    }
}

impl From<&Path> for CliPath {
    fn from(path: &Path) -> Self {
        CliPath {
            path: path.to_path_buf(),
        }
    }
}

impl FromStr for CliPath {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<CliPath> {
        match Path::new(s).canonicalize() {
            Ok(p) => Ok(CliPath::from(p)),
            Err(e) => bail!(io_error_to_string(&e)),
        }
    }
}

impl Display for CliPath {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        Display::fmt(&self.path.display(), f)
    }
}

impl Debug for CliPath {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        Debug::fmt(&self.path, f)
    }
}

impl From<CliPath> for PathBuf {
    fn from(obj: CliPath) -> PathBuf {
        obj.path
    }
}

impl AsRef<Path> for CliPath {
    fn as_ref(&self) -> &Path {
        &self.path
    }
}

impl Deref for CliPath {
    type Target = Path;

    fn deref(&self) -> &Path {
        &self.path
    }
}

/// A subcommand that just passes through to an arbitrary command on the path
#[derive(Parser, Debug)]
#[clap(
    allow_hyphen_values = true,
    disable_help_flag = true,
    trailing_var_arg = true
)]
pub struct AllArgs {
    pub input: Vec<String>,
}

/// Run an arbitrary command from the PATH
/// Pass the arguments after the subcommand to a subprocess `cmd`.
pub fn external_subcommand(cmd: &str, args: &AllArgs) -> Result<ExitCode> {
    let mut command = Command::new(cmd);

    for r in &args.input {
        command.arg(r);
    }

    Ok(command
        .status()
        .with_context(|| format!("Running {cmd}"))?
        .into_exit_code())
}

/// Parse and validate an identifier, for use with Clap's value_parser.
/// A valid indentifier contains only letters, digits, underscores, and dashes.
pub fn validate_ascii_identifier(id: &str) -> Result<String> {
    ensure!(
        id.chars()
            .all(|c| matches!(c, '0'..='9' | 'A'..='Z' | 'a'..='z' | '_' | '-')),
        "must contain only letters, digits, underscores, and dashes."
    );
    Ok(String::from(id))
}

/// Max allowed length of the --id argument
const MAX_ID_LEN: usize = 64;

/// Parse and validate the --id argument, for use with Clap's value_parser.
pub fn validate_id(id: &str) -> Result<String> {
    ensure!(
        id.len() <= MAX_ID_LEN,
        "The --id parameter must be {MAX_ID_LEN} characters or less, please use a shorter string. \
         You can use the --description flag to store additional information.",
    );
    validate_ascii_identifier(id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_id() {
        for s in ["A-Z", "SAMPLE-123", "__SAMPLE-3", "__---00"] {
            assert!(validate_id(s).is_ok(), "{s} should be a valid ID");
        }
        for s in ["_A*S", "_SAMPLE/123", "(?:)", "ΔδΔ"] {
            assert!(validate_id(s).is_err(), "{s} should be an invalid ID");
        }
    }
}
