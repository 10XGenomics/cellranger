// Possibly useful in the future -- these are currently unneccesary.

use crate::telemetry::{CollectTelemetry, TelemetryCollector};
use crate::utils::AllArgs;
use crate::{check_deprecated_os, IntoExitCode};
use anyhow::{Context, Result};
use clap::parser::ValueSource;
use clap::{CommandFactory, Parser};
use cr_types::constants::COMMAND_LINE_ENV_VARIABLE_NAME;
use serde::{Deserialize, Serialize};
use shell_escape::escape;
use std::collections::HashMap;
use std::env;
use std::ffi::OsString;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

/// Specfiy how commands initialize the environment.
pub trait Initialize: Parser + CommandFactory + CollectTelemetry {
    /// Return the path to the env JSON override file, if this command provides it.
    fn override_env_json(&self) -> Option<&Path>;

    fn initialize() -> Result<(Self, PkgEnv, TelemetryCollector)> {
        // Setup environment variables.
        set_env_cmdline();
        set_env_columns();

        let bare_cmd = Self::command();

        // Canonicalize the CLI args and use the matches for parsing.
        let (canonical_cli_args, subcommand, cmd) = {
            let matches = Self::command().get_matches();
            let (subcommand, submatches) = matches.subcommand().unwrap();

            let argmap: HashMap<_, _> = bare_cmd
                .get_subcommands()
                .find(|sc| sc.get_name() == subcommand)
                .unwrap()
                .get_arguments()
                .filter_map(|arg| arg.get_long().map(|f| (arg.get_id(), f)))
                .collect();
            let args = submatches
                .ids()
                .filter(|id| argmap.contains_key(id))
                .filter(|id| submatches.value_source(id.as_str()) == Some(ValueSource::CommandLine))
                .filter_map(|id| {
                    submatches
                        .get_raw(id.as_str())
                        .map(|rv| rv.map(|val| (argmap[id], val)))
                })
                .flatten()
                .map(|(flag, val)| {
                    let mut arg: OsString = format!("--{flag}=").into();
                    arg.push(val);
                    arg
                })
                .collect();

            (
                args,
                subcommand.to_string(),
                Self::from_arg_matches(&matches)?,
            )
        };

        // Initialize pkg environment.
        let pkg_env = setup_env(
            cmd.override_env_json(),
            bare_cmd.get_name().to_string(),
            bare_cmd.get_version().unwrap().to_string(),
        )?;
        // Initialize telemetry collection, and collect initial telemetry.
        let mut telemetry =
            TelemetryCollector::new_for_command(&cmd, subcommand, canonical_cli_args, &pkg_env);
        telemetry.collect(None, None, None);

        check_deprecated_os()?;
        Ok((cmd, pkg_env, telemetry))
    }
}

/// Set CMDLINE to the command line arguments, needed to create the pipeline output file _cmdline.
fn set_env_cmdline() {
    let cmdline: Vec<_> = std::env::args().map(|x| escape(x.into())).collect();
    std::env::set_var(COMMAND_LINE_ENV_VARIABLE_NAME, cmdline.join(" "));
}

// Set COLUMNS to 80 when the terminal size is unknown.
// Wrap the output of --help to 80 columns when the terminal size is unknown.
// The default value of clap is 100.
fn set_env_columns() {
    if terminal_size::terminal_size().is_none() && std::env::var_os("COLUMNS").is_none() {
        std::env::set_var("COLUMNS", "80");
    }
}

pub(crate) const TENX_COPYRIGHT: &str = "Copyright 2023 10x Genomics, Inc. All rights reserved.";

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "lowercase")]
pub enum EnvModType {
    /// Specifies that the environment variable should be replaced with the
    /// given value, e.g. `export KEY=value`
    #[serde(alias = "string")]
    Set,
    //// Specifies that the environment variable should have the given value
    //// prepended, e.g. `export KEY=value:$KEY`
    #[serde(alias = "path_prepend")]
    Prepend,
    /// Specifies that the environment variable should have the given value
    /// appended, e.g. `export KEY=$KEY:value`
    #[serde(alias = "path_append")]
    Append,
    /// Specifies that the given environment variable should be unset.
    /// The value is ignored in this case.
    #[serde(alias = "unset")]
    Unset,
    /// Specifies that the given environment variable should be unset.
    /// The value of a variable VAR should be stored in  _TENX_VAR.
    #[serde(alias = "setaside")]
    Setaside,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct EnvironmentMod {
    pub key: String,

    pub value: String,

    #[serde(rename = "type")]
    pub mod_type: EnvModType,
}

fn get_paths(var: &str) -> Vec<PathBuf> {
    let current_paths: Vec<PathBuf> = env::var_os(var)
        .map(|v| env::split_paths(&v).collect())
        .unwrap_or_default();

    current_paths
}

fn new_path(base: &Path, rest: &str) -> PathBuf {
    let mut path = PathBuf::from(base);
    path.push(rest);
    path
}

impl EnvironmentMod {
    pub fn apply(&self, base_path: &Path) -> Result<()> {
        match &self.mod_type {
            EnvModType::Set => {
                env::set_var(&self.key, &self.value);
            }

            EnvModType::Prepend => {
                let mut paths: Vec<PathBuf> = vec![new_path(base_path, &self.value)];
                paths.extend(get_paths(&self.key));
                env::set_var(&self.key, env::join_paths(paths)?);
            }

            EnvModType::Append => {
                let mut paths = get_paths(&self.key);
                paths.push(new_path(base_path, &self.value));
                env::set_var(&self.key, env::join_paths(paths)?);
            }

            EnvModType::Setaside => {
                let setaside_key = format!("_TENX_{}", self.key);
                let old_val = env::var_os(&self.key).unwrap_or_default();
                env::set_var(setaside_key, old_val);
                if self.value.is_empty() {
                    env::remove_var(&self.key);
                } else {
                    env::set_var(&self.key, &self.value);
                }
            }

            EnvModType::Unset => {
                env::remove_var(&self.key);
            }
        }

        Ok(())
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct EnvSpec {
    /// The path to the mrp executable for this package.
    mro_path: String,

    // Name of base pipeline
    name: String,

    envs: Vec<EnvironmentMod>,
}

impl EnvSpec {
    fn apply(&self, build_path: &Path) -> Result<()> {
        // appy env modifications
        for env_mod in &self.envs {
            env_mod.apply(build_path)?;
        }

        // set MROPATH
        let mro_paths: Vec<_> = env::split_paths(&self.mro_path).collect();

        for p in mro_paths {
            let m = EnvironmentMod {
                key: "MROPATH".to_string(),
                value: p.to_str().unwrap().to_string(),
                mod_type: EnvModType::Append,
            };

            m.apply(build_path)?;
        }

        Ok(())
    }
}

pub struct PkgEnv {
    pub build_path: PathBuf,
    pub tenx_product: String,
    pub tenx_version: String,
}

impl PkgEnv {
    /// Return the path to the binary `cmd`, specified relative to the installation root.
    pub fn subcmd_path(&self, cmd: &str) -> PathBuf {
        self.build_path.join(cmd)
    }

    /// Run the binary `cmd`, specified relative to the installation root.
    pub fn run_subcmd(&self, cmd: &str, args: &AllArgs) -> Result<ExitCode> {
        let path = &self.subcmd_path(cmd);

        let mut scriptdir = self.build_path.clone();
        scriptdir.push("bin");
        env::set_var("TENX_SCRIPTDIR", scriptdir);

        let subcmd = path.file_name().unwrap();
        env::set_var("TENX_SUBCMD", subcmd);

        let mut command = Command::new(path);
        for r in &args.input {
            command.arg(r);
        }

        Ok(command
            .status()
            .with_context(|| format!("Running {cmd}"))?
            .into_exit_code())
    }
}

/// Prepare the environment for running a pipestance in MRP.
/// Find the root of the build dir, load `env.json` and update the environment accordingly.
fn setup_env(
    override_env_json: Option<&Path>,
    tenx_product: String,
    tenx_version: String,
) -> Result<PkgEnv> {
    // Path of installation is one above the binary
    let exe = std::env::current_exe()?;
    let exe_name = exe.file_name().unwrap().to_str().unwrap();
    let exe_dir = exe.parent().unwrap();

    let runfiles_env_json = exe_dir.join(format!("{exe_name}.runfiles/cellranger/env.json"));
    let exe_env_json = exe_dir.parent().unwrap().join("env.json");

    let env_json = if let Some(override_env_json) = &override_env_json {
        override_env_json
    } else if runfiles_env_json.exists() {
        runfiles_env_json.as_path()
    } else {
        exe_env_json.as_path()
    };

    let build_path = env_json.parent().unwrap();
    let env_spec: EnvSpec = serde_json::from_reader(
        File::open(env_json).with_context(|| env_json.display().to_string())?,
    )?;
    env_spec
        .apply(build_path)
        .with_context(|| "error setting up pipeline environment")?;

    // additional environment needed by some subcommands
    env::set_var("TENX_PRODUCT", &tenx_product);
    env::set_var("TENX_VERSION", &tenx_version);
    env::set_var("TENX_COPYRIGHT", TENX_COPYRIGHT);

    Ok(PkgEnv {
        build_path: build_path.to_path_buf(),
        tenx_product,
        tenx_version,
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_cs_env_json() -> Result<()> {
        let mut rdr = File::open("test/test_cs_env.json")?;
        let env_spec: EnvSpec = serde_json::from_reader(&mut rdr)?;

        println!("{env_spec:#?}");
        Ok(())
    }
}
