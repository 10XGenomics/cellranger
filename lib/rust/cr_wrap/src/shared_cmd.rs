use crate::deprecated_os::oscheck;
use crate::env::PkgEnv;
use crate::mkref::Mkref;
use crate::utils::{external_subcommand, AllArgs};
use anyhow::Result;
use clap::{self, Parser};
use std::process::ExitCode;

// Shared between spaceranger/cellranger
#[derive(Parser, Debug)]
#[allow(clippy::large_enum_variant)]
pub enum RnaSharedCmd {
    /// Convert a feature-barcode matrix to CSV format
    #[clap(name = "mat2csv")]
    Mat2csv(AllArgs),

    /// Prepare a reference for use with 10x analysis software. Requires a GTF and FASTA.
    #[clap(name = "mkref")]
    Mkref(Mkref),

    /// Filter a GTF file by attribute prior to creating a 10x reference
    #[clap(name = "mkgtf")]
    Mkgtf(AllArgs),
}

// Shared between cellranger/spaceranger/cellranger-arc/cellranger-atac
#[derive(Parser, Debug)]
pub enum SharedCmd {
    /// Upload analysis logs to 10x Genomics support.
    #[clap(name = "upload")]
    Upload(AllArgs),

    /// Collect Linux system configuration information.
    #[clap(name = "sitecheck")]
    Sitecheck(AllArgs),
}

#[derive(Parser, Debug)]
pub enum HiddenCmd {
    /// Invoke the python interpreter inside the 10x software installation environment
    #[clap(name = "python", hide = true)]
    Python(AllArgs),

    /// Invoke redstone inside the 10x software installation environment
    #[clap(name = "redstone", hide = true)]
    Redstone(AllArgs),

    /// Invoke mrc inside the 10x software installation environment
    #[clap(name = "mro", hide = true)]
    Mro(AllArgs),

    /// Invoke mrp inside the 10x software installation environment
    #[clap(name = "mrp", hide = true)]
    Mrp(AllArgs),

    /// Invoke bamtofastq inside the 10x software installation environment
    #[clap(name = "bamtofastq", hide = true)]
    Bamtofastq(AllArgs),

    /// Invoke tarmri to bundle analysis pipeline output files for transmission to 10x Genomics.
    #[clap(name = "tarmri", hide = true)]
    Tarmri(AllArgs),

    /// Invoke any arbitrary command in the 10x software installation environment
    #[clap(name = "pass", hide = true)]
    Pass(AllArgs),

    /// Deprecated OS checks from within a sanitized environment
    #[clap(name = "oscheck", hide = true)]
    OsCheck(AllArgs),
}

pub fn run_rna_shared(pkg_env: &PkgEnv, args: RnaSharedCmd) -> Result<ExitCode> {
    match args {
        RnaSharedCmd::Mat2csv(args) => pkg_env.run_subcmd("bin/rna/mat2csv", &args),
        RnaSharedCmd::Mkref(mut args) => {
            args.shared.populate_version(pkg_env.tenx_version);
            args.execute()
        }
        RnaSharedCmd::Mkgtf(args) => pkg_env.run_subcmd("bin/rna/mkgtf", &args),
    }
}

pub fn run_shared(pkg_env: &PkgEnv, args: SharedCmd) -> Result<ExitCode> {
    match args {
        SharedCmd::Upload(args) => pkg_env.run_subcmd("bin/tenkit/upload", &args),
        SharedCmd::Sitecheck(args) => pkg_env.run_subcmd("bin/tenkit/sitecheck", &args),
    }
}

pub fn run_hidden(args: HiddenCmd) -> Result<ExitCode> {
    match args {
        // pass-thru commands
        HiddenCmd::Python(args) => external_subcommand("python", &args),
        HiddenCmd::Redstone(args) => external_subcommand("redstone", &args),
        HiddenCmd::Mro(args) => external_subcommand("mro", &args),
        HiddenCmd::Mrp(args) => external_subcommand("mrp", &args),
        HiddenCmd::Bamtofastq(args) => external_subcommand("bamtofastq", &args),
        HiddenCmd::Tarmri(args) => external_subcommand("tarmri", &args),
        HiddenCmd::OsCheck(_) => oscheck().map(|()| ExitCode::SUCCESS),
        // pass through into the general CR env
        HiddenCmd::Pass(args) => {
            if !args.input.is_empty() {
                let cmd = args.input[0].clone();
                let args = AllArgs {
                    input: args.input.into_iter().skip(1).collect(),
                };
                external_subcommand(&cmd, &args)
            } else {
                Ok(ExitCode::SUCCESS)
            }
        }
    }
}
