use crate::env::PkgEnv;
use crate::mrp_args::MrpArgs;
use crate::telemetry::TelemetryCollector;
use crate::utils::{validate_ascii_identifier, CliPath};
use crate::{execute_to_status, make_mro, IntoExitCode};
use anyhow::{bail, ensure, Result};
use clap::{self, value_parser, Parser};
use cr_types::reference::reference_info::MULTI_GENOME_SEPARATOR;
use serde::Serialize;
use std::fs;
use std::path::{self, Path};
use std::process::ExitCode;

#[derive(Parser, Debug, Clone, Serialize)]
pub struct MkrefShared {
    /// Maximum memory (GB) used.
    #[clap(long = "memgb", default_value_t = 16, value_parser = value_parser!(u32).range(1..))]
    pub mem_gb: u32,

    /// Optional reference version string to include with reference.
    #[clap(long = "ref-version")]
    pub ref_version: Option<String>,

    /// Version of the mkref program. Populated after arg parsing.
    #[clap(skip)]
    pub mkref_version: String,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[serde(skip)]
    #[clap(long)]
    pub dry: bool,

    #[serde(skip)]
    #[clap(flatten)]
    pub mrp: MrpArgs,
}

impl MkrefShared {
    /// Set mkref_version to "x.y.z". `version_description` is of the form
    /// product-x.y.z or yyyy.mmdd.z or yyyy.mmdd.z-n-hash.
    pub fn populate_version(&mut self, pkg_env: &PkgEnv) {
        let version = if pkg_env.tenx_version.matches('-').count() == 1 {
            // product-x.y.z
            pkg_env.tenx_version.split_once('-').unwrap().1
        } else {
            // yyyy.mmdd.z or yyyy.mmdd.z-n-hash
            &pkg_env.tenx_version
        };
        self.mkref_version = format!("{0}-{1}", pkg_env.tenx_product, version);
    }
}

#[derive(Parser, Debug, Clone, Serialize)]
pub struct Mkref {
    /// Unique genome name, used to name output folder [a-zA-Z0-9_-]+.
    /// Specify multiple genomes by specifying this argument multiple times; the
    /// output folder will be <name1>_and_<name2>.
    #[clap(long = "genome", required = true, value_parser = validate_ascii_identifier)]
    pub genome_names: Vec<String>,

    /// Path to FASTA file containing your genome reference.
    /// Specify multiple genomes by specifying this argument multiple times.
    #[clap(long = "fasta", required = true)]
    pub fasta_files: Vec<CliPath>,

    /// Path to genes GTF file containing annotated genes for your genome reference.
    /// Specify multiple genomes by specifying this argument multiple times.
    #[clap(long = "genes", required = true)]
    pub gtf_files: Vec<CliPath>,

    /// Number of threads used during STAR genome index generation. Defaults to 1.
    #[clap(long = "nthreads", default_value_t = 1, value_parser = value_parser!(u32).range(1..))]
    pub num_threads: u32,

    #[clap(flatten)]
    #[serde(flatten)]
    pub shared: MkrefShared,
}

impl Mkref {
    /// Execute the mkref tool as a martian pipeline.
    /// If successful, rename the outs folder into the current directory and delete the rest of
    /// the pipestance folder.
    pub fn execute(&self, telemetry: &mut TelemetryCollector) -> Result<ExitCode> {
        ensure!(
            self.genome_names.len() == self.fasta_files.len(),
            "provided {} genome names but {} FASTA files",
            self.genome_names.len(),
            self.fasta_files.len(),
        );

        ensure!(
            self.genome_names.len() == self.gtf_files.len(),
            "provided {} genome names but {} GTF files",
            self.genome_names.len(),
            self.gtf_files.len(),
        );

        let output_dir_name = self.genome_names.join(MULTI_GENOME_SEPARATOR);
        let pipestance_name = format!("mkref_{output_dir_name}");

        // Check to make sure the final output path is clear, but only if the
        // user didn't specify a custom output directory.
        if self.shared.mrp.output_dir.is_none() {
            ensure!(
                !Path::new(&output_dir_name).try_exists()?,
                "output path '{output_dir_name}' already exists"
            );
        }

        let mro = make_mro("MAKE_REFERENCE", self, "rna/make_reference.mro")?;
        let exit_status = execute_to_status(
            &pipestance_name,
            &mro,
            &self.shared.mrp,
            self.shared.dry,
            telemetry,
        )?;
        if !exit_status.success() {
            return Ok(exit_status.into_exit_code());
        }

        move_outputs(&pipestance_name, &self.shared.mrp, &output_dir_name)?;

        let abs_path = path::absolute(&output_dir_name)?;
        let abs_path_str = abs_path.to_string_lossy();
        println!(
            ">>> Reference successfully created! <<<\n\
            \n\
            Reference location:\n\
            {abs_path_str}\n\
            \n\
             You can now specify this reference on the command line:\n\
             cellranger count --transcriptome={output_dir_name} ..."
        );
        Ok(ExitCode::SUCCESS)
    }
}

#[derive(Parser, Debug, Clone, Serialize)]
pub struct Mkvdjref {
    /// Unique genome name, used to name output folder [a-zA-Z0-9_-]+.
    #[clap(long = "genome", required = true, value_parser = validate_ascii_identifier)]
    pub genome_name: String,

    /// Path to FASTA file containing your genome reference.
    #[clap(long = "fasta")]
    pub fasta_file: Option<CliPath>,

    /// Path to genes GTF file containing annotated genes for your genome reference.
    /// Specify multiple genomes by specifying this argument multiple times.
    #[clap(long = "genes")]
    pub gtf_files: Vec<CliPath>,

    /// Path to a FASTA file that directly specifies V(D)J sequences.
    /// This is mutually exclusive with the "fasta" and "genes" args.
    #[clap(long = "seqs")]
    pub seq_file: Option<CliPath>,

    /// Path to text file with transcript IDs to ignore. This
    /// file should have one transcript ID per line where
    /// the IDs correspond to the "transcript_id" key in the
    /// GTF info column.
    #[clap(long = "rm-transcripts")]
    pub remove_transcripts_file: Option<CliPath>,

    #[clap(flatten)]
    #[serde(flatten)]
    pub shared: MkrefShared,
}

impl Mkvdjref {
    /// Execute the mkref tool as a martian pipeline.
    /// If successful, rename the outs folder into the current directory and delete the rest of
    /// the pipestance folder.
    pub fn execute(&self, telemetry: &mut TelemetryCollector) -> Result<ExitCode> {
        let output_dir_name = &self.genome_name;
        let pipestance_name = format!("mkvdjref_{output_dir_name}");

        let output_dir = Path::new(output_dir_name);

        // Check to make sure the final output path is clear, but only if the
        // user didn't specify a custom output directory.
        if self.shared.mrp.output_dir.is_none() {
            ensure!(
                !output_dir.try_exists()?,
                "output path '{output_dir_name}' already exists"
            );
        }

        // Check for mutually-exclusive args.
        match (self.fasta_file.is_some(), self.seq_file.is_some()) {
            (true, false) => {
                ensure!(
                    !self.gtf_files.is_empty(),
                    "please provide at least one genes file"
                );
            }
            (false, true) => {
                ensure!(
                    self.gtf_files.is_empty(),
                    "--genes is mututally exclusive with --seqs"
                );
                ensure!(
                    self.remove_transcripts_file.is_none(),
                    "the --rm-transcripts option cannot be used with --seqs"
                );
            }
            (true, true) => {
                bail!("--fasta cannot be specified if --seqs is specified.");
            }
            (false, false) => {
                bail!("must specify either (--fasta and --genes) or --seqs");
            }
        }

        let mro = make_mro("MAKE_VDJ_REFERENCE", self, "rna/make_vdj_reference.mro")?;
        let exit_status = execute_to_status(
            &pipestance_name,
            &mro,
            &self.shared.mrp,
            self.shared.dry,
            telemetry,
        )?;
        if !exit_status.success() {
            return Ok(exit_status.into_exit_code());
        }

        // If the build succeeded, extract the contents of the outs
        // folder and move it to the final path.
        move_outputs(&pipestance_name, &self.shared.mrp, output_dir_name)?;

        println!(">>> Reference successfully created! <<<");
        Ok(ExitCode::SUCCESS)
    }
}

fn move_outputs(pipestance_name: &str, args: &MrpArgs, output_ref_dir_name: &str) -> Result<()> {
    let pipestance_path = Path::new(args.output_dir.as_deref().unwrap_or(pipestance_name));
    let outs_path = pipestance_path.join("outs");
    let output_reference_path = outs_path.join("reference");

    if args.output_dir.is_some() {
        // If the user specified an output directory, the pipestance will have
        // run in it. Delete everything except the output reference, then move
        // it into the specified directory.
        // Delete everything besides outs.
        for item in fs::read_dir(pipestance_path)? {
            let item = item?;
            if item.file_name() == "outs" {
                continue;
            }
            if item.file_type()?.is_dir() {
                fs::remove_dir_all(item.path())?;
            } else {
                fs::remove_file(item.path())?;
            }
        }
        // Now lift the contents of the reference folder up into the user-
        // specified one, and delete the outs folder.
        for item in fs::read_dir(output_reference_path)? {
            let item = item?;
            fs::rename(item.path(), pipestance_path.join(item.file_name()))?;
        }
        // Don't fail if we didn't successfully delete.
        let _ = fs::remove_dir_all(outs_path.clone());
    } else {
        // We ran in the current working dir, so move the outputs into it
        // and delete the pipestance directory.
        fs::rename(output_reference_path, output_ref_dir_name)?;
        // Don't fail if we didn't successfully delete.
        let _ = fs::remove_dir_all(pipestance_path);
    }
    Ok(())
}
