use crate::utils::CliPath;
use anyhow::{anyhow, bail, Result};
use clap::builder::NonEmptyStringValueParser;
use clap::{self, Parser};
use cr_types::sample_def::{FastqMode, SampleDef, COUNT_HELP};
use cr_types::LibraryType;
use fastq_set::filenames::fastq_dir::FastqChecker;
use itertools::Itertools;
use metric::TxHashSet;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

pub fn get_target_set_name(target_set_csv: &Option<CliPath>) -> Option<String> {
    target_set_csv
        .as_ref()
        .map(|t| t.file_stem().unwrap().to_string_lossy().into_owned())
}

macro_rules! impl_get_sample_defs {
    ($name:ident) => {
        impl $name {
            /// Convert the command-line FASTQ finding arguments into SampleDefs
            pub fn get_sample_defs(&self, library_type: LibraryType) -> Result<Vec<SampleDef>> {
                let mut sample_defs = Vec::new();

                // Get the FASTQs paths that are compatible with the given bcl2fastq project
                let fastq_paths = get_fastq_paths(&self.fastqs, self.project.as_deref())?;

                let mut samples_found = TxHashSet::default();

                for path in fastq_paths {
                    let sample_names = FastqChecker::bcl2fastq_check_and_infer_sample_names(
                        &path,
                        &self.sample,
                        &self.lanes,
                        COUNT_HELP,
                    )?;
                    samples_found.extend(sample_names.clone());

                    sample_defs.push(SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        gem_group: None,
                        lanes: self.lanes.clone(),
                        library_type: Some(library_type),
                        r1_length: None,
                        r2_length: None,
                        read_path: path,
                        sample_indices: Some(vec!["any".to_string()]),
                        // NOTE: subsampling parameters could be added here
                        subsample_rate: None,
                        sample_names: Some(sample_names.into_iter().collect()),
                        fastq_id: None,
                    });
                }

                // We should have found all the requested samples in at least one sample def
                if let Some(ref samples) = self.sample {
                    let missing: Vec<_> = samples
                        .iter()
                        .filter(|&s| !samples_found.contains(s))
                        .collect();
                    if !missing.is_empty() {
                        bail!(
                            "The following sample(s) were not found in any of the \
                            fastq paths specified:\n{}",
                            missing.iter().join("\n")
                        );
                    }
                }

                Ok(sample_defs)
            }
        }
    };
}

#[derive(Parser, Debug, Clone)]
pub struct FastqArgs {
    /// Path to input FASTQ data
    #[clap(
        long,
        required_unless_present = "libraries",
        conflicts_with = "libraries",
        value_delimiter = ',',
        value_name = "PATH"
    )]
    pub fastqs: Vec<CliPath>,

    /// Name of the project folder within a mkfastq or
    /// bcl2fastq-generated folder from which to pick FASTQs
    #[clap(long, conflicts_with = "libraries", value_name = "TEXT")]
    pub project: Option<String>,

    /// Prefix of the filenames of FASTQs to select
    #[clap(
        long,
        conflicts_with = "libraries",
        value_delimiter = ',',
        value_name = "PREFIX",
        value_parser = NonEmptyStringValueParser::new(),
    )]
    pub sample: Option<Vec<String>>,

    /// Only use FASTQs from selected lanes
    #[clap(
        long,
        conflicts_with = "libraries",
        value_delimiter = ',',
        value_name = "NUMS"
    )]
    pub lanes: Option<Vec<usize>>,
}
impl_get_sample_defs!(FastqArgs);

// reimplement FastqArgs because we don't want
// "libraries" in cellranger-atac
#[derive(Parser, Debug, Clone)]
pub struct FastqArgsNoLibraries {
    /// Path to input FASTQ data
    #[clap(long, required = true, value_delimiter = ',', value_name = "PATH")]
    pub fastqs: Vec<CliPath>,

    /// Name of the project folder within a mkfastq or
    /// bcl2fastq-generated folder to pick FASTQs from.
    #[clap(
        long,
        value_name = "TEXT",
        value_parser = NonEmptyStringValueParser::new(),
    )]
    pub project: Option<String>,

    /// Prefix of the filenames of FASTQs to select.
    #[clap(
        long,
        value_delimiter = ',',
        value_name = "PREFIX",
        value_parser = NonEmptyStringValueParser::new(),
    )]
    pub sample: Option<Vec<String>>,

    /// Only use FASTQs from selected lanes
    #[clap(long, value_delimiter = ',', value_name = "NUMS")]
    pub lanes: Option<Vec<usize>>,
}
impl_get_sample_defs!(FastqArgsNoLibraries);

/// Convert an input FASTQ paths to a 'project' fastq paths.
pub fn get_fastq_paths(input_paths: &[CliPath], project: Option<&str>) -> Result<Vec<PathBuf>> {
    let outputs: Vec<_> = input_paths
        .iter()
        .filter_map(|path| get_output_path(path, project).transpose())
        .try_collect()?;
    match project {
        Some(project) if outputs.is_empty() => {
            bail!("Could not find any paths that matched --project value: {project}")
        }
        _ => Ok(outputs),
    }
}

///If the path contains a pipestance or bcl2fastq output in a
///known location, return the path to the bcl2fastq output
///root.
/// - path: The FASTQ input path.
/// - project: The supplied project field.
///
/// Returns the FASTQ path that is legal for a downstream run
/// Returns Err if the path/project path is illegal
fn get_bcl2fastq_output_folder(path: &Path) -> Option<PathBuf> {
    let link = path.join("outs/fastq_path");
    if link.exists() {
        return Some(link);
    }

    let link = path.join("fastq_path");
    if link.exists() {
        return Some(link);
    }

    if path.join("Reports").exists() || path.join("Stats").exists() {
        return Some(path.to_path_buf());
    }

    None
}

/// Get the bcl2fastq output dirs corresponding to bcl2fastq projects
fn get_projects(bcl2fastq_path: &Path) -> Result<Vec<String>> {
    let mut projects = Vec::new();

    for entry in std::fs::read_dir(bcl2fastq_path)? {
        let entry = entry?;
        let path = entry.path();
        // JAG-1314: added `Logs` to the deny list of project folders
        if path.is_dir()
            && path.file_name() != Some(OsStr::new("Reports"))
            && path.file_name() != Some(OsStr::new("Stats"))
            && path.file_name() != Some(OsStr::new("Logs"))
        {
            let project = path.file_name().unwrap().to_str().unwrap();
            projects.push(project.to_string());
        }
    }

    Ok(projects)
}

fn get_output_path(fastq_path: &Path, project: Option<&str>) -> Result<Option<PathBuf>> {
    let path = &expand_tilde(fastq_path)
        .ok_or_else(|| anyhow!("Can't expand '~' in fastq path. No home directory is known"))?;
    if !fastq_path.exists() {
        bail!("FASTQ path doesn't exist: {}", path.display());
    }

    let bcl_output_folder = get_bcl2fastq_output_folder(&path.canonicalize()?);
    if bcl_output_folder.is_none() {
        return Ok(Some(fastq_path.to_path_buf()));
    }

    let mut bcl_output_folder = bcl_output_folder.unwrap();

    let projects = get_projects(&bcl_output_folder)?;

    if projects.is_empty() {
        return Ok(Some(bcl_output_folder));
    }

    if let Some(p) = project {
        if projects.iter().any(|x| x == p) {
            bcl_output_folder.push(p);
            Ok(Some(bcl_output_folder))
        } else {
            Ok(None)
        }
    } else if projects.len() > 1 {
        bail!(
            "The --project argument must be specified if BCLs were demultiplexed into multiple project folders. Options:\n{}",
            projects.join("\n")
        )
    } else {
        assert!(projects.len() == 1);
        bcl_output_folder.push(&projects[0]);
        Ok(Some(bcl_output_folder))
    }
}

/// Expand a `~` to the current users home dir.
fn expand_tilde(path: &Path) -> Option<PathBuf> {
    if path.starts_with("~") {
        if path == Path::new("~") {
            dirs::home_dir()
        } else {
            dirs::home_dir().map(|mut h| {
                if h == Path::new("/") {
                    // Corner case: `h` root directory;
                    // don't prepend extra `/`, just drop the tilde.
                    path.strip_prefix("~").unwrap().to_path_buf()
                } else {
                    h.push(path.strip_prefix("~/").unwrap());
                    h
                }
            })
        }
    } else {
        Some(path.to_path_buf())
    }
}
