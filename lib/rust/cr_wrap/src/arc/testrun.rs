use crate::arc::count::CountMro;
use crate::env;
use crate::fastqs::FastqArgs;
use crate::mrp_args::MrpArgs;
use crate::utils::{validate_id, CliPath};
use anyhow::Result;
use clap::{self, Parser};
use cr_types::LibraryType;
use ordered_float::NotNan;
use std::path::PathBuf;

/// Executes the count pipeline on a tiny dataset that is bundled with the tarball
#[derive(Parser, Debug, Clone)]
pub struct TestrunArgs {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+ of maximum length
    /// 64 characters
    #[clap(long, value_name = "ID", required = true, value_parser = validate_id)]
    pub id: String,

    /// Sample description to embed in output files.
    #[clap(long, value_name = "TEXT")]
    description: Option<String>,

    #[clap(flatten)]
    pub mrp: MrpArgs,
}

impl TestrunArgs {
    pub fn to_mro_args(&self, pkg_env: &env::PkgEnv) -> Result<CountMro> {
        let t = self.clone();
        // determine paths relative to build
        let mut files_path: PathBuf = pkg_env.build_path.clone();
        files_path.push("external");
        files_path.push("arc_testrun_files");
        let mut reference_path = files_path.clone();
        reference_path.push("reference");
        let mut fastq_path = files_path;
        fastq_path.push("fastqs");

        // sample defs
        let mut sample_defs = Vec::new();
        sample_defs.extend(
            FastqArgs {
                fastqs: vec![CliPath::from(fastq_path.clone())],
                project: None,
                sample: Some(vec![String::from("tiny_arc_atac")]),
                lanes: None,
            }
            .get_sample_defs(LibraryType::Atac)?,
        );
        sample_defs.extend(
            FastqArgs {
                fastqs: vec![CliPath::from(fastq_path)],
                project: None,
                sample: Some(vec![String::from("tiny_arc_gex")]),
                lanes: None,
            }
            .get_sample_defs(LibraryType::Gex)?,
        );
        Ok(CountMro {
            sample_id: t.id.clone(),
            sample_desc: t.description.clone().unwrap_or(t.id),
            reference_path: CliPath::from(reference_path),
            sample_def: sample_defs,
            force_cells: None,
            check_library_compatibility: true,
            feature_linkage_max_dist_mb: NotNan::try_from(5.0).ok(),
            rna_include_introns: true,
            custom_peaks: None,
            peak_qval: None,
            k_means_max_clusters: None,
            no_bam: false,
        })
    }
}
