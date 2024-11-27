use crate::arc::types::{
    validate_distance, validate_strict_fraction, ForceCellsArgs, MinCounts, MAX_CLUSTERS_RANGE,
};
use crate::create_bam_arg::CreateBam;
use crate::fastqs::FastqArgs;
use crate::mrp_args::MrpArgs;
use crate::utils::{validate_id, CliPath};
use anyhow::{bail, Result};
use clap::{self, value_parser, Parser};
use cr_types::sample_def::SampleDef;
use cr_types::LibraryType;
use ordered_float::NotNan;
use serde::{self, Serialize};
use std::collections::HashMap;

/// A subcommand for controlling testing
#[derive(Parser, Debug, Clone)]
pub struct CountArgs {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+ of maximum length
    /// 64 characters
    #[clap(long, value_name = "ID", required = true, value_parser = validate_id)]
    pub id: String,

    /// Sample description to embed in output files.
    #[clap(long, value_name = "TEXT")]
    description: Option<String>,

    /// Path to folder containing cellranger-arc-compatible reference. Reference
    /// packages can be downloaded from support.10xgenomics.com or constructed using
    /// the `cellranger-arc mkref` command.
    #[clap(long, value_name = "PATH", required = true)]
    reference: CliPath,

    /// Path to 3-column CSV file defining the paths to ATAC and gene expression FASTQ
    /// data generated with the Chromium Single Cell Multiome ATAC + Gene Expression solution.
    /// A template CSV would look as follows (blank lines are ignored):
    ///
    /// fastqs,sample,library_type
    ///
    /// /data/HAWT7ADXX/outs/fastq_path,myATAC,Chromatin Accessibility
    ///
    /// /data/H35KCBCXY/outs/fastq_path,myGEX,Gene Expression
    ///
    #[clap(long, value_name = "CSV", required = true)]
    libraries: CliPath,

    #[clap(flatten)]
    force_cells: ForceCellsArgs,

    /// Override peak caller: specify peaks to use in downstream analyses from
    /// supplied 3-column BED file.
    /// The supplied peaks file must be sorted by position and not contain overlapping peaks;
    /// comment lines beginning with `#` are allowed.
    #[clap(long, value_name = "BED")]
    peaks: Option<CliPath>,

    /// Disable counting of intronic reads. In this mode, only reads that are
    /// exonic and compatible with annotated splice junctions in the reference
    /// are counted. Note: using this mode will reduce the UMI counts in the
    /// feature-barcode matrix.
    #[clap(long)]
    gex_exclude_introns: bool,

    /// Hidden option: Whether to check for barcode compatibility between libraries. [default: true]
    #[clap(
        long,
        hide = true,
        value_name = "true|false",
        conflicts_with = "skip_compatibility_check"
    )]
    check_library_compatibility: Option<bool>,

    /// Hidden option: Equivalent to --check-library-compatibility=false. Provided for backward compatibility.
    #[clap(long, hide = true, conflicts_with = "check_library_compatibility")]
    skip_compatibility_check: bool,

    /// Hidden option: override peak caller with custom q-value
    #[clap(long, hide = true, value_parser = validate_strict_fraction)]
    peak_qval: Option<NotNan<f64>>,

    /// Hidden option: change feature linkage window size
    #[clap(long, hide = true, value_parser = validate_distance)]
    feature_linkage_max_dist_mb: Option<NotNan<f64>>,

    /// Hidden option: change the max k of kmeans clustering
    #[clap(long, hide = true, value_parser = value_parser!(u64).range(MAX_CLUSTERS_RANGE))]
    k_means_max_clusters: Option<u64>,

    #[clap(flatten)]
    create_bam: CreateBam,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    pub dry: bool,

    #[clap(flatten)]
    pub mrp: MrpArgs,
}

impl CountArgs {
    pub fn to_mro_args(&self) -> Result<CountMro> {
        let c = self.clone();
        let mut sample_defs = Vec::new();

        // parse the libraries.csv & convert it into a set of SampleDefs.
        let libraries = cr_types::parse_legacy_libraries_csv(&c.libraries)?;

        let mut seen_atac = false;
        let mut seen_gex = false;

        for l in libraries {
            let fq = FastqArgs {
                fastqs: vec![CliPath::from(l.fastqs)],
                project: l.project,
                sample: Some(vec![l.sample]),
                lanes: None,
            };
            match l.library_type {
                LibraryType::Atac => seen_atac = true,
                LibraryType::Gex => seen_gex = true,
                v => {
                    bail!(
                        "Invalid value in library_type column {}: cellranger-arc \
                        count is only compatible with `{}` and \
                        `{}` libraries.",
                        v,
                        LibraryType::Gex,
                        LibraryType::Atac,
                    );
                }
            };
            sample_defs.extend(fq.get_sample_defs(l.library_type)?);
        }
        if !seen_atac {
            bail!("Invalid libraries file: missing `Chromatin Accessibility` FASTQ files");
        }
        if !seen_gex {
            bail!("Invalid libraries file: missing `Gene Expression` FASTQ files");
        }

        Ok(CountMro {
            sample_id: c.id.clone(),
            sample_desc: c.description.unwrap_or(c.id),
            reference_path: c.reference,
            sample_def: sample_defs,
            force_cells: c.force_cells.to_mro_arg()?,
            check_library_compatibility: c
                .check_library_compatibility
                .unwrap_or(!c.skip_compatibility_check),
            feature_linkage_max_dist_mb: c.feature_linkage_max_dist_mb,
            rna_include_introns: !c.gex_exclude_introns,
            custom_peaks: c.peaks,
            peak_qval: c.peak_qval,
            k_means_max_clusters: c.k_means_max_clusters,
            no_bam: !self.create_bam.validated()?,
        })
    }
}

#[derive(Serialize)]
pub struct CountMro {
    pub sample_id: String,
    pub sample_desc: String,
    pub reference_path: CliPath,
    pub sample_def: Vec<SampleDef>,
    pub force_cells: Option<HashMap<String, MinCounts>>,
    pub check_library_compatibility: bool,
    pub rna_include_introns: bool,
    pub custom_peaks: Option<CliPath>,
    pub peak_qval: Option<NotNan<f64>>,
    pub feature_linkage_max_dist_mb: Option<NotNan<f64>>,
    pub k_means_max_clusters: Option<u64>,
    pub no_bam: bool,
}
