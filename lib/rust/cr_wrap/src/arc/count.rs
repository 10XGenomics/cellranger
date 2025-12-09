#![expect(missing_docs)]
use crate::arc::types::{
    ForceCellsArgs, MAX_CLUSTERS_RANGE, MinCounts, validate_distance, validate_strict_fraction,
};
use crate::create_bam_arg::CreateBam;
use crate::fastqs::FastqArgs;
use crate::mrp_args::MrpArgs;
use crate::utils::{CliPath, validate_id};
use anyhow::{Result, bail, ensure};
use clap::{self, Parser, value_parser};
use cr_types::LibraryType;
use cr_types::sample_def::SampleDef;
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

    /// The path to the 10x Cloud Analysis user token used to enable cell annotation. If not provided, will default to the location stored through cellranger-arc cloud auth setup
    #[clap(long, value_name = "PATH")]
    tenx_cloud_token_path: Option<CliPath>,

    /// Cell annotation model to use. Valid model names can be viewed by running `cellranger-arc cloud annotation models` or on the 10x Genomics Support site. If "auto", uses the default model for the species. If not provided, does not run cell annotation
    #[clap(long, value_name = "MODEL")]
    cell_annotation_model: Option<String>,

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

    /// Disable secondary analysis, e.g. clustering
    #[clap(long = "nosecondary")]
    no_secondary_analysis: bool,

    /// Trim the input Read 1 for GEX data to this length before analysis
    #[clap(long, value_name = "NUM")]
    rna_r1_length: Option<usize>,

    /// Trim the input Read 2 for GEX data to this length before analysis
    #[clap(long, value_name = "NUM")]
    rna_r2_length: Option<usize>,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    pub dry: bool,

    #[clap(flatten)]
    pub mrp: MrpArgs,
}

impl CountArgs {
    pub fn to_mro_args(&self) -> Result<CountMro> {
        use LibraryType::{Atac, GeneExpression};

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
                Atac => seen_atac = true,
                GeneExpression => seen_gex = true,
                library_type => {
                    bail!(
                        "Invalid value in library_type column {library_type}: cellranger-arc count \
                         is only compatible with `{GeneExpression}` and `{Atac}` libraries.",
                    );
                }
            };
            sample_defs.extend(fq.get_sample_defs(l.library_type)?);
        }

        ensure!(
            seen_atac,
            "Invalid libraries file: missing `Chromatin Accessibility` FASTQ files"
        );
        ensure!(
            seen_gex,
            "Invalid libraries file: missing `Gene Expression` FASTQ files"
        );

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
            no_secondary_analysis: c.no_secondary_analysis,
            rna_r1_length: c.rna_r1_length,
            rna_r2_length: c.rna_r2_length,
            tenx_cloud_token_path: c.tenx_cloud_token_path.clone(),
            cell_annotation_model: c.cell_annotation_model.clone(),
            skip_cell_annotation: c.cell_annotation_model.is_none(),
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
    pub no_secondary_analysis: bool,
    pub rna_r1_length: Option<usize>,
    pub rna_r2_length: Option<usize>,
    pub tenx_cloud_token_path: Option<CliPath>,
    pub cell_annotation_model: Option<String>,
    pub skip_cell_annotation: bool,
}
