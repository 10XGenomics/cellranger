//! cellranger
#![deny(missing_docs)]

use anyhow::{Context, Result, bail, ensure};
use clap::{self, Parser};
use cr_types::cell_annotation::CellAnnotationModel;
use cr_types::sample_def::SampleDef;
use cr_types::types::FileOrBytes;
use cr_types::{LibraryType, TargetingMethod};
use cr_wrap::chemistry_arg::validate_chemistry;
use cr_wrap::create_bam_arg::CreateBam;
use cr_wrap::env::{Initialize, get_version};
use cr_wrap::fastqs::{FastqArgs, FastqArgsNoLibraries};
use cr_wrap::mkref::Mkvdjref;
use cr_wrap::mrp_args::MrpArgs;
use cr_wrap::shared_cmd::{self, HiddenCmd, RnaSharedCmd, SharedCmd};
use cr_wrap::telemetry::CollectTelemetry;
use cr_wrap::utils::{AllArgs, CliPath, validate_id};
use cr_wrap::{cloud, env, execute, make_mro, make_mro_with_comment, mkfastq};
use multi::config::{ChemistryParam, ChemistrySet, MultiConfigCsv};
use serde::{self, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{File, read_to_string};
use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;
use std::str::FromStr;
const CMD: &str = "cellranger";

#[derive(Debug, Serialize, Clone, Copy)]
struct ForceCells(usize);

// NOTE: This is duplicated in cr_lib/src/parse_multi_config.rs and the value 10 is written out in
// the help text explicitly
const MIN_FORCE_CELLS: usize = 10;

impl FromStr for ForceCells {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<ForceCells> {
        let Ok(val) = s.parse() else {
            bail!("Expecting a positive integer here.");
        };
        if val < MIN_FORCE_CELLS {
            bail!(
                "Needs to be at least {MIN_FORCE_CELLS}. The value you have specified is {val} \
                 which is too low."
            );
        }
        Ok(ForceCells(val))
    }
}

/// Process 10x Genomics Gene Expression, Feature Barcode, and Immune Profiling data
#[derive(Parser, Debug)]
#[clap(name = CMD, version = get_version(), before_help = format!("{CMD} {}", get_version()))]
struct CellRanger {
    #[clap(subcommand)]
    subcmd: SubCommand,

    /// Provide a path to the environment definition json for the build
    #[clap(long, hide = true)]
    env_json: Option<PathBuf>,
}

impl Initialize for CellRanger {
    fn override_env_json(&self) -> Option<&std::path::Path> {
        self.env_json.as_deref()
    }
}

impl CollectTelemetry for CellRanger {
    fn should_collect_telemetry(&self) -> bool {
        !matches!(
            self.subcmd,
            SubCommand::Shared(SharedCmd::Telemetry(_)) | SubCommand::Hidden(HiddenCmd::OsCheck(_))
        )
    }
}

#[derive(Parser, Debug)]
enum SubCommand {
    /// Count gene expression and/or feature barcode reads from a single sample and GEM well.
    #[clap(name = "count")]
    Count(Count),

    /// Analyze multiplexed data or combined gene expression/immune profiling/feature barcode data.
    #[clap(name = "multi")]
    Multi(Multi),

    /// Output a multi config CSV template.
    #[clap(name = "multi-template")]
    MultiTemplate(MultiTemplate),

    /// Assembles single-cell VDJ receptor sequences from 10x Immune Profiling libraries.
    #[clap(name = "vdj")]
    Vdj(Vdj),

    /// Aggregate data from multiple Cell Ranger runs.
    #[clap(name = "aggr")]
    Aggr(Aggr),

    /// Annotate cell-types from outputs of a CellRanger run.
    #[clap(name = "annotate")]
    Annotate(Annotate),

    /// Re-run secondary analysis (dimensionality reduction, clustering, etc).
    #[clap(name = "reanalyze")]
    Reanalyze(Reanalyze),

    /// Prepare a reference for use with CellRanger VDJ.
    ///
    /// Build a Cell Ranger V(D)J-compatible reference folder from:
    /// 1. A user-supplied genome FASTA and gene GTF files.
    ///    For example, using files from ENSEMBL.
    ///
    /// OR
    ///
    /// 2. A FASTA file containing V(D)J segments as per the mkvdjref spec.
    ///    For example, using files from IMGT.
    ///
    /// Creates a new folder named after the genome.
    #[clap(name = "mkvdjref")]
    Mkvdjref(Mkvdjref),

    /// Run Illumina demultiplexer on sample sheets that contain 10x-specific
    /// sample index sets.
    #[clap(name = "mkfastq")]
    Mkfastq(AllArgs),

    /// Execute the 'count' pipeline on a small test dataset
    #[clap(name = "testrun")]
    Testrun(Testrun),

    /// Invoke cloud commands
    #[clap(name = "cloud")]
    Cloud(AllArgs),

    /// Auxiliary commands
    #[clap(flatten)]
    RnaShared(RnaSharedCmd),

    /// Commands to transmit metadata to 10x Genomics
    #[clap(flatten)]
    Shared(SharedCmd),

    /// Not user facing
    #[clap(flatten)]
    Hidden(HiddenCmd),
}

/// A subcommand for controlling testing
#[derive(Parser, Debug, Clone)]
struct Count {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long, value_name = "ID", value_parser = validate_id, required = true)]
    id: String,

    /// Sample description to embed in output files.
    #[clap(long, default_value = "", value_name = "TEXT")]
    description: String,

    /// Path of folder containing 10x-compatible transcriptome reference.
    #[clap(long, value_name = "PATH", required_unless_present = "feature_ref")]
    transcriptome: Option<CliPath>,

    // Arguments for specify fastq data -- should be shared across most pipelines.
    #[clap(flatten)]
    fastqs: FastqArgs,

    /// CSV file declaring input library data sources.
    #[clap(long, value_name = "CSV", required_unless_present = "fastqs")]
    libraries: Option<CliPath>,

    /// Feature reference CSV file, declaring Feature Barcode
    /// constructs and associated barcodes.
    #[clap(long, value_name = "CSV")]
    feature_ref: Option<CliPath>,

    /// CSV file specifying the probe set used, if any.
    #[clap(hide = true, long, value_name = "CSV")]
    probe_set: Option<String>, // don't want to do validation of path, just print error message if specified

    /// Expected number of recovered cells, used as input to cell calling algorithm.
    #[clap(long, value_name = "NUM")]
    expect_cells: Option<usize>,

    /// Force pipeline to use this number of cells, bypassing cell calling algorithm. [MINIMUM: 10]
    #[clap(long, value_name = "NUM")]
    force_cells: Option<ForceCells>,

    #[clap(flatten)]
    create_bam: CreateBam,

    /// Disable secondary analysis, e.g. clustering. Optional.
    #[clap(long = "nosecondary")]
    no_secondary_analysis: bool,

    ///  Hard trim the input Read 1 to this length before
    /// analysis.
    #[clap(long, value_name = "NUM")]
    r1_length: Option<usize>,

    ///  Hard trim the input Read 2 to this length before
    /// analysis.
    #[clap(long, value_name = "NUM")]
    r2_length: Option<usize>,

    /// Include intronic reads in count.
    #[clap(long, value_name = "true|false", default_value = "true")]
    include_introns: Option<bool>,

    /// Assay configuration. NOTE: by default the assay
    /// configuration is detected automatically, which is
    /// the recommended mode. You usually will not need to
    /// specify a chemistry. Options are: 'auto' for
    /// autodetection, 'threeprime' for Single Cell 3',
    /// 'fiveprime' for  Single Cell 5', 'SC3Pv1' or
    /// 'SC3Pv2' or 'SC3Pv3' or 'SC3Pv4' for Single Cell 3' v1/v2/v3/v4,
    /// 'SC3Pv3LT' for Single Cell 3' v3 LT,
    /// 'SC3Pv3HT' for Single Cell 3' v3 HT,
    /// 'SC5P-PE' or 'SC5P-PE-v3' or 'SC5P-R2' or 'SC5P-R2-v3' for Single Cell 5',
    /// paired-end/R2-only, 'SC-FB' for Single Cell
    /// Antibody-only 3' v2 or 5'. To analyze the GEX portion
    /// of multiome data, chemistry must be set to 'ARC-v1'.
    #[clap(long, value_name = "CHEM", default_value = "auto", value_parser=validate_chemistry)]
    chemistry: ChemistryParam,

    /// Proceed with processing using a --feature-ref but no
    /// Feature Barcode libraries specified with the
    /// 'libraries' flag.
    #[clap(long)]
    no_libraries: bool,

    /// Whether to check for barcode compatibility between libraries. [default: true]
    #[clap(long, value_name = "true|false")]
    check_library_compatibility: Option<bool>,

    /// The path to the 10x Cloud Analysis user token used to enable cell annotation.
    /// If not provided, will default to the location stored through cellranger cloud auth setup.
    #[clap(long = "tenx-cloud-token-path", value_name = "PATH")]
    tenx_cloud_token_path: Option<String>,

    /// Cell annotation model to use. Valid model names can be viewed by running
    /// `cellranger cloud annotation models` or on the 10x Genomics Support site
    /// (https://www.10xgenomics.com/support).
    /// If "auto", uses the default model for the species.
    /// If not provided, does not run cell annotation.
    #[clap(long = "cell-annotation-model", value_name = "MODEL")]
    cell_annotation_model: Option<CellAnnotationModel>,

    /// Enable or disable antibody and antigen aggregate filtering during cell calling.
    #[clap(long, value_name = "true|false", default_value = "true", hide = true)]
    filter_aggregates: Option<bool>,

    /// Minimum CRISPR UMI threshold
    #[clap(
        long = "min-crispr-umi",
        default_value = "3",
        value_name = "NUM",
        value_parser=clap::value_parser!(usize)
    )]
    min_crispr_umi_threshold: usize,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    dry: bool,

    #[clap(flatten)]
    mrp: MrpArgs,
}

impl Count {
    pub fn to_mro_args(&self) -> Result<CountCsMro> {
        let c = self.clone();
        let mut tenx_cloud_token_path = c.tenx_cloud_token_path;

        // probe-set argument allowed in PD but bail with error message for command-line
        ensure!(
            c.probe_set.is_none(),
            "--probe-set argument not supported by cellranger count. \
             Please use cellranger multi for all Flex samples."
        );

        let mut sample_defs = Vec::new();
        if let Some(ref libraries) = c.libraries {
            // parse the libraries.csv & convert it into a set of SampleDefs.
            let libraries = cr_types::parse_legacy_libraries_csv(libraries)?;

            for l in libraries {
                let fq = FastqArgs {
                    fastqs: vec![CliPath::from(l.fastqs)],
                    project: l.project,
                    sample: Some(vec![l.sample]),
                    lanes: None,
                };
                sample_defs.extend(fq.get_sample_defs(l.library_type)?);
            }
        } else {
            // convert the fastq cmd-line args to a SampleDef
            sample_defs.extend(c.fastqs.get_sample_defs(LibraryType::Gex)?);
        }

        // if the token path is not supplied, default to the value expected by
        // cellranger cloud auth setup, which is $HOME/.config/txg/credential
        // cellranger cloud [args], which is $HOME/.config/txg/credentials
        if tenx_cloud_token_path.is_none() && c.cell_annotation_model.is_some() {
            if let Ok(default_token_path) = cloud_utils::default_token_path() {
                tenx_cloud_token_path = Some(default_token_path);
            } else {
                eprintln!("{}", cloud_utils::CELL_ANNOTATION_HOMEDIR_MSG);
            }
        }

        Ok(CountCsMro {
            sample_id: c.id,
            sample_def: sample_defs,
            target_set_name: None,
            target_set: None,
            sample_desc: c.description,
            reference_path: c.transcriptome,
            recovered_cells: c.expect_cells,
            no_bam: !c.create_bam.validated()?,
            no_secondary_analysis: c.no_secondary_analysis,
            force_cells: c.force_cells.map(|fc| fc.0),
            chemistry: c.chemistry,
            r1_length: c.r1_length,
            r2_length: c.r2_length,
            targeting_method: None,
            aligner: None,
            trim_polya_min_score: Some(20),
            trim_tso_min_score: Some(20),
            feature_reference: c.feature_ref,
            include_introns: c.include_introns.unwrap(),
            check_library_compatibility: c.check_library_compatibility.unwrap_or(true),
            disable_ab_aggregate_detection: !c.filter_aggregates.unwrap_or(true),
            min_crispr_umi_threshold: c.min_crispr_umi_threshold,
            tenx_cloud_token_path,
            skip_cell_annotation: c.cell_annotation_model.is_none(),
            cell_annotation_model: c
                .cell_annotation_model
                .as_ref()
                .and_then(CellAnnotationModel::to_pipeline_inputs),
        })
    }
}
#[derive(Serialize)]
struct CountCsMro {
    sample_id: String,
    sample_def: Vec<SampleDef>,
    target_set: Option<CliPath>,
    target_set_name: Option<String>,
    sample_desc: String,
    reference_path: Option<CliPath>,
    recovered_cells: Option<usize>,
    no_bam: bool,
    no_secondary_analysis: bool,
    force_cells: Option<usize>,
    chemistry: ChemistryParam,
    r1_length: Option<usize>,
    r2_length: Option<usize>,
    targeting_method: Option<TargetingMethod>,
    aligner: Option<String>,
    trim_polya_min_score: Option<usize>,
    trim_tso_min_score: Option<usize>,
    feature_reference: Option<CliPath>,
    include_introns: bool,
    check_library_compatibility: bool,
    disable_ab_aggregate_detection: bool,
    min_crispr_umi_threshold: usize,
    tenx_cloud_token_path: Option<String>,
    cell_annotation_model: Option<String>,
    skip_cell_annotation: bool,
}

/// A subcommand for controlling testing
#[derive(Parser, Debug, Clone)]
struct Multi {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long, value_name = "ID", value_parser = validate_id, required = true)]
    id: String,

    /// Sample description to embed in output files.
    #[clap(long, default_value = "", value_name = "TEXT")]
    description: String,

    /// Path of CSV file enumerating input libraries and analysis parameters.
    #[clap(long, value_name = "CSV")]
    csv: CliPath,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    dry: bool,

    #[clap(flatten)]
    mrp: MrpArgs,
}
#[derive(Serialize)]
struct MultiCsMro {
    sample_id: String,
    sample_desc: String,
    config: FileOrBytes,
    config_hash: String,
    no_preflight: bool,
}

impl Multi {
    pub fn to_mro_args(&self) -> Result<MultiCsMro> {
        let config_hash = {
            let data = std::fs::read(&self.csv).with_context(|| self.csv.to_string())?;
            let mut hasher = Sha256::new();
            hasher.update(&data);
            hex::encode(hasher.finalize())
        };
        Ok(MultiCsMro {
            sample_id: self.id.clone(),
            sample_desc: self.description.clone(),
            config: FileOrBytes {
                file: Some(self.csv.clone().into()),
                bytes: None,
            },
            config_hash,
            no_preflight: self.mrp.nopreflight,
        })
    }

    fn validate_config(&self) -> Result<()> {
        let _config = MultiConfigCsv::from_csv(self.csv.clone())?;
        Ok(())
    }
}

#[derive(Parser, Debug, Clone)]
struct MultiTemplate {
    /// Optional output file path.
    #[clap(short, long, value_name = "CSV")]
    output: Option<String>,

    /// Print multi config parameter options and descriptions
    #[clap(short, long)]
    parameters: bool,
}

impl MultiTemplate {
    fn print(&self) -> Result<ExitCode> {
        let message = if self.parameters {
            include_str!("parameter-descriptions.txt")
        } else {
            include_str!("multi-config-template.csv")
        };
        if let Some(filename) = &self.output {
            write!(File::create(filename)?, "{message}")?;
        } else {
            print!("{message}");
        }
        Ok(ExitCode::SUCCESS)
    }
}

/// Aggregates the feature/cell count data
/// generated from multiple runs of the 'cellranger count' pipeline.
// To run this pipeline, supply a CSV that enumerates the paths to the
// molecule_info.h5 files produced by 'cellranger count'.
#[derive(Parser, Debug, Clone, Serialize)]
struct Aggr {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long = "id", value_name = "ID", value_parser = validate_id, required = true)]
    sample_id: String,

    /// Sample description to embed in output files.
    #[clap(long = "description", default_value = "", value_name = "TEXT")]
    sample_desc: String,

    /// Path of CSV file enumerating 'cellranger count/vdj/multi' outputs.
    #[clap(long = "csv", value_name = "CSV")]
    aggregation_csv: CliPath,

    /// Library depth normalization mode.
    #[clap(
        long = "normalize",
        default_value = "mapped",
        value_name = "MODE",
        value_parser = ["mapped", "none"],
    )]
    normalization_mode: String,

    /// Disable secondary analysis, e.g. clustering.
    #[clap(long = "nosecondary")]
    no_secondary_analysis: bool,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[serde(skip)]
    #[clap(long)]
    dry: bool,

    /// Minimum CRISPR UMI threshold
    #[clap(
        long = "min-crispr-umi",
        default_value = "3",
        value_name = "NUM",
        value_parser=clap::value_parser!(usize)
    )]
    min_crispr_umi_threshold: usize,

    /// Enable t-SNE projection
    #[clap(
        long = "enable-tsne",
        default_value = "false",
        value_name = "true|false"
    )]
    enable_tsne: Option<bool>,

    // not a cmd-line arg -- should be filled in with the working dir
    #[clap(hide = true, default_value = ".")]
    pipestance_root: PathBuf,

    #[serde(skip)]
    #[clap(flatten)]
    mrp: MrpArgs,
}

/// Performs cell type annotation on the feature-barcode matrix
/// generated from a run of the Cell Ranger 'count' or the
/// 'multi' pipeline. To run this pipeline, filtered_matrix.h5
/// files produced by 'count' or 'multi', a cell annotation token,
/// and a cloupe file (optional) are required.
#[derive(Parser, Debug, Clone, Serialize)]
struct Annotate {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long = "id", value_name = "ID", value_parser = validate_id, required = true)]
    sample_id: String,

    /// Sample description to embed in output files.
    #[clap(long = "description", default_value = "", value_name = "TEXT")]
    sample_desc: String,

    /// Path to 10x Cloud cell annotation token file.
    #[clap(long = "tenx-cloud-token-path", value_name = "PATH")]
    tenx_cloud_token_path: Option<String>,

    /// Cell annotation model to use.
    #[clap(long = "cell-annotation-model")]
    cell_annotation_model: Option<CellAnnotationModel>,

    /// Cloupe file to use.
    #[clap(long = "cloupe", value_name = "CLOUPE")]
    sample_cloupe: Option<CliPath>,

    /// Name of track to add onto the cloupe file.
    #[clap(
        long = "cloupe-group-name",
        requires = "sample_cloupe",
        value_name = "TEXT"
    )]
    cas_track_name: Option<String>,

    /// A feature-barcode matrix containing data for one genome.
    /// Should be the filtered version.
    #[clap(long = "matrix", value_name = "MATRIX_H5", required = true)]
    filtered_matrix: CliPath,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[serde(skip)]
    #[clap(long)]
    dry: bool,

    #[serde(skip)]
    #[clap(flatten)]
    mrp: MrpArgs,
}

#[derive(Serialize)]
struct AnnotateMro {
    sample_id: String,
    sample_desc: String,
    tenx_cloud_token_path: Option<String>,
    cell_annotation_model: Option<String>,
    sample_cloupe: Option<CliPath>,
    cas_track_name: Option<String>,
    filtered_matrix: CliPath,
}

impl Annotate {
    // set the default token path if necessary
    fn with_default_cas_token(&self) -> Result<AnnotateMro> {
        let mut a = self.clone();
        if self.tenx_cloud_token_path.is_none() {
            if let Ok(tenx_cloud_token_path) = cloud_utils::default_token_path() {
                a.tenx_cloud_token_path = Some(tenx_cloud_token_path);
            } else {
                eprintln!("{}", cloud_utils::CELL_ANNOTATION_HOMEDIR_MSG);
            }
        }
        Ok(AnnotateMro {
            sample_id: a.sample_id,
            sample_desc: a.sample_desc,
            tenx_cloud_token_path: a.tenx_cloud_token_path,
            cell_annotation_model: a
                .cell_annotation_model
                .as_ref()
                .and_then(CellAnnotationModel::to_pipeline_inputs),
            sample_cloupe: a.sample_cloupe,
            cas_track_name: a.cas_track_name,
            filtered_matrix: a.filtered_matrix,
        })
    }
}

#[derive(Parser, Debug, Clone, Serialize)]
struct Reanalyze {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long = "id", value_name = "ID", value_parser = validate_id, required = true)]
    sample_id: String,

    /// Sample description to embed in output files.
    #[clap(long = "description", default_value = "", value_name = "TEXT")]
    sample_desc: String,

    /// A feature-barcode matrix containing data for one genome.
    /// Should be the filtered version, unless using --force-cells.
    #[clap(long = "matrix", value_name = "MATRIX_H5")]
    filtered_matrices_h5: CliPath,

    /// A CSV file specifying analysis parameters. Optional.
    #[clap(long = "params", value_name = "PARAMS_CSV")]
    params_csv: Option<CliPath>,

    /// A CSV file containing a list of cell barcodes to use for
    /// reanalysis, e.g. barcodes exported from Loupe Browser. Optional.
    #[clap(long = "barcodes", value_name = "BARCODES_CSV")]
    barcodes_csv: Option<CliPath>,

    /// A CSV file containing a list of feature IDs to use for
    /// reanalysis. For gene expression, this should
    /// correspond to the gene_id field in the reference GTF
    /// should be \(e.g. ENSG... for ENSEMBL-based
    /// references\). Optional.
    #[clap(long = "genes", value_name = "GENES_CSV")]
    genes_csv: Option<CliPath>,

    /// A CSV file containing a list of feature IDs to exclude
    /// from reanalysis. For gene expression, this should
    /// correspond to the gene_id field in the reference GTF
    /// \(e.g., ENSG... for ENSEMBL-based references\). The
    /// exclusion is applied after --genes. Optional.
    #[clap(long = "exclude-genes", value_name = "GENES_CSV")]
    exclude_genes_csv: Option<CliPath>,

    /// If the input matrix was produced by 'aggr',
    /// you may pass the same aggregation CSV in order to
    /// retain per-library tag information in the resulting
    /// .cloupe file.  This argument is required to enable
    /// chemistry batch correction. Optional.
    #[clap(long = "agg", value_name = "AGGREGATION_CSV")]
    aggregation_csv: Option<CliPath>,

    /// Force pipeline to use this number of cells, bypassing cell calling algorithm. [MINIMUM: 10]
    #[clap(long = "force-cells", value_name = "NUM")]
    force_cells: Option<ForceCells>,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[serde(skip)]
    #[clap(long)]
    dry: bool,

    #[serde(skip)]
    #[clap(flatten)]
    mrp: MrpArgs,
}

// TODO(CELLRANGER-7889) collapse this with the other VDJ chain type types
#[derive(Serialize, Debug, Clone, Copy)]
enum VdjChainType {
    TR,
    IG,
    #[serde(rename = "auto")]
    Auto,
}

impl FromStr for VdjChainType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<VdjChainType> {
        Ok(match s {
            "TR" => VdjChainType::TR,
            "IG" => VdjChainType::IG,
            "auto" => VdjChainType::Auto,
            "TR_GD" => bail!("Please use Cellranger multi to run Gamma/Delta TCR datasets."),
            _ => bail!("unknown chain type \"{s}\""),
        })
    }
}

/// A subcommand for controlling testing
#[derive(Serialize, Parser, Debug)]
struct Vdj {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long = "id", value_name = "ID", value_parser = validate_id, required = true)]
    sample_id: String,

    /// Sample description to embed in output files.
    #[clap(long = "description", default_value = "", value_name = "TEXT")]
    sample_desc: String,

    /// Path of folder containing 10x-compatible VDJ reference.
    /// Optional if '--denovo' is specified.
    #[clap(long = "reference", value_name = "PATH")]
    vdj_reference_path: Option<CliPath>,

    // Arguments for specify fastq data -- should be shared across most pipelines.
    #[clap(flatten)]
    #[serde(skip)]
    fastqs: FastqArgsNoLibraries,

    // Sample defs -- this will be filled in based on fastqs field above
    #[clap(skip)]
    sample_def: Vec<SampleDef>,

    /// Run in reference-free mode (do not use annotations).
    #[clap(long)]
    denovo: bool,

    /// Disable clonotyping.
    #[clap(long, hide = true)]
    skip_clonotyping: bool,

    /// Chain type to display metrics for: 'TR' for T cell receptors,
    /// 'IG' for B cell receptors, or 'auto' to autodetect.
    #[clap(long = "chain", default_value = "auto", value_name = "CHAIN_SPEC")]
    chain_type: VdjChainType,

    /// If inner enrichment primers other than those provided in
    /// the 10x kits are used, they need to be specified
    /// here as a textfile with one primer per line.
    /// Disable secondary analysis, e.g. clustering.
    #[clap(long, value_name = "PATH")]
    inner_enrichment_primers: Option<CliPath>,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    dry: bool,

    #[serde(skip)]
    #[clap(flatten)]
    mrp: MrpArgs,
}

#[derive(Parser, Debug, Clone)]
struct Testrun {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long, value_name = "ID", required = true, value_parser = validate_id)]
    id: String,

    /// Sample description to embed in output files.
    #[clap(long, value_name = "TEXT")]
    description: Option<String>,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[clap(long)]
    dry: bool,

    #[clap(flatten)]
    mrp: MrpArgs,
}

impl Testrun {
    fn to_mro_args(&self, pkg_env: &env::PkgEnv) -> Result<CountCsMro> {
        // determine paths relative to build
        let files_path = pkg_env.build_path.join("external");
        let fastqs = vec![CliPath::from(files_path.join("cellranger_tiny_fastq"))];
        let reference_path = Some(CliPath::from(files_path.join("cellranger_tiny_ref")));

        // sample defs
        let mut sample_def = Vec::new();
        sample_def.extend(
            FastqArgs {
                fastqs,
                project: None,
                sample: Some(vec![String::from("tinygex")]),
                lanes: None,
            }
            .get_sample_defs(LibraryType::Gex)?,
        );

        let t = self.clone();
        Ok(CountCsMro {
            sample_id: t.id,
            sample_def,
            target_set: None,
            target_set_name: None,
            sample_desc: t.description.unwrap_or_default(),
            reference_path,
            recovered_cells: None,
            no_bam: false,
            no_secondary_analysis: false,
            force_cells: None,
            chemistry: ChemistryParam::Set(ChemistrySet::ThreePrimeV3),
            r1_length: None,
            r2_length: None,
            targeting_method: None,
            aligner: None,
            trim_polya_min_score: None,
            trim_tso_min_score: None,
            feature_reference: None,
            include_introns: false,
            check_library_compatibility: true,
            disable_ab_aggregate_detection: false,
            min_crispr_umi_threshold: 3,
            tenx_cloud_token_path: None,
            cell_annotation_model: None,
            skip_cell_annotation: true,
        })
    }
}

fn inner_main() -> Result<ExitCode> {
    let (opts, pkg_env, mut telemetry) = CellRanger::initialize()?;

    match opts.subcmd {
        SubCommand::Count(c) => {
            let has_gex_library = if !c.fastqs.fastqs.is_empty() {
                true
            } else if let Some(libraries_csv) = &c.libraries {
                cr_types::parse_legacy_libraries_csv(libraries_csv)?
                    .iter()
                    .any(|x| x.library_type == LibraryType::GeneExpression)
            } else {
                // CLAP ensures that either --fastqs or --libraries is present.
                unreachable!();
            };
            if has_gex_library {
                ensure!(
                    c.transcriptome.is_some(),
                    "A reference transcriptome is required to analyze a Gene Expression library."
                );
            }

            // If --feature-ref is provided, require --libraries or --no-libraries.
            //
            // --no-libraries historically was used to allow a feature reference to be added to
            // pre-feature-barcoding GEX data in order to allow aggr of that data w/ GEX+FB
            // libraries by specifying an (unused, but compatible) feature reference for the
            // GEX-only run.
            //
            if c.feature_ref.is_some() && c.libraries.is_none() && !c.no_libraries {
                bail!(
"You specified --feature-ref, but not --libraries. Did you mean to input feature barcode libraries?
If you have 1 or more feature barcode libraries:
    Use --libraries to specify your input FASTQs and the associated library types.
If you want to proceed with a feature barcode reference, but no feature barcode data:
    Add the --no-libraries flag.");
            }
            let mro = if let Some(libraries) = &c.libraries {
                let libraries_lines =
                    read_to_string(libraries).with_context(|| libraries.to_string())?;
                make_mro_with_comment(
                    "SC_RNA_COUNTER_CS",
                    &c.to_mro_args()?,
                    "rna/sc_rna_counter_cs.mro",
                    &libraries_lines,
                )?
            } else {
                make_mro(
                    "SC_RNA_COUNTER_CS",
                    &c.to_mro_args()?,
                    "rna/sc_rna_counter_cs.mro",
                )?
            };

            execute(&c.id, &mro, &c.mrp, c.dry, &mut telemetry)
        }

        SubCommand::Multi(m) => {
            m.validate_config()?;
            let mro = make_mro_with_comment(
                "SC_MULTI_CS",
                &m.to_mro_args()?,
                "rna/sc_multi_cs.mro",
                &read_to_string(&m.csv).with_context(|| m.csv.to_string())?,
            )?;
            execute(&m.id, &mro, &m.mrp, m.dry, &mut telemetry)
        }

        SubCommand::MultiTemplate(mt) => mt.print(),

        SubCommand::Aggr(mut aggr) => {
            // Custom validation

            // fill in pipestance_root
            aggr.pipestance_root = std::env::current_dir()?;
            let mro = make_mro_with_comment(
                "SC_RNA_AGGREGATOR_CS",
                &aggr,
                "rna/sc_rna_aggregator_cs.mro",
                &read_to_string(&aggr.aggregation_csv)
                    .with_context(|| aggr.aggregation_csv.to_string())?,
            )?;
            execute(&aggr.sample_id, &mro, &aggr.mrp, aggr.dry, &mut telemetry)
        }

        SubCommand::Annotate(annotate) => {
            let acopy = annotate.with_default_cas_token()?;
            let mro = make_mro(
                "CELLRANGER_ANNOTATE_CS",
                &acopy,
                "rna/cell_annotation_service.mro",
            )?;
            execute(
                &annotate.sample_id,
                &mro,
                &annotate.mrp,
                annotate.dry,
                &mut telemetry,
            )
        }

        SubCommand::Reanalyze(ra) => {
            // Custom validation

            let mro = make_mro("SC_RNA_REANALYZER_CS", &ra, "rna/sc_rna_reanalyzer_cs.mro")?;
            execute(&ra.sample_id, &mro, &ra.mrp, ra.dry, &mut telemetry)
        }

        SubCommand::Vdj(mut vdj) => {
            // convert fastq args to sample defs
            vdj.sample_def = vdj.fastqs.get_sample_defs(LibraryType::VdjAuto)?;

            if let Some(ref ref_folder) = vdj.vdj_reference_path {
                vdj_reference::VdjReference::check(ref_folder).with_context(|| {
                    format!(
                        "Failed while verifying the contents of the V(D)J reference at \"{}\"",
                        ref_folder.as_ref().display()
                    )
                })?;
            }

            let mro = make_mro("SC_VDJ_ASSEMBLER_CS", &vdj, "rna/sc_vdj_assembler_cs.mro")?;
            execute(&vdj.sample_id, &mro, &vdj.mrp, vdj.dry, &mut telemetry)
        }

        SubCommand::Mkfastq(m) => mkfastq::run_mkfastq(&m, "_cellranger_internal"),
        SubCommand::Mkvdjref(mut args) => {
            args.shared.populate_version(&pkg_env);
            args.execute(&mut telemetry)
        }
        SubCommand::Testrun(t) => {
            let mro_args = t.to_mro_args(&pkg_env)?;
            let mro = make_mro("SC_RNA_COUNTER_CS", &mro_args, "rna/sc_rna_counter_cs.mro")?;
            execute(&t.id, &mro, &t.mrp, t.dry, &mut telemetry)
        }

        SubCommand::Cloud(args) => cloud::run_cloud(&pkg_env, &args),
        SubCommand::RnaShared(args) => shared_cmd::run_rna_shared(&pkg_env, args, &mut telemetry),
        SubCommand::Shared(args) => shared_cmd::run_shared(&pkg_env, args),
        SubCommand::Hidden(args) => shared_cmd::run_hidden(args),
    }
}

fn main() -> ExitCode {
    match inner_main() {
        Ok(exit_code) => exit_code,
        Err(err) => {
            cr_wrap::utils::print_error_chain(&err);
            ExitCode::FAILURE
        }
    }
}
