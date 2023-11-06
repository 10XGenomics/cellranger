use crate::mrp_args::MrpArgs;
use crate::utils::{validate_id, CliPath};
use clap::{self, Parser};
use serde::Serialize;

/// Analyze targeted enrichment performance by comparing
/// a targeted sample to its cognate parent WTA sample.
#[derive(Parser, Debug, Clone, Serialize)]
pub struct TargetedCompare {
    /// A unique run id and output folder name [a-zA-Z0-9_-]+.
    #[clap(long = "id", display_order = 0, value_name = "ID", value_parser = validate_id, required = true)]
    pub sample_id: String,

    /// Sample description to embed in output files.
    #[clap(long = "description", default_value = "", value_name = "TEXT")]
    pub sample_desc: String,

    /// Path to the targeted molecule_info.h5 from a 'count --target-panel'
    /// analysis run (targeted gene expression run)
    #[clap(long = "targeted", value_name = "MOL_INFO_H5")]
    pub targeted_molecule_info: CliPath,

    /// Path to the parent molecule_info.h5 from a 'count'
    /// analysis run (parent unbiased gene expression run)
    #[clap(long = "parent", value_name = "MOL_INFO_H5")]
    pub parent_molecule_info: CliPath,

    /// A CSV file declaring the target gene panel used in the targeted
    /// experiment. Must be the same target panel CSV file specified
    /// in the 'count --target-panel' analysis run.
    #[clap(long = "target-panel", value_name = "CSV")]
    pub target_set: CliPath,

    /// Do not execute the pipeline.
    /// Generate a pipeline invocation (.mro) file and stop.
    #[serde(skip)]
    #[clap(long)]
    pub dry: bool,

    #[serde(skip)]
    #[clap(flatten)]
    pub mrp: MrpArgs,
}
