//! Common code for annotate subcommand of cellranger and spaceranger

use crate::utils::CliPath;
use serde::Serialize;

#[derive(Serialize)]
/// The MRO arguments for the annotate CS pipeline.
pub struct AnnotateMro {
    /// Sample ID
    pub sample_id: String,
    /// Sample description
    pub sample_desc: String,
    /// Path to the 10x Cloud token
    pub tenx_cloud_token_path: Option<String>,
    /// Path to the cell annotation model. Needed to run 10x Cloud annotation,
    /// but not for local annotation
    pub cell_annotation_model: Option<String>,
    /// Path to the sample cloupe file
    pub sample_cloupe: Option<CliPath>,
    /// Name of the CAS track
    pub cas_track_name: Option<String>,
    /// Path to the filtered matrix
    pub filtered_matrix: CliPath,
    /// Type of the pipestance
    pub pipestance_type: String,
    /// Whether to run aggregate and annotate or just annotate
    pub aggregate_and_annotate: bool,
}
