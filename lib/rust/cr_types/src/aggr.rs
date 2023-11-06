//! Types used in aggr pipelines

use martian_derive::MartianType;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Pointer to one logical set of FASTQ from a unique (library, gem_group) tuple
#[derive(Serialize, Deserialize, Clone, PartialOrd, PartialEq, Eq, Debug, MartianType)]
pub struct SampleDef {
    pub library_id: String,
    pub molecule_h5: PathBuf,
}

#[derive(Serialize, Deserialize, Clone, PartialOrd, PartialEq, Debug, MartianType, Ord, Eq)]
pub struct LibraryInfo {
    pub aggr_id: String,
    pub batch_id: u16,
    pub batch_name: String,
    pub gem_group: u16,
    pub library_id: u16,
    pub library_type: Option<String>,
    pub old_gem_group: u16,
    pub old_library_index: u16,
    pub target_set_name: Option<String>,
}
