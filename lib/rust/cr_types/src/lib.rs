//! cr_types
#![expect(missing_docs)]

use anyhow::{Context, Result, ensure};
use csv_parser::CsvParser;
use std::path::{Path, PathBuf};

pub mod aggr;
pub mod barcode_index;
pub mod cell_annotation;
pub mod chemistry;
pub mod clonotype;
pub mod constants;
pub mod csv_parser;
pub mod filtered_barcodes;
pub mod mempool;
mod metrics_file;
pub mod probe_set;
pub mod reference;
pub mod rna_read;
pub mod sample_def;
mod serde_helpers;
pub mod spill_vec;
pub mod target_panel_summary;
pub mod types;
pub mod utils;
pub mod websummary;

pub use barcode_index::*;
pub use fastq_set::ERROR_CODE_INFO;
pub use metrics_file::*;
pub use types::*;

#[derive(Debug)]
pub struct LegacyLibrariesEntry {
    pub fastqs: PathBuf,
    pub library_type: LibraryType,
    pub sample: String,
    pub project: Option<String>,
}

pub fn parse_legacy_libraries_csv(file: &Path) -> Result<Vec<LegacyLibrariesEntry>> {
    let required_fields = ["fastqs", "library_type", "sample"];

    let mut parser = CsvParser::new(file, required_fields, "Libraries CSV")
        .with_context(|| format!("Error opening --libraries file: {}", file.display()))?;

    let mut libraries: Vec<LegacyLibrariesEntry> = Vec::new();

    for i in 0..parser.len() {
        parser.set_line(i);

        let fastqs = parser.require_string("fastqs")?;
        let sample = parser.require_string("sample")?;
        let library_type: LibraryType = parser.parse_field("library_type", "valid library_type")?;
        // CELLRANGER-7889: this previously only parsed "auto VDJ", so maintain this invariant.
        if let Some(chain_type) = library_type.vdj_chain_type() {
            ensure!(
                chain_type == VdjChainType::Auto,
                "invalid library_type {library_type}"
            );
        }

        let project = parser.get_extra_data().get("project").cloned();

        let def = LegacyLibrariesEntry {
            fastqs: fastqs.into(),
            library_type,
            sample,
            project,
        };

        libraries.push(def);
    }

    Ok(libraries)
}

// initialize insta test harness
#[cfg(test)]
#[ctor::ctor]
fn init() {
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    unsafe { std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root) }
}
