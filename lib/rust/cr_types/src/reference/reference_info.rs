#![deny(missing_docs)]
use super::probe_set_reference::TargetSetFile;
use crate::GenomeName;
use crate::probe_set::ProbeSetReferenceMetadata;
use anyhow::{Context, Result, ensure};
use itertools::Itertools;
use martian::MartianFileType;
use martian_derive::MartianStruct;
use martian_filetypes::FileTypeRead;
use martian_filetypes::json_file::JsonFile;
use metric::{JsonReport, JsonReporter};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha1::{Digest, Sha1};
use std::fs::File;
use std::path::{Path, PathBuf};

/// Use this string to join multiple genome names in a multi-genome reference.
pub const MULTI_GENOME_SEPARATOR: &str = "_and_";

/// Metadata associated with the reference from the `reference.json` file in the reference directory
/// This is primarily used by a number of functions to query the `genomes` present in the reference.
///
/// NOTE: We use `#[serde(default)]` for certain fields because some older versions of references,
/// for example `hg19-1.0.0`, did not populate these fields in the `reference.json`.
#[derive(Clone, Default, Deserialize)]
struct TranscriptomeReferenceJson {
    #[serde(default)]
    pub fasta_hash: String,
    #[serde(default)]
    pub gtf_hash: String,
    #[serde(default, rename = "gtf_hash.gz")]
    pub gtf_hash_gz: String,
    #[serde(default)]
    pub input_fasta_files: Vec<String>,
    #[serde(default)]
    pub input_gtf_files: Vec<String>,
    pub mem_gb: f64,
    #[serde(default)]
    pub mkref_version: String,
    #[serde(default)]
    #[allow(dead_code)]
    pub path: PathBuf,
    pub version: Option<String>,
    pub genomes: Vec<GenomeName>,
}

#[derive(Clone, Default, Serialize, Deserialize, MartianStruct, Debug)]
/// Stores information that is extracted from the reference.json file
/// but not extractable from the probe set CSV file.
pub struct TranscriptomeFileInfo {
    /// SHA1 hash of the reference FASTA file.
    pub fasta_hash: String,
    /// SHA1 hash of the reference GTF file.
    pub gtf_hash: String,
    /// SHA1 hash of the gzipped reference GTF file.
    pub gtf_hash_gz: String,
    /// List of input FASTA files used to build the reference.
    pub input_fasta_files: Vec<String>,
    /// List of input GTF files used to build the reference.
    pub input_gtf_files: Vec<String>,
    /// Memory requirement for the reference.
    pub mem_gb: f64,
    /// mkref version used to build the reference.
    pub mkref_version: String,
    /// Path to the reference directory.
    pub reference_path: PathBuf,
}

/// Reference information extracted from a reference.json (if reference_path
/// is provided) or from a probe set CSV file (if probe_set_csv is provided).
#[derive(Clone, Default, Serialize, Deserialize, MartianStruct, Debug)]
pub struct ReferenceInfo {
    /// List of genomes in the reference.
    pub genomes: Vec<GenomeName>,
    /// Reference version.
    pub version: Option<String>,
    /// Transcriptome files information extracted from <reference_path>/reference.json
    pub transcriptome_info: Option<TranscriptomeFileInfo>,
}

impl ReferenceInfo {
    /// Create a `ReferenceInfo` from a reference.json file
    /// that sits inside reference_path directory.
    pub fn from_reference_path(reference_path: &Path) -> Result<ReferenceInfo> {
        let transcriptome_json: TranscriptomeReferenceJson =
            JsonFile::new(reference_path, "reference.json").read()?;
        Ok(ReferenceInfo {
            genomes: transcriptome_json.genomes,
            version: transcriptome_json.version,
            transcriptome_info: Some(TranscriptomeFileInfo {
                reference_path: reference_path.to_path_buf(),
                fasta_hash: transcriptome_json.fasta_hash,
                gtf_hash: transcriptome_json.gtf_hash,
                gtf_hash_gz: transcriptome_json.gtf_hash_gz,
                input_fasta_files: transcriptome_json.input_fasta_files,
                input_gtf_files: transcriptome_json.input_gtf_files,
                mem_gb: transcriptome_json.mem_gb,
                mkref_version: transcriptome_json.mkref_version,
            }),
        })
    }

    /// Create a `ReferenceInfo` from a probe set CSV file.
    pub fn from_probe_set_csv(probe_set: &TargetSetFile) -> Result<ReferenceInfo> {
        let metadata = ProbeSetReferenceMetadata::load_from(probe_set)?;
        Ok(ReferenceInfo {
            genomes: metadata.reference_genomes().into_iter().sorted().collect(),
            version: Some(metadata.reference_version().to_string()),
            transcriptome_info: None,
        })
    }

    /// Validate that associated data complies with expectations.
    /// Check that FASTA and GTF files exist.
    /// Check that FASTA and GTF files match checksums, if present.
    pub fn validate(&self) -> Result<()> {
        // Checksum the FASTA file.
        let Some(file_info) = &self.transcriptome_info else {
            return Ok(());
        };
        if !file_info.fasta_hash.is_empty() {
            sha1_checksum_file(
                &file_info.reference_path.join("fasta/genome.fa"),
                &file_info.fasta_hash,
            )
            .context("failed to validate FASTA file for reference")?;
        }
        // Checksum the GTF file.
        if !file_info.gtf_hash.is_empty() {
            sha1_checksum_file(
                &file_info.reference_path.join("genes/genes.gtf"),
                &file_info.gtf_hash,
            )
            .context("failed to validate GTF file for reference")?;
        }
        Ok(())
    }

    /// Get the path to the reference directory.
    pub fn get_reference_path(&self) -> Option<&Path> {
        self.transcriptome_info
            .as_ref()
            .map(|info| info.reference_path.as_path())
    }
}

impl JsonReport for ReferenceInfo {
    /// Convert the `ReferenceInfo` into a JSON object.
    fn to_json_reporter(&self) -> JsonReporter {
        // Mandatory fields.
        let mut key_values = Vec::from([
            (
                "reference_genomes",
                self.genomes.iter().join(MULTI_GENOME_SEPARATOR),
            ),
            (
                "reference_version",
                self.version.clone().unwrap_or_default(),
            ),
        ]);

        // Optional fields.
        if let Some(file_info) = self.transcriptome_info.clone() {
            key_values.extend(Vec::from([
                ("reference_fasta_hash", file_info.fasta_hash),
                ("reference_gtf_hash", file_info.gtf_hash),
                ("reference_gtf_hash.gz", file_info.gtf_hash_gz),
                (
                    "reference_input_fasta_files",
                    file_info.input_fasta_files.join(", "),
                ),
                (
                    "reference_input_gtf_files",
                    file_info.input_gtf_files.join(", "),
                ),
                ("reference_mkref_version", file_info.mkref_version),
                ("reference_type", "Transcriptome".to_string()),
                (
                    "reference_path",
                    file_info.reference_path.display().to_string(),
                ),
            ]));
        }
        key_values
            .into_iter()
            .map(|(k, v)| (k.to_string(), Value::from(v)))
            .collect()
    }
}

/// Compute the SHA1 checksum of the file at the provided path.
/// Return an error if the computed checksum doesn't match the expected value.
fn sha1_checksum_file(path: &Path, expected: &str) -> Result<()> {
    let mut hasher = Sha1::new();
    std::io::copy(
        &mut File::open(path).with_context(|| path.display().to_string())?,
        &mut hasher,
    )?;
    let computed = format!("{:x}", hasher.finalize());
    ensure!(
        expected == computed,
        ChecksumMismatch {
            expected: expected.to_string(),
            computed,
            path: path.to_owned(),
        }
    );
    Ok(())
}
#[derive(Debug, thiserror::Error)]
#[error("checksum \"{computed}\" of file \"{}\" does not match expected value \"{expected}\"", .path.display())]
struct ChecksumMismatch {
    expected: String,
    computed: String,
    path: PathBuf,
}

#[cfg(test)]
mod test {
    use super::*;
    use insta::assert_json_snapshot;

    #[test]
    fn test_reference_info_single_report() -> Result<()> {
        let ref_info =
            ReferenceInfo::from_reference_path(Path::new("test/reference/GRCh38_ref_tiny"))?;
        assert_json_snapshot!(ref_info.to_json_reporter());
        Ok(())
    }

    #[test]
    fn test_reference_info_multi_report() -> Result<()> {
        let ref_info = ReferenceInfo::from_reference_path(Path::new(
            "test/reference/GRCh38-and-mm10_ref_tiny",
        ))?;
        assert_json_snapshot!(ref_info.to_json_reporter());
        Ok(())
    }

    #[test]
    fn test_reference_info_from_probe_set_report() -> Result<()> {
        let ref_info = ReferenceInfo::from_probe_set_csv(&TargetSetFile::from(
            "test/probe_set_merge/GRCh38-fmt3-refv24_a.csv",
        ))?;
        assert_json_snapshot!(ref_info.to_json_reporter());
        Ok(())
    }

    #[test]
    fn test_validate_reference() -> Result<()> {
        let mut ref_info = ReferenceInfo::from_reference_path(Path::new(
            "../dui_tests/test_resources/reference/GRCh38-2020-A-chrM",
        ))?;
        ref_info.validate()?;

        // Break the checksums and ensure we get failures.
        ref_info.transcriptome_info.as_mut().unwrap().fasta_hash = "foobarbaz".to_string();

        let err = ref_info
            .validate()
            .expect_err("expected ref validation to fail");
        assert!(
            err.root_cause()
                .downcast_ref::<ChecksumMismatch>()
                .is_some(),
            "expected a checksum mismatch, but error was {err:?}"
        );
        Ok(())
    }
}
