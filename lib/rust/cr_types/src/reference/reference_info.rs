use crate::GenomeName;
use anyhow::{ensure, Context, Result};
use itertools::Itertools;
use martian::MartianFileType;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
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
pub struct ReferenceInfo {
    #[serde(default)]
    pub fasta_hash: String,
    pub genomes: Vec<GenomeName>,
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
    pub path: PathBuf,
    pub version: Option<String>,
}

impl ReferenceInfo {
    pub fn from_reference_path(reference_path: &Path) -> Result<ReferenceInfo> {
        Ok(ReferenceInfo {
            path: reference_path.to_path_buf(),
            ..JsonFile::new(reference_path, "reference.json").read()?
        })
    }

    /// Validate that associated data complies with expectations.
    /// Check that FASTA and GTF files exist.
    /// Check that FASTA and GTF files match checksums, if present.
    pub fn validate(&self) -> Result<()> {
        // Checksum the FASTA file.
        if !self.fasta_hash.is_empty() {
            sha1_checksum_file(&self.path.join("fasta/genome.fa"), &self.fasta_hash)
                .context("failed to validate FASTA file for reference")?;
        }
        // Checksum the GTF file.
        if !self.gtf_hash.is_empty() {
            sha1_checksum_file(&self.path.join("genes/genes.gtf"), &self.gtf_hash)
                .context("failed to validate GTF file for reference")?;
        }
        Ok(())
    }

    pub fn into_report(self) -> ReferenceInfoReport {
        ReferenceInfoReport {
            reference_fasta_hash: self.fasta_hash,
            reference_genomes: self.genomes.into_iter().join(MULTI_GENOME_SEPARATOR),
            reference_gtf_hash: self.gtf_hash,
            reference_gtf_hash_gz: self.gtf_hash_gz,
            reference_input_fasta_files: self.input_fasta_files.join(", "),
            reference_input_gtf_files: self.input_gtf_files.join(", "),
            reference_mkref_version: self.mkref_version,
            reference_path: self.path,
            reference_type: "Transcriptome".to_string(),
            reference_version: self.version.unwrap_or_default(),
        }
    }

    pub fn into_json_report(self) -> Map<String, Value> {
        if let Value::Object(map) = serde_json::to_value(self.into_report()).unwrap() {
            map
        } else {
            unreachable!()
        }
    }
}

/// Reference info metrics
#[derive(Serialize)]
pub struct ReferenceInfoReport {
    pub reference_fasta_hash: String,
    pub reference_genomes: String,
    pub reference_gtf_hash: String,
    #[serde(rename = "reference_gtf_hash.gz")]
    pub reference_gtf_hash_gz: String,
    pub reference_input_fasta_files: String,
    pub reference_input_gtf_files: String,
    pub reference_mkref_version: String,
    pub reference_version: String,
    pub reference_type: String,
    pub reference_path: PathBuf,
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
        let ref_info: ReferenceInfo =
            JsonFile::from("test/reference/single_genome_reference.json").read()?;
        assert_json_snapshot!(ref_info.into_json_report());
        Ok(())
    }

    #[test]
    fn test_reference_info_multi_report() -> Result<()> {
        let ref_info: ReferenceInfo =
            JsonFile::from("test/reference/multi_genome_reference.json").read()?;
        assert_json_snapshot!(ref_info.into_json_report());
        Ok(())
    }

    #[test]
    fn test_validate_reference() -> Result<()> {
        let mut ref_info = ReferenceInfo::from_reference_path(Path::new(
            "../dui_tests/test_resources/reference/GRCh38-2020-A-chrM",
        ))?;
        ref_info.validate()?;

        // Break the checksums and ensure we get failures.
        ref_info.fasta_hash = "foobarbaz".to_string();

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
