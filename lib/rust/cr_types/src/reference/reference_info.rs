use anyhow::{ensure, Context, Result};
use martian::MartianFileType;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use metric::{JsonReport, JsonReporter, Metric};
use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};
use std::fs::File;
use std::path::{Path, PathBuf};

pub type GenomeName = String;

/// Use this string to join multiple genome names in a multi-genome reference.
pub const MULTI_GENOME_SEPARATOR: &str = "_and_";

/// Metadata associated with the reference from the `reference.json` file in the reference directory
///
/// This is primarily used by a number of functions to query the `genomes` present in the reference.
/// This also appears in the metrics summary json with `reference_` prefix. See the `JsonReport`
/// impl for details. Some of these fields are also added to the molecule info metrics json.
///
/// NOTE: We use `#[serde(default)]` for certain fields because some older versions of references,
/// for example `hg19-1.0.0`, did not populate these fields in the `reference.json`.
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
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
    pub version: Option<String>,
    #[serde(default)]
    pub path: PathBuf,
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
}

impl JsonReport for ReferenceInfo {
    fn to_json_reporter(&self) -> JsonReporter {
        let mut reporter = JsonReporter::new();
        reporter.insert("fasta_hash", self.fasta_hash.clone());
        reporter.insert("genomes", self.genomes.join(MULTI_GENOME_SEPARATOR));
        reporter.insert("gtf_hash", self.gtf_hash.clone());
        reporter.insert("gtf_hash.gz", self.gtf_hash_gz.clone());
        reporter.insert("input_fasta_files", self.input_fasta_files.join(", "));
        reporter.insert("input_gtf_files", self.input_gtf_files.join(", "));
        reporter.insert("mkref_version", self.mkref_version.clone());
        match self.version {
            Some(ref version) => reporter.insert("version", version.clone()),
            None => reporter.insert("version", ""),
        }
        reporter.insert("type", "Transcriptome");
        reporter.insert("path", self.path.to_str());
        reporter.add_prefix("reference");
        reporter
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
        let ref_info: ReferenceInfo =
            JsonFile::from("test/reference/single_genome_reference.json").read()?;
        assert_json_snapshot!(ref_info.to_json_reporter());
        Ok(())
    }

    #[test]
    fn test_reference_info_multi_report() -> Result<()> {
        let ref_info: ReferenceInfo =
            JsonFile::from("test/reference/multi_genome_reference.json").read()?;
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
