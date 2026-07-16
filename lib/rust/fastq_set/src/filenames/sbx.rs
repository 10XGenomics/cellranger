//! Gathering FASTQ files from Roche SBX sequencers.
use crate::read_pair_iter::FastqSet;
use anyhow::{Context, Result, anyhow};
use regex::bytes::Regex;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::sync::LazyLock;

/// Gather FASTQs from a Roche SBX output directory.
/// TODO: this currently just grabs every fastq it finds in all subdirectories.
pub fn find_fastqs(read_path: &Path) -> Result<Vec<FastqSet>> {
    let mut fastqs = vec![];
    for entry in read_path
        .read_dir()
        .with_context(|| read_path.display().to_string())?
    {
        let path = entry
            .with_context(|| read_path.to_string_lossy().to_string())?
            .path();
        if path.is_dir() {
            fastqs.extend(find_fastqs(&path)?);
        } else if FASTQ_REGEX.is_match(path.as_os_str().as_bytes()) {
            fastqs.push(FastqSet {
                r1: path
                    .to_str()
                    .ok_or_else(|| {
                        anyhow!(
                            "FASTQ path contains non-unicode characters: {}",
                            path.display()
                        )
                    })?
                    .to_string(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: false,
            });
        }
    }
    Ok(fastqs)
}

static FASTQ_REGEX: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\.fastq(\.gz)?$").unwrap());
