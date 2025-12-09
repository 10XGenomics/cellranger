//! I/O functionality
#![deny(missing_docs)]

use anyhow::{Context, Result};
use itertools::Itertools;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// Count the number of lines in a file.
/// The number of newline characters is counted.
/// A final line missing a terminal newline character is not counted.
pub(crate) fn count_lines(path: &Path) -> Result<usize> {
    BufReader::new(File::open(path).with_context(|| path.display().to_string())?)
        .bytes()
        .fold_ok(0, |acc, c| acc + if c == b'\n' { 1 } else { 0 })
        .with_context(|| path.display().to_string())
}
