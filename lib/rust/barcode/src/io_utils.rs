use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

const GZ_BUF_SIZE: usize = 1 << 20;

/// Open a (possibly gzipped) file into a BufReader.
pub(crate) fn open_with_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| path.display().to_string())?;
    match path.extension().unwrap_or_default().to_str().unwrap() {
        "gz" => Ok(Box::new(BufReader::with_capacity(
            GZ_BUF_SIZE,
            MultiGzDecoder::new(f),
        ))),
        _ => Ok(Box::new(BufReader::with_capacity(32 * 1024, f))),
    }
}
