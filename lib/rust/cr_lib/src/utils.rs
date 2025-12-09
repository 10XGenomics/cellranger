#![deny(missing_docs)]
use anyhow::Result;
use martian::prelude::{MartianFileType, MartianRover};
use std::path::Path;

/// Attempts to create a hard link a file and creates a copy if hard linking fails
///
/// # Arguments
///
/// * `src` - The source file path
/// * `dst` - The destination link path
pub fn hard_link_file<P: AsRef<Path>>(src: P, dst: P) -> Result<()> {
    std::fs::hard_link(&src, &dst).or_else(|_| std::fs::copy(src, &dst).map(|_| ()))?;
    Ok(())
}

/// Creates a hard link (falls back to copy) to an external file in the directory of a running Martian stage
///
/// # Arguments
///
/// * `p` - The path to the external file
/// * `rover` - Martian rover for the running stage
///
/// # Returns
///
/// * The path of the newly created link
pub fn hard_link_martianfile<P: MartianFileType + AsRef<Path>>(
    p: P,
    rover: &MartianRover,
) -> Result<P> {
    let stem: &str = p
        .as_ref()
        .file_stem()
        .expect("Could not get file stem.")
        .to_str()
        .expect("Not valid unicode.");
    let new_link: P = rover.make_path(stem);
    hard_link_file(&p, &new_link)?;
    Ok(new_link)
}

pub mod estimate_mem {
    use serde_json::Value;
    use std::collections::HashMap;

    /// Return the total number of valid barcodes from a metrics JSON.
    /// This is computed by the BARCODE_CORRECTION stage as the size of the union
    /// set of the barcodes of all library types including barcodes_under_tissue.*
    pub fn get_total_barcodes_detected(metrics: &HashMap<String, Value>) -> usize {
        metrics
            .get("total_barcodes_detected")
            .and_then(Value::as_u64)
            .unwrap() as usize
    }

    /// Return the memory requirement for a data structure given barcode_count,
    /// bytes_per_barcode, and offset.
    pub fn barcode_mem_gib(
        barcodes_count: usize,
        bytes_per_barcode: usize,
        offset_gib: usize,
    ) -> isize {
        let mem_bytes = bytes_per_barcode * barcodes_count;

        (offset_gib + mem_bytes / 1024 / 1024 / 1024) as isize
    }
}
