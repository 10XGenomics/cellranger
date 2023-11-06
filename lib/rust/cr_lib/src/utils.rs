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
