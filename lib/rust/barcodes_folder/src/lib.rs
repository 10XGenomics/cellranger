//! Crate for resolving paths from lib/python/{assay}/barcodes folder.
//!

use anyhow::{bail, Context, Result};
use itertools::Itertools;
use std::path::PathBuf;

const PD_CANONICAL_SLIDE_NAME: &str = "visium_hd_rc1";
const CS_CANONICAL_SLIDE_NAME: &str = "visium_hd_v1";

fn search_for_whitelist(
    whitelist_name: &str,
    translation: bool,
    mut p: PathBuf,
) -> Option<PathBuf> {
    if translation {
        return search_for_whitelist(whitelist_name, false, p.join("translation"));
    }
    if p.is_dir() {
        p.push(whitelist_name);

        let p_gz = p.with_extension("txt.gz");
        if p_gz.is_file() {
            return Some(p_gz);
        }

        p.set_extension("txt");
        if p.is_file() {
            return Some(p);
        }

        // Slide design file
        p.set_extension("slide");
        if p.is_file() {
            return Some(p);
        }
    }

    None
}

pub fn find_slide_design(slide_name: &str) -> Result<PathBuf> {
    let fname = if slide_name == PD_CANONICAL_SLIDE_NAME || slide_name == CS_CANONICAL_SLIDE_NAME {
        find_whitelist(PD_CANONICAL_SLIDE_NAME, false)
            .or_else(|_|find_whitelist(CS_CANONICAL_SLIDE_NAME, false))
            .context( format!("Could not find slide files for slide designs {CS_CANONICAL_SLIDE_NAME} or {PD_CANONICAL_SLIDE_NAME}"))
    } else {
        find_whitelist(slide_name, false)
    }?;
    if fname.extension().is_some_and(|ext| ext == "slide") {
        Ok(fname)
    } else {
        bail!("Could not find slide design file for slide {}", slide_name)
    }
}

/// Find the path to a barcode whitelist file, given the whitelist name.
/// This searches the known relative paths between built rust executables
/// and the barcode whitelist folder.
pub fn find_whitelist(whitelist_name: &str, translation: bool) -> Result<PathBuf> {
    find_whitelist_in_folder(whitelist_name, translation, "cellranger")
}

/// Find the path to a barcode whitelist file for ATAC.
/// This is shim code to support this deprecated pattern for legacy ATAC code.
pub fn find_atac_whitelist(whitelist_name: &str) -> Result<PathBuf> {
    find_whitelist_in_folder(whitelist_name, false, "atac")
}

/// Find the path to a barcode whitelist file, given the whitelist name.
/// This searches the known relative paths between built rust executables
/// and the barcode whitelist folder.
/// The barcode folder is provided as the final argument.
fn find_whitelist_in_folder(
    whitelist_name: &str,
    translation: bool,
    folder: &str,
) -> Result<PathBuf> {
    let exe = bazel_utils::current_exe()?;
    let exe_path = exe.parent().unwrap().to_path_buf();

    let mut barcode_search_paths = Vec::new();

    // Option 1: find relative to the PYTHONPATH and default to ../python, since in a Cell
    // Ranger build the exe will be at lib/bin/<exe>, and the barcodes are at
    // lib/python/barcodes. So we try "<exe path>/../python/cellranger/barcodes"
    for ppath in std::env::var("PYTHONPATH")
        .unwrap_or_else(|_| String::from("../python"))
        .split(':')
        .filter(|path| path.ends_with("python"))
    {
        barcode_search_paths.push(format!("{ppath}/{folder}/barcodes"));
    }

    // Option 2: in a checkout on disk and a cargo build, the exe will be at
    // lib/rust/target/release/<exe>, and the barcodes are at lib/python/barcodes.
    // So we try "<exe path>/../../../python/cellranger/barcodes"
    barcode_search_paths.push(format!("../../../python/{folder}/barcodes"));

    // Option 3: in a checkout on disk, the cargo test build exe will be at
    // lib/rust/target/release/deps/<exe>, and the barcodes are at lib/python/barcodes.
    // So we try "<exe path>/../../../../python/cellranger/barcodes"
    barcode_search_paths.push(format!("../../../../python/{folder}/barcodes"));

    // Option 4: in a checkout on disk using the --target flag, the cargo test build exe will be at
    // lib/rust/target/release/deps/<target>/<exe>, and the barcodes are at lib/python/barcodes.
    // So we try "<exe path>/../../../../../python/cellranger/barcodes"
    barcode_search_paths.push(format!("../../../../../python/{folder}/barcodes"));

    // Option 5: find relative to the current working directory, in case we're
    // running under cargo test. When running tests under bazel test the cwd is
    // lib/rust/cr_[atac|lib] => we want to go ../../python/cellranger/barcodes
    barcode_search_paths.push(format!("../../python/{folder}/barcodes"));

    for rel_path in &barcode_search_paths {
        let wl_path = search_for_whitelist(whitelist_name, translation, exe_path.join(rel_path));
        if let Some(p) = wl_path {
            return Ok(std::fs::canonicalize(p)?);
        }
    }

    // cr_lib can be used as a library, use debug_assertions as proxy for dev builds
    // Option 6: find relative to the CARGO_MANIFEST_DIR, which may happen
    // if you've set your CARGO_HOME to a non-default value.
    #[cfg(debug_assertions)]
    {
        let path = format!(
            "{}/../../python/{}/barcodes",
            env!("CARGO_MANIFEST_DIR"),
            folder,
        );
        let wl_path = search_for_whitelist(whitelist_name, translation, path.into());
        if let Some(p) = wl_path {
            return Ok(std::fs::canonicalize(p)?);
        }
    }

    bail!(
        "Couldn't find the barcode whitelist {} (translation: {})in the following paths:\n{}",
        whitelist_name,
        translation,
        barcode_search_paths
            .iter()
            .enumerate()
            .map(|(i, p)| format!("{i} {:?}", exe_path.join(p)))
            .join("\n")
    )
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_find_whitelist() -> Result<()> {
        let _ = find_whitelist("737K-august-2016", false)?;
        assert!(find_whitelist("737K-august-2016", true).is_err());
        let _ = find_whitelist("3M-february-2018", false)?;
        let _ = find_whitelist("3M-february-2018", true)?;
        let _ = find_whitelist("9K-LT-march-2021", false)?;
        let _ = find_whitelist("9K-LT-march-2021", true)?;
        let _ = find_atac_whitelist("737K-arc-v1")?;
        Ok(())
    }

    #[test]
    fn test_find_slide_design() {
        assert!(find_slide_design("737K-august-2016").is_err());
        assert!(find_slide_design("visium_hd_rc1").is_ok());
    }
}
