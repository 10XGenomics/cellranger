use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::path::Path;

pub fn ensure_no_diff(actual: &Path, expected: &Path) {
    println!("Comparing {} and {}", actual.display(), expected.display());
    assert!(
        file_diff::diff(actual.to_str().unwrap(), expected.to_str().unwrap()),
        "Files {:?} and {:?} differ",
        actual.display(),
        expected.display()
    );
    println!(
        "Files {} and {} are identical",
        actual.display(),
        expected.display()
    );
}

/// Copy from src to dest. Panic if dest exists.
pub fn safe_copy(src: &Path, dest: &Path) -> Result<()> {
    assert!(
        !dest.exists(),
        "Error in safe_copy: destination file exists: {}",
        dest.display()
    );
    std::fs::copy(src, dest).with_context(|| {
        format!(
            "Error: unable to copy {} to {}",
            src.display(),
            dest.display()
        )
    })?;
    set_permissions(dest)
}

/// Set group write permission on coped files (linux only)
pub fn set_permissions(path: &Path) -> Result<()> {
    assert!(
        path.exists(),
        "Path {} does not exist to set any permission",
        path.display()
    );
    #[cfg(target_os = "linux")]
    {
        use std::os::unix::fs::PermissionsExt;
        let permissions = std::fs::Permissions::from_mode(0o775);
        std::fs::set_permissions(path, permissions)?;
    }
    Ok(())
}

pub fn ensure_identical_set_of_lines(actual: &Path, expected: &Path) {
    let actual_contents = std::fs::read_to_string(actual).unwrap();
    let actual_lines: HashSet<_> = actual_contents.lines().collect();
    let expected_contents = std::fs::read_to_string(expected).unwrap();
    let expected_lines: HashSet<_> = expected_contents.lines().collect();
    assert_eq!(
        actual_lines,
        expected_lines,
        "Files {} and {} differ",
        actual.display(),
        expected.display()
    );
}

pub fn diff_metrics(
    actual_metrics: HashMap<String, serde_json::value::Value>,
    expected_metrics: HashMap<String, serde_json::value::Value>,
) {
    let no_changed_metrics = expected_metrics
        .iter()
        .map(|(k, e)| match actual_metrics.get(k) {
            Some(_) if k == "cellranger_version" => true,
            // These metrics are subsampled in CR3, but not in CR4.
            Some(_)
                if k.ends_with("top_raw_barcodes")
                    || k.ends_with("top_raw_umis")
                    || k.ends_with("top_read_prefixes") =>
            {
                true
            }
            // these metrics include secondary alignments in CR3, but not in CR4
            Some(_)
                if k.ends_with("_intergenic_mapped_reads_frac")
                    || k.ends_with("_transcriptome_mapped_reads_frac")
                    || k.ends_with("_intronic_mapped_reads_frac")
                    || k.ends_with("_exonic_mapped_reads_frac")
                    || k.ends_with("genome_mapped_reads_frac")
                    || k == "multi_insert_size_histogram" =>
            {
                true
            }
            // PE metrics now always present in CR4
            Some(_)
                if k.ends_with("discordant_pairs_frac") || k.ends_with("improper_pairs_frac") =>
            {
                true
            }
            Some(a) if a != e => {
                println!("Key: {k} Expected: {e:?} Actual: {a:?}");
                false
            }
            None => {
                println!("Key: {k} Expected: {e:?} Actual: missing");
                false
            }
            _ => true,
        })
        .collect::<Vec<_>>()
        .iter()
        .all(|x| *x);

    let no_new_metrics = actual_metrics
        .iter()
        .map(|(k, a)| {
            // PE metrics now always present in CR4
            if expected_metrics.contains_key(k)
                || k == "cellranger_version"
                || k.ends_with("discordant_pairs_frac")
                || k.ends_with("improper_pairs_frac")
            {
                true
            } else {
                println!("Key: {k} Expected: missing Actual: {a:?}");
                false
            }
        })
        .collect::<Vec<_>>()
        .iter()
        .all(|x| *x);

    if no_changed_metrics && no_new_metrics {
        println!("Metrics are identical.");
    } else {
        panic!("Error: Metrics differ.");
    }
}
