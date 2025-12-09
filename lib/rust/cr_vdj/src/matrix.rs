#![deny(missing_docs)]

use anyhow::Result;
use hdf5::Group;
use hdf5::types::FixedAscii;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

martian_derive::martian_filetype! {H5File, "h5"}

pub fn open_matrix_group(matrix_path: &Path) -> Result<Group> {
    Ok(hdf5::File::open(matrix_path)?.group("matrix")?)
}

fn load_barcodes_from_group(matrix_group: &Group) -> Result<Vec<String>> {
    const MAX_BARCODE_AND_GEM_GROUP_LEN: usize = 44 + 2;
    let dataset = matrix_group.dataset("barcodes")?;
    if dataset.ndim() == 0 {
        Ok(vec![])
    } else {
        Ok(dataset
            .read_1d()?
            .into_iter()
            .map(|v: FixedAscii<MAX_BARCODE_AND_GEM_GROUP_LEN>| v.as_str().to_owned())
            .collect::<Vec<String>>())
    }
}

pub fn load_barcodes_from_matrix(matrix_path: &H5File) -> Result<Vec<String>> {
    load_barcodes_from_group(&open_matrix_group(matrix_path)?)
}

// Would be better to reuse scan-rs, but not currently possible due to dependency conflicts
pub fn gex_umi_counts_per_barcode(raw_matrix: &Path) -> Result<HashMap<String, u32>> {
    let matrix_group = open_matrix_group(raw_matrix)?;

    // Zero-based index into data / indices of the start of each column
    let csc_pointers: Vec<usize> = matrix_group.dataset("indptr")?.read_1d()?.to_vec();

    let is_gex_feature: Vec<_> = matrix_group
        .group("features")?
        .dataset("feature_type")?
        .read_1d()?
        .into_iter()
        .map(|feat: FixedAscii<256>| feat.as_str() == "Gene Expression")
        .collect();

    let barcodes = load_barcodes_from_group(&matrix_group)?;

    // FIXME: Whole matrix in memory. Load in chunks
    let counts = matrix_group.dataset("data")?.read_1d::<u32>()?.to_vec();
    let feat_indices: Vec<usize> = matrix_group.dataset("indices")?.read_1d()?.to_vec();

    let mut per_barcode_umi_counts = HashMap::new();

    for (barcode, (start, end)) in barcodes
        .into_iter()
        .zip_eq(csc_pointers.into_iter().tuple_windows::<(usize, usize)>())
    {
        per_barcode_umi_counts.insert(
            barcode,
            (start..end)
                .filter_map(|i| {
                    if is_gex_feature[feat_indices[i]] {
                        Some(counts[i])
                    } else {
                        None
                    }
                })
                .sum(),
        );
    }

    println!(
        "UMI Counts checksum: {}",
        per_barcode_umi_counts.values().sum::<u32>()
    );

    Ok(per_barcode_umi_counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_barcodes_empty_matrix() {
        assert!(
            load_barcodes_from_matrix(&H5File::from("test/count_matrix_no_bcs.h5"))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn test_load_barcodes() {
        assert_eq!(
            load_barcodes_from_matrix(&H5File::from("test/pbmc4k_tiny_97_bcs.h5"))
                .unwrap()
                .len(),
            97
        );
    }
}
