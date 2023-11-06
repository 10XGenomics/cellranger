// Read a molecule info file and compute usable reads per each library it contains
use cr_h5::molecule_info::{MoleculeInfoIterator, MoleculeInfoReader};
use ndarray::Axis;
use pyo3::prelude::*;
use std::collections::HashSet;
use std::path::Path;

/// Method to count reads for each given library from a molecule info file.
/// Takes in a molecule info file as a path and a list of library IDs.
/// iterates over the records, checks if a given molecule needs to come from a cell barcode,
/// or a certain feature (like a targeted gene), and increments read count accordingly.
#[pyfunction]
pub(crate) fn count_reads_per_library(
    _py: Python<'_>,
    path: String,
    lib_indices: Vec<usize>,
    filter_feature_idx: Vec<u64>,
    exclude_noncells: bool,
    custom_barcode_idx: HashSet<(u64, u16)>,
) -> PyResult<Vec<u64>> {
    let cur_path = Path::new(&path);
    let iterator = MoleculeInfoIterator::new(cur_path).unwrap();

    // get the set of cell barcodes
    let cell_barcodes: HashSet<(u64, u16)> = match !custom_barcode_idx.is_empty() & exclude_noncells
    {
        true => custom_barcode_idx,
        false => {
            let bcs = MoleculeInfoReader::read_barcode_info(cur_path).unwrap().0;
            let mut cell_barcodes_local = HashSet::with_capacity(bcs.nrows());
            for cell_bc in bcs.axis_iter(Axis(0)) {
                cell_barcodes_local.insert((cell_bc[0], cell_bc[1] as u16));
            }
            cell_barcodes_local
        }
    };

    // get the set of features to filter (typically targeted genes)
    let filter_feature_idx_set: HashSet<_> = filter_feature_idx.iter().copied().collect();

    // iterate over mol info records and increment read counts for each library,
    // provided the molecule passes certain filters (if given)
    let mut result = vec![0; lib_indices.len()];
    iterator.for_each(|x| {
        let library_idx = x.umi_data.library_idx;
        let feature_idx = x.umi_data.feature_idx as u64;
        let barcode_idx = x.barcode_idx;

        let mut count = x.umi_data.read_count as u64;
        if exclude_noncells && !cell_barcodes.contains(&(barcode_idx, library_idx)) {
            count = 0;
        }
        if !filter_feature_idx_set.contains(&feature_idx) {
            count = 0;
        }

        result[library_idx as usize] += count;
    });
    Ok(result)
}
