mod barcode_counter;
pub mod barcode_index;
mod compute_extra_multiplexing_metrics;
mod diff_exp;
mod filtered_barcodes;
mod gdna_analysis;
mod library_read_counter;
mod molecule_info;
mod multi_graph;

use crate::barcode_counter::counts_per_barcode;
use crate::compute_extra_multiplexing_metrics::tag_read_counts;
use crate::diff_exp::{compute_sseq_params_o3, sseq_differential_expression_o3};
use crate::gdna_analysis::count_umis_per_probe;
use crate::library_read_counter::count_reads_per_library;
use crate::molecule_info::{
    concatenate_molecule_infos, count_reads_and_reads_in_cells, count_usable_reads,
    downsample_molinfo, get_num_umis_per_barcode,
};
use crate::multi_graph::MultiGraph;
use barcode::binned::SquareBinIndex;
use barcode::Barcode;
use barcode_index::MatrixBarcodeIndex;
use cr_h5::molecule_info::{MoleculeInfoIterator, MoleculeInfoReader};
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::{
    BcCountDataType, BcCountFormat, LibraryType, TotalBcCountDataType, TotalBcCountFormat,
};
use filtered_barcodes::FilteredBarcodes;
use martian_filetypes::FileTypeRead;
use metric::CountMetric;
use multi_graph::PyFingerprint;
use ndarray::prelude::*;
use ndarray::Array1;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use pyo3::wrap_pyfunction;
use std::collections::HashSet;
use std::path::Path;

// Method counts the total number of reads for each feature
#[allow(clippy::type_complexity)]
#[pyfunction]
fn reads_per_feature(
    _py: Python<'_>,
    path: String,
    separate_by_cell_call: bool,
) -> PyResult<(Vec<u32>, Vec<String>, Vec<String>, Option<Vec<u32>>)> {
    // Sum up all the read counts for a list of features, returning a tuple of
    // cnts, feature_id, feature_type
    // method is inefficient because it parses all data columns in mol_info but only uses 2,
    // also is efficient because it doesn't use a lot of memory of once and streams the mol_info
    let cur_path = Path::new(&path);
    let iterator = MoleculeInfoIterator::new(cur_path).unwrap();
    let mut cell_bcs_cnts = match separate_by_cell_call {
        true => {
            let bcs = MoleculeInfoReader::read_barcode_info(cur_path).unwrap().0;
            let mut cell_barcodes = HashSet::with_capacity(bcs.nrows());
            for cell_bc in bcs.axis_iter(Axis(0)) {
                cell_barcodes.insert((cell_bc[0], cell_bc[1] as u16));
            }
            let cell_cnts: Array1<u32> = Array1::zeros(iterator.feature_ref.feature_defs.len());
            Some((cell_barcodes, cell_cnts))
        }
        false => None,
    };
    let mut cnts: Array1<u32> = Array1::zeros(iterator.feature_ref.feature_defs.len());
    // Note assumes (as should be the case) feature_defs are in order
    let feature_names: Vec<String> = iterator
        .feature_ref
        .feature_defs
        .iter()
        .map(|feat| feat.id.to_string())
        .collect();
    let feature_types: Vec<String> = iterator
        .feature_ref
        .feature_defs
        .iter()
        .map(|feat| feat.feature_type.to_string())
        .collect();
    // Counts up reads per feature
    let mut last_valid_bc_lib: (u64, u16) = (u64::MAX, u16::MAX);
    let mut last_invalid_bc_lib: (u64, u16) = (u64::MAX, u16::MAX);

    iterator.for_each(|x| {
        cnts[[x.umi_data.feature_idx as usize]] += x.umi_data.read_count;
        if let Some((ref valid, ref mut cell_cnts)) = cell_bcs_cnts {
            let bidx = x.barcode_idx;
            let lidx = x.umi_data.library_idx;
            let is_cell: bool;
            let bidx_lidx = (bidx, lidx);
            // Cache results for the same barcode/library
            if bidx_lidx == last_valid_bc_lib {
                is_cell = true;
            } else if bidx_lidx == last_invalid_bc_lib {
                is_cell = false;
            } else {
                is_cell = valid.contains(&bidx_lidx);
                if is_cell {
                    last_valid_bc_lib = bidx_lidx;
                } else {
                    last_invalid_bc_lib = bidx_lidx;
                }
            }
            if is_cell {
                cell_cnts[[x.umi_data.feature_idx as usize]] += x.umi_data.read_count;
            }
        }
    });
    let cell_cnts = cell_bcs_cnts.map(|(_, y)| y.to_vec());
    Ok((cnts.to_vec(), feature_names, feature_types, cell_cnts))
}

#[pyfunction]
fn get_gex_barcode_cnts(py: Python<'_>, file_path: String) -> Vec<(&PyBytes, i64)> {
    // Gets the total_barcode_counts file output by BARCODE_CORRECTION stage and
    // allows us to easily do analysis on this.
    // Extracts the Gene Expression section.
    let data: BcCountFormat = BcCountFormat::from(file_path);
    let mut data_by_feature: BcCountDataType = data.read().unwrap();
    data_by_feature
        .remove(&LibraryType::Gex)
        .unwrap()
        .drain()
        .map(|(k, v): (Barcode, CountMetric)| (PyBytes::new(py, k.sequence_bytes()), v.count()))
        .collect()
}

#[pyfunction]
fn get_total_barcode_cnts(py: Python<'_>, file_path: String) -> PyResult<&PyDict> {
    let reader: TotalBcCountFormat = TotalBcCountFormat::from(file_path);
    let mut data: TotalBcCountDataType = reader.read().unwrap();
    let dict_to_return = PyDict::new(py);
    for (k, v) in data.drain() {
        dict_to_return.set_item(PyBytes::new(py, k.sequence_bytes()), v.count())?;
    }
    Ok(dict_to_return)
}

#[pyfunction]
fn validate_reference(_py: Python<'_>, reference_path: &str) -> PyResult<()> {
    ReferenceInfo::from_reference_path(Path::new(&reference_path))
        .and_then(|ref_info| ref_info.validate())
        .map_err(|err| PyErr::new::<PyException, _>(format!("{err:#}")))
}

/// The python module in our extension module that will be imported.
#[pymodule]
fn fast_utils(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reads_per_feature, m)?)?;
    m.add_function(wrap_pyfunction!(tag_read_counts, m)?)?;
    m.add_function(wrap_pyfunction!(counts_per_barcode, m)?)?;
    m.add_function(wrap_pyfunction!(count_reads_per_library, m)?)?;
    m.add_function(wrap_pyfunction!(downsample_molinfo, m)?)?;
    m.add_function(wrap_pyfunction!(concatenate_molecule_infos, m)?)?;
    m.add_function(wrap_pyfunction!(get_gex_barcode_cnts, m)?)?;
    m.add_function(wrap_pyfunction!(get_total_barcode_cnts, m)?)?;
    m.add_function(wrap_pyfunction!(count_umis_per_probe, m)?)?;
    m.add_function(wrap_pyfunction!(count_reads_and_reads_in_cells, m)?)?;
    m.add_function(wrap_pyfunction!(count_usable_reads, m)?)?;
    m.add_function(wrap_pyfunction!(sseq_differential_expression_o3, m)?)?;
    m.add_function(wrap_pyfunction!(compute_sseq_params_o3, m)?)?;
    m.add_function(wrap_pyfunction!(validate_reference, m)?)?;
    m.add_function(wrap_pyfunction!(get_num_umis_per_barcode, m)?)?;
    m.add_class::<MultiGraph>()?;
    m.add_class::<PyFingerprint>()?;
    m.add_class::<FilteredBarcodes>()?;
    m.add_class::<MatrixBarcodeIndex>()?;
    m.add_function(wrap_pyfunction!(
        filtered_barcodes::save_filtered_bcs_groups,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        filtered_barcodes::load_filtered_bcs_groups,
        m
    )?)?;
    m.add_class::<SquareBinIndex>()?;

    Ok(())
}
