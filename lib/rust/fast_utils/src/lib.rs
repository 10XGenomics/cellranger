//! fast_utils
#![deny(missing_docs)]
mod barcode_counter;
mod barcode_index;
mod compute_extra_multiplexing_metrics;
mod diff_exp;
mod feature_collection_geojson_validator;
mod filtered_barcodes;
mod gdna_analysis;
mod matrix;
mod molecule_info;
mod multi_graph;

use crate::barcode_counter::counts_per_barcode;
use crate::compute_extra_multiplexing_metrics::tag_read_counts;
use crate::diff_exp::{compute_sseq_params_o3, sseq_differential_expression_o3};
use crate::feature_collection_geojson_validator::FeatureCollectionGeoJson;
use crate::gdna_analysis::count_umis_per_probe;
use crate::matrix::load_matrix_by_filtered_feature_type;
use crate::molecule_info::{
    count_reads_and_reads_in_cells, count_usable_reads, downsample_molinfo,
    get_num_umis_per_barcode,
};
use crate::multi_graph::MultiGraph;
use barcode::binned::SquareBinIndex;
use barcode::cell_name::CellId;
use barcode_index::MatrixBarcodeIndex;
use cr_h5::count_matrix::CountMatrixFile;
use cr_h5::molecule_info::{MoleculeInfoIterator, MoleculeInfoReader};
use cr_types::SortedBarcodeCountFile;
use cr_types::probe_set::merge_probe_set_csvs;
use cr_types::reference::feature_reference::FeatureType;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::ReferenceInfo;
use filtered_barcodes::FilteredBarcodes;
use geo::{LineString, SimplifyIdx};
use itertools::Itertools;
use martian::MartianFileType;
use martian_filetypes::LazyFileTypeIO;
use ndarray::Array1;
use ndarray::prelude::*;
use pyanyhow::PyAnyhowError;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::wrap_pyfunction;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::str::FromStr;

// Method counts the total number of reads for each feature
#[allow(clippy::type_complexity)]
#[pyfunction]
fn reads_per_feature(
    _py: Python<'_>,
    path: String,
    separate_by_cell_call: bool,
) -> PyResult<(Vec<u32>, Vec<String>, Vec<&'static str>, Option<Vec<u32>>)> {
    // Sum up all the read counts for a list of features, returning a tuple of
    // cnts, feature_id, feature_type
    // method is inefficient because it parses all data columns in mol_info but only uses 2,
    // also is efficient because it doesn't use a lot of memory of once and streams the mol_info
    let cur_path = Path::new(&path);
    let iterator = MoleculeInfoIterator::new(cur_path).unwrap();
    let mut cell_bcs_cnts = match separate_by_cell_call {
        true => {
            let bcs = MoleculeInfoReader::read_barcode_info_pass_filter(cur_path).unwrap();
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
    let feature_types = iterator
        .feature_ref
        .feature_defs
        .iter()
        .map(|feat| feat.feature_type.as_str())
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
fn get_gex_barcode_cnts(py: Python<'_>, file_path: String) -> Vec<(Bound<'_, PyBytes>, usize)> {
    // The file_path should point to a SortedBarcodeCountFile containing the
    // GEX barcode counts.
    SortedBarcodeCountFile::from(file_path)
        .lazy_reader()
        .unwrap()
        .process_results(|bc_counts| {
            bc_counts
                .map(|(bc, count)| (PyBytes::new(py, bc.sequence_bytes()), count))
                .collect()
        })
        .unwrap()
}

/// Given path to a Feature Barcode Matrix and a library type, gets the  UMIs in that
/// library per barcode. Does not load the entire matrix into memory
#[pyfunction]
fn per_barcode_umis_in_library(
    _py: Python<'_>,
    matrix_path: String,
    library: String,
) -> PyResult<HashMap<String, i64>> {
    let fbm_path = Path::new(&matrix_path);
    let library_type = FeatureType::from_str(&library).map_err(PyAnyhowError::from)?;
    let fbm = CountMatrixFile::from_path(fbm_path)
        .read_streaming()
        .map_err(PyAnyhowError::from)?;
    let per_bc_counts: Result<HashMap<_, _>, _> = fbm
        .barcode_counts_for_feature_type(library_type)
        .map_err(PyAnyhowError::from)?
        .map(|x| {
            let y = x.map_err(PyAnyhowError::from)?;
            Ok((String::from(y.0.as_str()), y.1))
        })
        .collect();
    per_bc_counts
}

#[pyfunction]
fn get_streaming_fbc_mem_gb_estimate(_py: Python<'_>, matrix_path: String) -> PyResult<f64> {
    let fbm_path = Path::new(&matrix_path);
    Ok(CountMatrixFile::from_path(fbm_path)
        .estimate_mem_gib_without_counts()
        .map_err(PyAnyhowError::from)?)
}

#[pyfunction]
fn validate_feature_collection_geojson(_py: Python<'_>, geojson_path: &str) -> PyResult<()> {
    FeatureCollectionGeoJson::from_geojson_path(Path::new(&geojson_path).to_path_buf())
        .map(|_| ())
        .map_err(|err| PyErr::new::<PyException, _>(format!("{err:#}")))
}

#[pyfunction]
fn validate_reference(_py: Python<'_>, reference_path: &str) -> PyResult<()> {
    ReferenceInfo::from_reference_path(Path::new(&reference_path))
        .and_then(|ref_info| ref_info.validate())
        .map_err(|err| PyErr::new::<PyException, _>(format!("{err:#}")))
}

#[pyfunction]
fn merge_probe_sets(
    _py: Python<'_>,
    csv_paths: Vec<PathBuf>,
    output_path: PathBuf,
) -> pyanyhow::Result<()> {
    let inputs: Vec<TargetSetFile> = csv_paths.iter().map(Into::into).collect();
    let mut output = File::create(&output_path)?;
    merge_probe_set_csvs(&inputs, &mut output, None)?;
    Ok(())
}

/// Simplify the provided line using Ramer-Douglas-Peucker.
///
/// Return the indices of the line that should be retained.
#[pyfunction]
pub fn simplify_line(points: Vec<(f64, f64)>, epsilon: f64) -> Vec<usize> {
    LineString::from(points).simplify_idx(epsilon)
}

/// The python module in our extension module that will be imported.
#[pymodule]
fn fast_utils(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reads_per_feature, m)?)?;
    m.add_function(wrap_pyfunction!(tag_read_counts, m)?)?;
    m.add_function(wrap_pyfunction!(counts_per_barcode, m)?)?;
    m.add_function(wrap_pyfunction!(downsample_molinfo, m)?)?;
    m.add_function(wrap_pyfunction!(get_gex_barcode_cnts, m)?)?;
    m.add_function(wrap_pyfunction!(count_umis_per_probe, m)?)?;
    m.add_function(wrap_pyfunction!(per_barcode_umis_in_library, m)?)?;
    m.add_function(wrap_pyfunction!(get_streaming_fbc_mem_gb_estimate, m)?)?;
    m.add_function(wrap_pyfunction!(count_reads_and_reads_in_cells, m)?)?;
    m.add_function(wrap_pyfunction!(count_usable_reads, m)?)?;
    m.add_function(wrap_pyfunction!(sseq_differential_expression_o3, m)?)?;
    m.add_function(wrap_pyfunction!(compute_sseq_params_o3, m)?)?;
    m.add_function(wrap_pyfunction!(validate_feature_collection_geojson, m)?)?;
    m.add_function(wrap_pyfunction!(validate_reference, m)?)?;
    m.add_function(wrap_pyfunction!(merge_probe_sets, m)?)?;
    m.add_function(wrap_pyfunction!(get_num_umis_per_barcode, m)?)?;
    m.add_function(wrap_pyfunction!(simplify_line, m)?)?;
    m.add_function(wrap_pyfunction!(load_matrix_by_filtered_feature_type, m)?)?;
    m.add_class::<MultiGraph>()?;
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
    m.add_class::<CellId>()?;

    Ok(())
}
