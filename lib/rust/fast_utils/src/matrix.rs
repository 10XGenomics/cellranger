use std::path::Path;
use std::str::FromStr;

use cr_h5::count_matrix::{
    BarcodeCountOffset, Count, CountMatrixFile, CountMatrixStreaming, FeatureIdx, MatrixDimensions,
};
use cr_types::reference::feature_reference::FeatureType;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// A function that loads an H5 Matrix while filtering out certain feature types during loading.
///
/// We stored matrices as column sparse so to filter to certain features/rows is very inefficient.
/// In python land, we would usually accomplish this by loading the whole matrix, and then making a
/// subset copy, which can use a lot of memory and take a lot of time.  If only one feature type
/// is needed, this function allows you to do the feature subsetting during loading and uses an
/// async reader with ~2 threads to load the file as well, leading to much better performance.
#[allow(clippy::type_complexity)]
#[pyfunction]
pub(super) fn load_matrix_by_filtered_feature_type<'py>(
    py: Python<'py>,
    matrix_path: String,
    feature_type: String,
) -> PyResult<(
    Bound<'py, PyArray1<BarcodeCountOffset>>,
    Bound<'py, PyArray1<FeatureIdx>>,
    Bound<'py, PyArray1<Count>>,
)> {
    let feature_type_resolved = FeatureType::from_str(&feature_type);
    match feature_type_resolved {
        Ok(ft) => {
            // Make our vectors, with capacity up to the max size
            let matrix_file = CountMatrixFile::from(Path::new(&matrix_path));
            let mat_dims: MatrixDimensions = matrix_file
                .load_dimensions()
                .expect("Could not read matrix dimension sfrom file.");

            let mut row_indices = Vec::with_capacity(mat_dims.num_non_zeros);
            let mut data = Vec::with_capacity(mat_dims.num_non_zeros);
            let mut indptr = Vec::with_capacity(mat_dims.num_barcodes + 1);
            indptr.push(0);
            {
                let matrix_streaming: CountMatrixStreaming = matrix_file
                    .read_streaming_async()
                    .expect("Could not open matrix file.");

                let mut cur_cnt = 0;
                let mut cur_bc: usize = 0;
                let is_selected_feature: Vec<bool> = matrix_streaming
                    .bcs_and_features
                    .feature_reference
                    .feature_defs
                    .iter()
                    .map(|x| x.feature_type == ft)
                    .collect();
                for raw_cnt_res in matrix_streaming
                    .raw_counts()
                    .expect("Could not build iterator over matrix")
                {
                    let raw_cnt = raw_cnt_res.expect("Could not get raw count");
                    let bc = raw_cnt.barcode_idx;
                    if bc != cur_bc {
                        indptr.push(cur_cnt);
                        cur_bc = bc;
                    }
                    if is_selected_feature[raw_cnt.feature_idx] {
                        row_indices.push(raw_cnt.feature_idx as FeatureIdx);
                        data.push(raw_cnt.count);
                        cur_cnt += 1;
                    }
                }
                indptr.push(cur_cnt);
            }
            indptr.shrink_to_fit();
            row_indices.shrink_to_fit();
            data.shrink_to_fit();
            Ok((
                indptr.into_pyarray(py),
                row_indices.into_pyarray(py),
                data.into_pyarray(py),
            ))
        }
        Err(e) => Err(PyErr::new::<PyTypeError, _>(format!(
            "Failed to convert FeatureType {feature_type}: {e}"
        ))),
    }
}
