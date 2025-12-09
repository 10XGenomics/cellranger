#![deny(missing_docs)]
use anyhow::{Context, Result};
use diff_exp::diff_exp::{DiffExpResult, SSeqParams};
use diff_exp::{compute_sseq_params, sseq_differential_expression};
use numpy::{Element, IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBool, PyDict, PyList};
use rayon::ThreadPoolBuilder;
use sprs::CsMatBase;
use sqz::mat::AdaptiveMatOwned;
use sqz::matrix_map::MatrixIntoMap;
use std::sync::Once;

static INIT_THREAD_POOL: Once = Once::new();

fn init_global_thread_pool(num_threads: usize) {
    INIT_THREAD_POOL.call_once(|| {
        ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .expect("Failed to build global thread pool");
    });
}

/// Extract a python array of a desired type.
fn extract_pyarray<'a, T: Element + 'a>(
    py: &Bound<'a, PyAny>,
    attr: &'static str,
) -> Result<PyReadonlyArray1<'a, T>> {
    py.getattr(attr).expect(attr).extract().context(attr)
}

/// Transmute a &[i32] to a &[u32] and validate that all elements are positive.
fn transmute_i32_to_u32(xs: &[i32]) -> &[u32] {
    assert!(xs.iter().all(|&x| x >= 0));
    unsafe { std::slice::from_raw_parts(xs.as_ptr().cast(), xs.len()) }
}

/// Convert a scipy sparse matrix to a rust sparse matrix.
fn convert_matrix(
    _py: Python<'_>,
    matrix: &Bound<'_, PyAny>,
) -> AdaptiveMatOwned<u32, MatrixIntoMap<u32, u32>> {
    let shape: (usize, usize) = matrix.getattr("shape").unwrap().extract().unwrap();
    let data_py = extract_pyarray(matrix, "data").unwrap();
    let data = transmute_i32_to_u32(data_py.as_slice().unwrap());

    // indices and indptr may be either i32 or i64. They have the same type.
    if let Ok(indptr) = extract_pyarray::<i32>(matrix, "indptr") {
        let indptr = indptr.as_slice().unwrap();
        let indices = extract_pyarray::<i32>(matrix, "indices").unwrap();
        let indices = indices.as_slice().unwrap();
        AdaptiveMatOwned::from_csmat(&CsMatBase::new_csc(shape, indptr, indices, data))
    } else {
        let indptr = extract_pyarray::<i64>(matrix, "indptr")
            .context("indptr is neither i32 nor i64")
            .unwrap();
        let indptr = indptr.as_slice().unwrap();
        let indices = extract_pyarray::<i64>(matrix, "indices").unwrap();
        let indices = indices.as_slice().unwrap();
        AdaptiveMatOwned::from_csmat(&CsMatBase::new_csc(shape, indptr, indices, data))
    }
}

#[pyfunction]
pub(super) fn compute_sseq_params_o3<'a>(
    py: Python<'a>,
    matrix: &Bound<'a, PyAny>,
    zeta_quartile: f64,
) -> PyResult<Bound<'a, PyDict>> {
    let converted_matrix = convert_matrix(py, matrix);
    let sseq_params = compute_sseq_params(&converted_matrix, Some(zeta_quartile), None, None);
    let params = PyDict::new(py);
    params.set_item("N", converted_matrix.cols())?;
    params.set_item("G", converted_matrix.rows())?;
    params.set_item("size_factors", sseq_params.size_factors.into_pyarray(py))?;
    params.set_item("mean_g", sseq_params.gene_means.into_pyarray(py))?;
    params.set_item("var_g", sseq_params.gene_variances.into_pyarray(py))?;
    params.set_item("use_g", sseq_params.use_genes.into_pyarray(py))?;
    params.set_item("phi_mm_g", sseq_params.gene_moment_phi.into_pyarray(py))?;
    params.set_item("zeta_hat", sseq_params.zeta_hat)?;
    params.set_item("delta", sseq_params.delta)?;
    params.set_item("phi_g", sseq_params.gene_phi)?;
    Ok(params)
}

// Note this originally returned a typed #[pyclass] class, but that led to a copy on the returned value
// that I couldn't figure out how to avoid, so am returning a dictionary to avoid the
#[pyfunction]
#[pyo3(signature = (matrix, cond_a, cond_b, sseq_params, big_count=None, threads=1))]
pub(super) fn sseq_differential_expression_o3<'a>(
    py: Python<'a>,
    matrix: &Bound<'_, PyAny>,
    cond_a: Vec<usize>,
    cond_b: Vec<usize>,
    sseq_params: &Bound<'_, PyDict>,
    big_count: Option<u64>,
    threads: usize,
) -> PyResult<Bound<'a, PyDict>> {
    init_global_thread_pool(threads);
    let converted_matrix = convert_matrix(py, matrix);
    // converting _bool doesn't seem to work automatically,
    let use_g = sseq_params.get_item("use_g")?.expect("No use_g found.");
    let new_use_g: Vec<bool> = use_g.downcast::<PyArray1<bool>>()?.to_vec()?;
    let params = SSeqParams {
        num_cells: sseq_params
            .get_item("N")
            .unwrap()
            .expect("N not found in sseq params")
            .extract()
            .unwrap(),
        num_genes: sseq_params
            .get_item("G")
            .unwrap()
            .expect("G not found in sseq params")
            .extract()
            .unwrap(),
        size_factors: sseq_params
            .get_item("size_factors")
            .unwrap()
            .expect("size_factors not found in sseq params")
            .extract()
            .unwrap(),
        gene_means: sseq_params
            .get_item("mean_g")
            .unwrap()
            .expect("mean_g not found in sseq params")
            .extract()
            .unwrap(),
        gene_variances: sseq_params
            .get_item("var_g")
            .expect("var_g not found in sseq params")
            .unwrap()
            .extract()
            .unwrap(),
        use_genes: new_use_g,
        gene_moment_phi: sseq_params
            .get_item("phi_mm_g")
            .expect("phi_mm_g not found in sseq params")
            .unwrap()
            .extract()
            .unwrap(),
        zeta_hat: sseq_params
            .get_item("zeta_hat")
            .expect("zeta_hat not found in sseq params")
            .unwrap()
            .extract()
            .unwrap(),
        delta: sseq_params
            .get_item("delta")
            .unwrap()
            .expect("delta not found in sseq params")
            .extract()
            .unwrap(),
        gene_phi: sseq_params
            .get_item("phi_g")
            .unwrap()
            .expect("phi_g not found in sseq params")
            .extract()
            .unwrap(),
    };
    let result: DiffExpResult =
        sseq_differential_expression(&converted_matrix, &cond_a, &cond_b, &params, big_count);
    drop(converted_matrix);
    let converted_result = PyDict::new(py);
    converted_result.set_item::<&str, _>(
        "genes_tested",
        PyList::new(py, result.genes_tested.iter().map(|x| PyBool::new(py, *x)))?,
    )?;
    converted_result.set_item::<&str, _>("sums_in", result.sums_in.into_pyarray(py))?;
    converted_result.set_item::<&str, _>("sums_out", result.sums_out.into_pyarray(py))?;
    converted_result.set_item::<&str, _>("common_mean", result.common_mean.into_pyarray(py))?;
    converted_result.set_item::<&str, _>(
        "common_dispersion",
        result.common_dispersion.into_pyarray(py),
    )?;
    converted_result.set_item::<&str, _>(
        "normalized_mean_in",
        result.normalized_mean_in.into_pyarray(py),
    )?;
    converted_result.set_item::<&str, _>(
        "normalized_mean_out",
        result.normalized_mean_out.into_pyarray(py),
    )?;
    converted_result.set_item::<&str, _>("p_values", result.p_values.into_pyarray(py))?;
    converted_result.set_item::<&str, _>(
        "adjusted_p_values",
        result.adjusted_p_values.into_pyarray(py),
    )?;
    converted_result
        .set_item::<&str, _>("log2_fold_change", result.log2_fold_change.into_pyarray(py))?;
    Ok(converted_result)
}
