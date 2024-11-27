use anyhow::{bail, Context, Result};
use diff_exp::diff_exp::{DiffExpResult, SSeqParams};
use diff_exp::{compute_sseq_params, sseq_differential_expression};
use numpy::{Element, IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBool, PyDict, PyList};
use sprs::CsMatBase;
use sqz::mat::AdaptiveMatOwned;
use sqz::matrix_map::MatrixIntoMap;

/// Extract a python array of a desired type.
fn extract_pyarray<'a, T: Element>(
    py: &'a PyAny,
    attr: &'static str,
) -> Result<PyReadonlyArray1<'a, T>> {
    py.getattr(attr)
        .context(attr)
        .unwrap()
        .extract()
        .context(attr)
        .map(|x: &PyArray1<T>| x.readonly())
}

/// Transmute a &[i32] to a &[u32] and validate that all elements are positive.
fn transmute_i32_to_u32(xs: &[i32]) -> &[u32] {
    assert!(xs.iter().all(|&x| x >= 0));
    unsafe { std::slice::from_raw_parts(xs.as_ptr() as *const u32, xs.len()) }
}

/// Convert a scipy sparse matrix to a rust sparse matrix.
fn convert_matrix(
    _py: Python<'_>,
    matrix: &PyAny,
) -> Result<AdaptiveMatOwned<u32, MatrixIntoMap<u32, u32>>> {
    let shape: (usize, usize) = matrix.getattr("shape").unwrap().extract().unwrap();
    let data_py = extract_pyarray(matrix, "data").unwrap();
    let data = transmute_i32_to_u32(data_py.as_slice().unwrap());

    // indices and indptr may be either i32 or i64. They have the same type.
    let csc_matrix = if let Ok(indptr) = extract_pyarray::<i32>(matrix, "indptr") {
        let indices = extract_pyarray::<i32>(matrix, "indices").unwrap();
        AdaptiveMatOwned::from_csmat(&CsMatBase::new_csc(
            shape,
            indptr.as_slice().unwrap(),
            indices.as_slice().unwrap(),
            data,
        ))
    } else if let Ok(indptr) = extract_pyarray::<i64>(matrix, "indptr") {
        let indices = extract_pyarray::<i64>(matrix, "indices").unwrap();
        AdaptiveMatOwned::from_csmat(&CsMatBase::new_csc(
            shape,
            indptr.as_slice().unwrap(),
            indices.as_slice().unwrap(),
            data,
        ))
    } else {
        bail!("indptr is neither i32 nor i64")
    };
    Ok(csc_matrix)
}

#[pyfunction]
pub(crate) fn compute_sseq_params_o3<'a>(
    py: Python<'a>,
    matrix: &'a PyAny,
    zeta_quartile: f64,
) -> PyResult<&'a PyDict> {
    let converted_matrix = convert_matrix(py, matrix).unwrap();
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
pub(crate) fn sseq_differential_expression_o3<'a>(
    py: Python<'a>,
    matrix: &PyAny,
    cond_a: Vec<usize>,
    cond_b: Vec<usize>,
    sseq_params: &PyDict,
    big_count: Option<u64>,
) -> PyResult<&'a PyDict> {
    let converted_matrix = convert_matrix(py, matrix).unwrap();
    // converting _bool doesn't seem to work automatically,
    let use_g: &PyArray1<bool> = sseq_params.get_item("use_g").unwrap().downcast()?;
    let new_use_g: Vec<bool> = use_g.to_vec()?;
    let params = SSeqParams {
        num_cells: sseq_params.get_item("N").unwrap().extract().unwrap(),
        num_genes: sseq_params.get_item("G").unwrap().extract().unwrap(),
        size_factors: sseq_params
            .get_item("size_factors")
            .unwrap()
            .extract()
            .unwrap(),
        gene_means: sseq_params.get_item("mean_g").unwrap().extract().unwrap(),
        gene_variances: sseq_params.get_item("var_g").unwrap().extract().unwrap(),
        use_genes: new_use_g,
        gene_moment_phi: sseq_params.get_item("phi_mm_g").unwrap().extract().unwrap(),
        zeta_hat: sseq_params.get_item("zeta_hat").unwrap().extract().unwrap(),
        delta: sseq_params.get_item("delta").unwrap().extract().unwrap(),
        gene_phi: sseq_params.get_item("phi_g").unwrap().extract().unwrap(),
    };
    let result: DiffExpResult =
        sseq_differential_expression(&converted_matrix, &cond_a, &cond_b, &params, big_count);
    std::mem::drop(converted_matrix);
    let converted_result: &PyDict = PyDict::new(py);
    converted_result.set_item::<&str, &PyList>(
        "genes_tested",
        PyList::new(py, result.genes_tested.iter().map(|x| PyBool::new(py, *x))),
    )?;
    converted_result
        .set_item::<&str, &PyArray1<u64>>("sums_in", result.sums_in.into_pyarray(py))?;
    converted_result
        .set_item::<&str, &PyArray1<u64>>("sums_out", result.sums_out.into_pyarray(py))?;
    converted_result
        .set_item::<&str, &PyArray1<f64>>("common_mean", result.common_mean.into_pyarray(py))?;
    converted_result.set_item::<&str, &PyArray1<f64>>(
        "common_dispersion",
        result.common_dispersion.into_pyarray(py),
    )?;
    converted_result.set_item::<&str, &PyArray1<f64>>(
        "normalized_mean_in",
        result.normalized_mean_in.into_pyarray(py),
    )?;
    converted_result.set_item::<&str, &PyArray1<f64>>(
        "normalized_mean_out",
        result.normalized_mean_out.into_pyarray(py),
    )?;
    converted_result
        .set_item::<&str, &PyArray1<f64>>("p_values", result.p_values.into_pyarray(py))?;
    converted_result.set_item::<&str, &PyArray1<f64>>(
        "adjusted_p_values",
        result.adjusted_p_values.into_pyarray(py),
    )?;
    converted_result.set_item::<&str, &PyArray1<f64>>(
        "log2_fold_change",
        result.log2_fold_change.into_pyarray(py),
    )?;
    Ok(converted_result)
}
