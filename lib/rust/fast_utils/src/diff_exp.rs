use anyhow::Result;
use diff_exp::diff_exp::{DiffExpResult, SSeqParams};
use diff_exp::{compute_sseq_params, sseq_differential_expression};
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBool, PyDict, PyList};
//use scan_types::label_class::make_labelclass_from_feature_type_vector;
use sprs::CsMatI;
use sqz::mat::AdaptiveMatOwned;
use sqz::matrix_map::MatrixIntoMap;

type MatrixType = CsMatI<u32, u32>;

// fn push_attr_to_list(list: &mut Vec<String>, obj: &PyAny, attr: &str) {
//     list.push(obj.getattr(attr).unwrap().extract().unwrap());
// }
// Converts from a Python CountMatrix class to a Rust FeatureBarcodeMatrix struct
// fn convert_matrix(  _py: Python<'_>,
//     count_matrix: &PyAny) -> Result<AdaptiveFeatureBarcodeMatrix> {
//     let matrix = count_matrix.getattr("m").unwrap();
//     let indices: Vec<u32> = matrix.getattr("indices").unwrap().extract().unwrap();
//     let indptr: Vec<u32> = matrix.getattr("indptr").unwrap().extract().unwrap();
//     let data: Vec<u32> = matrix.getattr("data").unwrap().extract().unwrap();
//     let barcodes: Vec<Vec<u8>> = count_matrix.getattr("bcs").unwrap().extract().unwrap();
//     let barcodes : Vec<String> = barcodes.iter().map(|x| {String::from_utf8(x.to_vec()).unwrap()}).collect();
//     let feature_defs: &PyList = count_matrix.getattr("feature_ref").unwrap().getattr("feature_defs").unwrap().extract().unwrap();
//     let n_feats = feature_defs.len();
//     let mut feature_ids = Vec::<String>::with_capacity(n_feats);
//     let mut feature_names = Vec::<String>::with_capacity(n_feats);
//     let mut feature_types = Vec::<String>::with_capacity(n_feats);
//     for feature_def in feature_defs {
//         push_attr_to_list(&mut feature_ids, feature_def, "id");
//         push_attr_to_list(&mut feature_names, feature_def, "name");
//         push_attr_to_list(&mut feature_types, feature_defs, "feature_type");
//     }
//     let feature_type_label = make_labelclass_from_feature_type_vector(&feature_types)?;
//     let csc_mat = MatrixType::new_csc((feature_ids.len(), barcodes.len()), indptr, indices, data);
//     let mat = AdaptiveMatOwned::from_csmat(&csc_mat);
//     let final_matrix = AdaptiveFeatureBarcodeMatrix {
//         name: "Dummy".to_string(),
//         barcodes,
//         feature_ids,
//         feature_names,
//         feature_types: feature_type_label,
//         matrix: mat,
//     };
//     Ok(final_matrix)
// }

fn get_attr(matrix: &PyAny, attr: &str) -> Vec<u32> {
    let pyvalues: &PyArray1<i32> = matrix.getattr(attr).unwrap().extract().unwrap();
    pyvalues
        .to_vec()
        .unwrap()
        .into_iter()
        .map(|x| x as u32)
        .collect()
}

// Convert a scipy sparse matrix to a rust sparse matrix
fn convert_matrix(
    _py: Python<'_>,
    matrix: &PyAny,
) -> Result<AdaptiveMatOwned<u32, MatrixIntoMap<u32, u32>>> {
    let shape: (usize, usize) = matrix.getattr("shape").unwrap().extract().unwrap();
    // Avoid a simple cast from PyAny to vec which iteratively grabs a python object for each array
    // element
    let indices: Vec<u32> = get_attr(matrix, "indices");
    let indptr: Vec<u32> = get_attr(matrix, "indptr");
    let data: Vec<u32> = get_attr(matrix, "data");

    let csc_mat = MatrixType::new_csc(shape, indptr, indices, data);
    let final_matrix = AdaptiveMatOwned::from_csmat(&csc_mat);
    Ok(final_matrix)
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
