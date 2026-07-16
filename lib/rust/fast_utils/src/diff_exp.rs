#![deny(missing_docs)]

use anyhow::{Context, Result};
use diff_exp::diff_exp::{DiffExpResult, SSeqParams, sseq_de_from_sums, sseq_params_from_moments};
use diff_exp::stat::Statistics;
use diff_exp::{NbExactBackend, compute_sseq_params, sseq_differential_expression_backend};
use ndarray::Array1;
use numpy::{Element, IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBool, PyDict, PyList};
use rayon::ThreadPoolBuilder;
use sprs::CsMatBase;
use sqz::mat::AdaptiveMatOwned;
use sqz::matrix_map::{MatrixIntoMap, MatrixMap, TransposeMap};
use std::sync::{Arc, Once};

/// Default big_count threshold (mirrors diff_exp::diff_exp::BIG_COUNT_DEFAULT).
static BIG_COUNT_DEFAULT: u64 = 900;

/// MatrixMap that normalizes a raw count `v` by the per-cell total `count_total[c]`,
/// yielding the intrinsic per-cell-normalized value w = v / count_total[c] (independent of any
/// per-perturbation size-factor median). Mirrors diff_exp's private `SizeNormalized` but divides
/// by precomputed per-cell totals, so the control-side aggregates can be computed once and reused.
#[derive(Clone)]
struct CountNormalized {
    count_total: Arc<Vec<f64>>,
}

impl PartialEq for CountNormalized {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.count_total, &other.count_total)
    }
}
impl Eq for CountNormalized {}

impl MatrixMap<u32, f64> for CountNormalized {
    type T = TransposeMap<u32, f64, CountNormalized>;

    #[inline(always)]
    fn map(&self, v: u32, _r: usize, c: usize) -> f64 {
        let d = self.count_total[c];
        if d == 0.0 { 0.0 } else { v as f64 / d }
    }

    fn t(&self) -> Self::T {
        TransposeMap::new(self.clone())
    }
}

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
    let result: DiffExpResult = sseq_differential_expression_backend(
        &converted_matrix,
        &cond_a,
        &cond_b,
        &params,
        big_count,
        NbExactBackend::Ratio,
    );
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

/// Per-cell total counts for the cells `batch_perturbation_de` touches (control ∪ all
/// perturbations), indexed by full-matrix column. Cells not referenced stay 0. Extracted so the
/// reduce / `nb_de_from_stats` decomposition produces the same `count_total` as the monolithic path.
fn compute_count_total(
    mat: &sqz::AdaptiveMat,
    control: &[usize],
    perturbations: &[Vec<usize>],
) -> Arc<Vec<f64>> {
    let mut needed: Vec<usize> = control.to_vec();
    for a in perturbations {
        needed.extend_from_slice(a);
    }
    needed.sort_unstable();
    needed.dedup();
    let needed_totals = mat.sum_cols::<f64>(&needed);
    let mut ct = vec![0.0f64; mat.cols()];
    for (i, &c) in needed.iter().enumerate() {
        ct[c] = needed_totals[i];
    }
    Arc::new(ct)
}

/// Per-(cell group) sufficient statistics for the sSeq batched DE — everything
/// `nb_de_from_stats` needs about one group (the control, or one perturbation), reduced from the
/// matrix. Every field is an additive sum over the group's cells, so the same statistics can be
/// produced by a single streaming pass over the matrix (the Part-2 SPLIT redesign) instead of an
/// in-memory matrix.
struct GroupStats {
    /// Σ raw count per gene over the group's cells.
    feature_sums: Array1<u64>,
    /// Σ w per gene, with w = count/count_total. Computed as `mean·n` to match the matrix path.
    sum_w: Vec<f64>,
    /// Σ w² per gene. Computed as `(var + mean²)·n` to match the matrix path.
    sum_wsq: Vec<f64>,
    /// Σ count_total over the group's cells.
    sum_ct: f64,
    /// Σ 1/count_total over the group's cells (cells with zero total skipped).
    sum_inv: f64,
    /// Number of cells in the group.
    n: usize,
    /// Per-cell count_total in the group's cell order (for the `median_ab` over a ∪ control).
    count_totals: Vec<f64>,
}

/// Reduce one cell group to its sufficient statistics. Mirrors the per-group aggregate computation
/// in `batch_perturbation_de` exactly (`sum_w`/`sum_wsq` via `mean_var_rows`, the same `sum_ct` /
/// `sum_inv` folds), so `nb_de_from_stats` over these reproduces `batch_perturbation_de` bit-for-bit.
fn reduce_group(
    mat: &sqz::AdaptiveMat,
    count_total: &Arc<Vec<f64>>,
    cells: &[usize],
) -> GroupStats {
    let w_view = mat.view().compose_map(CountNormalized {
        count_total: count_total.clone(),
    });
    let nf = cells.len() as f64;
    let (mw, vw) = w_view.mean_var_rows(cells);
    let sum_w: Vec<f64> = mw.iter().map(|&m| m * nf).collect();
    let sum_wsq: Vec<f64> = mw
        .iter()
        .zip(vw.iter())
        .map(|(&m, &v)| (v + m * m) * nf)
        .collect();
    let feature_sums: Array1<u64> = mat.sum_rows::<u64>(cells);
    let sum_ct: f64 = cells.iter().map(|&c| count_total[c]).sum();
    let sum_inv: f64 = cells
        .iter()
        .filter(|&&c| count_total[c] != 0.0)
        .map(|&c| 1.0 / count_total[c])
        .sum();
    let count_totals: Vec<f64> = cells.iter().map(|&c| count_total[c]).collect();
    GroupStats {
        feature_sums,
        sum_w,
        sum_wsq,
        sum_ct,
        sum_inv,
        n: cells.len(),
        count_totals,
    }
}

/// sSeq DE params + exact/asymptotic test from precomputed per-group sufficient statistics — the
/// CRISPR shared-control batching half of `batch_perturbation_de` (the matrix is never touched
/// here). CR owns only the w-space combine, the `median_ab`-as-scalar factoring and the
/// control-aggregate reuse; the generic sSeq math (`sseq_params_from_moments` +
/// `sseq_de_from_sums`) lives in scan-rs's `diff-exp` crate and is shared with the matrix DE path.
/// `median_ab` factors out as a per-perturbation scalar, so the control aggregates are shared across
/// all perturbations. The streaming SPLIT reduction may reassociate `sum_w²`, which is the
/// (TOLERANT) Part-2 drift source vs the in-memory `reduce_group` fold.
fn nb_de_from_stats(
    control: &GroupStats,
    perts: &[GroupStats],
    n_genes: usize,
    zeta_quintile: f64,
    big_count: u64,
) -> Vec<(DiffExpResult, f64)> {
    let nb = control.n;
    perts
        .iter()
        .map(|a| {
            let na = a.n;
            let m = (na + nb) as f64;

            // median of count_total over a ∪ control (matches size_factors()'s median over cells).
            let mut cts: Vec<f64> = Vec::with_capacity(na + nb);
            cts.extend_from_slice(&a.count_totals);
            cts.extend_from_slice(&control.count_totals);
            let median_ab = cts.as_slice().median();
            let sum_size_factors = median_ab * (a.sum_inv + control.sum_inv);

            // Combine control + a per gene -> mean_g, var_g (norm = w * median_ab).
            let mut mean_g = Array1::<f64>::zeros(n_genes);
            let mut var_g = Array1::<f64>::zeros(n_genes);
            for g in 0..n_genes {
                let sum_w = control.sum_w[g] + a.sum_w[g];
                let sum_wsq = control.sum_wsq[g] + a.sum_wsq[g];
                let mean_w = sum_w / m;
                let var_w = sum_wsq / m - mean_w * mean_w;
                mean_g[g] = median_ab * mean_w;
                var_g[g] = median_ab * median_ab * var_w;
            }

            // Method-of-moments dispersion + shrinkage, then the per-gene NB test + BH + log2FC +
            // normalized means — both owned by scan-rs (`diff-exp` crate) so the math lives in one
            // place shared with the matrix DE path. CR keeps only the CRISPR shared-control
            // batching above (the w-space combine, median_ab scaling, control-aggregate reuse).
            let params = sseq_params_from_moments(
                mean_g.as_slice().unwrap(),
                var_g.as_slice().unwrap(),
                sum_size_factors,
                m,
                n_genes as f64,
                zeta_quintile,
            );
            // The batched path uses the transcendental-free ratio-recurrence exact-test backend
            // (scan-rs PR #788) — same conditional NB p-value as the log-space kernel (max|Δp|≈1e-9)
            // with a log-space deep-tail fallback; this is the MEASURE_PERTURBATIONS hotspot. Passing
            // `NbExactBackend::Ratio` preserves the pre-move behavior (the matrix API defaults to
            // `LogSpace`).
            let size_factor_a = a.sum_ct / median_ab;
            let size_factor_b = control.sum_ct / median_ab;
            let result = sseq_de_from_sums(
                a.feature_sums.as_slice().unwrap(),
                control.feature_sums.as_slice().unwrap(),
                size_factor_a,
                size_factor_b,
                &params,
                NbExactBackend::Ratio,
                big_count,
            );
            (result, median_ab)
        })
        .collect()
}

/// Batched sSeq DE (the production path): reduce the control + each perturbation to sufficient
/// statistics, then run the params + test over them. This is the in-memory composition of the
/// Part-2 split — `reduce_group` here, `nb_de_from_stats` there. Validated against scan-rs's
/// independent per-pair matrix DE path by `test_batched_matches_per_perturbation` (tolerant).
fn batch_perturbation_de(
    mat: &sqz::AdaptiveMat,
    control: &[usize],
    perturbations: &[Vec<usize>],
    zeta_quintile: f64,
    big_count: u64,
) -> Vec<(DiffExpResult, f64)> {
    let n_genes = mat.rows();
    let count_total = compute_count_total(mat, control, perturbations);
    let control_stats = reduce_group(mat, &count_total, control);
    let pert_stats: Vec<GroupStats> = perturbations
        .iter()
        .map(|a| reduce_group(mat, &count_total, a))
        .collect();
    nb_de_from_stats(
        &control_stats,
        &pert_stats,
        n_genes,
        zeta_quintile,
        big_count,
    )
}

/// Pack one perturbation's `(DiffExpResult, median_ab)` into the Python dict the MEASURE chunk
/// consumes. Shared by `batch_perturbation_de_o3` (matrix path) and `nb_de_from_stats_o3`
/// (sufficient-statistics path) so the two entries return byte-identical dict shapes.
fn de_result_to_pydict<'a>(
    py: Python<'a>,
    result: DiffExpResult,
    median_ab: f64,
) -> PyResult<Bound<'a, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("median_ab", median_ab)?;
    d.set_item(
        "genes_tested",
        PyList::new(py, result.genes_tested.iter().map(|x| PyBool::new(py, *x)))?,
    )?;
    d.set_item("sums_in", result.sums_in.into_pyarray(py))?;
    d.set_item("sums_out", result.sums_out.into_pyarray(py))?;
    d.set_item("common_mean", result.common_mean.into_pyarray(py))?;
    d.set_item(
        "common_dispersion",
        result.common_dispersion.into_pyarray(py),
    )?;
    d.set_item(
        "normalized_mean_in",
        result.normalized_mean_in.into_pyarray(py),
    )?;
    d.set_item(
        "normalized_mean_out",
        result.normalized_mean_out.into_pyarray(py),
    )?;
    d.set_item("p_values", result.p_values.into_pyarray(py))?;
    d.set_item(
        "adjusted_p_values",
        result.adjusted_p_values.into_pyarray(py),
    )?;
    d.set_item("log2_fold_change", result.log2_fold_change.into_pyarray(py))?;
    Ok(d)
}

/// PyO3 wrapper for `batch_perturbation_de`. Returns a list (one per perturbation, in input order)
/// of dicts with the same keys as `sseq_differential_expression_o3`.
#[pyfunction]
#[pyo3(signature = (matrix, control, perturbations, zeta_quintile, big_count=None, threads=1))]
pub(super) fn batch_perturbation_de_o3<'a>(
    py: Python<'a>,
    matrix: &Bound<'_, PyAny>,
    control: Vec<usize>,
    perturbations: Vec<Vec<usize>>,
    zeta_quintile: f64,
    big_count: Option<u64>,
    threads: usize,
) -> PyResult<Bound<'a, PyList>> {
    init_global_thread_pool(threads);
    let converted_matrix = convert_matrix(py, matrix);
    let big_count = big_count.unwrap_or(BIG_COUNT_DEFAULT);
    let results = batch_perturbation_de(
        &converted_matrix,
        &control,
        &perturbations,
        zeta_quintile,
        big_count,
    );
    drop(converted_matrix);
    let out = PyList::empty(py);
    for (result, median_ab) in results {
        out.append(de_result_to_pydict(py, result, median_ab)?)?;
    }
    Ok(out)
}

/// Build a `GroupStats` from the per-group arrays the MEASURE chunk reads out of SPLIT's
/// sufficient-statistics H5. `feature_sums`/`sum_w`/`sum_wsq` are `n_genes`-long.
fn group_stats_from_arrays(
    feature_sums: &[u64],
    sum_w: &[f64],
    sum_wsq: &[f64],
    sum_ct: f64,
    sum_inv: f64,
    n: usize,
    count_totals: Vec<f64>,
) -> GroupStats {
    // GroupStats owns its per-gene vectors (they outlive the borrowed numpy arrays and are consumed
    // by the parallel DE), so one copy out of the input is unavoidable here.
    GroupStats {
        feature_sums: Array1::from(feature_sums.to_vec()),
        sum_w: sum_w.to_vec(),
        sum_wsq: sum_wsq.to_vec(),
        sum_ct,
        sum_inv,
        n,
        count_totals,
    }
}

/// PyO3 wrapper for `nb_de_from_stats`: the sufficient-statistics DE path. Takes the control group's
/// statistics (computed once by SPLIT) and the chunk's perturbation statistics (stacked
/// `[n_perts, n_genes]` for the per-gene fields), and returns the same list of dicts as
/// `batch_perturbation_de_o3`. This replaces the per-chunk matrix `hstack` + `batch_perturbation_de_o3`
/// in the pre-split MEASURE chunk. Round-tripping `reduce_group`'s output through the array packing
/// `group_stats_from_arrays` performs is bit-for-bit (`test_group_stats_array_roundtrip_bit_identical`);
/// end-to-end against the matrix path the `feature_sums`, log2 fold change and bootstrap CIs stay
/// bit-identical, while the DE p-values are TOLERANT (~1e-11, well inside the rtol=1e-6 gate) — the
/// streaming reducer's `count_total` sums reassociate vs the old hstacked-matrix order.
#[pyfunction]
#[pyo3(signature = (
    n_genes,
    control_feature_sums, control_sum_w, control_sum_wsq, control_sum_ct, control_sum_inv,
    control_n, control_count_totals,
    pert_feature_sums, pert_sum_w, pert_sum_wsq, pert_sum_ct, pert_sum_inv, pert_n,
    pert_count_totals,
    zeta_quintile, big_count=None, threads=1
))]
#[allow(clippy::too_many_arguments)]
pub(super) fn nb_de_from_stats_o3<'a>(
    py: Python<'a>,
    n_genes: usize,
    control_feature_sums: PyReadonlyArray1<'a, u64>,
    control_sum_w: PyReadonlyArray1<'a, f64>,
    control_sum_wsq: PyReadonlyArray1<'a, f64>,
    control_sum_ct: f64,
    control_sum_inv: f64,
    control_n: usize,
    control_count_totals: PyReadonlyArray1<'a, f64>,
    pert_feature_sums: PyReadonlyArray2<'a, u64>,
    pert_sum_w: PyReadonlyArray2<'a, f64>,
    pert_sum_wsq: PyReadonlyArray2<'a, f64>,
    pert_sum_ct: Vec<f64>,
    pert_sum_inv: Vec<f64>,
    pert_n: Vec<usize>,
    pert_count_totals: Vec<Vec<f64>>,
    zeta_quintile: f64,
    big_count: Option<u64>,
    threads: usize,
) -> PyResult<Bound<'a, PyList>> {
    init_global_thread_pool(threads);
    let big_count = big_count.unwrap_or(BIG_COUNT_DEFAULT);

    let control = group_stats_from_arrays(
        control_feature_sums.as_slice()?,
        control_sum_w.as_slice()?,
        control_sum_wsq.as_slice()?,
        control_sum_ct,
        control_sum_inv,
        control_n,
        control_count_totals.as_slice()?.to_vec(),
    );

    let n_perts = pert_n.len();
    // The stacked `[n_perts, n_genes]` arrays are C-contiguous (fresh fancy-indexed numpy arrays), so
    // take each as one flat slice and sub-slice per row — no per-row view/copy. `as_slice()?` returns
    // a proper error instead of panicking if a caller ever passes a non-contiguous array.
    let fs = pert_feature_sums.as_slice()?;
    let sw = pert_sum_w.as_slice()?;
    let swsq = pert_sum_wsq.as_slice()?;
    let perts: Vec<GroupStats> = (0..n_perts)
        .map(|p| {
            let row = p * n_genes..(p + 1) * n_genes;
            group_stats_from_arrays(
                &fs[row.clone()],
                &sw[row.clone()],
                &swsq[row],
                pert_sum_ct[p],
                pert_sum_inv[p],
                pert_n[p],
                pert_count_totals[p].clone(),
            )
        })
        .collect();

    let results = nb_de_from_stats(&control, &perts, n_genes, zeta_quintile, big_count);
    let out = PyList::empty(py);
    for (result, median_ab) in results {
        out.append(de_result_to_pydict(py, result, median_ab)?)?;
    }
    Ok(out)
}

#[cfg(test)]
mod batch_tests {
    use super::*;
    use diff_exp::{
        compute_sseq_params as ref_compute_sseq_params, sseq_differential_expression as ref_sseq_de,
    };
    use sqz::mat::AdaptiveMatOwned;
    use sqz::matrix_map::MatrixIntoMap;

    /// Build a deterministic dense CSC matrix (genes × cells) with varied per-gene counts,
    /// mirroring `convert_matrix`'s construction (no rng dependency).
    fn build_mat(n_genes: usize, n_cells: usize) -> AdaptiveMatOwned<u32, MatrixIntoMap<u32, u32>> {
        let mut indptr: Vec<i32> = Vec::with_capacity(n_cells + 1);
        let mut indices: Vec<i32> = Vec::with_capacity(n_cells * n_genes);
        let mut data: Vec<u32> = Vec::with_capacity(n_cells * n_genes);
        indptr.push(0);
        for c in 0..n_cells {
            for g in 0..n_genes {
                indices.push(g as i32);
                // varied, deterministic, non-negative counts with non-zero per-gene variance
                data.push(((c * 5 + g * 7 + 3) % 11) as u32);
            }
            indptr.push((indices.len()) as i32);
        }
        AdaptiveMatOwned::from_csmat(&CsMatBase::new_csc(
            (n_genes, n_cells),
            indptr.as_slice(),
            indices.as_slice(),
            data.as_slice(),
        ))
    }

    /// Batched DE must match the per-perturbation compute_sseq_params + sseq_differential_expression
    /// (on the full matrix with cell_indices = a ∪ control) to within a tight float tolerance.
    #[test]
    fn test_batched_matches_per_perturbation() {
        let n_genes = 40usize;
        let n_cells = 200usize;
        let mat = build_mat(n_genes, n_cells);

        // control = first 120 cells; two perturbations of differing sizes.
        let control: Vec<usize> = (0..120).collect();
        let pert0: Vec<usize> = (120..160).collect();
        let pert1: Vec<usize> = (160..200).collect();
        let perturbations = vec![pert0.clone(), pert1.clone()];

        let batched = batch_perturbation_de(&mat, &control, &perturbations, 0.995, 900);

        for (a, (got, _median_ab)) in perturbations.iter().zip(batched.iter()) {
            let cells: Vec<usize> = a.iter().chain(control.iter()).copied().collect();
            let params = ref_compute_sseq_params(&mat, Some(0.995), Some(&cells), None);
            let want = ref_sseq_de(&mat, a, &control, &params, Some(900));

            for g in 0..n_genes {
                assert_eq!(got.sums_in[g], want.sums_in[g], "sums_in g={g}");
                assert_eq!(got.sums_out[g], want.sums_out[g], "sums_out g={g}");
                let close = |x: f64, y: f64, ctx: &str| {
                    if x.is_nan() && y.is_nan() {
                        return;
                    }
                    assert!(
                        (x - y).abs() <= 1e-9 + 1e-6 * y.abs(),
                        "{ctx} g={g}: {x} vs {y}"
                    );
                };
                close(got.common_mean[g], want.common_mean[g], "mean");
                close(got.common_dispersion[g], want.common_dispersion[g], "phi");
                close(got.log2_fold_change[g], want.log2_fold_change[g], "log2fc");
                close(got.p_values[g], want.p_values[g], "pval");
                close(
                    got.adjusted_p_values[g],
                    want.adjusted_p_values[g],
                    "adj_pval",
                );
            }
        }
    }

    /// `nb_de_from_stats_o3` rebuilds each `GroupStats` from flat arrays (the H5 roundtrip
    /// `group_stats_from_arrays` performs) before calling `nb_de_from_stats`. Round-tripping
    /// `reduce_group`'s output through that array packing must leave the DE result BIT-IDENTICAL —
    /// guards against a field mis-wiring (e.g. swapping `sum_w`/`sum_wsq`) in the PyO3 marshalling.
    #[test]
    fn test_group_stats_array_roundtrip_bit_identical() {
        let n_genes = 40usize;
        let n_cells = 200usize;
        let mat = build_mat(n_genes, n_cells);
        let control: Vec<usize> = (0..120).collect();
        let perturbations = vec![
            (120..160).collect::<Vec<usize>>(),
            (160..200).collect::<Vec<usize>>(),
        ];
        let count_total = compute_count_total(&mat, &control, &perturbations);

        let pack = |gs: &GroupStats| -> GroupStats {
            group_stats_from_arrays(
                gs.feature_sums.as_slice().unwrap(),
                &gs.sum_w,
                &gs.sum_wsq,
                gs.sum_ct,
                gs.sum_inv,
                gs.n,
                gs.count_totals.clone(),
            )
        };
        let control_stats = reduce_group(&mat, &count_total, &control);
        let pert_stats: Vec<GroupStats> = perturbations
            .iter()
            .map(|a| reduce_group(&mat, &count_total, a))
            .collect();

        let direct = nb_de_from_stats(&control_stats, &pert_stats, n_genes, 0.995, 900);
        let packed_control = pack(&control_stats);
        let packed_perts: Vec<GroupStats> = pert_stats.iter().map(pack).collect();
        let roundtrip = nb_de_from_stats(&packed_control, &packed_perts, n_genes, 0.995, 900);

        assert_eq!(direct.len(), roundtrip.len());
        for ((d, dm), (r, rm)) in direct.iter().zip(roundtrip.iter()) {
            assert_eq!(dm, rm, "median_ab");
            for g in 0..n_genes {
                assert_eq!(d.p_values[g], r.p_values[g], "pval g={g}");
                assert_eq!(d.sums_in[g], r.sums_in[g], "sums_in g={g}");
                assert_eq!(d.log2_fold_change[g], r.log2_fold_change[g], "log2fc g={g}");
            }
        }
    }
}
