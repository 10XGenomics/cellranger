//! PCA analysis code
#![deny(missing_docs)]

use crate::types::PcaResult;
use anyhow::{Error, Result};
use cr_types::FeatureBarcodeType;
use cr_types::reference::feature_reference::FeatureType;
use log::warn;
use ndarray::linalg::Dot;
use ndarray::{Array, Array1, Array2, Axis, s};
use ndarray_stats::Quantile1dExt;
use ndarray_stats::interpolate::Linear;
use noisy_float::types::{N64, n64};
use scan_rs::dim_red::Pca;
use scan_rs::dim_red::bk_svd::BkSvd;
use scan_rs::normalization::{
    FixedPointFormat, LogBase, Normalization, log1p_normalize_fixed_point,
    normalize_with_size_factor,
};
use scan_rs::stats::median_mut;
use sqz::{AdaptiveMat, AdaptiveVec, LowRankOffset, MatrixMap};
use std::cmp::Ordering;
use std::iter::zip;
use std::ops::Deref;

const DISPERSION_BINS: usize = 20;

fn binned_median(
    binned: &Array1<f64>,
    values: &Array1<f64>,
    bin_edges: &[N64],
) -> Result<(Array1<f64>, Vec<usize>)> {
    let max_index = bin_edges.len() - 1;
    let mut bin_indices = Vec::with_capacity(values.len());
    let mut binned_values = vec![Vec::new(); max_index];
    for (i, &v) in binned.iter().enumerate() {
        let idx = match bin_edges[1..max_index].binary_search(&n64(v)) {
            Ok(i) => i + 1,
            Err(i) => i,
        };
        bin_indices.push(idx);
        binned_values[idx].push(n64(values[i]));
    }
    let mut medians = Vec::with_capacity(max_index);
    for vec in binned_values {
        let mut arr = Array1::from(vec);
        let median = median_mut(&mut arr).map_or(f64::NAN, N64::raw);
        medians.push(median);
    }
    Ok((Array1::from(medians), bin_indices))
}

fn unselect(len: usize, vals: &Array1<f64>, idxs: &[usize]) -> Array1<f64> {
    let mut ret = Array1::from_elem(len, f64::NAN);
    for (&idx, &val) in zip(idxs, vals) {
        ret[idx] = val;
    }
    ret
}

fn compute_normalized_dispersion(
    matrix: &AdaptiveMat,
    pca_threshold: f64,
) -> Result<(Array1<f64>, AdaptiveMat, Vec<usize>)> {
    let (submatrix, _, features, _) = matrix.view().partition_on_threshold(pca_threshold);
    let inv_scale_factors = {
        let counts_per_bc = submatrix.sum_axis::<f64>(Axis(0));
        let median_count = median_mut(&mut counts_per_bc.mapv(n64)).map_or(f64::NAN, N64::raw);
        Some(counts_per_bc / median_count)
    };
    let (mean, var) = submatrix
        .view()
        .scale(Axis(0), inv_scale_factors)
        .mean_var_axis(Axis(1));
    let dispersion = (&var - &mean) / mean.mapv(|v| v * v);
    let qs = Array::linspace(0.0, 1.0, DISPERSION_BINS + 1).mapv(n64);
    let mut quantiles = mean
        .mapv(n64)
        .quantiles_mut(&qs, &Linear)
        .map(|xs| {
            let (xs, offset) = xs.into_raw_vec_and_offset();
            assert_eq!(offset, Some(0));
            xs
        })
        .unwrap_or_default();
    quantiles.dedup();
    if quantiles.len() <= 1 {
        return Ok((
            unselect(matrix.rows(), &dispersion, &features),
            submatrix,
            features,
        ));
    }
    let (medians, bin_indices) = binned_median(&mean, &dispersion, &quantiles)?;
    let medians = Array1::from(
        bin_indices
            .into_iter()
            .map(|i| medians[i])
            .collect::<Vec<_>>(),
    );
    let deviations = (&dispersion - &medians).mapv_into(f64::abs);
    let (mads, bin_indices) = binned_median(&mean, &deviations, &quantiles)?;
    let mads = Array1::from(bin_indices.into_iter().map(|i| mads[i]).collect::<Vec<_>>());
    let dispersion = (&dispersion - &medians) / &mads;
    Ok((
        unselect(matrix.rows(), &dispersion, &features),
        submatrix,
        features,
    ))
}

fn select_features(
    matrix: &AdaptiveMat,
    max_features: usize,
    pca_threshold: f64,
) -> Result<(AdaptiveMat, Array1<f64>, Vec<usize>)> {
    let (dispersion, submatrix, features) = compute_normalized_dispersion(matrix, pca_threshold)?;
    // sort in reverse, such than NANs are at the end
    let mut feature_indices = features.into_iter().enumerate().collect::<Vec<_>>();
    feature_indices.sort_by(|&(_, a), &(_, b)| {
        let da = dispersion[a];
        let db = dispersion[b];
        if da.is_nan() && db.is_nan() {
            Ordering::Equal
        } else if da.is_nan() {
            Ordering::Less
        } else if db.is_nan() {
            Ordering::Greater
        } else {
            da.partial_cmp(&db).unwrap()
        }
        .reverse()
    });
    feature_indices.truncate(max_features);
    let mut submatrix_features = Vec::with_capacity(feature_indices.len());
    let mut features = Vec::with_capacity(feature_indices.len());
    for &(j, i) in &feature_indices {
        submatrix_features.push(j);
        features.push(i);
    }
    let filtered_matrix = submatrix.select_rows(&submatrix_features);
    Ok((filtered_matrix, dispersion, features))
}

type BkSvdPcaResult = (Array2<f64>, Array1<f64>, Array2<f64>);

fn run_pca_on_filtered_norm_matrix<D>(
    matrix: LowRankOffset<D, impl MatrixMap<u32, f64>>,
    num_pcs: usize,
) -> Result<(BkSvdPcaResult, usize), Error>
where
    D: Deref<Target = [AdaptiveVec]>,
{
    let [num_features, num_cells] = matrix.shape();
    let min_dim = std::cmp::min(num_features, num_cells);
    let num_pcs = if min_dim < num_pcs {
        warn!(
            "matrix shape {:?} < requested PCs {}, reducing to {}",
            matrix.shape(),
            num_pcs,
            min_dim
        );
        min_dim
    } else {
        num_pcs
    };
    match BkSvd::new().run_pca(&matrix, num_pcs) {
        Ok(outs) => Ok((outs, num_pcs)),
        Err(e) => Err(e),
    }
}

fn get_components_and_transform_matrix<D>(
    matrix: LowRankOffset<D, impl MatrixMap<u32, f64>>,
    u: Array2<f64>,
    num_pcs: usize,
    selected_features: &[usize],
) -> (Array2<f64>, Array2<f64>)
where
    D: Deref<Target = [AdaptiveVec]>,
{
    let mut components = Array2::from_elem((num_pcs, matrix.rows()), 0.0);
    for (j, &i) in selected_features.iter().enumerate() {
        components.slice_mut(s![.., i]).assign(&u.slice(s![j, ..]));
    }
    let transformed_pca_matrix = matrix.t().dot(&components.t()); // (bcs, feat) x (feat, pcs)
    (components, transformed_pca_matrix)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_pca<'a>(
    matrix: &AdaptiveMat,
    feature_ids: &'a [String],
    feature_type: FeatureType,
    max_features: usize,
    num_pcs: usize,
    threshold: f64,
    is_spatial: bool,
    fixed_point: Option<FixedPointFormat>,
) -> Result<PcaResult<'a>> {
    let (filtered_matrix, dispersion, selected_features) =
        select_features(matrix, max_features, threshold)?;

    let selected_feature_ids = selected_features
        .iter()
        .rev()
        .map(|&i| feature_ids[i].as_str())
        .collect::<Vec<_>>();

    let ((u, sigma, _), adjusted_num_pcs) = match (feature_type, is_spatial) {
        (FeatureType::Barcode(FeatureBarcodeType::Antibody), true) => {
            run_pca_on_filtered_norm_matrix(
                normalize_with_size_factor(filtered_matrix, Normalization::LogTransform, None),
                num_pcs,
            )?
        }
        (FeatureType::ProteinExpression, false) | (FeatureType::ProteinExpression, true) => {
            match fixed_point {
                Some(fp) => run_pca_on_filtered_norm_matrix(
                    log1p_normalize_fixed_point(filtered_matrix, LogBase::Two, fp),
                    num_pcs,
                )?,
                None => panic!("missing fixed point details for protein expression"),
            }
        }
        _ => run_pca_on_filtered_norm_matrix(
            normalize_with_size_factor(filtered_matrix, Normalization::CellRanger, None),
            num_pcs,
        )?,
    };

    let (components, transformed_pca_matrix) = match (feature_type, is_spatial) {
        (FeatureType::Barcode(FeatureBarcodeType::Antibody), true) => {
            get_components_and_transform_matrix(
                normalize_with_size_factor(matrix.view(), Normalization::LogTransform, None),
                u,
                adjusted_num_pcs,
                &selected_features,
            )
        }
        (FeatureType::ProteinExpression, false) | (FeatureType::ProteinExpression, true) => {
            match fixed_point {
                Some(fp) => get_components_and_transform_matrix(
                    log1p_normalize_fixed_point(matrix.view(), LogBase::Two, fp),
                    u,
                    adjusted_num_pcs,
                    &selected_features,
                ),
                None => panic!("missing fixed point details for protein expression"),
            }
        }
        _ => get_components_and_transform_matrix(
            normalize_with_size_factor(matrix.view(), Normalization::CellRanger, None),
            u,
            adjusted_num_pcs,
            &selected_features,
        ),
    };
    let num_bcs = transformed_pca_matrix.shape()[0];
    let num_features = selected_features.len();
    let variance_explained =
        sigma.mapv_into(|v| v * v / ((num_bcs - 1) as f64 * num_features as f64));

    Ok(PcaResult::new(
        components,
        dispersion,
        feature_type,
        selected_feature_ids,
        transformed_pca_matrix,
        variance_explained,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use sqz::AdaptiveMatOwned;

    #[test]
    fn test_get_normalized_matrix_antibody_spatial() {
        // # Python code to reproduce test
        // mat = np.array([[136, 936, 0, 0, 264],
        //         [134, 682, 417, 8, 391],
        //         [0, 133, 780, 0, 0],
        //         [396, 76, 96, 198, 0],
        //  ])
        // almost_processed_mat = np.log2(1 + mat)
        // centering_factor = almost_processed_mat.mean(axis = 1).reshape((4,1))
        // scaling_factor = 1/np.std(almost_processed_mat, axis=1)
        // norm_mat = np.diag(scaling_factor).dot(almost_processed_mat - centering_factor)
        // norm_mat
        // array([[ 0.50075509,  1.16407001, -1.1965938 , -1.1965938 ,  0.72836249],
        //     [-0.14245194,  0.89844192,  0.58318993, -1.88113806,  0.54195815],
        //     [-0.80111703,  0.89623633,  1.50711477, -0.80111703, -0.80111703],
        //     [ 0.92609909,  0.14507504,  0.25503138,  0.59722303, -1.92342854]])

        let dense: Array2<u32> = array![
            [136, 936, 0, 0, 264],
            [134, 682, 417, 8, 391],
            [0, 133, 780, 0, 0],
            [396, 76, 96, 198, 0],
        ];
        let expected_out_dense = array![
            [0.50075509, 1.16407001, -1.1965938, -1.1965938, 0.72836249],
            [-0.14245194, 0.89844192, 0.58318993, -1.88113806, 0.54195815],
            [
                -0.80111703,
                0.89623633,
                1.50711477,
                -0.80111703,
                -0.80111703
            ],
            [0.92609909, 0.14507504, 0.25503138, 0.59722303, -1.92342854]
        ];
        let mtx = AdaptiveMatOwned::<u32>::from_dense(dense.view());
        let norm_mat = normalize_with_size_factor(mtx, Normalization::LogTransform, None);

        assert!(expected_out_dense.abs_diff_eq(&norm_mat.to_dense(), 1e-6));
    }

    #[test]
    fn test_get_normalized_matrix_gex_spatial() {
        // # Python code to reconstruct this test
        // mat = np.array([[136, 936, 0, 0, 264],
        //     [134, 682, 417, 8, 391],
        //     [0, 133, 780, 0, 0],
        //     [396, 76, 96, 198, 0],
        //         ])
        // scale_factor = mat.sum(axis=0)
        // target_umi_count = np.median(mat.sum(axis=0))
        // half_processed_mat = mat.dot(np.diag(target_umi_count/scale_factor))
        // almost_processed_mat = np.log2(1 + half_processed_mat)
        // centering_factor = almost_processed_mat.mean(axis = 1).reshape((4,1))
        // scaling_factor = 1/np.std(almost_processed_mat, axis=1)
        // norm_mat = np.diag(scaling_factor).dot(almost_processed_mat - centering_factor)
        // norm_mat
        // array([[ 0.61392149,  0.95459951, -1.21707302, -1.21707302,  0.86562504],
        //    [-0.11878431,  0.54279925,  0.38607315, -1.85660965,  1.04652156],
        //    [-0.78758751,  0.76437149,  1.59839105, -0.78758751, -0.78758751],
        //    [ 0.88718256, -0.25584717, -0.01048423,  1.09574143, -1.71659259]])

        let dense: Array2<u32> = array![
            [136, 936, 0, 0, 264],
            [134, 682, 417, 8, 391],
            [0, 133, 780, 0, 0],
            [396, 76, 96, 198, 0],
        ];
        let expected_out_dense = array![
            [0.61392149, 0.95459951, -1.21707302, -1.21707302, 0.86562504],
            [-0.11878431, 0.54279925, 0.38607315, -1.85660965, 1.04652156],
            [
                -0.78758751,
                0.76437149,
                1.59839105,
                -0.78758751,
                -0.78758751
            ],
            [
                0.88718256,
                -0.25584717,
                -0.01048423,
                1.09574143,
                -1.71659259
            ]
        ];
        let mtx = AdaptiveMatOwned::<u32>::from_dense(dense.view());
        let norm_mat = normalize_with_size_factor(mtx, Normalization::CellRanger, None);

        assert!(expected_out_dense.abs_diff_eq(&norm_mat.to_dense(), 1e-6));
    }
}
