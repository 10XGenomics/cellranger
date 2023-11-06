#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Functions to compute various statistical operations on sparse matrices."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.stats
import sklearn.utils.sparsefuncs as sparsefuncs

if TYPE_CHECKING:
    from scipy.sparse import csc_matrix, csr_matrix

    from cellranger.matrix import CountMatrix


def normalize_by_umi(matrix: CountMatrix) -> csc_matrix | csr_matrix:
    counts_per_bc = matrix.get_counts_per_bc()
    median_counts_per_bc = max(1.0, np.median(counts_per_bc))
    scaling_factors = median_counts_per_bc / counts_per_bc

    # Normalize each barcode's total count by median total count
    m = matrix.m.copy().astype(np.float64)
    sparsefuncs.inplace_column_scale(m, scaling_factors)

    return m


def normalize_by_idf(matrix: CountMatrix) -> csc_matrix | csr_matrix:
    """Perform feature normalization."""
    numbcs_per_feature = matrix.get_numbcs_per_feature()
    scaling_factors_row = np.log(matrix.bcs_dim + 1) - np.log(1 + numbcs_per_feature)

    m = matrix.m.copy().astype(np.float64)
    sparsefuncs.inplace_row_scale(m, scaling_factors_row)

    # Extremely Rare Case (1 out of 1000s of samples tested):
    # Either the scaling or the count may be zero for all features for some barcode
    # This would lead to zero-ing out entire barcode upon normalization, which leads to a null
    # projection as well. This is harmful to analysis code that depends on at least a non-zero norm
    # for each barcode (e.g. spherical clustering and normalized tsne). We sprinkle in a small
    # value that ensures an nnz for the all-zero barcode, after finding such barcodes.

    # find zeroed barcodes and assign nnz to first feature (these barcodes are indistinguishable
    # anyway). We run the very small risk of making it similar to another barcode that is also nnz
    # in the first feature only
    zeroed = np.where(np.squeeze(np.asarray(m.sum(axis=0))) == 0)
    for bc_ix in zeroed:
        m[0, bc_ix] = 1e-15

    return m


def summarize_columns(matrix: csc_matrix | csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    """Calculate mean and variance of each column, in a sparsity-preserving way."""
    mu, var = sparsefuncs.mean_variance_axis(matrix, axis=0)
    return np.array([mu]), np.array([var])


def get_normalized_dispersion(
    mat_mean: np.ndarray, mat_var: np.ndarray, nbins: int = 20
) -> np.ndarray:
    """Calculates the normalized dispersion.

    The dispersion is calculated for each feature
    and then normalized to see how its dispersion compares to samples that had a
    similar mean value.
    """
    # See equation in https://academic.oup.com/nar/article/40/10/4288/2411520
    # If a negative binomial is parameterized with mean m, and variance = m + d * m^2
    # then this d = dispersion as calculated below
    mat_disp: np.ndarray = (mat_var - mat_mean) / np.square(mat_mean)

    quantiles = np.percentile(mat_mean, np.arange(0, 100, 100 // nbins))
    quantiles = np.append(quantiles, mat_mean.max())

    # merge bins with no difference in value
    quantiles = np.unique(quantiles)

    if len(quantiles) <= 1:
        # pathological case: the means are all identical. just return raw dispersion.
        return mat_disp

    # calc median dispersion per bin
    (disp_meds, _, disp_bins) = scipy.stats.binned_statistic(
        mat_mean, mat_disp, statistic="median", bins=quantiles
    )

    # calc median absolute deviation of dispersion per bin
    disp_meds_arr = disp_meds[disp_bins - 1]  # 0th bin is empty since our quantiles start from 0
    disp_abs_dev = abs(mat_disp - disp_meds_arr)
    (disp_mads, _, disp_bins) = scipy.stats.binned_statistic(
        mat_mean, disp_abs_dev, statistic="median", bins=quantiles
    )

    # calculate normalized dispersion
    disp_mads_arr = disp_mads[disp_bins - 1]
    disp_norm = (mat_disp - disp_meds_arr) / disp_mads_arr
    return disp_norm
