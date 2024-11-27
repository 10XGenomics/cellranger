#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#


from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING

import numpy as np
from six import ensure_str
from sklearn.utils import sparsefuncs

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
import cellranger.analysis.stats as analysis_stats
import cellranger.h5_constants as h5_constants
from cellranger.analysis.analysis_types import PCA
from cellranger.analysis.irlb import irlb

if TYPE_CHECKING:
    from cellranger.feature_ref import FeatureDef
    from cellranger.matrix import CountMatrix


# The RUNPCA stage attempts to run the PCA at this threshold, and if that
# fails it reruns at zero.  In the event thresholding prevents us from
# returning the requested number of components and we are at this threshold
# value, we throw an exception.
DEFAULT_RUNPCA_THRESHOLD = 2


class MatrixRankTooSmallException(Exception):
    pass


def get_original_columns_used(cols_not_removed, cols_used_after_removal):
    """If a matrix is subset down to only have columns indexed by cols_not_removed, and then is further subset to.

    only contain cols_used_after removal, in that order, than this method returns the index of which columns in the old
    matrix correspond the the columns in the new matrix.
    """
    return [cols_not_removed[x] for x in cols_used_after_removal]


def run_pca(
    matrix: CountMatrix,
    pca_features: int | None = None,
    pca_bcs: int | None = None,
    n_pca_components: int | None = None,
    random_state: int | None = None,
    min_count_threshold: int = 0,
):
    """Run a PCA on the matrix using the IRLBA matrix factorization algorithm.

    Prior to the PCA analysis, the
    matrix is modified so that all barcodes/columns have the same counts, and then the counts are transformed
    by a log2(1+X) operation.

    If desired, only a subset of features (e.g. sample rows) can be selected for PCA analysis.  Each feature is ranked
    by its dispersion relative to other features that have a similar mean count.  The top `pca_features` as ranked by
    this method will then be used for the PCA.

    One can also select to subset number of barcodes to use (e.g. sample columns), but in this case they are simply
    randomly sampled.

    Args:
        matrix (CountMatrix): The matrix to perform PCA on.
        pca_features (int): Number of features to subset from matrix and use in PCA. The top pca_features ranked by
                            dispersion are used
        pca_bcs (int): Number of barcodes to randomly sample for the matrix.
        n_pca_components (int): How many PCA components should be used.
        random_state (int): The seed for the RNG
        min_count_threshold (int): The minimum sum of each row/column for that row/column to be passed to PCA
                                   (this filter is prior to any subsetting that occurs).

    Returns:
        A PCA object
    """
    if random_state is None:
        random_state = analysis_constants.RANDOM_STATE
    np.random.seed(0)

    # Threshold the rows/columns of matrix, will throw error if an empty matrix results.
    thresholded_matrix, _, thresholded_features = matrix.select_axes_above_threshold(
        min_count_threshold
    )

    # If requested, we can subsample some of the barcodes to get a smaller matrix for PCA
    pca_bc_indices = np.arange(thresholded_matrix.bcs_dim)
    if pca_bcs is None:
        pca_bcs = thresholded_matrix.bcs_dim
        pca_bc_indices = np.arange(thresholded_matrix.bcs_dim)
    elif pca_bcs < thresholded_matrix.bcs_dim:
        pca_bc_indices = np.sort(
            np.random.choice(np.arange(thresholded_matrix.bcs_dim), size=pca_bcs, replace=False)
        )
    elif pca_bcs > thresholded_matrix.bcs_dim:
        msg = (
            f"You requested {pca_bcs} barcodes but the matrix after thresholding only "
            f"included {thresholded_matrix.bcs_dim}, so the smaller amount is being used."
        )
        print(msg)
        pca_bcs = thresholded_matrix.bcs_dim
        pca_bc_indices = np.arange(thresholded_matrix.bcs_dim)

    # If requested, select fewer features to use by selecting the features with highest normalized dispersion
    if pca_features is None:
        pca_features = thresholded_matrix.features_dim
    elif pca_features > thresholded_matrix.features_dim:
        msg = (
            f"You requested {pca_features} features but the matrix after thresholding only included {thresholded_matrix.features_dim} features,"
            "so the smaller amount is being used."
        )
        print(msg)
        pca_features = thresholded_matrix.features_dim
    # Calc mean and variance of counts after normalizing
    # But don't transform to log space, in order to preserve the mean-variance relationship
    m = analysis_stats.normalize_by_umi(thresholded_matrix)
    # Get mean and variance of rows
    (mu, var) = analysis_stats.summarize_columns(m.T)
    dispersion = analysis_stats.get_normalized_dispersion(
        mu.squeeze(), var.squeeze()
    )  # TODO set number of bins?
    # pylint: disable=invalid-unary-operand-type
    pca_feature_indices = np.argsort(dispersion, kind="stable")[-pca_features:]

    # Now determine how many components.
    if n_pca_components is None:
        n_pca_components = analysis_constants.PCA_N_COMPONENTS_DEFAULT
    likely_matrix_rank = min(pca_features, pca_bcs)
    if likely_matrix_rank < n_pca_components:
        if min_count_threshold == DEFAULT_RUNPCA_THRESHOLD:
            # Kick back to run_pca stage so it can retry with no threshold, this is for historical reasons
            raise MatrixRankTooSmallException("Matrix rank is too small")
        else:
            print(
                f"There are fewer nonzero features or barcodes ({likely_matrix_rank}) than requested "
                f"PCA components ({n_pca_components}); reducing the number of components."
            )
            n_pca_components = likely_matrix_rank

    if (likely_matrix_rank * 0.5) <= float(n_pca_components):
        print(
            "Requested number of PCA components is large relative to the matrix "
            "size, an exact approach to matrix factorization may be faster."
        )

    # Note, after subsetting it is possible some rows/cols in pca_mat have counts below the threshold.
    # However, we are not performing a second thresholding as in practice subsetting is not used and we explain
    # that thresholding occurs prior to subsetting in the doc string.
    pca_mat = thresholded_matrix.select_barcodes(pca_bc_indices).select_features(
        pca_feature_indices
    )
    (pca_norm_mat, pca_center, pca_scale) = normalize_and_transpose(pca_mat)
    # we scale here, instead of using irlb's scaling, which is occasionally broken
    sparsefuncs.inplace_column_scale(pca_norm_mat, 1.0 / np.atleast_1d(pca_scale.squeeze()))
    pca_center /= pca_scale
    (_, d, v, _, _) = irlb(
        pca_norm_mat,
        n_pca_components,
        center=np.atleast_1d(pca_center.squeeze()),
        random_state=random_state,
    )

    # make sure to project the matrix before centering, to avoid densification
    (full_norm_mat, full_center, full_scale) = normalize_and_transpose(matrix)
    sparsefuncs.inplace_column_scale(
        full_norm_mat, 1.0 / np.atleast_1d(full_scale.squeeze())
    )  # can have some zeros here
    # Get a coordinate map so we know which columns in the old matrix correspond to columns in the new
    org_cols_used = get_original_columns_used(thresholded_features, pca_feature_indices)
    transformed_irlba_matrix = full_norm_mat[:, org_cols_used].dot(v) - (full_center / full_scale)[
        :, org_cols_used
    ].dot(v)
    irlba_components = np.zeros((n_pca_components, matrix.features_dim))
    irlba_components[:, org_cols_used] = v.T

    # calc proportion of variance explained
    variance_sum = len(
        pca_feature_indices
    )  # each feature has variance=1, mean=0 after normalization
    variance_explained: np.ndarray = np.square(d) / ((len(pca_bc_indices) - 1) * variance_sum)
    features_selected: np.ndarray = np.array([f.id for f in matrix.feature_ref.feature_defs])[
        org_cols_used
    ]

    # Now project back up the dispersion to return.
    full_dispersion = np.empty(matrix.features_dim)
    full_dispersion[:] = np.nan
    full_dispersion[thresholded_features] = dispersion

    # sanity check dimensions
    assert transformed_irlba_matrix.shape == (matrix.bcs_dim, n_pca_components)
    assert irlba_components.shape == (n_pca_components, matrix.features_dim)
    assert variance_explained.shape == (n_pca_components,)

    return PCA(
        transformed_irlba_matrix,
        irlba_components,
        variance_explained,
        full_dispersion,
        features_selected,
    )


def normalize_and_transpose(matrix: CountMatrix):
    matrix.tocsc()

    m = analysis_stats.normalize_by_umi(matrix)

    # Use log counts
    m.data = np.log2(1 + m.data)

    # Transpose
    m = m.T

    # compute centering (mean) and scaling (stdev)
    (c, v) = analysis_stats.summarize_columns(m)
    # TODO: Inputs to this function shouldn't have zero variance columns
    v[np.where(v == 0.0)] = 1.0

    s = np.sqrt(v)
    return (m, c, s)


# IRLBA adds 20 'temporary' PCs during processing
EXTRA_IRLBA_PCS = 20
# Overhead of the matrix per feature and per barcode
BYTES_PER_FEATURE = 1e3
BYTES_PER_BC = 1e3


def get_irlb_mem_gb_from_matrix_dim(nonzero_entries, features, bcs, pcs):
    """An approximate model of the memory consumption of PCA preprocessing plus IRLBA.

    The key factor is analysis_constants.NUM_IRLB_MATRIX_ENTRIES_PER_MEM_GB
    which is set empirically from looking at the memory consumption of many jobs.

    See lib/python/cellranger/test/test_pca.py for details & a test.
    """
    irlba_mem_gb = (
        nonzero_entries / analysis_constants.NUM_IRLB_MATRIX_ENTRIES_PER_MEM_GB
        + features * (pcs + EXTRA_IRLBA_PCS) * 8 / 1e9
        + features * BYTES_PER_FEATURE / 1e9
        + bcs * BYTES_PER_BC / 1e9
        + bcs * (pcs + 20) * 8 / 1e9
        + analysis_constants.IRLB_BASE_MEM_GB
    )

    return max(h5_constants.MIN_MEM_GB, np.ceil(irlba_mem_gb))


def save_pca_csv(pca_map: Mapping[int, PCA], matrix: CountMatrix, base_dir: str):
    save_pca_csv_with_bc_feature(pca_map, matrix.bcs, matrix.feature_ref.feature_defs, base_dir)


def save_pca2_csv(pca_map: Mapping[int, PCA], barcodes, base_dir: str, library_type: str):
    """Used only for saving the results of PCA2 in the batch correction pipeline."""
    for n_components, pca in pca_map.items():
        n_components_dir = os.path.join(base_dir, f"{library_type}_{n_components}_components")
        os.makedirs(n_components_dir, exist_ok=True)

        matrix_fn = os.path.join(n_components_dir, "projection.csv")
        n_columns = pca.transformed_pca_matrix.shape[1]
        assert n_columns <= n_components
        matrix_header = ["Barcode"] + ["PC-%d" % (i + 1) for i in range(n_columns)]
        analysis_io.save_matrix_csv(matrix_fn, pca.transformed_pca_matrix, matrix_header, barcodes)


def save_pca_csv_with_bc_feature(
    pca_map: Mapping[int, PCA],
    barcodes: Iterable[bytes],
    features: Iterable[FeatureDef],
    base_dir: str,
):
    """Used only for saving the results of python-based PCA in the ATAC pipeline."""
    for n_components, pca in pca_map.items():
        n_components_dir = os.path.join(base_dir, f"{n_components}_components")
        os.makedirs(n_components_dir, exist_ok=True)

        matrix_fn = os.path.join(n_components_dir, "projection.csv")
        n_columns = pca.transformed_pca_matrix.shape[1]
        assert n_columns <= n_components
        matrix_header = ["Barcode"] + ["PC-%d" % (i + 1) for i in range(n_columns)]
        analysis_io.save_matrix_csv(matrix_fn, pca.transformed_pca_matrix, matrix_header, barcodes)

        # FBPCA presently provides 0-sized entries for the following PCA() member variables.
        #   This allows us to distinguish FBPCA from IRLBA, and also avoids weird empty files.
        if pca.components.size > 0:
            components_fn = os.path.join(n_components_dir, "components.csv")
            components_header = ["PC"] + [ensure_str(f.id) for f in features]
            analysis_io.save_matrix_csv(
                components_fn, pca.components, components_header, range(1, n_components + 1)
            )

        if pca.variance_explained.size > 0:
            variance_fn = os.path.join(n_components_dir, "variance.csv")
            variance_header = ["PC", "Proportion.Variance.Explained"]
            analysis_io.save_matrix_csv(
                variance_fn, pca.variance_explained, variance_header, range(1, n_components + 1)
            )

        if pca.dispersion.size > 0:
            dispersion_fn = os.path.join(n_components_dir, "dispersion.csv")
            dispersion_header = ["Feature", "Normalized.Dispersion"]
            analysis_io.save_matrix_csv(
                dispersion_fn, pca.dispersion, dispersion_header, [f.id for f in features]
            )

        if pca.features_selected.size > 0:
            features_fn = os.path.join(n_components_dir, "features_selected.csv")
            # TODO: there are two columns here, but only 1 entry in the header...BAD
            features_header = ["Feature"]
            analysis_io.save_matrix_csv(
                features_fn,
                pca.features_selected,
                features_header,
                range(1, len(pca.features_selected) + 1),
            )


def save_pca_h5(pca_map: Mapping[int, PCA], fname: str):
    analysis_io.save_dimension_reduction_h5(
        pca_map, fname, analysis_constants.ANALYSIS_H5_PCA_GROUP
    )


def save_pca2_h5(pca_map: Mapping[int, PCA], fname: str, library_type: str):
    analysis_io.save_pca2_dimension_reduction_h5(
        pca_map,
        fname,
        analysis_constants.ANALYSIS_H5_PCA_GROUP,
        library_type,
    )


def load_pca_from_h5(filename: str):
    """Load just the PCA info from an analysis h5."""
    return analysis_io.load_dimension_reduction_from_h5(
        filename, analysis_constants.ANALYSIS_H5_PCA_GROUP, PCA
    )
