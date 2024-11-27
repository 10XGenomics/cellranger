#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
"""Differential expression analysis for single-cell RNA-seq."""

from __future__ import annotations

import os
import sys
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scipy.sparse

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
from cellranger.analysis.analysis_types import DifferentialExpression
from cellranger.fast_utils import (  # pylint: disable=no-name-in-module, invalid-name
    compute_sseq_params_o3,
    sseq_differential_expression_o3,
)

if TYPE_CHECKING:
    import scipy.stats

SSEQ_ZETA_QUANTILE = 0.995


def estimate_size_factors(x: scipy.sparse.csc_matrix) -> np.ndarray[int, np.dtype[np.float64]]:
    """Estimate size factors (related to cell RNA content and GEM-to-GEM technical variance).

    Args:
      x: Sparse matrix (csc) of counts (feature x cell)

    Returns:
      Array of floats, one per cell.
    """
    counts_per_cell = np.squeeze(np.asarray(x.sum(axis=0)))
    size_factors = counts_per_cell.astype(np.float64) / np.median(counts_per_cell)
    return size_factors


def compute_sseq_params(x: scipy.sparse.csc_matrix, zeta_quantile=SSEQ_ZETA_QUANTILE):
    """Compute global parameters for the sSeq differential expression method.

    The key parameters are the shrunken feature-wise dispersions.

    This method was published in:
    Yu D, et al. (2013) Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq experiments with small sample size.
    Bioinformatics. 29: 1275-1282. doi: 10.1093/bioinformatics/btt143

    Args:
      x: Sparse matrix (csc) of counts (feature x cell)
      zeta_quantile (float): Quantile of method-of-moments dispersion estimates to
                             use as the shrinkage target zeta.

    Returns:
      A dictionary containing the sSeq parameters and some diagnostic info.
    """
    # pylint: disable=invalid-name
    if not x.has_sorted_indices:
        x = x.sorted_indices()
    return compute_sseq_params_o3(x, zeta_quantile)


def get_local_sseq_params(x: scipy.sparse.csc_matrix, group_a, group_b):
    """Compute locally-distinguishing sseq parameters.

    For perturbation vs control analysis for CRISPR
    and meta sample comparison within an aggregated matrix.
    """
    print("...Computing params for this comparison...")
    sys.stdout.flush()
    both_conditions = np.concatenate([group_a, group_b])
    matrix_groups: scipy.sparse.csc_matrix = x[:, both_conditions]
    if not matrix_groups.has_sorted_indices:
        matrix_groups = matrix_groups.sorted_indices()
    nbc = len(both_conditions)

    new_group_a = np.fromiter(range(len(group_a)), dtype=int, count=len(group_a))
    new_group_b = np.fromiter(range(len(group_a), nbc), dtype=int, count=nbc - len(group_a))
    return compute_sseq_params(matrix_groups), new_group_a, new_group_b, matrix_groups


def adjust_pvalue_bh(p):
    """Multiple testing correction of p-values using the Benjamini-Hochberg procedure."""
    descending = np.argsort(p)[::-1]
    # q = p * N / k where p = p-value, N = # tests, k = p-value rank
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    # pylint: disable=no-member
    qual = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

    # Return to original order
    return qual[np.argsort(descending)]


def sseq_differential_expression(x, cond_a, cond_b, sseq_params, big_count=900):
    """Run sSeq pairwise differential expression test.

    Args:
      x: Sparse matrix (csc) of counts (feature x cell)
      cond_a (np.array(int)): Indices of cells in group A
      cond_b (np.array(int)): Indices of cells in group B
      sseq_params (dict): Precomputed global parameters
      big_count (int): Use asymptotic approximation if both counts > this

    Returns:
      pd.DataFrame: DE results for group A relative to group B.
    """
    # Scipy is not guaranteed to have sorted indices, make a copy if they are not
    if not x.has_sorted_indices:
        x = x.sorted_indices()
    diff_exp_results = sseq_differential_expression_o3(
        x, cond_a, cond_b, sseq_params, big_count
    )  # pylint: disable=invalid-name
    de_result = pd.DataFrame(
        {
            "tested": sseq_params["use_g"],
            "sum_a": diff_exp_results["sums_in"],
            "sum_b": diff_exp_results["sums_out"],
            "common_mean": diff_exp_results["common_mean"],
            "common_dispersion": diff_exp_results["common_dispersion"],
            "norm_mean_a": diff_exp_results["normalized_mean_in"],
            "norm_mean_b": diff_exp_results["normalized_mean_out"],
            "p_value": diff_exp_results["p_values"],
            "adjusted_p_value": diff_exp_results["adjusted_p_values"],
            # Introduce a pseudocount into log2(fold_change)
            "log2_fold_change": diff_exp_results["log2_fold_change"],
        }
    )
    return de_result


def run_differential_expression(
    matrix,
    clusters: np.ndarray[Any, np.dtype[np.integer]],
    sseq_params: dict[str, Any] | None = None,
):
    """Compute differential expression for each cluster vs all other cells.

    Args:
        matrix (GeneBCMatrix):  feature expression data
        clusters (np.ndarray[int, int]):  1-based cluster labels
        sseq_params (dict):  params from compute_sseq_params
    """
    n_clusters = np.max(clusters)

    if sseq_params is None:
        print("Computing params...")
        sys.stdout.flush()
        sseq_params = compute_sseq_params(matrix.m)

    # Create a numpy array with 3*K columns;
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((matrix.features_dim, 3 * n_clusters))

    for cluster in range(1, 1 + n_clusters):
        in_cluster = clusters == cluster
        group_a = np.flatnonzero(in_cluster)
        group_b = np.flatnonzero(np.logical_not(in_cluster))
        print("Computing DE for cluster %d..." % cluster)
        sys.stdout.flush()

        de_result = sseq_differential_expression(matrix.m, group_a, group_b, sseq_params)
        all_de_results[:, 0 + 3 * (cluster - 1)] = de_result["norm_mean_a"]
        all_de_results[:, 1 + 3 * (cluster - 1)] = de_result["log2_fold_change"]
        all_de_results[:, 2 + 3 * (cluster - 1)] = de_result["adjusted_p_value"]

    return DifferentialExpression(all_de_results)


def save_differential_expression_csv(
    clustering_key,
    de,  # pylint: disable=invalid-name
    matrix,
    base_dir,
    cluster_names=None,
    file_name="differential_expression",
    cell_types=None,
):
    """Write diffexp results to CSV."""
    out_dir = base_dir
    if clustering_key is not None:
        out_dir = os.path.join(base_dir, clustering_key)
    os.makedirs(out_dir, exist_ok=True)

    diff_expression_fn = os.path.join(out_dir, file_name + ".csv")
    diff_expression_header = ["Feature ID", "Feature Name"]

    n_clusters = de.data.shape[1] // 3
    for i in range(n_clusters):
        if cell_types:
            diff_expression_header += [
                f"{cell_types[i]}, Mean Counts",
                f"{cell_types[i]}, Log2 fold change",
                f"{cell_types[i]}, Adjusted p value",
            ]
        elif cluster_names is None:
            diff_expression_header += [
                "Cluster %d Mean Counts" % (i + 1),
                "Cluster %d Log2 fold change" % (i + 1),
                "Cluster %d Adjusted p value" % (i + 1),
            ]
        else:
            diff_expression_header += [
                f"Perturbation {cluster_names[i]}, Mean Counts",
                f"Perturbation {cluster_names[i]}, Log2 fold change",
                f"Perturbation {cluster_names[i]}, Adjusted p value",
            ]

    diff_expression_prefixes = [(f.id, f.name) for f in matrix.feature_ref.feature_defs]
    analysis_io.save_matrix_csv(
        diff_expression_fn, de.data, diff_expression_header, diff_expression_prefixes
    )


def save_differential_expression_csv_from_features(
    clustering_key, de, diff_expression_prefixes, base_dir  # pylint: disable=invalid-name
):
    """Write diffexp results to CSV."""
    out_dir = os.path.join(base_dir, clustering_key)
    os.makedirs(out_dir, exist_ok=True)

    diff_expression_fn = os.path.join(out_dir, "differential_expression.csv")
    diff_expression_header = ["Feature ID", "Feature Name"]

    n_clusters = de.data.shape[1] // 3
    for i in range(n_clusters):
        diff_expression_header += [
            "Cluster %d Mean Counts" % (i + 1),
            "Cluster %d Log2 fold change" % (i + 1),
            "Cluster %d Adjusted p value" % (i + 1),
        ]

    analysis_io.save_matrix_csv(
        diff_expression_fn, de.data, diff_expression_header, diff_expression_prefixes
    )


def save_differential_expression_h5(f, clustering_key, de):
    """Write diffexp results to H5File `f`."""
    # pylint: disable=invalid-name
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP)

    analysis_io.save_h5(f, group, clustering_key, de)

    cr_clustering.create_legacy_kmeans_nodes(
        f,
        analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP,
        analysis_constants.ANALYSIS_H5_KMEANS_DIFFERENTIAL_EXPRESSION_GROUP,
        DifferentialExpression,
        clustering_key,
    )
