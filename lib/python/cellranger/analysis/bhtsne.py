#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
import os.path

import numpy as np
import tsne as tsne_bh

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
from cellranger.analysis.analysis_types import TSNE


def run_tsne(
    transformed_pca_matrix,
    name: str = "TSNE",
    key: str = "TSNE",
    tsne_dims=None,
    input_pcs=None,
    perplexity=None,
    theta=None,
    max_iter=None,
    stop_lying_iter=None,
    mom_switch_iter=None,
    copy_data=False,
    random_state=None,
):
    if tsne_dims is None:
        tsne_dims = analysis_constants.TSNE_N_COMPONENTS

    if perplexity is None:
        perplexity = analysis_constants.TSNE_DEFAULT_PERPLEXITY

    if theta is None:
        theta = analysis_constants.TSNE_THETA

    if random_state is None:
        random_state = analysis_constants.RANDOM_STATE

    if max_iter is None:
        max_iter = analysis_constants.TSNE_MAX_ITER

    if stop_lying_iter is None:
        stop_lying_iter = analysis_constants.TSNE_STOP_LYING_ITER

    if mom_switch_iter is None:
        mom_switch_iter = analysis_constants.TSNE_MOM_SWITCH_ITER

    if input_pcs is not None:
        transformed_pca_matrix = transformed_pca_matrix[:, :input_pcs]

    # Make sure perplexity satisfies 'tsne' requirements
    N = transformed_pca_matrix.shape[0]
    perplexity = min(perplexity, max(1, -1 + float(N - 1) / 3))

    # At sufficiently low cell count, we cannot satisfy perplexity >= 1 and this condition below,
    # so the projection is defined as all zeros
    if N - 1 < 3 * perplexity:
        transformed_tsne_matrix = np.zeros((N, tsne_dims))
    else:
        transformed_tsne_matrix = tsne_bh.bh_sne(
            transformed_pca_matrix,
            d=tsne_dims,
            theta=theta,
            perplexity=perplexity,
            max_iter=max_iter,
            stop_lying_iter=stop_lying_iter,
            mom_switch_iter=mom_switch_iter,
            copy_data=copy_data,
            random_state=np.random.RandomState(random_state),
        )

    return TSNE(transformed_tsne_matrix, name=name, key=key)


def save_tsne_csv(tsne, barcodes, base_dir):
    """Save a TSNE object to CSV."""
    # Preserve backward compatibility with pre-3.0 CSV files
    #   where the CSV directory was named "2_components" and the HDF5 dataset was named "_2"
    key = tsne.key + "_components"

    tsne_dir = os.path.join(base_dir, key)
    os.makedirs(tsne_dir, exist_ok=True)

    matrix_fn = os.path.join(tsne_dir, "projection.csv")
    n_tsne_components = tsne.transformed_tsne_matrix.shape[1]
    matrix_header = ["Barcode"] + ["TSNE-%d" % (i + 1) for i in range(n_tsne_components)]
    analysis_io.save_matrix_csv(matrix_fn, tsne.transformed_tsne_matrix, matrix_header, barcodes)


def save_tsne_h5(tsne, fname):
    """Save a TSNE object to HDF5."""
    t_map = {tsne.key: tsne}
    analysis_io.save_dimension_reduction_h5(t_map, fname, analysis_constants.ANALYSIS_H5_TSNE_GROUP)
