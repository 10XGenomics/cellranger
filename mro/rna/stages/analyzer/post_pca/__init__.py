#!/usr/bin/env python
#
# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#
"""Write pickle files for downstream consumption from the dimension reduced and full matrix."""

import pickle

import martian
import numpy as np

from cellranger.matrix import CountMatrix

__MRO__ = """
stage POST_PCA(
    in  h5     matrix_h5,
    in  h5     filt_matrix_h5,
    in  npy    dimred_matrix_npy,
    out pickle dimred_matrix,
    out pickle matrix_barcode_feature_info,
    src py     "stages/analyzer/post_pca",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    (num_features, num_barcodes, nnz) = CountMatrix.load_dims_from_h5(args.matrix_h5)
    matrix_mem_gib = CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=1)
    mem_gib = 3 + 1.1 * matrix_mem_gib
    print(f"{num_features=},{num_barcodes=},{nnz=},{matrix_mem_gib=},{mem_gib=}")
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(args, outs, chunk_defs, chunk_outs):
    all_bcs_matrix = CountMatrix.load_h5_file(args.matrix_h5)
    matrix = CountMatrix.load_h5_file(args.filt_matrix_h5)
    dimred_matrix = np.load(args.dimred_matrix_npy)

    # restore the zero count entries to the dimred matrix
    assert matrix.get_shape()[1] == dimred_matrix.shape[0]  # number of barcodes
    restored_matrix = np.zeros((all_bcs_matrix.get_shape()[1], dimred_matrix.shape[1]), dtype=float)

    for i, bc in enumerate(matrix.bcs):
        restored_matrix[all_bcs_matrix.bc_to_int(bc), :] = dimred_matrix[i, :]

    bc_feature_info = {
        "barcodes": all_bcs_matrix.bcs,
        "features": all_bcs_matrix.feature_ref.feature_defs,
    }

    outs.dimred_matrix = martian.make_path("dimred_matrix.pickle")
    with open(outs.dimred_matrix, "wb") as f_p:
        pickle.dump(restored_matrix, f_p, pickle.HIGHEST_PROTOCOL)

    outs.matrix_barcode_feature_info = martian.make_path("matrix_barcode_feature_info.pickle")
    with open(outs.matrix_barcode_feature_info, "wb") as f_p:
        pickle.dump(bc_feature_info, f_p, pickle.HIGHEST_PROTOCOL)
