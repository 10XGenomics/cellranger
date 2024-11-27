#!/usr/bin/env python
#
# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#
"""Write pickle files for downstream consumption from the dimension reduced and full matrix.

files
"""
import math
import pickle

import martian
import numpy as np

import cellranger.matrix as cr_matrix

__MRO__ = """
stage POST_PCA(
    in  h5     matrix_h5,
    in  h5     filt_matrix_h5,
    in  npy    dimred_matrix_npy,
    in  float  mem_gb,
    out pickle dimred_matrix,
    out pickle matrix_barcode_feature_info,
    src py     "stages/analyzer/post_pca",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    mem_needed = math.ceil(3 + args.mem_gb * 0.14)
    return {
        "chunks": [],
        "join": {
            "__mem_gb": mem_needed,
            "__vmem_gb": mem_needed + 10,
        },
    }


def join(args, outs, chunk_defs, chunk_outs):
    all_bcs_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
    matrix = cr_matrix.CountMatrix.load_h5_file(args.filt_matrix_h5)
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
