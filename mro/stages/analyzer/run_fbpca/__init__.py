#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
import numpy as np
import cPickle
import martian

import sklearn.preprocessing as sk_preprocessing
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_util
import cellranger.h5_constants as h5_constants
from cellranger.library_constants import GENE_EXPRESSION_LIBRARY_TYPE
import cellranger.analysis.constants as analysis_constants

from fbpca.fbpca import pca

__MRO__  = """
stage RUN_FBPCA(
    in  h5     matrix_h5,
    in  map[]  library_info,
    in  int    num_pcs,
    in  bool   skip,
    out pickle dimred_matrix,
    out pickle matrix_barcode_feature_info,
    src py     "stages/analyzer/run_fbpca",
) split (
)
"""

MAX_MEM_GB = 64

def fbpca_reduce_dimension(matrix, dimred):
    """ Fast Randomized PCA """
    X = matrix.m.T
    k = min((dimred, X.shape[0], X.shape[1]))
    U, s, Vt = pca(X, k=k) # Automatically centers.
    return U[:, range(k)] * s[range(k)]

def split(args):
    if args.skip:
        return {'chunks': []}

    # memory usage to load a h5 matrix
    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.matrix_h5)
    (features_dim, bcs_dim, nonzero_entries) = matrix_dims
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_dim(bcs_dim, nonzero_entries)

    mem_gb = min(MAX_MEM_GB, max(matrix_mem_gb, h5_constants.MIN_MEM_GB)) 
    return {'chunks': [], 'join': {'__mem_gb': mem_gb}}

def main(args, outs):
    return

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        return

    gg_id_to_batch_id, batch_id_to_name = {}, {}

    for lib in args.library_info:
        gg_id_to_batch_id[lib['gem_group']] = lib['batch_id']
        batch_id_to_name[lib['batch_id']] = lib['batch_name']

    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
    matrix = matrix.select_features_by_type(GENE_EXPRESSION_LIBRARY_TYPE)

    batch_ids = np.array([gg_id_to_batch_id[cr_util.split_barcode_seq(bc)[1]] for bc in matrix.bcs])

    # select intersect of non-zero feature in each batch
    feature_mask = np.ones(matrix.features_dim)
    for b_id in batch_id_to_name:
        batch_bc_indices = np.where(batch_ids == b_id)[0]
        matrix_view = cr_matrix.CountMatrixView(matrix, bc_indices=batch_bc_indices)
        feature_mask = np.logical_and(feature_mask, matrix_view.sum(axis=1))

    matrix = matrix.select_features(np.flatnonzero(feature_mask))

    # filter barcodes with zero count
    bc_indices = np.flatnonzero(matrix.get_counts_per_bc())
    matrix = matrix.select_barcodes(bc_indices)

    # l2 norm
    matrix.m = sk_preprocessing.normalize(matrix.m, axis=0)

    n_pcs = args.num_pcs if args.num_pcs is not None else analysis_constants.CBC_N_COMPONENTS_DEFAULT
    dimred_matrix = fbpca_reduce_dimension(matrix, n_pcs)

    outs.dimred_matrix = martian.make_path('dimred_matrix.pickle')
    with open(outs.dimred_matrix, 'wb') as fp:
        cPickle.dump(dimred_matrix, fp, cPickle.HIGHEST_PROTOCOL)

    bc_feature_info = {
        'barcodes' : matrix.bcs,
        'features' : matrix.feature_ref.feature_defs,
    }
    outs.matrix_barcode_feature_info = martian.make_path('matrix_barcode_feature_info.pickle')
    with open(outs.matrix_barcode_feature_info, 'wb') as fp:
        cPickle.dump(bc_feature_info, fp, cPickle.HIGHEST_PROTOCOL)