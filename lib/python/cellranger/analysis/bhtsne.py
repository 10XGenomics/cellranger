#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

import collections
import numpy as np
import os
import tsne as tsne_bh

TSNE = collections.namedtuple('TSNE', ['transformed_tsne_matrix'])

def run_tsne(transformed_pca_matrix, tsne_dims=None, input_pcs=None, perplexity=None, theta=None,
             max_iter=None, stop_lying_iter=None, mom_switch_iter=None, copy_data=False, random_state=None):

    if tsne_dims is None:
        tsne_dims = cr_constants.TSNE_N_COMPONENTS

    if perplexity is None:
        perplexity = cr_constants.TSNE_DEFAULT_PERPLEXITY

    if theta is None:
        theta = cr_constants.TSNE_THETA

    if random_state is None:
        random_state = cr_constants.RANDOM_STATE

    if max_iter is None:
        max_iter = cr_constants.TSNE_MAX_ITER

    if stop_lying_iter is None:
        stop_lying_iter = cr_constants.TSNE_STOP_LYING_ITER

    if mom_switch_iter is None:
        mom_switch_iter = cr_constants.TSNE_MOM_SWITCH_ITER

    if input_pcs is not None:
        transformed_pca_matrix = transformed_pca_matrix[:, :input_pcs]

    # Make sure perplexity satisfies 'tsne' requirements
    N = transformed_pca_matrix.shape[0]
    perplexity = min(perplexity, max(1, -1 + float((N-1))/3))

    transformed_tsne_matrix = tsne_bh.bh_sne(transformed_pca_matrix,
                                             d = tsne_dims,
                                             theta = theta,
                                             perplexity = perplexity,
                                             max_iter = max_iter,
                                             stop_lying_iter = stop_lying_iter,
                                             mom_switch_iter = mom_switch_iter,
                                             copy_data = copy_data,
                                             random_state = np.random.RandomState(random_state))

    return TSNE(transformed_tsne_matrix)

def save_tsne_csv(tsne_map, matrix, base_dir):
    for n_tsne_components, tsne in tsne_map.iteritems():
        n_tsne_components_dir = os.path.join(base_dir, '%d_components' % n_tsne_components)
        cr_utils.makedirs(n_tsne_components_dir, allow_existing=True)

        matrix_fn = os.path.join(n_tsne_components_dir, 'projection.csv')
        matrix_header = ['Barcode'] + ['TSNE-%d' % (i+1) for i in xrange(n_tsne_components)]
        cr_io.save_matrix_csv(matrix_fn, tsne.transformed_tsne_matrix, matrix_header, matrix.bcs)

def save_tsne_h5(tsne_map, f):
    group = f.create_group(f.root, cr_constants.ANALYSIS_H5_TSNE_GROUP)
    for n_tsne_components, tsne in tsne_map.iteritems():
        cr_io.save_h5(f, group, str(n_tsne_components), tsne)
