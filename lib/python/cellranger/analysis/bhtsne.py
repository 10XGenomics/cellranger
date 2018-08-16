#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.io as cr_io

import collections
import numpy as np
import os
import tsne as tsne_bh

TSNE = collections.namedtuple('TSNE', ['transformed_tsne_matrix',
                                       'name', # Human readable form
                                       'key', # Machine queryable form, must be unique
])

def run_tsne(transformed_pca_matrix, name='TSNE', key='TSNE', tsne_dims=None, input_pcs=None, perplexity=None, theta=None,
             max_iter=None, stop_lying_iter=None, mom_switch_iter=None, copy_data=False, random_state=None):

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

    return TSNE(transformed_tsne_matrix, name=name, key=key)

def save_tsne_csv(tsne, matrix, base_dir):
    """Save a TSNE object to CSV"""
    # Preserve backward compatibility with pre-3.0 CSV files
    #   where the CSV directory was named "2_components" and the HDF5 dataset was named "_2"
    key = tsne.key + '_components'

    tsne_dir = os.path.join(base_dir, key)
    cr_io.makedirs(tsne_dir, allow_existing=True)

    matrix_fn = os.path.join(tsne_dir, 'projection.csv')
    n_tsne_components = tsne.transformed_tsne_matrix.shape[1]
    matrix_header = ['Barcode'] + ['TSNE-%d' % (i+1) for i in xrange(n_tsne_components)]
    analysis_io.save_matrix_csv(matrix_fn, tsne.transformed_tsne_matrix, matrix_header, matrix.bcs)

def save_tsne_h5(tsne, f):
    """Save a TSNE object to HDF5"""
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_TSNE_GROUP)
    analysis_io.save_h5(f, group, tsne.key, tsne)
