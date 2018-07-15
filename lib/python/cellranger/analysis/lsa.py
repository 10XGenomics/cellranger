#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.analysis.stats as analysis_stats

import collections
from irlb import irlb
import numpy as np
import os
import tables

LSA = collections.namedtuple('LSA', ['transformed_lsa_matrix', 'components', 'variance_explained', 'dispersion', 'features_selected'])

def run_lsa(matrix, lsa_features=None, lsa_bcs=None, n_lsa_components=None, random_state=None):
    if lsa_features is None:
        lsa_features = matrix.features_dim
    if lsa_bcs is None:
        lsa_bcs = matrix.bcs_dim
    if n_lsa_components is None:
        n_lsa_components = analysis_constants.LSA_N_COMPONENTS_DEFAULT
        if n_lsa_components > lsa_features:
            print "There are fewer nonzero features than LSA components; reducing the number of components."
            n_lsa_components = lsa_features
    if random_state is None:
        random_state=analysis_constants.RANDOM_STATE

    np.random.seed(0)

    # perform idf transform, which is suited for lsa
    full_norm_mat  = normalize_and_transpose(matrix)

    # initialize LSA subsets
    lsa_bc_indices = np.arange(matrix.bcs_dim)
    lsa_feature_indices = np.arange(matrix.features_dim)

    # Calc mean and variance of counts after normalizing
    # Don't transform to log space in LSA
    # Dispersion is not exactly meaningful after idf transform. This is retained simply to follow PCA code
    m = analysis_stats.normalize_by_idf(matrix)
    (mu, var) = analysis_stats.summarize_columns(m.T)
    dispersion = analysis_stats.get_normalized_dispersion(mu.squeeze(), var.squeeze()) # TODO set number of bins?

    lsa_feature_indices = np.argsort(dispersion)[-lsa_features:]

    if lsa_bcs < matrix.bcs_dim:
        lsa_bc_indices = np.sort(np.random.choice(np.arange(matrix.bcs_dim), size=lsa_bcs, replace=False))

    lsa_mat, _, lsa_features_nonzero = matrix.select_barcodes(lsa_bc_indices).select_features(lsa_feature_indices).select_nonzero_axes()
    lsa_feature_nonzero_indices = lsa_feature_indices[lsa_features_nonzero]

    if lsa_mat.features_dim < 2 or lsa_mat.bcs_dim < 2:
        print "Matrix is too small for further downsampling - num_lsa_bcs and num_lsa_features will be ignored."
        lsa_mat, _, lsa_features_nonzero = matrix.select_nonzero_axes()
        lsa_feature_nonzero_indices = lsa_features_nonzero

    lsa_norm_mat = normalize_and_transpose(lsa_mat)

    (u, d, v, _, _) = irlb(lsa_norm_mat, n_lsa_components, random_state=random_state)

    # project the matrix to complete the transform: X --> X*v = u*d
    transformed_irlba_matrix = full_norm_mat[:,lsa_feature_nonzero_indices].dot(v) 
    irlba_components = np.zeros((n_lsa_components, matrix.features_dim))
    irlba_components[:,lsa_feature_nonzero_indices] = v.T

    # calc proportion of variance explained
    variance_explained = np.square(d) / np.sum(lsa_norm_mat.data**2)

    features_selected = np.array([f.id for f in matrix.feature_ref.feature_defs])[lsa_feature_nonzero_indices]

    # sanity check dimensions
    assert transformed_irlba_matrix.shape == (matrix.bcs_dim, n_lsa_components)
    assert irlba_components.shape == (n_lsa_components, matrix.features_dim)
    assert variance_explained.shape == (n_lsa_components,)

    return LSA(transformed_irlba_matrix, irlba_components, variance_explained, dispersion, features_selected)

def normalize_and_transpose(matrix):
    matrix.tocsc()

    m = analysis_stats.normalize_by_idf(matrix)

    # Transpose
    m = m.T

    return m

def get_irlb_mem_gb_from_matrix_dim(nonzero_entries):
    irlba_mem_gb = round(np.ceil(1.0 * nonzero_entries / analysis_constants.NUM_IRLB_MATRIX_ENTRIES_PER_MEM_GB)) + analysis_constants.IRLB_BASE_MEM_GB
    return h5_constants.MATRIX_MEM_GB_MULTIPLIER * max(h5_constants.MIN_MEM_GB, irlba_mem_gb)

def save_lsa_csv(lsa_map, matrix, base_dir):
    for n_components, lsa in lsa_map.iteritems():
        n_components_dir = os.path.join(base_dir, '%d_components' % n_components)
        cr_io.makedirs(n_components_dir, allow_existing=True)

        matrix_fn = os.path.join(n_components_dir, 'projection.csv')
        matrix_header = ['Barcode'] + ['PC-%d' % (i+1) for i in xrange(n_components)]
        analysis_io.save_matrix_csv(matrix_fn, lsa.transformed_lsa_matrix, matrix_header,
                              matrix.bcs)

        components_fn = os.path.join(n_components_dir, 'components.csv')
        components_header = ['PC'] + [f.id for f in matrix.feature_ref.feature_defs]
        analysis_io.save_matrix_csv(components_fn, lsa.components, components_header,
                              range(1, n_components+1))

        variance_fn = os.path.join(n_components_dir, 'variance.csv')
        variance_header = ['PC','Proportion.Variance.Explained']
        analysis_io.save_matrix_csv(variance_fn, lsa.variance_explained, variance_header,
                              range(1, n_components+1))

        dispersion_fn = os.path.join(n_components_dir, 'dispersion.csv')
        dispersion_header = ['Feature','Normalized.Dispersion']
        analysis_io.save_matrix_csv(dispersion_fn, lsa.dispersion, dispersion_header,
                              [f.id for f in matrix.feature_ref.feature_defs])

        features_fn = os.path.join(n_components_dir, 'features_selected.csv')
        features_header = ['Feature']
        analysis_io.save_matrix_csv(features_fn, lsa.features_selected, features_header, range(1, len(lsa.features_selected)+1))

def save_lsa_h5(lsa_map, f):
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_LSA_GROUP)
    for n_components, lsa in lsa_map.iteritems():
        analysis_io.save_h5(f, group, str(n_components), lsa)

def load_lsa_from_h5(filename):
    """ Load just the LSA info from an analysis h5 """
    with tables.open_file(filename, 'r') as f:
        group = f.root._v_groups[analysis_constants.ANALYSIS_H5_LSA_GROUP]
        # Just take the first LSA object, assuming we never have multiple
        for _, lsa in analysis_io.load_h5_iter(group, LSA):
            return lsa
