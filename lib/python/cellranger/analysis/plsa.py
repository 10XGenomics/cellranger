#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.analysis.stats as analysis_stats

import tenkit.log_subprocess as tk_subproc

import subprocess
import collections
import numpy as np
import os
import tables

import scipy.io as sp_io

PLSA = collections.namedtuple('PLSA', ['transformed_plsa_matrix', 'components', 'variance_explained', 'dispersion', 'features_selected'])
PLSA_BINPATH='plsa'

def run_plsa(matrix, temp_dir, plsa_features=None, plsa_bcs=None, n_plsa_components=None, random_state=None, threads=1):
    if not os.path.exists(temp_dir):
        raise Exception('Temporary directory does not exist. Need it to run plsa binary. Aborting..')

    if plsa_features is None:
        plsa_features = matrix.features_dim
    if plsa_bcs is None:
        plsa_bcs = matrix.bcs_dim
    if n_plsa_components is None:
        n_plsa_components = analysis_constants.PLSA_N_COMPONENTS_DEFAULT
        if n_plsa_components > plsa_features:
            print "There are fewer nonzero features than PLSA components; reducing the number of components."
            n_plsa_components = plsa_features
    if random_state is None:
        random_state=analysis_constants.RANDOM_STATE

    np.random.seed(random_state)

    # initialize PLSA subsets
    plsa_bc_indices = np.arange(matrix.bcs_dim)
    plsa_feature_indices = np.arange(matrix.features_dim)

    # NOTE: This is retained simply to follow PCA code
    # Calc mean and variance of counts after normalizing
    # Don't transform to log space in PLSA
    # Dispersion is not exactly meaningful after idf transform.
    m = analysis_stats.normalize_by_idf(matrix)
    (mu, var) = analysis_stats.summarize_columns(m.T)
    dispersion = analysis_stats.get_normalized_dispersion(mu.squeeze(), var.squeeze()) # TODO set number of bins?

    plsa_feature_indices = np.argsort(dispersion)[-plsa_features:]

    if plsa_bcs < matrix.bcs_dim:
        plsa_bc_indices = np.sort(np.random.choice(np.arange(matrix.bcs_dim), size=plsa_bcs, replace=False))

    plsa_mat, _, plsa_features_nonzero = matrix.select_barcodes(plsa_bc_indices).select_features(plsa_feature_indices).select_nonzero_axes()
    plsa_feature_nonzero_indices = plsa_feature_indices[plsa_features_nonzero]

    if plsa_mat.features_dim < 2 or plsa_mat.bcs_dim < 2:
        print "Matrix is too small for further downsampling - num_plsa_bcs and num_plsa_features will be ignored."
        plsa_mat, _, plsa_features_nonzero = matrix.select_nonzero_axes()
        plsa_feature_nonzero_indices = plsa_features_nonzero

    ### Write out sparse matrix without transforms
    plsa_mat.tocoo()
    out_matrix_fn = os.path.join(temp_dir, 'matrix.mtx')
    sp_io.mmwrite(out_matrix_fn, plsa_mat.m, field='integer', symmetry='general')

    ### Run plsa module, reading in sparse matrix
    proc = tk_subproc.Popen([PLSA_BINPATH,
                            out_matrix_fn,
                            temp_dir,
                            '--topics', str(n_plsa_components),
                            '--nt', str(threads),
                            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_data, stderr_data = proc.communicate()
    if proc.returncode != 0:
        print stdout_data
        raise Exception("%s returned error code while running plsa binary %d: %s" % (proc, proc.returncode, stderr_data))

    ### Read back data
    transformed_plsa_em_matrix_file = os.path.join(temp_dir, "transformed_matrix.csv")
    n_components_file = os.path.join(temp_dir, "components.csv")
    variance_explained_file = os.path.join(temp_dir, "topic_relevance.csv")
    transformed_plsa_em_matrix = np.genfromtxt(transformed_plsa_em_matrix_file, delimiter=",").astype('float64')
    plsa_em_components = np.zeros((n_plsa_components, matrix.features_dim))
    plsa_em_components[:, plsa_feature_nonzero_indices] = np.genfromtxt(n_components_file, delimiter=",").astype('float64')
    variance_explained = np.genfromtxt(variance_explained_file, delimiter=",").astype('float64')

    ### reorder components by variance explained as PLSA binary gives arbitrary order
    new_order = range(n_plsa_components)
    variance_explained, new_order = zip(*sorted(zip(variance_explained, new_order), reverse=True))
    variance_explained = np.array(variance_explained)
    plsa_em_components = plsa_em_components[new_order, :]
    transformed_plsa_em_matrix = transformed_plsa_em_matrix[:, new_order]

    ### delete files
    cr_io.remove(transformed_plsa_em_matrix_file, allow_nonexisting=True)
    cr_io.remove(n_components_file, allow_nonexisting=True)
    cr_io.remove(variance_explained_file, allow_nonexisting=True)
    cr_io.remove(out_matrix_fn, allow_nonexisting=True)

    features_selected = np.array([f.id for f in matrix.feature_ref.feature_defs])[plsa_feature_nonzero_indices]

    # sanity check dimensions
    assert transformed_plsa_em_matrix.shape == (matrix.bcs_dim, n_plsa_components)
    assert plsa_em_components.shape == (n_plsa_components, matrix.features_dim)
    assert variance_explained.shape == (n_plsa_components,)

    return PLSA(transformed_plsa_em_matrix, plsa_em_components, variance_explained, dispersion, features_selected)

def get_plsa_mem_gb_from_matrix_dim(nonzero_entries):
    em_mem_gb = round(np.ceil(1.0 * nonzero_entries / analysis_constants.NUM_PLSA_EM_MATRIX_ENTRIES_PER_MEM_GB)) + analysis_constants.PLSA_EM_BASE_MEM_GB
    return h5_constants.MATRIX_MEM_GB_MULTIPLIER * max(h5_constants.MIN_MEM_GB, em_mem_gb)

def save_plsa_csv(plsa_map, matrix, base_dir):
    for n_components, plsa in plsa_map.iteritems():
        n_components_dir = os.path.join(base_dir, '%d_components' % n_components)
        cr_io.makedirs(n_components_dir, allow_existing=True)

        matrix_fn = os.path.join(n_components_dir, 'projection.csv')
        n_columns = plsa.transformed_plsa_matrix.shape[1]
        assert n_columns <= n_components
        matrix_header = ['Barcode'] + ['PC-%d' % (i+1) for i in xrange(n_columns)]
        analysis_io.save_matrix_csv(matrix_fn, plsa.transformed_plsa_matrix, matrix_header,
                              matrix.bcs)

        components_fn = os.path.join(n_components_dir, 'components.csv')
        components_header = ['PC'] + [f.id for f in matrix.feature_ref.feature_defs]
        analysis_io.save_matrix_csv(components_fn, plsa.components, components_header,
                              range(1, n_components+1))

        variance_fn = os.path.join(n_components_dir, 'variance.csv')
        variance_header = ['PC','Proportion.Variance.Explained']
        analysis_io.save_matrix_csv(variance_fn, plsa.variance_explained, variance_header,
                              range(1, n_components+1))

        dispersion_fn = os.path.join(n_components_dir, 'dispersion.csv')
        dispersion_header = ['Feature','Normalized.Dispersion']
        analysis_io.save_matrix_csv(dispersion_fn, plsa.dispersion, dispersion_header,
                              [f.id for f in matrix.feature_ref.feature_defs])

        features_fn = os.path.join(n_components_dir, 'features_selected.csv')
        # TODO: there are two columns here, but only 1 entry in the header...BAD
        features_header = ['Feature']
        analysis_io.save_matrix_csv(features_fn, plsa.features_selected, features_header, range(1, len(plsa.features_selected)+1))

def save_plsa_h5(plsa_map, f):
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_PLSA_GROUP)
    for n_components, plsa in plsa_map.iteritems():
        analysis_io.save_h5(f, group, str(n_components), plsa)

def load_plsa_from_h5(filename):
    """ Load just the PLSA info from an analysis h5 """
    with tables.open_file(filename, 'r') as f:
        group = f.root._v_groups[analysis_constants.ANALYSIS_H5_PLSA_GROUP]
        # Just take the first PLSA object, assuming we never have multiple
        for _, plsa in analysis_io.load_h5_iter(group, PLSA):
            return plsa
