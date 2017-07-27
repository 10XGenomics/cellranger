#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

import collections
from irlb import irlb
import numpy as np
import os
import tables

import scipy
from sklearn.utils import sparsefuncs

PCA = collections.namedtuple('PCA', ['transformed_pca_matrix', 'components', 'variance_explained', 'dispersion', 'genes_selected'])

def run_pca(matrix, pca_genes=None, pca_bcs=None, n_pca_components=None, random_state=None):
    if pca_genes is None:
        pca_genes = matrix.genes_dim
    if pca_bcs is None:
        pca_bcs = matrix.bcs_dim
    if n_pca_components is None:
        n_pca_components = cr_constants.PCA_N_COMPONENTS_DEFAULT
        if n_pca_components > pca_genes:
            print "There are fewer nonzero genes than PCA components; reducing the number of components."
            n_pca_components = pca_genes
    if random_state is None:
        random_state=cr_constants.RANDOM_STATE

    np.random.seed(0)

    (full_norm_mat, full_center, full_scale) = normalize_and_transpose(matrix)

    # initialize PCA subsets
    pca_bc_indices = np.arange(matrix.bcs_dim)
    pca_gene_indices = np.arange(matrix.genes_dim)

    # Calc mean and variance of counts after normalizing
    # But don't transform to log space, in order to preserve the mean-variance relationship
    m = normalize_by_umi(matrix)
    (mu, var) = summarize_columns(m.T)
    dispersion = get_normalized_dispersion(mu.squeeze(), var.squeeze()) # TODO set number of bins?

    pca_gene_indices = np.argsort(dispersion)[-pca_genes:]

    if pca_bcs < matrix.bcs_dim:
        pca_bc_indices = np.sort(np.random.choice(np.arange(matrix.bcs_dim), size=pca_bcs, replace=False))

    pca_mat, _, pca_genes_nonzero = matrix.select_barcodes(pca_bc_indices).select_genes(pca_gene_indices).select_nonzero_axes()
    pca_gene_nonzero_indices = pca_gene_indices[pca_genes_nonzero]

    if pca_mat.genes_dim < 2 or pca_mat.bcs_dim < 2:
        print "Matrix is too small for further downsampling - num_pca_bcs and num_pca_genes will be ignored."
        pca_mat, _, pca_genes_nonzero = matrix.select_nonzero_axes()
        pca_gene_nonzero_indices = pca_genes_nonzero

    (pca_norm_mat, pca_center, pca_scale) = normalize_and_transpose(pca_mat)

    (u, d, v, _, _) = irlb(pca_norm_mat, n_pca_components, center=pca_center.squeeze(), scale=pca_scale.squeeze(), random_state=random_state)

    # make sure to project the matrix before centering, to avoid densification
    sparsefuncs.inplace_column_scale(full_norm_mat, 1 / full_scale.squeeze())
    transformed_irlba_matrix = full_norm_mat[:,pca_gene_nonzero_indices].dot(v) - (full_center / full_scale)[:,pca_gene_nonzero_indices].dot(v)
    irlba_components = np.zeros((n_pca_components, matrix.genes_dim))
    irlba_components[:,pca_gene_nonzero_indices] = v.T

    # calc proportion of variance explained
    variance_sum = len(pca_gene_indices) # each gene has variance=1, mean=0 after normalization
    variance_explained = np.square(d)/((len(pca_bc_indices)-1) * variance_sum)

    genes_selected = np.array([gene.id for gene in matrix.genes])[pca_gene_nonzero_indices]

    # sanity check dimensions
    assert transformed_irlba_matrix.shape == (matrix.bcs_dim, n_pca_components)
    assert irlba_components.shape == (n_pca_components, matrix.genes_dim)
    assert variance_explained.shape == (n_pca_components,)

    return PCA(transformed_irlba_matrix, irlba_components, variance_explained, dispersion, genes_selected)

def normalize_by_umi(matrix):
    reads_per_bc = matrix.get_reads_per_bc()
    median_reads_per_bc = np.median(reads_per_bc)
    scaling_factors = median_reads_per_bc / reads_per_bc

    # Normalize each barcode's total count by median total count
    m = matrix.m.copy().astype(np.float64)
    sparsefuncs.inplace_column_scale(m, scaling_factors)

    return m

def summarize_columns(matrix):
    ''' Calculate mean and variance of each column, in a sparsity-preserving way.'''
    mu = matrix.mean(axis=0).A

    # sparse variance = E(col^2) - E(col)^2
    mu2 = matrix.multiply(matrix).mean(axis=0).A
    var = mu2 - mu**2

    return mu, var

def normalize_and_transpose(matrix):
    matrix.tocsc()

    m = normalize_by_umi(matrix)

    # Use log counts
    m.data = np.log2(1 + m.data)

    # Transpose
    m = m.T

    # compute centering (mean) and scaling (stdev)
    (c,v) = summarize_columns(m)
    s = np.sqrt(v)

    return (m, c, s)

def get_normalized_dispersion(mat_mean, mat_var, nbins=20):
    mat_disp = (mat_var - mat_mean) / np.square(mat_mean)

    quantiles = np.percentile(mat_mean, np.arange(0, 100, 100 / nbins))
    quantiles = np.append(quantiles, mat_mean.max())

    # merge bins with no difference in value
    quantiles = np.unique(quantiles)

    if len(quantiles) <= 1:
        # pathological case: the means are all identical. just return raw dispersion.
        return mat_disp

    # calc median dispersion per bin
    (disp_meds, _, disp_bins) = scipy.stats.binned_statistic(mat_mean, mat_disp, statistic='median', bins=quantiles)

    # calc median absolute deviation of dispersion per bin
    disp_meds_arr = disp_meds[disp_bins-1] # 0th bin is empty since our quantiles start from 0
    disp_abs_dev = abs(mat_disp - disp_meds_arr)
    (disp_mads, _, disp_bins) = scipy.stats.binned_statistic(mat_mean, disp_abs_dev, statistic='median', bins=quantiles)

    # calculate normalized dispersion
    disp_mads_arr = disp_mads[disp_bins-1]
    disp_norm = (mat_disp - disp_meds_arr) / disp_mads_arr
    return disp_norm

def get_irlb_mem_gb_from_matrix_dim(nonzero_entries):
    irlba_mem_gb = round(np.ceil(1.0 * nonzero_entries / cr_constants.NUM_IRLB_MATRIX_ENTRIES_PER_MEM_GB)) + cr_constants.IRLB_BASE_MEM_GB
    return cr_constants.MATRIX_MEM_GB_MULTIPLIER * max(cr_constants.MIN_MEM_GB, irlba_mem_gb)

def save_pca_csv(pca_map, matrix, base_dir):
    for n_components, pca in pca_map.iteritems():
        n_components_dir = os.path.join(base_dir, '%d_components' % n_components)
        cr_utils.makedirs(n_components_dir, allow_existing=True)

        matrix_fn = os.path.join(n_components_dir, 'projection.csv')
        matrix_header = ['Barcode'] + ['PC-%d' % (i+1) for i in xrange(n_components)]
        cr_io.save_matrix_csv(matrix_fn, pca.transformed_pca_matrix, matrix_header,
                              matrix.bcs)

        components_fn = os.path.join(n_components_dir, 'components.csv')
        components_header = ['PC'] + [gene.id for gene in matrix.genes]
        cr_io.save_matrix_csv(components_fn, pca.components, components_header,
                              range(1, n_components+1))

        variance_fn = os.path.join(n_components_dir, 'variance.csv')
        variance_header = ['PC','Proportion.Variance.Explained']
        cr_io.save_matrix_csv(variance_fn, pca.variance_explained, variance_header,
                              range(1, n_components+1))

        dispersion_fn = os.path.join(n_components_dir, 'dispersion.csv')
        dispersion_header = ['Gene','Normalized.Dispersion']
        cr_io.save_matrix_csv(dispersion_fn, pca.dispersion, dispersion_header,
                              [gene.id for gene in matrix.genes])

        genes_fn = os.path.join(n_components_dir, 'genes_selected.csv')
        genes_header = ['Gene']
        cr_io.save_matrix_csv(genes_fn, pca.genes_selected, genes_header, range(1, len(pca.genes_selected)+1))

def save_pca_h5(pca_map, f):
    group = f.create_group(f.root, cr_constants.ANALYSIS_H5_PCA_GROUP)
    for n_components, pca in pca_map.iteritems():
        cr_io.save_h5(f, group, str(n_components), pca)

def load_pca_from_h5(filename):
    """ Load just the PCA info from an analysis h5 """
    with tables.open_file(filename, 'r') as f:
        group = f.root._v_groups[cr_constants.ANALYSIS_H5_PCA_GROUP]
        # Just take the first PCA object, assuming we never have multiple
        for _, pca in cr_io.load_h5_iter(group, PCA):
            return pca
