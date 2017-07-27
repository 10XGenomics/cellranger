#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix

from cellranger.analysis.bhtsne import TSNE
import cellranger.analysis.clustering as cr_clustering
from cellranger.analysis.diffexp import DIFFERENTIAL_EXPRESSION
from cellranger.analysis.pca import PCA

import numpy as np
import os
import tables

class SingleGenomeAnalysis:
    def __init__(self, matrix):
        self.matrix = matrix

        # parameters
        self.random_state = cr_constants.RANDOM_STATE
        self.n_pca_components = cr_constants.PCA_N_COMPONENTS_DEFAULT
        self.min_n_clusters = cr_constants.MIN_N_CLUSTERS
        self.max_n_clusters = cr_constants.MAX_N_CLUSTERS_DEFAULT
        self.tsne_input_dims = cr_constants.PCA_N_COMPONENTS_DEFAULT
        self.n_tsne_components = cr_constants.TSNE_N_COMPONENTS
        self.perplexity = cr_constants.TSNE_DEFAULT_PERPLEXITY
        self.theta = cr_constants.TSNE_THETA
        self.pca_bcs = matrix.bcs_dim
        self.pca_genes = matrix.genes_dim

        # PCA: n_components -> PCA
        self.pca = {}

        # clustering_key -> CLUSTERING
        self.clusterings = {}

        # DE: clustering_key -> DIFFERENTIAL_EXPRESSION
        self.differential_expression = {}

        # t-SNE: n_tsne_components -> TSNE
        self.tsne = {}

        # Fix random seed
        np.random.seed(0)

    def is_zero_matrix(self):
        return self.matrix.bcs_dim == 0 or self.matrix.genes_dim == 0

    def get_cluster_sizes(self):
        return {k: cr_clustering.get_cluster_sizes(v) for k,v in self.clusterings.iteritems() }

    def _select_bc_indices(self, cell_bc_indices):
        self.matrix = self.matrix.select_barcodes(cell_bc_indices)

        self.pca_bcs = min(self.pca_bcs, len(cell_bc_indices))
        for n_components, pca in self.pca.iteritems():
            self.pca[n_components] = PCA(pca.transformed_pca_matrix[cell_bc_indices, :], pca.components, pca.variance_explained, pca.dispersion, pca.genes_selected)

        # Subset all the cluster label vectors
        for key, clustering in self.clusterings.iteritems():
            self.clusterings[key] = cr_clustering.subselect_barcodes(clustering, cell_bc_indices)

        for n_tsne_components, tsne in self.tsne.iteritems():
            self.tsne[n_tsne_components] = TSNE(tsne.transformed_tsne_matrix[cell_bc_indices, :])

    def subsample_bcs(self, num_bcs):
        """ Subsample barcodes across entire analysis (matrix, PCA, etc) """
        if num_bcs >= self.matrix.bcs_dim:
            return

        cell_bc_indices = np.sort(np.random.choice(np.arange(self.matrix.bcs_dim), size=num_bcs, replace=False))
        self._select_bc_indices(cell_bc_indices)

    def get_pca(self, n_components=cr_constants.PCA_N_COMPONENTS_DEFAULT):
        return self.pca[n_components]

    def get_clustering(self, cluster_key):
        return self.clusterings[cluster_key]

    def get_tsne(self, n_tsne_components=cr_constants.TSNE_N_COMPONENTS):
        return self.tsne[n_tsne_components]

    @staticmethod
    def load_h5(filename):
        with tables.open_file(filename, 'r') as f:
            group = f.root._v_groups[cr_constants.ANALYSIS_H5_MATRIX_GROUP]
            matrix = cr_matrix.GeneBCMatrix.load(group)

            analysis = SingleGenomeAnalysis(matrix)
            group = f.root._v_groups[cr_constants.ANALYSIS_H5_PCA_GROUP]
            analysis._load_pca_h5(group)

            group = f.root._v_groups[cr_constants.ANALYSIS_H5_CLUSTERING_GROUP]
            analysis._load_clustering_h5(group)

            group = f.root._v_groups[cr_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP]
            analysis._load_differential_expression_h5(group)

            group = f.root._v_groups[cr_constants.ANALYSIS_H5_TSNE_GROUP]
            analysis._load_tsne_h5(group)

        return analysis

    def _load_pca_h5(self, group):
        for n_components, pca in cr_io.load_h5_iter(group, PCA):
            self.pca[int(n_components)] = pca

    def _load_clustering_h5(self, group):
        for clustering_key, clustering in cr_io.load_h5_iter(group, cr_clustering.CLUSTERING):
            self.clusterings[clustering_key] = clustering

    def _load_differential_expression_h5(self, group):
        for clustering_key, de in cr_io.load_h5_iter(group, DIFFERENTIAL_EXPRESSION):
            self.differential_expression[clustering_key] = de

    def _load_tsne_h5(self, group):
        for n_tsne_components, tsne in cr_io.load_h5_iter(group, TSNE):
            self.tsne[int(n_tsne_components)] = tsne

    @staticmethod
    def load_pca_from_h5(filename):
        """ Load just the PCA info from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = f.root._v_groups[cr_constants.ANALYSIS_H5_PCA_GROUP]
            # Just take the first PCA object, assuming we never have multiple
            for _, pca in cr_io.load_h5_iter(group, PCA):
                return pca

    @staticmethod
    def load_clustering_keys_from_h5(filename):
        """ Load just the clustering keys from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = getattr(f.root, cr_constants.ANALYSIS_H5_CLUSTERING_GROUP)
            return [node._v_name[1:] for node in group]

    @staticmethod
    def load_clustering_from_h5(filename, clustering_key):
        """ Load a single clustering from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = getattr(f.root, cr_constants.ANALYSIS_H5_CLUSTERING_GROUP)
            for subgroup in group:
                if subgroup._v_name == '_' + clustering_key:
                    return cr_io.load_h5_namedtuple(subgroup, cr_clustering.CLUSTERING)
            raise ValueError("Could not find clustering key: %s in HDF5 file %s" % (clustering_key, filename))

    @staticmethod
    def load_bcs_from_matrix_h5(filename):
        """ Load just the barcodes from a matrix h5 """
        with tables.open_file(filename, 'r') as f:
            # Take the first group, assuming a single-genome matrix
            group = list(f.list_nodes(f.root))[0]
            return cr_matrix.GeneBCMatrix.load_bcs_from_h5_group(group)

    @staticmethod
    def load_default_format(base_dir):
        h5_file_path = cr_io.h5_path(base_dir)
        if os.path.exists(h5_file_path):
            return SingleGenomeAnalysis.load_h5(h5_file_path)
        else:
            return None
