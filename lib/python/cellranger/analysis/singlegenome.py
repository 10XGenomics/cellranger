#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.matrix as cr_matrix

from cellranger.analysis.bhtsne import TSNE
import cellranger.analysis.clustering as cr_clustering
from cellranger.analysis.diffexp import DIFFERENTIAL_EXPRESSION
from cellranger.analysis.pca import PCA
from cellranger.analysis.lsa import LSA
from cellranger.analysis.plsa import PLSA

import h5py as h5
import numpy as np
import os
import tables

COMPONENTS = {'pca': analysis_constants.PCA_N_COMPONENTS_DEFAULT,
              'lsa': analysis_constants.LSA_N_COMPONENTS_DEFAULT,
              'plsa': analysis_constants.PLSA_N_COMPONENTS_DEFAULT}

class SingleGenomeAnalysis:
    def __init__(self, matrix, method):
        self.matrix = matrix

        # parameters
        self.method = method
        self.random_state = analysis_constants.RANDOM_STATE
        self.min_n_clusters = analysis_constants.MIN_N_CLUSTERS
        self.max_n_clusters = analysis_constants.MAX_N_CLUSTERS_DEFAULT
        self.n_dimensionality_reduction_components = COMPONENTS[method]
        self.tsne_input_dims = COMPONENTS[method]
        self.n_tsne_components = analysis_constants.TSNE_N_COMPONENTS
        self.perplexity = analysis_constants.TSNE_DEFAULT_PERPLEXITY
        self.theta = analysis_constants.TSNE_THETA
        self.dr_bcs = matrix.bcs_dim
        self.dr_features = matrix.features_dim

        # DR: n_components -> DR
        self.dimensionality_reducted_matrix = {}

        # clustering_key -> CLUSTERING
        self.clusterings = {}

        # DE: clustering_key -> DIFFERENTIAL_EXPRESSION
        self.differential_expression = {}

        # t-SNE: n_tsne_components -> TSNE
        self.tsne = {}

        # Fix random seed
        np.random.seed(0)

    def is_zero_matrix(self):
        return self.matrix.bcs_dim == 0 or self.matrix.features_dim == 0

    def get_cluster_sizes(self):
        return {k: cr_clustering.get_cluster_sizes(v) for k, v in self.clusterings.iteritems()}

    def _select_bc_indices(self, cell_bc_indices):
        self.matrix = self.matrix.select_barcodes(cell_bc_indices)

        self.dr_bcs = min(self.dr_bcs, len(cell_bc_indices))
        if self.method == 'pca':
            for n_components, pca in self.dimensionality_reducted_matrix.iteritems():
                self.dimensionality_reducted_matrix[n_components] = PCA(pca.transformed_pca_matrix[cell_bc_indices, :],
                                                                        pca.components, pca.variance_explained,
                                                                        pca.dispersion, pca.features_selected)
        elif self.method == 'lsa':
            for n_components, lsa in self.dimensionality_reducted_matrix.iteritems():
                self.dimensionality_reducted_matrix[n_components] = LSA(lsa.transformed_lsa_matrix[cell_bc_indices, :],
                                                                        lsa.components, lsa.variance_explained,
                                                                        lsa.dispersion, lsa.features_selected)
        elif self.method == 'plsa':
            for n_components, plsa in self.dimensionality_reducted_matrix.iteritems():
                self.dimensionality_reducted_matrix[n_components] = PLSA(plsa.transformed_plsa_matrix[cell_bc_indices, :],
                                                                         plsa.components, plsa.variance_explained,
                                                                         plsa.dispersion, plsa.features_selected)


        # Subset all the cluster label vectors
        for key, clustering in self.clusterings.iteritems():
            self.clusterings[key] = cr_clustering.subselect_barcodes(clustering, cell_bc_indices)

        for name, tsne in self.tsne.iteritems():
            self.tsne[name] = TSNE(tsne.transformed_tsne_matrix[cell_bc_indices, :],
                                   name=name)

    def subsample_bcs(self, num_bcs):
        """ Subsample barcodes across entire analysis (matrix, DR, etc) """
        if num_bcs >= self.matrix.bcs_dim:
            return

        cell_bc_indices = np.sort(np.random.choice(np.arange(self.matrix.bcs_dim), size=num_bcs, replace=False))
        self._select_bc_indices(cell_bc_indices)

    def get_dimensionality_reduced_matrix(self, n_components=None):
        if not n_components:
            if self.method in COMPONENTS.keys():
                n_components = COMPONENTS[self.method]
            else:
                raise ValueError('method not found')
        return self.dimensionality_reducted_matrix[n_components]

    def get_clustering(self, cluster_key):
        return self.clusterings[cluster_key]

    def get_tsne(self, key=analysis_constants.TSNE_DEFAULT_KEY):
        return self.tsne[key]

    @staticmethod
    def load_h5(filename, method):

        # CountMatrix uses h5py, not pytables
        with h5.File(filename, 'r') as f:
            group = f[analysis_constants.ANALYSIS_H5_MATRIX_GROUP]
            matrix = cr_matrix.CountMatrix.load(group)

        analysis = SingleGenomeAnalysis(matrix, method)

        with tables.open_file(filename, 'r') as f:
            grp = None
            if method == 'pca':
                grp = analysis_constants.ANALYSIS_H5_PCA_GROUP
            elif method == 'lsa':
                grp = analysis_constants.ANALYSIS_H5_LSA_GROUP
            elif method == 'plsa':
                grp = analysis_constants.ANALYSIS_H5_PLSA_GROUP
            else:
                raise ValueError('method invalid')
            group = f.root._v_groups[grp]
            analysis._load_dimensionality_reduced_matrix_h5(group)

            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP]
            analysis._load_clustering_h5(group)

            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP]
            analysis._load_differential_expression_h5(group)

            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_TSNE_GROUP]
            analysis._load_tsne_h5(group)

        return analysis

    def _load_dimensionality_reduced_matrix_h5(self, group):
        if self.method == 'pca':
            for n_components, pca in analysis_io.load_h5_iter(group, PCA):
                self.dimensionality_reducted_matrix[int(n_components)] = pca
            return
        elif self.method == 'lsa':
            for n_components, lsa in analysis_io.load_h5_iter(group, LSA):
                self.dimensionality_reducted_matrix[int(n_components)] = lsa
            return
        elif self.method == 'plsa':
            for n_components, plsa in analysis_io.load_h5_iter(group, PLSA):
                self.dimensionality_reducted_matrix[int(n_components)] = plsa
            return
        raise ValueError('method {} not allowed'.format(self.method))

    def _load_clustering_h5(self, group):
        for clustering_key, clustering in analysis_io.load_h5_iter(group, cr_clustering.CLUSTERING):
            self.clusterings[clustering_key] = clustering

    def _load_differential_expression_h5(self, group):
        for clustering_key, de in analysis_io.load_h5_iter(group, DIFFERENTIAL_EXPRESSION):
            self.differential_expression[clustering_key] = de

    def _load_tsne_h5(self, group):
        for _, tsne in analysis_io.load_h5_iter(group, TSNE):
            self.tsne[tsne.key] = tsne

    @staticmethod
    def load_pca_from_h5(filename):
        """ Load just the PCA info from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_PCA_GROUP]
            # Just take the first PCA object, assuming we never have multiple
            for _, pca in analysis_io.load_h5_iter(group, PCA):
                return pca

    def load_lsa_from_h5(filename):
        """ Load just the LSA info from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_LSA_GROUP]
            # Just take the first LSA object, assuming we never have multiple
            for _, lsa in analysis_io.load_h5_iter(group, LSA):
                return lsa

    def load_plsa_from_h5(filename):
        """ Load just the PLSA info from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_PLSA_GROUP]
            # Just take the first PLSA object, assuming we never have multiple
            for _, plsa in analysis_io.load_h5_iter(group, PLSA):
                return plsa

    @staticmethod
    def load_clustering_keys_from_h5(filename):
        """ Load just the clustering keys from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = getattr(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
            return [node._v_name[1:] for node in group]

    @staticmethod
    def load_clustering_from_h5(filename, clustering_key):
        """ Load a single clustering from an analysis h5 """
        with tables.open_file(filename, 'r') as f:
            group = getattr(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
            for subgroup in group:
                if subgroup._v_name == '_' + clustering_key:
                    return analysis_io.load_h5_namedtuple(subgroup, cr_clustering.CLUSTERING)
            raise ValueError("Could not find clustering key: %s in HDF5 file %s" % (clustering_key, filename))

    @staticmethod
    def load_bcs_from_matrix_h5(filename):
        """ Load just the barcodes from a matrix h5 """
        with h5.File(filename, 'r') as f:
            # Take the first group, assuming a single-genome matrix
            # TODO: fixme when we have a single matrix group
            group_name = f.keys()[0]
            return cr_matrix.CountMatrix.load_bcs_from_h5_group(f[group_name])

    @staticmethod
    def load_default_format(base_dir, method):
        h5_file_path = analysis_io.h5_path(base_dir)
        if os.path.exists(h5_file_path):
            return SingleGenomeAnalysis.load_h5(h5_file_path, method)
        else:
            return None
