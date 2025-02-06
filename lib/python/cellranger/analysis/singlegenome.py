#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os.path
import re
from collections.abc import Sequence
from typing import TypeAlias

import h5py as h5
import numpy as np
from six import ensure_binary

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
import cellranger.matrix as cr_matrix
from cellranger.analysis.analysis_types import (
    LSA,
    PCA,
    PLSA,
    TSNE,
    UMAP,
    DifferentialExpression,
)
from cellranger.wrapped_tables import tables

COMPONENTS = {
    "pca": analysis_constants.PCA_N_COMPONENTS_DEFAULT,
    "lsa": analysis_constants.LSA_N_COMPONENTS_DEFAULT,
    "plsa": analysis_constants.PLSA_N_COMPONENTS_DEFAULT,
}

TSNE_NAME = "tsne"
UMAP_NAME = "umap"
PROJECTION_TITLE = {TSNE_NAME: "t-SNE", UMAP_NAME: "UMAP"}

Projection: TypeAlias = str


class SingleGenomeAnalysis:
    def __init__(
        self, matrix: cr_matrix.CountMatrix, method: str, projections: Sequence[str] = (TSNE_NAME,)
    ):
        """Initialize secondary analysis results for a single genome analysis.

        Args:
            matrix (cr_matrix.CountMatrix): count matrix
            method (str): dimensionality reduction method pca|lsa|plsa
            projections (Sequence[str], optional): List of manifold projections
                to load. Defaults to ("tsne",).
        """
        assert projections, "Must specify projection(s) to load"
        for proj in projections:
            assert proj in (TSNE_NAME, UMAP_NAME)
        self.projections = projections
        self.matrix = matrix

        # parameters
        self.method = method
        self.random_state = analysis_constants.RANDOM_STATE
        self.min_n_clusters = analysis_constants.MIN_N_CLUSTERS
        self.max_n_clusters = analysis_constants.MAX_N_CLUSTERS_DEFAULT
        self.n_dimensionality_reduction_components = []
        # TSNE
        self.tsne_input_dims = COMPONENTS[method]
        self.n_tsne_components = analysis_constants.TSNE_N_COMPONENTS
        self.perplexity = analysis_constants.TSNE_DEFAULT_PERPLEXITY
        self.theta = analysis_constants.TSNE_THETA
        # UMAP
        self.umap_input_dims = COMPONENTS[method]
        self.n_umap_components = analysis_constants.UMAP_N_COMPONENTS
        self.umap_min_dist = analysis_constants.UMAP_MIN_DIST
        self.umap_n_neighbors = analysis_constants.UMAP_DEFAULT_N_NEIGHBORS
        self.dr_bcs = matrix.bcs_dim
        self.dr_features = matrix.features_dim

        # DR: n_components -> DR
        self.dimensionality_reduced_matrix: dict[int, PCA | PLSA | LSA] = {}

        # clustering_key -> CLUSTERING
        self.clusterings: dict[str, cr_clustering.CLUSTERING] = {}

        # DE: clustering_key -> DifferentialExpression
        self.differential_expression: dict[str, DifferentialExpression] = {}

        # t-SNE: n_tsne_components -> TSNE
        self.tsne: dict[str | bytes | None, TSNE] = {}

        # UMAP: n_umap components -> UMAP
        self.umap: dict[str | bytes | None, UMAP] = {}

        # Fix random seed
        np.random.seed(0)

    def is_zero_matrix(self):
        return self.matrix.bcs_dim == 0 or self.matrix.features_dim == 0

    def get_cluster_sizes(self):
        return {k: cr_clustering.get_cluster_sizes(v) for k, v in self.clusterings.items()}

    def _select_bc_indices(self, cell_bc_indices):
        self.matrix = self.matrix.select_barcodes(cell_bc_indices)

        self.dr_bcs = min(self.dr_bcs, len(cell_bc_indices))
        if self.method == "pca":
            for n_components, pca in self.dimensionality_reduced_matrix.items():
                assert isinstance(pca, PCA)
                self.dimensionality_reduced_matrix[n_components] = PCA(
                    pca.transformed_pca_matrix[cell_bc_indices, :],
                    pca.components,
                    pca.variance_explained,
                    pca.dispersion,
                    pca.features_selected,
                )
        elif self.method == "lsa":
            for n_components, lsa in self.dimensionality_reduced_matrix.items():
                assert isinstance(lsa, LSA)
                self.dimensionality_reduced_matrix[n_components] = LSA(
                    lsa.transformed_lsa_matrix[cell_bc_indices, :],
                    lsa.components,
                    lsa.variance_explained,
                    lsa.dispersion,
                    lsa.features_selected,
                )
        elif self.method == "plsa":
            for n_components, plsa in self.dimensionality_reduced_matrix.items():
                assert isinstance(plsa, PLSA)
                self.dimensionality_reduced_matrix[n_components] = PLSA(
                    plsa.transformed_plsa_matrix[cell_bc_indices, :],
                    plsa.components,
                    plsa.variance_explained,
                    plsa.dispersion,
                    plsa.features_selected,
                )

        # Subset all the cluster label vectors
        for key, clustering in self.clusterings.items():
            self.clusterings[key] = cr_clustering.subselect_barcodes(clustering, cell_bc_indices)

        for name, tsne in self.tsne.items():
            self.tsne[name] = TSNE(
                tsne.transformed_tsne_matrix[cell_bc_indices, :], name=tsne.name, key=tsne.key
            )

        for name, umap in self.umap.items():
            self.umap[name] = UMAP(
                umap.transformed_umap_matrix[cell_bc_indices, :], name=umap.name, key=umap.key
            )

    def subsample_bcs(self, num_bcs):
        """Subsample barcodes across entire analysis (matrix, DR, etc)."""
        if num_bcs >= self.matrix.bcs_dim:
            return

        cell_bc_indices = np.sort(
            np.random.choice(np.arange(self.matrix.bcs_dim), size=num_bcs, replace=False)
        )
        self._select_bc_indices(cell_bc_indices)

    def get_dimensionality_reduced_matrix(self, n_components=None):
        if not n_components:
            if self.method not in COMPONENTS:
                raise ValueError(f"method {self.method} not found")
            default_n_components = COMPONENTS[self.method]
            # if there is a reduced matrix matching the default n_components, choose that
            if default_n_components in self.n_dimensionality_reduction_components:
                n_components = default_n_components
            # otherwise, if there is only one reduced matrix of the correct type, choose that
            elif len(self.n_dimensionality_reduction_components) == 1:
                n_components = self.n_dimensionality_reduction_components[0]
            else:
                raise ValueError(
                    f"Analysis has multiple reduced matrices of type {self.method}, but the default"
                    f"n_components {default_n_components} is not available: "
                    f"{self.n_dimensionality_reduction_components}. Please specify n_components"
                )
        return self.dimensionality_reduced_matrix[n_components]

    def get_clustering(self, cluster_key):
        return self.clusterings[cluster_key]

    def get_tsne(self, key=analysis_constants.TSNE_DEFAULT_KEY):
        """Given analysis object, return the tsne with a given key while being compatible with older cellranger versions.

        Starting from CR7.0 and SR2.0, all libraries have their feature type prefixed on the tSNE output.
        Before, Gene Expression libraries only had the _2 prefix.
        """
        if key in self.tsne:
            tsne = self.tsne[key]
        elif b"2" in self.tsne:
            # Backward compatibility with pre-CR7.0/SR2.0 tSNE outputs,
            # where the key was not prefixed with gene_expression
            key = b"2"
            tsne = self.tsne[key]
        elif analysis_constants.TSNE_DEFAULT_KEY not in self.tsne and None in self.tsne:
            # Backward compatibility with older analysis HDF5 files
            # for which the key ends up being None because the 'key'
            # dataset was not generated in older versions.
            key = None
            tsne = self.tsne[key]
        else:
            raise KeyError(f"{key} not found in tsne analysis")

        return tsne

    def get_umap(self, key=analysis_constants.TSNE_DEFAULT_KEY):
        """Given analysis object, return the tsne with a given key while being compatible with older cellranger versions.

        Starting from CR7.0 and SR2.0, all libraries have their feature type prefixed on the tSNE output.
        Before, Gene Expression libraries only had the _2 prefix.
        """
        if key in self.umap:
            umap = self.umap[key]
        elif b"2" in self.umap:
            # Backward compatibility with pre-CR7.0/SR2.0 tSNE outputs,
            # where the key was not prefixed with gene_expression
            key = b"2"
            umap = self.umap[key]
        elif analysis_constants.TSNE_DEFAULT_KEY not in self.umap and None in self.umap:
            # Backward compatibility with older analysis HDF5 files
            # for which the key ends up being None because the 'key'
            # dataset was not generated in older versions.
            key = None
            umap = self.umap[key]
        else:
            raise KeyError(f"{key} not found in umap analysis")

        return umap

    @staticmethod
    def load_h5(filename, method, projections: Sequence[str] = (TSNE_NAME,)):
        assert projections, "Must specify projection(s) to load"

        # CountMatrix uses h5py, not pytables
        with h5.File(ensure_binary(filename), "r") as f:
            group: h5.Group = f[analysis_constants.ANALYSIS_H5_MATRIX_GROUP]
            matrix = cr_matrix.CountMatrix.load(group)

        analysis = SingleGenomeAnalysis(matrix, method, projections=projections)

        with tables.open_file(filename, "r") as f:
            grp = None
            if method == "pca":
                grp = analysis_constants.ANALYSIS_H5_PCA_GROUP
            elif method == "lsa":
                grp = analysis_constants.ANALYSIS_H5_LSA_GROUP
            elif method == "plsa":
                grp = analysis_constants.ANALYSIS_H5_PLSA_GROUP
            else:
                raise ValueError("method invalid")
            group = f.root._v_groups[grp]
            analysis._load_dimensionality_reduced_matrix_h5(group)

            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP]
            analysis._load_clustering_h5(group)

            group = f.root._v_groups[analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP]
            analysis._load_differential_expression_h5(group)

            if TSNE_NAME in projections:
                group = f.root._v_groups[analysis_constants.ANALYSIS_H5_TSNE_GROUP]
                analysis._load_tsne_h5(group)

            if UMAP_NAME in projections:
                group = f.root._v_groups[analysis_constants.ANALYSIS_H5_UMAP_GROUP]
                analysis._load_umap_h5(group)

        return analysis

    def _load_dimensionality_reduced_matrix_h5(self, group):
        regex = re.compile(r"^([^\d]*?)_?(\d+)$")
        method_to_class = {"pca": PCA, "lsa": LSA, "plsa": PLSA}
        if self.method in method_to_class:
            dim_red_class = method_to_class[self.method]
            for key, dim_red in analysis_io.load_h5_iter(group, dim_red_class):
                m = regex.match(key)
                if m is None:
                    print(f"key: {key}")
                n_components = int(regex.match(key).group(2))
                self.n_dimensionality_reduction_components.append(n_components)
                self.dimensionality_reduced_matrix[n_components] = dim_red
            return
        raise ValueError(f"method {self.method} not allowed")

    def _load_clustering_h5(self, group):
        for clustering_key, clustering in analysis_io.load_h5_iter(group, cr_clustering.CLUSTERING):
            self.clusterings[clustering_key] = clustering

    def _load_differential_expression_h5(self, group):
        for clustering_key, de in analysis_io.load_h5_iter(group, DifferentialExpression):
            self.differential_expression[clustering_key] = de

    def _load_tsne_h5(self, group):
        for _, tsne in analysis_io.load_h5_iter(group, TSNE):
            self.tsne[tsne.key] = tsne

    def _load_umap_h5(self, group):
        for _, umap in analysis_io.load_h5_iter(group, UMAP):
            self.umap[umap.key] = umap

    @staticmethod
    def load_clustering_keys_from_h5(filename):
        """Load just the clustering keys from an analysis h5."""
        with tables.open_file(filename, "r") as f:
            group = getattr(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
            return [node._v_name[1:] for node in group]

    @staticmethod
    def load_clustering_from_h5(filename, clustering_key):
        """Load a single clustering from an analysis h5."""
        with tables.open_file(filename, "r") as f:
            group = getattr(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
            for subgroup in group:
                if subgroup._v_name == "_" + clustering_key:
                    return analysis_io.load_h5_namedtuple(subgroup, cr_clustering.CLUSTERING)
            raise ValueError(
                f"Could not find clustering key: {clustering_key} in HDF5 file {filename}"
            )

    @staticmethod
    def load_default_format(base_dir, *, method, projections: Sequence[Projection]):
        h5_file_path = analysis_io.h5_path(base_dir)
        if os.path.exists(h5_file_path):
            return SingleGenomeAnalysis.load_h5(h5_file_path, method, projections)
        else:
            return None
