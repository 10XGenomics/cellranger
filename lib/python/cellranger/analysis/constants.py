#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

######################################################
# DO NOT add new items to this file.
#
# - If a constant is only used from a single module, put it in that module.
# - If a constant is only used in association with a particular module, put it
#   in that module.
# - If a constant is used in a bunch of places, create a module with a more
#   specifically descriptive name to put those constants.
######################################################

# tsne
TSNE_N_COMPONENTS = 2
TSNE_DEFAULT_KEY = b"gene_expression_2"
TSNE_DEFAULT_PERPLEXITY = 30
TSNE_THETA = 0.5
RANDOM_STATE = 0
TSNE_MAX_ITER = 1000
TSNE_STOP_LYING_ITER = 250
TSNE_MOM_SWITCH_ITER = 250
ANALYSIS_H5_TSNE_GROUP = "tsne"

# umap
UMAP_N_COMPONENTS = 2
UMAP_DEFAULT_N_NEIGHBORS = 30
UMAP_MIN_DIST = 0.3
UMAP_DEFAULT_METRIC = "correlation"
ANALYSIS_H5_UMAP_GROUP = "umap"

# clustering
ANALYSIS_H5_CLUSTERING_GROUP = "clustering"
ANALYSIS_H5_KMEANS_GROUP = "kmeans"
ANALYSIS_H5_KMEDOIDS_GROUP = "kmedoids"  # Deprecated
KMEDOIDS_SQMAT_ENTERIES_PER_MEM_GB = 12e7
GRAPHCLUST_NEIGHBORS_DEFAULT = 1
GRAPHCLUST_NEIGHBOR_A_DEFAULT = -230.0
GRAPHCLUST_NEIGHBOR_B_DEFAULT = 120.0

# multigenome
GEM_CLASS_MULTIPLET = b"Multiplet"
GEM_CLASS_GENOME0 = b"genome0"
GEM_CLASS_GENOME1 = b"genome1"
GEM_CLASSES = [GEM_CLASS_GENOME0, GEM_CLASS_GENOME1, GEM_CLASS_MULTIPLET]
DEFAULT_MULTIPLET_THRESHOLD = 10
MULTIPLET_PROB_THRESHOLD = 0.10
COUNT_PURITY_OUTLIER_PROB_THRESHOLD = 0.01

# pca
PCA_N_COMPONENTS_DEFAULT = 10
NUM_IRLB_MATRIX_ENTRIES_PER_MEM_GB = 10e6  # based on empirical testing
IRLB_BASE_MEM_GB = 1
ANALYSIS_H5_PCA_GROUP = "pca"

# chemistry batch correction
# this upper limit was determined via testing,
#   larger numbers of cells will require _at least_ memory reservation changes
#   in CORRECT_CHEMISTRY_BATCH.join, if not substantial algorithmic changes
CBC_MAX_NCELLS = 800000
CBC_N_COMPONENTS_DEFAULT = 100
CBC_KNN = 10
CBC_ALPHA = 0.1
CBC_SIGMA = 150
CBC_REALIGN_PANORAMA = False

# lsa
LSA_N_COMPONENTS_DEFAULT = 15
ANALYSIS_H5_LSA_GROUP = "lsa"

# plsa
PLSA_N_COMPONENTS_DEFAULT = 10
ANALYSIS_H5_PLSA_GROUP = "plsa"
NUM_PLSA_EM_MATRIX_ENTRIES_PER_MEM_GB = 25e6
PLSA_EM_BASE_MEM_GB = 1

# single genome
MIN_N_CLUSTERS = 2
MAX_N_CLUSTERS_DEFAULT = 10
MAX_N_CLUSTERS_ATAC = 5
ANALYSIS_H5_MATRIX_GROUP = "matrix"
ANALYSIS_H5_KMEANS_DIFFERENTIAL_EXPRESSION_GROUP = "differential_expression"  # Deprecated
ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP = "all_differential_expression"
ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_FEATURE_LIST = "diffexp_feature_indices"
