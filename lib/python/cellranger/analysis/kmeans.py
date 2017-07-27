#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import collections
import numpy as np
import scipy.spatial.distance as sp_dist
import sklearn.cluster as sk_cluster
import cellranger.analysis.io as cr_io
import cellranger.analysis.clustering as cr_clustering
import cellranger.constants as cr_constants


KMEANS = collections.namedtuple('KMEANS', ['clusters', 'cluster_score'])

# Clamp centroid distances to this value for DBI calc
MIN_CENTROID_DIST = 1e-3

def compute_db_index(matrix, kmeans):
    '''
    Compute Davies-Bouldin index, a measure of clustering quality.
    Faster and possibly more reliable than silhouette score.
    '''
    (n, m) = matrix.shape
    k = kmeans.n_clusters

    centers = kmeans.cluster_centers_
    labels = kmeans.labels_

    centroid_dists = sp_dist.squareform(sp_dist.pdist(centers))
    # Avoid divide-by-zero
    centroid_dists[np.abs(centroid_dists) < MIN_CENTROID_DIST] = MIN_CENTROID_DIST

    wss = np.zeros(k)
    counts = np.zeros(k)

    for i in xrange(n):
        label = labels[i]
        # note: this is 2x faster than scipy sqeuclidean
        sqdist = np.square(matrix[i,:] - centers[label,:]).sum()
        wss[label] += sqdist
        counts[label] += 1

    # Handle empty clusters
    counts[counts == 0] = 1

    scatter = np.sqrt(wss / counts)
    mixitude = (scatter + scatter[:, np.newaxis]) / centroid_dists
    np.fill_diagonal(mixitude, 0.0)

    worst_case_mixitude = np.max(mixitude, axis=1)
    db_score = worst_case_mixitude.sum() / k

    return db_score

def run_kmeans(transformed_pca_matrix, n_clusters, random_state=None):
    if random_state is None:
        random_state=cr_constants.RANDOM_STATE

    kmeans = sk_cluster.KMeans(n_clusters=n_clusters, random_state=random_state)
    clusters = kmeans.fit_predict(transformed_pca_matrix) + 1

    cluster_score = compute_db_index(transformed_pca_matrix, kmeans)

    clusters = cr_clustering.relabel_by_size(clusters)

    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_KMEANS, n_clusters)

    return cr_clustering.create_clustering(clusters=clusters,
                                           num_clusters=n_clusters,
                                           cluster_score=cluster_score,
                                           clustering_type=cr_clustering.CLUSTER_TYPE_KMEANS,
                                           global_sort_key=n_clusters,
                                           description=cr_clustering.humanify_clustering_key(clustering_key))

def save_kmeans_h5(f, n_clusters, kmeans):
    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_KMEANS, n_clusters)

    group = f.create_group(f.root, cr_constants.ANALYSIS_H5_CLUSTERING_GROUP)
    cr_io.save_h5(f, group, clustering_key, kmeans)

    cr_clustering.create_legacy_kmeans_nodes(f,
                                             cr_constants.ANALYSIS_H5_CLUSTERING_GROUP,
                                             cr_constants.ANALYSIS_H5_KMEANS_GROUP,
                                             cr_clustering.CLUSTERING,
                                             clustering_key)
