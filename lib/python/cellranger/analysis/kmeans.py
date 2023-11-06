#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""K-means clustering."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import scipy.spatial.distance as sp_dist
import sklearn.cluster as sk_cluster

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io


class KMEANS(NamedTuple):
    clusters: int
    cluster_score: float


# Clamp centroid distances to this value for DBI calc
MIN_CENTROID_DIST = 1e-3


def compute_db_index(matrix, kmeans) -> float:
    """Compute Davies-Bouldin index, a measure of clustering quality.

    Faster and possibly more reliable than silhouette score.
    """
    n = matrix.shape[0]
    k = kmeans.n_clusters

    centers = kmeans.cluster_centers_
    labels = kmeans.labels_

    centroid_dists = sp_dist.squareform(sp_dist.pdist(centers))
    # Avoid divide-by-zero
    centroid_dists[np.abs(centroid_dists) < MIN_CENTROID_DIST] = MIN_CENTROID_DIST

    wss = np.zeros(k)
    counts = np.zeros(k)

    for i in range(n):
        label = labels[i]
        # note: this is 2x faster than scipy sqeuclidean
        sqdist = np.square(matrix[i, :] - centers[label, :]).sum()
        wss[label] += sqdist
        counts[label] += 1

    # Handle empty clusters
    counts[counts == 0] = 1

    scatter: np.ndarray[int, np.dtype[np.float64]] = np.sqrt(wss / counts)
    mixitude = (scatter + scatter[:, np.newaxis]) / centroid_dists
    np.fill_diagonal(mixitude, 0.0)

    worst_case_mixitude = np.max(mixitude, axis=1)
    db_score = worst_case_mixitude.sum() / k

    return db_score


def run_kmeans(
    transformed_matrix,
    n_clusters: int,
    clustering_type=cr_clustering.CLUSTER_TYPE_KMEANS,
    random_state: int | None = None,
):
    """Run k-means clustering on the points in transformed_matrix.

    Find n_clusters clusters.
    """
    if random_state is None:
        random_state = analysis_constants.RANDOM_STATE

    kmeans = sk_cluster.KMeans(n_clusters=n_clusters, random_state=random_state)
    clusters = kmeans.fit_predict(transformed_matrix) + 1

    cluster_score = compute_db_index(transformed_matrix, kmeans)

    clusters = cr_clustering.relabel_by_size(clusters)
    clustering_key = cr_clustering.format_clustering_key(clustering_type, n_clusters)

    return cr_clustering.create_clustering(
        clusters=clusters,
        num_clusters=n_clusters,
        cluster_score=cluster_score,
        clustering_type=clustering_type,
        global_sort_key=n_clusters,
        description=cr_clustering.humanify_clustering_key(clustering_key),
    )


def save_kmeans_h5(f, n_clusters, kmeans, clustering_type):
    """Save k-means clustering results to HDF5 file object f.

    Args:
        f: The hdf5 object to save to.
        kmeans: Cluster centers.
        n_clusters: The number of clusters.
        clustering_type: whether GEX-based or Antibody-bases K-means
    """
    clustering_key = cr_clustering.format_clustering_key(clustering_type, n_clusters)

    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
    analysis_io.save_h5(f, group, clustering_key, kmeans)

    cr_clustering.create_legacy_kmeans_nodes(
        f,
        analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP,
        analysis_constants.ANALYSIS_H5_KMEANS_GROUP,
        cr_clustering.CLUSTERING,
        clustering_key,
    )
