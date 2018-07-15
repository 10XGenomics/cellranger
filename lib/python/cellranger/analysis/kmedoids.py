#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import collections
import numpy as np
import scipy.spatial.distance as sp_dist
import cellranger.analysis.io as analysis_io
import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
from sklearn.metrics import silhouette_score

KMEDOIDS = collections.namedtuple('KMEDOIDS', ['clusters', 'cluster_score'])

# Clamp centroid distances to this value for DBI calc
MIN_CENTROID_DIST = 1e-6

def cosine_dist(p, q, squared=False):
    """Cosine distance is calculated between two categorical distributions"""

    assert p.shape == q.shape
    coeff = np.sum(p * q) * 1.0 / (np.linalg.norm(p) + 1e-28) / (np.linalg.norm(q) + 1e-28)
    if squared:
        return (1.0 - coeff)**2
    else:
        return (1.0 - coeff)

class KMedoids:
    """
    Kmedoids implementation is provided here. This is an adaption of the implementation:
    https://github.com/letiantian/kmedoids/blob/master/kmedoids.py, which is an adaptation of
    Bauckhage C. Numpy/scipy Recipes for Data Science: k-Medoids Clustering[R]. Technical Report,
    University of Bonn, 2015. Our adaptation is immune to cases where a medoid is a true outlier
    and thus no points are assigned to the medoid during the iterations.
    """
    def __init__(self, n_clusters, random_state, tmax=300, force=True, metric='euclidean'):
        """Other options for the metric is cosine, which can be custom defined"""
        self.nmedoids = n_clusters
        self.tmax = tmax
        self.random_state = random_state
        self.labels_ = None
        self.medoids = None
        self.force = force
        self.metric = metric

    def fit_predict(self, matrix):
        distance_matrix = sp_dist.squareform(sp_dist.pdist(matrix, metric=self.metric))

        # determine dimensions of distance matrix D
        m, n = distance_matrix.shape

        if self.nmedoids > n:
            raise ValueError('too many medoids')

        np.random.seed(self.random_state)
        # randomly initialize an array of nmedoids medoid indices
        medoids = np.sort(np.random.choice(n, self.nmedoids))

        # create a copy of the array of medoid indices
        medoidsnew = np.copy(medoids)

        # initialize a dictionary to represent clusters
        clusters = {}

        Nrange = range(n)
        for t in xrange(self.tmax):
            # determine clusters, i.e. arrays of data indices
            J = np.argmin(distance_matrix[:, medoids], axis=1)
            for medoid_label in range(self.nmedoids):
                clusters[medoid_label] = np.where(J == medoid_label)[0]

            # update cluster medoids
            bad_medoid_label = []
            good_medoid_label = []
            for medoid_label in range(self.nmedoids):
                # NOTE: this can be memory intensive on a large bloc
                J = np.mean(distance_matrix[np.ix_(clusters[medoid_label], clusters[medoid_label])], axis=1)
                if J.size > 0:
                    j = np.argmin(J)
                    medoidsnew[medoid_label] = clusters[medoid_label][j]
                    good_medoid_label += [medoid_label]
                else:
                    bad_medoid_label += [medoid_label]
            # randomize bad medoid_labels
            if self.force:
                medoidsnew[bad_medoid_label] = np.random.choice([x for x in Nrange if x not in list(medoidsnew[good_medoid_label])], len(bad_medoid_label))
            np.sort(medoidsnew)

            # check for convergence, allowing for nans
            if ((medoids == medoidsnew) | (np.isnan(medoids) & np.isnan(medoidsnew))).all():
                break

            medoids = np.copy(medoidsnew)
        else:
            # final update of cluster memberships
            J = np.argmin(distance_matrix[:, medoids], axis=1)
            for medoid_label in range(self.nmedoids):
                clusters[medoid_label] = np.where(J == medoid_label)[0]

        # return results
        labels = np.full(distance_matrix.shape[0], 0)
        for n in clusters.keys():
            if len(clusters[n]) > 0:
                labels[clusters[n]] = n
        for n, medoid in enumerate(medoids):
            labels[medoid] = n

        self.medoids_ = matrix[medoids, :]
        self.labels_ = labels
        return labels

def compute_silhouette_index(matrix, kmedoids, metric='euclidean'):
    '''
    Compute Silhouette score, a measure of clustering quality.
    '''
    # TODO: potentially one could develop Davies-Bouldin index with custom metrics
    # Reasons:
    # DB index could be faster
    # A silhouette score close to 1 is equivalent to a DB index close to 0
    # A silhouette score close to -1 is equivalent to a DB index close to 1

    return silhouette_score(matrix, kmedoids.labels_, metric)

def run_kmedoids(transformed_matrix, n_clusters, random_state=None, metric='euclidean'):
    if random_state is None:
        random_state=analysis_constants.RANDOM_STATE

    kmedoids = KMedoids(n_clusters=n_clusters, random_state=random_state, metric=metric)
    clusters = kmedoids.fit_predict(transformed_matrix) + 1

    cluster_score = compute_silhouette_index(transformed_matrix, kmedoids, metric)

    clusters = cr_clustering.relabel_by_size(clusters)

    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_KMEDOIDS, n_clusters)

    return cr_clustering.create_clustering(clusters=clusters,
                                           num_clusters=n_clusters,
                                           cluster_score=cluster_score,
                                           clustering_type=cr_clustering.CLUSTER_TYPE_KMEDOIDS,
                                           global_sort_key=n_clusters,
                                           description=cr_clustering.humanify_clustering_key(clustering_key))

def save_kmedoids_h5(f, n_clusters, kmedoids):
    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_KMEDOIDS, n_clusters)

    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
    analysis_io.save_h5(f, group, clustering_key, kmedoids)
