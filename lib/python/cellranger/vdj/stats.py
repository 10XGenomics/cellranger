#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import sklearn.cluster as sk_cluster
import tenkit.stats as tk_stats

def compute_readpairs_per_umi_threshold(reads, subsample_rate):
    ''' Compute a threshold above which the UMIs are unlikely to be PCR off-products.
        reads (np.array(int)) - Read pairs for each UMI
        subsample_rate (float) - Subsample reads to this fraction.
        Returns threshold (int) - The RPPU threshold in the subsampled space '''

    if len(np.unique(reads)) < 2:
        print 'Skipping RPPU threshold calculation.'
        return 1

    print 'RPPU subsample rate: %0.4f' % subsample_rate

    reads = np.random.binomial(reads, subsample_rate)
    reads = reads[reads > 0]

    if len(np.unique(reads)) < 2:
        print 'Subsampling gave a degenerate distribution of RPPU. Skipping RPPU threshold calculation.'
        return 1

    new_n50 = tk_stats.NX(reads, 0.5)

    print 'New N50: %d:' % new_n50

    # Log-transform counts
    log_reads = np.log(reads)

    # Run K-Means. Reshape necessary because kmeans takes a matrix.
    kmeans = sk_cluster.KMeans(2).fit(log_reads.reshape((-1,1)))
    kmeans.predict(log_reads.reshape((-1,1)))

    # Take the cluster with the smallest mean
    min_cluster = np.argsort(np.ravel(kmeans.cluster_centers_))[0]

    print 'RPPU component means: ' + str(list(iter(np.exp(kmeans.cluster_centers_))))
    print 'RPPU component members: ' + str(np.bincount(kmeans.labels_))

    # Take the max element in the min-cluster
    threshold = np.max(reads[kmeans.labels_ == min_cluster])

    return threshold
