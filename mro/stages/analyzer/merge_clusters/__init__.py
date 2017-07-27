#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Merge similar clusters produced by graph-clustering

import numpy as np
import pandas as pd
import resource
import sys
from scipy.cluster.hierarchy import linkage
import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.io as cr_io
from cellranger.analysis.diffexp import compute_sseq_params, sseq_differential_expression
import cellranger.analysis.graphclust as cr_graphclust
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
from cellranger.matrix import GeneBCMatrices, GeneBCMatrix

__MRO__ = """
stage MERGE_CLUSTERS(
    in  h5   matrix_h5,
    in  h5   pca_h5,
    in  h5   clusters_h5,
    in  bool skip,
    in  bool is_multi_genome,
    out h5   clusters_h5,
    out path clusters_csv,
    src py   "stages/analyzer/merge_clusters",
) split using (
)
"""

# Adjusted p-value threshold for calling genes DE
MERGE_CLUSTERS_DE_ADJ_P_THRESHOLD = 0.05

# Do not merge cluster pairs if they have enough DE genes
MIN_DE_GENES = 1

def split(args):
    if args.skip or args.is_multi_genome:
        return {'chunks': [{}]}

    return {
        'chunks': [{}],
        'join': {
            '__mem_gb': 2 * GeneBCMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5),
        }
    }

def main(args, outs):
    pass

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or args.is_multi_genome:
        return

    np.random.seed(0)

    # Load the matrix
    matrices = GeneBCMatrices.load_h5(args.matrix_h5)
    if len(matrices.matrices) > 1:
        print "Multiple genomes detected. Skipping stage."
        return
    mat = matrices.matrices.values()[0]
    print mat.m.shape, mat.m.nnz

    barcodes = mat.bcs

    # Load graph-based clustering from analysis H5
    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_GRAPHCLUST, 0)
    clustering = SingleGenomeAnalysis.load_clustering_from_h5(args.clusters_h5, clustering_key)
    labels = clustering.clusters

    # Clusters that were 0 were unused in the clustering analysis (only relevant if the cluster stage was run by itself)
    total_bcs = len(labels)
    use_bcs = np.flatnonzero(labels > 0)
    expr_mat = mat.m[:,use_bcs]

    # Make cluster labels 0-based
    labels = labels[use_bcs] - 1

    # Convert PCA coords to dataframe
    pca = SingleGenomeAnalysis.load_pca_from_h5(args.pca_h5).transformed_pca_matrix[use_bcs,:]
    pca_df = pd.DataFrame(pca)
    print pca_df.shape

    # 1) Run hierarchical clustering on cluster medoids in PCA-space
    # 2) For each pair of clusters that are sibling leaves,
    #   3) Run a differential expression analysis
    #   4) Merge the clusters if not enough genes are differentially expressed
    #   5) If merged, stop considering cluster-pairs and goto 1)

    # Cache already-checked cluster-pairs
    # set of (frozenset, frozenset)
    checked_cluster_pairs = set()

    while True:
        print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        sys.stdout.flush()
        if len(np.bincount(labels)) == 1:
            # One cluster remains
            break

        # Compute medoids, perform hierarchical clustering
        pca_df['cluster'] = labels
        medoids = pca_df.groupby('cluster').apply(lambda x: x.median(axis=0)).as_matrix()[:,:-1]
        hc = linkage(medoids, 'complete')
        max_label = np.max(labels)

        print np.bincount(labels)
        print 'max=%d' % max_label

        any_merged = False
        for step in xrange(hc.shape[0]):
            if hc[step,0] <= max_label and hc[step,1] <= max_label:
                leaf0, leaf1 = hc[step,0], hc[step,1]

                group0 = np.flatnonzero(labels == leaf0)
                group1 = np.flatnonzero(labels == leaf1)

                # Skip this cluster pair if already checked
                set0 = frozenset(group0)
                set1 = frozenset(group1)
                cluster_pair = tuple(sorted([set0, set1]))
                if cluster_pair in checked_cluster_pairs:
                    continue
                checked_cluster_pairs.add(cluster_pair)

                print 'Comparing clusters (%d,%d)' % (1+leaf0,1+leaf1)
                submat = expr_mat[:,np.concatenate((group0,group1))]

                print '\tComputing params on (%d,%d) matrix' % submat.shape
                params = compute_sseq_params(submat)

                print '\tRunning DE on %d vs %d cells' % (len(group0), len(group1))
                group0_submat = np.arange(len(group0))
                group1_submat = np.arange(len(group0), len(group0)+len(group1))
                de_result = sseq_differential_expression(submat,
                                                         group0_submat, group1_submat,
                                                         params)

                n_de_genes = np.sum(de_result.adjusted_p_value < MERGE_CLUSTERS_DE_ADJ_P_THRESHOLD)
                if n_de_genes == 0:
                    print '\tFound %d DE genes. Merging clusters (%d,%d)' % (n_de_genes, 1+leaf0,1+leaf1)
                    # Relabel as the smaller-index cluster
                    labels[labels == leaf1] = leaf0

                    # Shift all labels above old label down
                    labels[labels > leaf1] = labels[labels > leaf1] - 1

                    any_merged = True
                    break

        sys.stdout.flush()

        if not any_merged:
            break

    # Convert back to one-based cluster labels
    labels += 1

    labels = cr_clustering.relabel_by_size(labels)

    # Convert back into original bc space, with 0s for unused bcs
    final_labels = np.zeros(total_bcs, dtype=int)
    final_labels[use_bcs] = labels

    # Save results
    with cr_io.open_h5_for_writing(outs.clusters_h5) as f:
        cr_graphclust.save_graphclust_h5(f, final_labels)

    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_GRAPHCLUST, 0)

    cr_clustering.save_clustering_csv(outs.clusters_csv, clustering_key, final_labels, barcodes)
