#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.kmeans as cr_kmeans
import cellranger.analysis.pca as cr_pca
import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils

__MRO__ = """
stage RUN_KMEANS(
    in  h5   matrix_h5,
    in  h5   pca_h5,
    in  bool is_multi_genome,
    in  bool skip,
    in  int  random_seed,
    in  int  max_clusters,
    in  int  num_bcs,
    in  int  num_pcs,
    out h5   kmeans_h5,
    out path kmeans_csv,
    src py   "stages/analyzer/run_kmeans",
) split using (
    in  int  n_clusters,
)
"""

MEM_FACTOR = 1.1

def split(args):
    if args.skip or args.is_multi_genome:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    chunks = []
    min_clusters = cr_constants.MIN_N_CLUSTERS
    max_clusters = args.max_clusters if args.max_clusters is not None else cr_constants.MAX_N_CLUSTERS_DEFAULT
    matrix_mem_gb = np.ceil(MEM_FACTOR * cr_matrix.GeneBCMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5))
    for n_clusters in xrange(min_clusters, max_clusters + 1):
        chunk_mem_gb = max(matrix_mem_gb, cr_constants.MIN_MEM_GB)
        chunks.append({
            'n_clusters': n_clusters,
            '__mem_gb': chunk_mem_gb,
        })

    return {'chunks': chunks}

def main(args, outs):
    np.random.seed(args.random_seed)

    if args.skip or args.is_multi_genome:
        return

    matrix = cr_matrix.GeneBCMatrix.load_h5(args.matrix_h5)
    pca = cr_pca.load_pca_from_h5(args.pca_h5)
    pca_mat = pca.transformed_pca_matrix

    # Subsample barcodes
    if args.num_bcs is not None:
        use_bcs = np.random.choice(pca_mat.shape[0], args.num_bcs, replace=False)
        matrix = matrix.select_barcodes(use_bcs)
        pca_mat = pca_mat[use_bcs,:]

    # Subset principal components
    if args.num_pcs is not None:
        pca_mat = pca_mat[:,np.arange(args.num_pcs)]

    kmeans = cr_kmeans.run_kmeans(pca_mat, args.n_clusters, random_state=args.random_seed)

    with cr_io.open_h5_for_writing(outs.kmeans_h5) as f:
        cr_kmeans.save_kmeans_h5(f, args.n_clusters, kmeans)

    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_KMEANS, args.n_clusters)

    cr_clustering.save_clustering_csv(outs.kmeans_csv, clustering_key, kmeans.clusters, matrix.bcs)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or args.is_multi_genome:
        return

    chunk_h5s = [chunk_out.kmeans_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.kmeans_csv for chunk_out in chunk_outs]

    cr_io.combine_h5_files(chunk_h5s, outs.kmeans_h5, [cr_constants.ANALYSIS_H5_CLUSTERING_GROUP,
                                                       cr_constants.ANALYSIS_H5_KMEANS_GROUP])

    for csv_dir in chunk_csv_dirs:
        cr_utils.copytree(csv_dir, outs.kmeans_csv, allow_existing=True)
