#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

import numpy as np

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
import cellranger.analysis.kmeans as cr_kmeans
import cellranger.analysis.pca as cr_pca
import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library

__MRO__ = """
stage RUN_KMEANS(
    in  h5   matrix_h5,
    in  h5   pca_h5,
    in  int  random_seed,
    in  int  max_clusters,
    in  int  num_bcs,
    in  int  num_pcs,
    out h5   kmeans_h5,
    out path kmeans_csv,
    src py   "stages/analyzer/run_kmeans",
) split using (
    in  int  n_clusters,
    in  string  library,
)
"""

MEM_FACTOR = 1.1


def split(args):
    chunks = []
    min_clusters = analysis_constants.MIN_N_CLUSTERS
    max_clusters = (
        args.max_clusters
        if args.max_clusters is not None
        else analysis_constants.MAX_N_CLUSTERS_DEFAULT
    )

    # edge case: very small number of barcodes
    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.matrix_h5)
    max_clusters = min(max_clusters, matrix_dims[1])

    matrix_mem_gb = np.ceil(
        MEM_FACTOR * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    )

    library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(args.matrix_h5)
    allowed_types = [rna_library.GENE_EXPRESSION_LIBRARY_TYPE, rna_library.ANTIBODY_LIBRARY_TYPE]

    for library in library_types:
        if library not in allowed_types:
            # K-means on other library types not supported
            continue
        pca = cr_pca.load_pca_from_h5(args.pca_h5)
        if library == rna_library.ANTIBODY_LIBRARY_TYPE and len(pca) > 1 and pca[1] is None:
            # For antigen and gex only case, antibody capture pca object will not be generated. Skip.
            continue
        for n_clusters in range(min_clusters, max_clusters + 1):
            chunk_mem_gb = max(matrix_mem_gb, h5_constants.MIN_MEM_GB)
            chunks.append(
                {
                    "n_clusters": n_clusters,
                    "library": library,
                    "__mem_gb": chunk_mem_gb,
                }
            )

    return {"chunks": chunks}


def main(args, outs):
    np.random.seed(args.random_seed)

    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)

    pca = cr_pca.load_pca_from_h5(args.pca_h5)
    if len(pca) > 1:  # only gex- and/or ab-based pca projections are present
        if args.library == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
            pca = pca[0]
        elif args.library == rna_library.ANTIBODY_LIBRARY_TYPE:
            pca = pca[1]
        elif args.library == rna_library.ANTIGEN_LIBRARY_TYPE:
            pca = pca[2]
    else:
        pca = pca[0]  # could be either gex- or ab- or ag-based
    pca_mat = pca.transformed_pca_matrix

    # Subsample barcodes
    if args.num_bcs is not None:
        use_bcs = np.random.choice(pca_mat.shape[0], args.num_bcs, replace=False)
        matrix = matrix.select_barcodes(use_bcs)
        pca_mat = pca_mat[use_bcs, :]

    # Subset principal components
    if args.num_pcs is not None:
        pca_mat = pca_mat[:, np.arange(args.num_pcs)]

    if args.library == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
        clustering_type = cr_clustering.CLUSTER_TYPE_KMEANS
    elif args.library == rna_library.ANTIBODY_LIBRARY_TYPE:
        clustering_type = cr_clustering.CLUSTER_TYPE_ANTIBODY_KMEANS

    kmeans = cr_kmeans.run_kmeans(
        pca_mat, args.n_clusters, clustering_type, random_state=args.random_seed
    )

    with analysis_io.open_h5_for_writing(outs.kmeans_h5) as f:
        cr_kmeans.save_kmeans_h5(f, args.n_clusters, kmeans, clustering_type)

    clustering_key = cr_clustering.format_clustering_key(clustering_type, args.n_clusters)
    cr_clustering.save_clustering_csv(outs.kmeans_csv, clustering_key, kmeans.clusters, matrix.bcs)


def join(args, outs, chunk_defs, chunk_outs):
    chunk_h5s = [chunk_out.kmeans_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.kmeans_csv for chunk_out in chunk_outs]

    analysis_io.combine_h5_files(
        chunk_h5s,
        outs.kmeans_h5,
        [
            analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP,
            cr_clustering.CLUSTER_TYPE_KMEANS,
            cr_clustering.CLUSTER_TYPE_ANTIBODY_KMEANS,
        ],
    )

    for csv_dir in chunk_csv_dirs:
        cr_io.hardlink_with_fallback(csv_dir, outs.kmeans_csv)
