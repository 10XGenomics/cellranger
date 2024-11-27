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
import cellranger.rna.library as rna_library
from cellranger.matrix import CountMatrix

__MRO__ = """
stage RUN_KMEANS(
    in  h5     matrix_h5,
    in  h5     pca_h5,
    in  int    random_seed,
    in  int    max_clusters,
    in  int    num_bcs,
    in  int    num_pcs,
    out h5     kmeans_h5,
    out path   kmeans_csv,
    src py     "stages/analyzer/run_kmeans",
) split (
    in  int    n_clusters,
    in  string library,
    in  int    which_pca,
) using (
    mem_gb   = 3,
    volatile = strict,
)
"""


def split(args):
    pcas = cr_pca.load_pca_from_h5(args.pca_h5)
    num_barcodes, num_pcs = pcas[0].transformed_pca_matrix.shape
    mem_gib = 1 + round(34 * num_pcs * num_barcodes / 1024**3, 1)
    print(f"{num_barcodes=},{num_pcs=},{mem_gib=}")

    max_clusters = min(
        (
            args.max_clusters
            if args.max_clusters is not None
            else analysis_constants.MAX_N_CLUSTERS_DEFAULT
        ),
        # edge case: very small number of barcodes
        num_barcodes,
    )

    library_types = CountMatrix.load_library_types_from_h5_file(args.matrix_h5)
    library_types = [
        x
        for x in [rna_library.GENE_EXPRESSION_LIBRARY_TYPE, rna_library.ANTIBODY_LIBRARY_TYPE]
        if x in library_types
    ]

    return {
        "chunks": [
            {
                "library": library_type,
                "which_pca": which_pca,
                "n_clusters": n_clusters,
                "__mem_gb": mem_gib,
            }
            for which_pca, (library_type, pca) in enumerate(zip(library_types, pcas))
            for n_clusters in range(analysis_constants.MIN_N_CLUSTERS, 1 + max_clusters)
            if pca is not None
        ],
        "join": {"__mem_gb": 1},
    }


def main(args, outs):
    np.random.seed(args.random_seed)

    pca_mat = cr_pca.load_pca_from_h5(args.pca_h5)[args.which_pca].transformed_pca_matrix

    barcodes = CountMatrix.load_bcs_from_h5(args.matrix_h5)
    if args.num_bcs is not None:
        # Subsample barcodes
        use_bcs = np.random.choice(pca_mat.shape[0], args.num_bcs, replace=False)
        use_bcs.sort()
        barcodes = barcodes[use_bcs]
        pca_mat = pca_mat[use_bcs, :]

    # Subset principal components
    if args.num_pcs is not None:
        pca_mat = pca_mat[:, np.arange(args.num_pcs)]

    clustering_type = {
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE: cr_clustering.CLUSTER_TYPE_KMEANS,
        rna_library.ANTIBODY_LIBRARY_TYPE: cr_clustering.CLUSTER_TYPE_ANTIBODY_KMEANS,
    }[args.library]
    kmeans = cr_kmeans.run_kmeans(
        pca_mat, args.n_clusters, clustering_type, random_state=args.random_seed
    )

    with analysis_io.open_h5_for_writing(outs.kmeans_h5) as f:
        cr_kmeans.save_kmeans_h5(f, args.n_clusters, kmeans, clustering_type)

    clustering_key = cr_clustering.format_clustering_key(clustering_type, args.n_clusters)
    cr_clustering.save_clustering_csv(outs.kmeans_csv, clustering_key, kmeans.clusters, barcodes)


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
