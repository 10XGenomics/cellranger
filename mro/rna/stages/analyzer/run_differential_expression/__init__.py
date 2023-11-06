#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.diffexp as cr_diffexp
import cellranger.analysis.io as analysis_io
import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants
import cellranger.hdf5 as cr_h5
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.analysis.singlegenome import SingleGenomeAnalysis

__MRO__ = """
stage RUN_DIFFERENTIAL_EXPRESSION(
    in  h5     matrix_h5,
    in  h5     clustering_h5,
    in  int    random_seed,
    in  int    max_clusters,
    in  bool   is_antibody_only,
    out h5     diffexp_h5,
    out path   diffexp_csv,
    src py     "stages/analyzer/run_differential_expression",
) split using (
    in  string clustering_key,
)
"""


def split(args):
    chunks = []
    matrix_mem_gb = 1.8 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    chunk_mem_gb = int(max(matrix_mem_gb, h5_constants.MIN_MEM_GB))

    threads = 1

    for key in SingleGenomeAnalysis.load_clustering_keys_from_h5(args.clustering_h5):
        chunks.append(
            {
                "clustering_key": key,
                "__mem_gb": chunk_mem_gb,
                "__threads": threads,
                "__vmem_gb": chunk_mem_gb + 8,
            }
        )

    return {"chunks": chunks, "join": {"__mem_gb": 1}}


def main(args, outs):
    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)

    # For now, only compute for gene expression features
    if args.is_antibody_only:
        matrix = matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
    else:
        matrix = matrix.select_features_by_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
    clustering = SingleGenomeAnalysis.load_clustering_from_h5(
        args.clustering_h5, args.clustering_key
    )

    diffexp = cr_diffexp.run_differential_expression(matrix, clustering.clusters)

    with analysis_io.open_h5_for_writing(outs.diffexp_h5) as f:
        cr_diffexp.save_differential_expression_h5(f, args.clustering_key, diffexp)

    cr_diffexp.save_differential_expression_csv(
        args.clustering_key, diffexp, matrix, outs.diffexp_csv
    )


def join(args, outs, chunk_defs, chunk_outs):
    if args.matrix_h5 is None:
        return

    chunk_h5s = [chunk_out.diffexp_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.diffexp_csv for chunk_out in chunk_outs]

    analysis_io.combine_h5_files(
        chunk_h5s,
        outs.diffexp_h5,
        [
            analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP,
            analysis_constants.ANALYSIS_H5_KMEANS_DIFFERENTIAL_EXPRESSION_GROUP,
        ],
    )

    feature_ref = cr_matrix.CountMatrix.load_feature_ref_from_h5_file(args.matrix_h5)

    # CELLRANGER-1792 annotate which feature indices are in the diffexp
    # This was originally stored as "feature_DE_map" in CR1.3, but that
    # structure has been subsumed by ATAC for motif z-indices (bad
    # move on Jeff's part).  Instead, add the feature list to the
    # diffexp struct explicitly.  The list should be the ordered list
    # of feature indices (in the larger FBM) used in the diffexp.
    feature_indices = feature_ref.get_indices_for_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
    if args.is_antibody_only:
        feature_indices = feature_ref.get_indices_for_type(rna_library.ANTIBODY_LIBRARY_TYPE)
    cr_h5.write_h5(
        outs.diffexp_h5,
        {analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_FEATURE_LIST: feature_indices},
        append=True,
    )

    for csv_dir in chunk_csv_dirs:
        cr_io.hardlink_with_fallback(csv_dir, outs.diffexp_csv)
