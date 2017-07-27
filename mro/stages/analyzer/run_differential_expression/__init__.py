#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.diffexp as cr_diffexp
import cellranger.analysis.io as cr_io
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils

__MRO__ = """
stage RUN_DIFFERENTIAL_EXPRESSION(
    in  h5     matrix_h5,
    in  h5     clustering_h5,
    in  bool   is_multi_genome,
    in  bool   skip,
    in  int    random_seed,
    in  int    max_clusters,
    out h5     diffexp_h5,
    out path   diffexp_csv,
    src py     "stages/analyzer/run_differential_expression",
) split using (
    in  string clustering_key,
)
"""

def split(args):
    if args.skip or args.is_multi_genome:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    chunks = []
    matrix_mem_gb = cr_matrix.GeneBCMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    chunk_mem_gb = max(matrix_mem_gb * 3, cr_constants.MIN_MEM_GB)

    # HACK - give big jobs more threads in order to avoid overloading a node
    threads = cr_utils.get_thread_request_from_mem_gb(chunk_mem_gb)

    for key in SingleGenomeAnalysis.load_clustering_keys_from_h5(args.clustering_h5):
        chunks.append({
            'clustering_key': key,
            '__mem_gb': chunk_mem_gb,
            '__threads': threads,
        })

    return {'chunks': chunks}

def main(args, outs):
    if args.skip or args.is_multi_genome:
        return

    matrix = cr_matrix.GeneBCMatrix.load_h5(args.matrix_h5)

    clustering = SingleGenomeAnalysis.load_clustering_from_h5(args.clustering_h5, args.clustering_key)

    diffexp = cr_diffexp.run_differential_expression(matrix, clustering.clusters)

    with cr_io.open_h5_for_writing(outs.diffexp_h5) as f:
        cr_diffexp.save_differential_expression_h5(f, args.clustering_key, diffexp)

    cr_diffexp.save_differential_expression_csv(args.clustering_key, diffexp, matrix, outs.diffexp_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or args.is_multi_genome:
        return

    chunk_h5s = [chunk_out.diffexp_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.diffexp_csv for chunk_out in chunk_outs]

    cr_io.combine_h5_files(chunk_h5s, outs.diffexp_h5, [cr_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP,
                                                        cr_constants.ANALYSIS_H5_KMEANS_DIFFERENTIAL_EXPRESSION_GROUP])

    for csv_dir in chunk_csv_dirs:
        cr_utils.copytree(csv_dir, outs.diffexp_csv, allow_existing=True)
