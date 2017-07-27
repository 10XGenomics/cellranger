#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.pca as cr_pca
import cellranger.analysis.bhtsne as cr_tsne
import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import tables

__MRO__ = """
stage RUN_TSNE(
    in  h5    matrix_h5,
    in  h5    pca_h5,
    in  bool  is_multi_genome,
    in  bool  skip,
    in  int   random_seed,
    in  int   perplexity,
    in  int   input_pcs,
    in  int   max_dims,
    in  int   max_iter,
    in  int   stop_lying_iter,
    in  int   mom_switch_iter,
    in  float theta,
    out h5    tsne_h5,
    out path  tsne_csv,
    src py    "stages/analyzer/run_tsne",
) split using (
    in  int   tsne_dims,
)
"""

def split(args):
    if args.skip or args.is_multi_genome:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    chunks = []
    matrix_mem_gb = cr_matrix.GeneBCMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    min_tsne_dims = cr_constants.TSNE_N_COMPONENTS
    max_tsne_dims = args.max_dims if args.max_dims is not None else min_tsne_dims
    for tsne_dims in xrange(min_tsne_dims, max_tsne_dims+1):
        chunks.append({
            'tsne_dims': tsne_dims,
            '__mem_gb': max(matrix_mem_gb, cr_constants.MIN_MEM_GB),
        })
    return {'chunks': chunks}

def main(args, outs):
    if args.skip or args.is_multi_genome:
        return

    tsne_dims = args.tsne_dims

    matrix = cr_matrix.GeneBCMatrix.load_h5(args.matrix_h5)
    pca = cr_pca.load_pca_from_h5(args.pca_h5)
    tsne = cr_tsne.run_tsne(pca.transformed_pca_matrix, input_pcs=args.input_pcs, perplexity=args.perplexity,
                     theta=args.theta, tsne_dims=tsne_dims, max_iter=args.max_iter, stop_lying_iter=args.stop_lying_iter,
                     mom_switch_iter=args.mom_switch_iter, random_state=args.random_seed)
    tsne_map = {tsne_dims: tsne}

    filters = tables.Filters(complevel = cr_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(outs.tsne_h5, 'w', filters = filters) as f:
        cr_tsne.save_tsne_h5(tsne_map, f)

    cr_tsne.save_tsne_csv(tsne_map, matrix, outs.tsne_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or args.is_multi_genome:
        return

    chunk_h5s = [chunk_out.tsne_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.tsne_csv for chunk_out in chunk_outs]
    cr_io.combine_h5_files(chunk_h5s, outs.tsne_h5, [cr_constants.ANALYSIS_H5_TSNE_GROUP])
    for csv_dir in chunk_csv_dirs:
        cr_utils.copytree(csv_dir, outs.tsne_csv, allow_existing=True)
