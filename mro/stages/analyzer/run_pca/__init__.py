#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.pca as cr_pca
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io
import tables

__MRO__ = """
stage RUN_PCA(
    in  h5   matrix_h5,
    in  bool skip,
    in  int  random_seed,
    in  int  num_bcs,
    in  int  num_genes,
    in  int  num_pcs,
    out h5   pca_h5,
    out path pca_csv,
    src py   "stages/analyzer/run_pca",
) split using (
)
"""

def split(args):
    if args.skip:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    # memory usage of just loading the full matrix

    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.matrix_h5)
    (features_dim, bcs_dim, nonzero_entries) = matrix_dims
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_dim(bcs_dim, nonzero_entries)

    # TODO adjust based on num PCs, etc?
    irlb_mem_gb = cr_pca.get_irlb_mem_gb_from_matrix_dim(nonzero_entries)

    mem_gb = max(irlb_mem_gb, matrix_mem_gb, h5_constants.MIN_MEM_GB)

    # HACK - give big jobs more threads in order to avoid overloading a node
    threads = cr_io.get_thread_request_from_mem_gb(mem_gb)

    chunks = [{
        '__mem_gb': mem_gb,
        '__threads': threads,
    }]
    return {'chunks': chunks, 'join': {'__mem_gb': 1}}

def main(args, outs):
    if args.skip:
        return

    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
    pca = cr_pca.run_pca(matrix, pca_features=args.num_genes, pca_bcs=args.num_bcs,
                         n_pca_components=args.num_pcs, random_state=args.random_seed)
    pca_key = args.num_pcs if args.num_pcs is not None else analysis_constants.PCA_N_COMPONENTS_DEFAULT
    pca_map = {pca_key: pca}

    filters = tables.Filters(complevel = h5_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(outs.pca_h5, 'w', filters = filters) as f:
        cr_pca.save_pca_h5(pca_map, f)

    cr_pca.save_pca_csv(pca_map, matrix, outs.pca_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        return

    chunk_out = chunk_outs[0]
    cr_io.copy(chunk_out.pca_h5, outs.pca_h5)
    cr_io.copytree(chunk_out.pca_csv, outs.pca_csv)
