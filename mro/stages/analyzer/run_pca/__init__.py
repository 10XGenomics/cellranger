#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.pca as cr_pca
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import tables

__MRO__ = """
stage RUN_PCA(
    in  h5   matrix_h5,
    in  bool is_multi_genome,
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
    if args.skip or args.is_multi_genome:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    # memory usage of just loading the full matrix
    matrix_dims = cr_matrix.GeneBCMatrices.load_dims_from_h5(args.matrix_h5)
    (genes_dim, bcs_dim, nonzero_entries) = matrix_dims.values()[0]
    matrix_mem_gb = cr_matrix.GeneBCMatrix.get_mem_gb_from_matrix_dim(nonzero_entries)

    # TODO adjust based on num PCs, etc?
    irlb_mem_gb = cr_pca.get_irlb_mem_gb_from_matrix_dim(nonzero_entries)

    mem_gb = max(irlb_mem_gb, matrix_mem_gb, cr_constants.MIN_MEM_GB)

    # HACK - give big jobs more threads in order to avoid overloading a node
    threads = cr_utils.get_thread_request_from_mem_gb(mem_gb)

    chunks = [{
        '__mem_gb': mem_gb,
        '__threads': threads,
    }]
    return {'chunks': chunks}

def main(args, outs):
    if args.skip or args.is_multi_genome:
        return

    matrix = cr_matrix.GeneBCMatrix.load_h5(args.matrix_h5)
    pca = cr_pca.run_pca(matrix, pca_genes=args.num_genes, pca_bcs=args.num_bcs,
                 n_pca_components=args.num_pcs, random_state=args.random_seed)
    pca_key = args.num_pcs if args.num_pcs is not None else cr_constants.PCA_N_COMPONENTS_DEFAULT
    pca_map = {pca_key: pca}

    filters = tables.Filters(complevel = cr_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(outs.pca_h5, 'w', filters = filters) as f:
        cr_pca.save_pca_h5(pca_map, f)

    cr_pca.save_pca_csv(pca_map, matrix, outs.pca_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or args.is_multi_genome:
        return

    chunk_out = chunk_outs[0]
    cr_utils.copy(chunk_out.pca_h5, outs.pca_h5)
    cr_utils.copytree(chunk_out.pca_csv, outs.pca_csv)
