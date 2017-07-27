#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import numpy as np

__MRO__ = """
stage PREPROCESS_MATRIX(
    in  h5   matrix_h5,
    in  bool skip,
    in  int  random_seed,
    in  csv  use_genes,
    in  csv  use_bcs,
    in  int  num_bcs,
    in  int  force_cells,
    out h5   preprocessed_matrix_h5,
    out bool is_multi_genome,
    src py   "stages/analyzer/preprocess_matrix",
) split using (
)
"""

def split(args):
    if args.skip:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    matrix_mem_gb = cr_matrix.GeneBCMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    chunks = [{
        '__mem_gb': max(matrix_mem_gb, cr_constants.MIN_MEM_GB)
    }]
    return {'chunks': chunks}

def main(args, outs):
    if args.skip:
        return

    if args.random_seed is not None:
        np.random.seed(args.random_seed)

    # detect barnyard
    genomes = cr_matrix.GeneBCMatrices.load_genomes_from_h5(args.matrix_h5)
    if len(genomes) > 1:
        outs.is_multi_genome = True
        cr_utils.copy(args.matrix_h5, outs.preprocessed_matrix_h5)
        return
    else:
        outs.is_multi_genome = False

    genome = genomes[0]
    matrix = cr_matrix.GeneBCMatrices.load_h5(args.matrix_h5).get_matrix(genome)
    matrix = cr_matrix.GeneBCMatrix.preprocess_matrix(matrix, num_bcs=args.num_bcs, use_bcs=args.use_bcs, use_genes=args.use_genes, force_cells=args.force_cells)

    gbm = cr_matrix.GeneBCMatrices()
    gbm.matrices[genome] = matrix
    matrix_attrs = cr_matrix.get_matrix_attrs(args.matrix_h5)
    gbm.save_h5(outs.preprocessed_matrix_h5, extra_attrs=matrix_attrs)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        return

    chunk_out = chunk_outs[0]
    cr_utils.copy(chunk_out.preprocessed_matrix_h5, outs.preprocessed_matrix_h5)
    outs.is_multi_genome = chunk_out.is_multi_genome
