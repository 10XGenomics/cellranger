#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.stats as cr_stats
import cellranger.io as cr_io
import cellranger.utils as cr_utils
import martian
import numpy as np
from cellranger.logperf import LogPerf

__MRO__ = """
stage PREPROCESS_MATRIX(
    in  h5   matrix_h5,
    in  bool skip,
    in  int  random_seed,
    in  csv  use_genes,
    in  csv  exclude_genes,
    in  csv  use_bcs,
    in  int  num_bcs,
    in  int  force_cells,
    out h5   cloupe_matrix_h5,
    out h5   preprocessed_matrix_h5,
    out bool is_multi_genome,
    src py   "stages/analyzer/preprocess_matrix",
) split using (
)
"""

def select_barcodes_and_features(matrix, num_bcs=None, use_bcs=None, use_genes=None,
                                 exclude_genes=None, force_cells=None):

    if force_cells is not None:
        bc_counts = matrix.get_counts_per_bc()
        bc_indices, _, _ = cr_stats.filter_cellular_barcodes_fixed_cutoff(bc_counts, force_cells)
        with LogPerf('f1'):
            matrix = matrix.select_barcodes(bc_indices)

    elif use_bcs is not None:
        bc_seqs = cr_io.load_csv_rownames(use_bcs)
        bc_indices = matrix.bcs_to_ints(bc_seqs)
        with LogPerf('f2'):
            matrix = matrix.select_barcodes(bc_indices)

    elif num_bcs is not None and num_bcs < matrix.bcs_dim:
        bc_indices = np.sort(np.random.choice(np.arange(matrix.bcs_dim), size=num_bcs, replace=False))
        with LogPerf('f3'):
            matrix = matrix.select_barcodes(bc_indices)

    include_indices = list(range(matrix.features_dim))
    if use_genes is not None:
        include_ids = cr_io.load_csv_rownames(use_genes)
        with LogPerf('f4'):
            try:
                include_indices = matrix.feature_ids_to_ints(include_ids)
            except KeyError as e:
                err_message = str(e).strip("'")
                err_message = err_message.replace("Specified feature ID",
                                                  "Feature ID specified in genes_csv")
                martian.exit(err_message)

    exclude_indices = []
    if exclude_genes is not None:
        exclude_ids = cr_io.load_csv_rownames(exclude_genes)
        with LogPerf('f5'):
            try:
                exclude_indices = matrix.feature_ids_to_ints(exclude_ids)
            except KeyError as e:
                err_message = str(e).strip("'")
                err_message = err_message.replace("Specified feature ID",
                                                  "Feature ID specified in exclude_genes_csv")
                martian.exit(err_message)

    gene_indices = np.array(sorted(list(set(include_indices) - set(exclude_indices))), dtype=int)
    with LogPerf('ff'):
        matrix = matrix.select_features(gene_indices)

    return matrix


def split(args):
    if args.skip:
        return {'chunks': []}

    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    return {
        'chunks': [], 
        'join': { '__mem_gb': max(matrix_mem_gb, h5_constants.MIN_MEM_GB) }
    }


def main(args, outs):
    pass

def join(args, outs, chunk_defs, chunk_outs):

    if args.skip:
        outs.cloupe_matrix_h5 = None
        outs.preprocessed_matrix_h5 = None
        outs.is_multi_genome = False
        return

    if args.random_seed is not None:
        np.random.seed(args.random_seed)

    # detect barnyard
    genomes = cr_matrix.CountMatrix.get_genomes_from_h5(args.matrix_h5)
    if len(genomes) > 1:
        outs.is_multi_genome = True
    else:
        outs.is_multi_genome = False

    with LogPerf('load'):
        matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)

    with LogPerf('select'):
        matrix = select_barcodes_and_features(matrix,
            num_bcs=args.num_bcs,
            use_bcs=args.use_bcs,
            use_genes=args.use_genes,
            exclude_genes=args.exclude_genes,
            force_cells=args.force_cells)


    # Preserve original matrix attributes
    matrix_attrs = cr_matrix.get_matrix_attrs(args.matrix_h5)

    # gem groups are needed for cloupe, and older versions of cellranger count
    # may not have added those to the matrix_attrs
    if cr_matrix.get_gem_group_index(args.matrix_h5) is None:
        chemistry = cr_matrix.CountMatrix.load_chemistry_from_h5(args.matrix_h5)

        bcs = cr_matrix.CountMatrix.load_bcs_from_h5(args.matrix_h5)
        gem_groups = list(set(cr_utils.split_barcode_seq(bc)[1] for bc in bcs))

        count_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, gem_groups, chemistry)
        matrix_attrs.update(count_attrs)

    # matrix h5 for cloupe (gene with zero count preserved)
    # this will only be used in reanalyzer, where user could include/exclude genes
    with LogPerf('w1'):
        matrix.save_h5_file(outs.cloupe_matrix_h5, extra_attrs=matrix_attrs)
    
    with LogPerf('selnz'):
        matrix, _, _ = matrix.select_nonzero_axes()

    with LogPerf('w2'):
        matrix.save_h5_file(outs.preprocessed_matrix_h5, extra_attrs=matrix_attrs)
