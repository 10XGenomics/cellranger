#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import martian
import os

import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils

__MRO__ = """
stage REANALYZER_PREFLIGHT(
    in h5 filtered_matrices_h5,
    src py "stages/analyzer/reanalyzer_preflight",
)
"""
def main(args, outs):
    if not (args.filtered_matrices_h5 and os.path.exists(args.filtered_matrices_h5)):
        martian.exit("Filtered matrices do not exist: %s" % args.filtered_matrices_h5)

    if not os.access(args.filtered_matrices_h5, os.R_OK):
        martian.exit("Filtered matrices file is not readable, please check file permissions: %s" % args.filtered_matrices_h5)

    h5_filetype = cr_utils.get_h5_filetype(args.filtered_matrices_h5)
    if h5_filetype and h5_filetype != cr_matrix.MATRIX_H5_FILETYPE:
        martian.exit("Input is a %s file, but a matrix file is required" % h5_filetype)

    flt_genomes = cr_matrix.GeneBCMatrices.load_genomes_from_h5(args.filtered_matrices_h5)

    if len(flt_genomes) != 1:
        martian.exit("Reanalyzer only supports matrices with one genome. This matrix has: %s" % flt_genomes)
