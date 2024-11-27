#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

import os

import martian

import cellranger.hdf5 as cr_h5
import cellranger.matrix as cr_matrix

__MRO__ = """
stage REANALYZER_PREFLIGHT(
    in  h5 filtered_matrices_h5,
    src py "stages/analyzer/reanalyzer_preflight",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if not (args.filtered_matrices_h5 and os.path.exists(args.filtered_matrices_h5)):
        martian.exit(f"Filtered matrices do not exist: {args.filtered_matrices_h5}")

    if not os.access(args.filtered_matrices_h5, os.R_OK):
        martian.exit(
            f"Filtered matrices file is not readable, please check file permissions: {args.filtered_matrices_h5}"
        )

    h5_filetype = cr_h5.get_h5_filetype(args.filtered_matrices_h5)
    if h5_filetype and h5_filetype != cr_matrix.MATRIX_H5_FILETYPE:
        martian.exit(f"Input is a {h5_filetype} file, but a matrix file is required")

    h5_version = cr_matrix.CountMatrix.get_format_version_from_h5(args.filtered_matrices_h5)
    if h5_version > cr_matrix.MATRIX_H5_VERSION:
        martian.exit(
            "Filtered matrices file format version (%d) "
            "is newer than this version of the software." % h5_version
        )

    if cr_matrix.get_gem_group_index(args.filtered_matrices_h5) is None:
        martian.exit(
            "Filtered matrices file was generated with an older version "
            "of cellranger that is incompatible with reanalyze. Please run "
            "cellranger count again to generate a new matrix."
        )

    flt_genomes = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrices_h5)

    if len(flt_genomes) > 1:
        martian.exit(
            f"Reanalyzer only supports matrices with one genome. This matrix has: {flt_genomes}"
        )
    if len(flt_genomes) == 0:
        martian.log_info("Only Antibody Capture library detected")
