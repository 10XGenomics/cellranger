#!/usr/bin/env python3
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Tool for converting feature-barcode matrices from sparse format to dense.

CSV format, for use by external programs.

The commands below should be preceded by '{cmd}':

Usage:
    mat2csv <input_path> <output_csv> [--genome=GENOME]
    mat2csv -h | --help | --version

Arguments:
    input_path          Path to a {product} feature-barcode matrix. Can be
                            either a feature-barcode h5 file (recommended) or a
                            path to a MEX {product} output folder.
    output_csv          Output CSV file.

Options:
    --genome=GENOME     Specify which genome to extract. This only applies to
                            multi-genome h5 input files.
    -h --help           Show this message.
    --version           Show version.
"""

from __future__ import annotations

import os
import pathlib
import sys

import docopt

import cellranger.cr_io as cr_io
from cellranger.matrix import CountMatrix
from cellranger.mtx_to_matrix_converter import load_mtx, save_dense_csv
from cellranger.products import get_cmd_names


def _parse_args(product_name):
    product, cmd = get_cmd_names(product_name)

    version = "{} {} {}\n{}".format(
        product_name,
        os.getenv("TENX_SUBCMD", ""),
        os.getenv("TENX_VERSION", ""),
        os.getenv("TENX_COPYRIGHT", ""),
    )
    return docopt.docopt(__doc__.format(cmd=cmd, product=product), version=version)


def main():
    args = _parse_args(os.getenv("TENX_PRODUCT", ""))

    output_csv = pathlib.Path(cr_io.get_output_path(args["<output_csv>"]))
    input_path = args["<input_path>"]
    genome = args["--genome"]

    if input_path.endswith(".h5"):
        input_path = cr_io.get_input_path(input_path)
        gbm = CountMatrix.load_h5_file(input_path)
    else:
        input_path = pathlib.Path(cr_io.get_input_path(input_path, is_dir=True))
        gbm = load_mtx(input_path)
        if genome is not None:
            sys.exit(
                "The '--genome' argument can only be use with .h5 input files, "
                "not with MEX directories"
            )

    if genome is None:
        matrix = gbm
    else:
        genomes = gbm.get_genomes()
        if genome not in genomes:
            sys.exit(f"Genome '{genome}' not found (genomes available: {genomes})")
        matrix = gbm.select_features_by_genome(genome)

    num_features, num_barcodes, num_entries = (
        matrix.features_dim,
        matrix.bcs_dim,
        matrix.get_num_nonzero(),
    )
    dense_size = num_features * num_barcodes
    zero_frac = float(dense_size - num_entries) * 100.0 / float(dense_size)
    print(
        """
    WARNING: this matrix has %d x %d (%d total) elements, %f%% of which are zero.
    Converting it to dense CSV format may be very slow and memory intensive.
    Moreover, other programs (e.g. Excel) may be unable to load it due to its size.
    To cancel this command, press <control key> + C.

    If you need to inspect the data, we recommend using Loupe Browser.
    """
        % (num_features, num_barcodes, dense_size, zero_frac)
    )
    sys.stdout.flush()

    try:
        save_dense_csv(matrix, output_csv)
    except KeyboardInterrupt:
        if output_csv.exists():
            output_csv.unlink()


if __name__ == "__main__":
    main()
