#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Stage checking if bcs in cloupe and feature-bc matrix are the same."""

import subprocess

import martian

import cellranger.matrix as cr_matrix
import tenkit.log_subprocess as tk_subproc

__MRO__ = """
stage CHECK_CLOUPE_MATRIX_CONSISTENT(
    in  cloupe sample_cloupe,
    in  h5     filtered_matrix,
    out string alert_string,
    src py     "stages/cas_cell_typing/check_cloupe_matrix_consistent",
)
"""


def main(args, outs):
    if not args.sample_cloupe or not args.filtered_matrix:
        martian.clear(outs)
        return

    tmp_bc_path = martian.make_path("bcs.txt").decode()

    call = [
        "getbarcodes",
        args.sample_cloupe,
        tmp_bc_path,
    ]

    unicode_call = [arg.encode("utf-8") for arg in call]

    martian.log_info("Running getbarcodes: {}".format(" ".join(call)))
    try:
        results = tk_subproc.check_output(unicode_call, stderr=subprocess.STDOUT)
        martian.log_info(f"getbarcodes output: {results}")
    except subprocess.CalledProcessError as err:
        martian.clear(outs)
        martian.throw(
            f"Could not run getbarcodes to get barcodes from cloupe file. Error: \n{err.output}"
        )

    with open(tmp_bc_path) as f:
        barcodes_from_cloupe = set(x.strip() for x in f.readlines())

    barcodes_from_matrix = set(
        x.decode() for x in cr_matrix.CountMatrix.load_bcs_from_h5(args.filtered_matrix)
    )

    if barcodes_from_matrix != barcodes_from_cloupe:
        common_bcs = len(barcodes_from_matrix.intersection(barcodes_from_cloupe))
        bcs_in_matrix_only = len(barcodes_from_matrix - barcodes_from_cloupe)
        bcs_in_cloupe_only = len(barcodes_from_cloupe - barcodes_from_matrix)
        total_bcs = len(barcodes_from_matrix.union(barcodes_from_cloupe))
        outs.alert_string = (
            f"{common_bcs} out of {total_bcs} are seen in both the feature barcode matrix and cloupe file. "
            f"{bcs_in_matrix_only} seen in the feature barcode matrix only. "
            f"{bcs_in_cloupe_only} seen in the cloupe file only. "
        )
    else:
        outs.alert_string = None
