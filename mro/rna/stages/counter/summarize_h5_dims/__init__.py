#!/usr/bin/env python
#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Summarize the dimensions of an H5 matrix."""


import json

from cellranger.matrix import CountMatrix

__MRO__ = """
stage SUMMARIZE_H5_DIMS(
    in  string prefix,
    in  h5     matrix,
    out json   summary,
    src py     "stages/counter/summarize_h5_dims",
) using (
    mem_gb   = 2,
    volatile = strict,
)
"""


def main(args, outs):
    """Summarize H5 matrix dimensions.

    Loads matrix dimensions from an H5 file and creates a JSON summary
    with prefixed keys for features_count, barcodes_count, and nonzero_entries_count.
    """
    features_count, barcodes_count, nonzero_entries_count = CountMatrix.load_dims_from_h5(
        args.matrix
    )
    summary = {
        f"{args.prefix}_features_count": features_count,
        f"{args.prefix}_barcodes_count": barcodes_count,
        f"{args.prefix}_nonzero_entries_count": nonzero_entries_count,
    }
    with open(outs.summary, "w") as f:
        json.dump(summary, f, indent=4, sort_keys=True)
