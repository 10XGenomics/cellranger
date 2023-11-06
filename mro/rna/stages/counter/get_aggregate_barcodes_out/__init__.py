#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Convert antibody_analysis folder out to a single file aggregate_barcodes.csv out."""

import os

from cellranger.cr_io import hard_link

__MRO__ = """
stage GET_AGGREGATE_BARCODES_OUT(
    in  path antibody_analysis,
    out csv  aggregate_barcodes,
    src py   "stages/counter/get_aggregate_barcodes_out",
)
"""


def main(args, outs):
    if args.antibody_analysis is not None and os.path.exists(
        os.path.join(args.antibody_analysis, "aggregate_barcodes.csv")
    ):
        outs.aggregate_barcodes = hard_link(
            os.path.join(args.antibody_analysis, "aggregate_barcodes.csv")
        )
    else:
        outs.aggregate_barcodes = None
