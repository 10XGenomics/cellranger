#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Convert antibody_analysis folder out to a single file aggregate_barcodes.csv out."""

import os

import cellranger.rna.library as rna_library
from cellranger.cr_io import hard_link

__MRO__ = """
stage GET_AGGREGATE_BARCODES_OUT(
    in  path   antibody_analysis,
    in  bool   is_multi,
    in  string multiplexing_method,
    out csv    aggregate_barcodes,
    src py     "stages/counter/get_aggregate_barcodes_out",
)
"""


def main(args, outs):
    make_aggregate_file = True
    if args.is_multi:
        if args.multiplexing_method is None:
            make_aggregate_file = False
        else:
            multiplexing_method = rna_library.BarcodeMultiplexingType(args.multiplexing_method)
            if multiplexing_method.is_read_multiplexed():
                make_aggregate_file = False

    if (
        args.antibody_analysis is not None
        and os.path.exists(os.path.join(args.antibody_analysis, "aggregate_barcodes.csv"))
        and make_aggregate_file
    ):
        outs.aggregate_barcodes = hard_link(
            os.path.join(args.antibody_analysis, "aggregate_barcodes.csv")
        )
    else:
        outs.aggregate_barcodes = None
