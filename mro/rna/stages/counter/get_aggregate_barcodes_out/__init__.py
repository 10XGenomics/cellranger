#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Convert antibody_analysis folder out to a single file aggregate_barcodes.csv out."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import martian

import cellranger.rna.library as rna_library
from cellranger.cr_io import hard_link

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage GET_AGGREGATE_BARCODES_OUT(
    in  path   antibody_analysis,
    in  bool   is_multi,
    in  string multiplexing_method,
    out csv    aggregate_barcodes,
    src py     "stages/counter/get_aggregate_barcodes_out",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """GET_AGGREGATE_BARCODES_OUT stage inputs."""

    antibody_analysis: None | mro_filetypes.path
    is_multi: bool
    multiplexing_method: None | str


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """GET_AGGREGATE_BARCODES_OUT stage outputs."""

    aggregate_barcodes: None | mro_filetypes.csv


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """GET_AGGREGATE_BARCODES_OUT chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """GET_AGGREGATE_BARCODES_OUT chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
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
