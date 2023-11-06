#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine method used for multiplexed data."""

import json

from cellranger.multi import config as multi_config

__MRO__ = """
stage MULTIPLEXING_METHOD(
    in  json multi_graph,
    out bool multiplexing_is_not_rtl,
    out bool multiplexing_is_not_cmo,
    out bool multiplexing_is_not_oh,
    out bool output_per_sample_raw_matrix,
    src py   "../rna/stages/multi/multiplexing_method",
) using (
    mem_gb   = 1,
    threads  = 1,
    volatile = strict,
)
"""


def main(args, outs):
    if args.multi_graph is None:
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo = True
        outs.multiplexing_is_not_oh = True
        return

    # Read in the multi graph
    with open(args.multi_graph) as in_file:
        config = multi_config.CrMultiGraph.from_json_val(json.load(in_file))
        multiplexing_type = config.get_cell_multiplexing_type()

        if multiplexing_type is None:
            outs.multiplexing_is_not_rtl = True
            outs.multiplexing_is_not_cmo = True
            outs.multiplexing_is_not_oh = True
            outs.output_per_sample_raw_matrix = False
        elif multiplexing_type == multi_config.CellMultiplexingType.CMO:
            outs.multiplexing_is_not_rtl = True
            outs.multiplexing_is_not_cmo = False
            outs.multiplexing_is_not_oh = True
            outs.output_per_sample_raw_matrix = False
        elif multiplexing_type == multi_config.CellMultiplexingType.RTL:
            outs.multiplexing_is_not_rtl = False
            outs.multiplexing_is_not_cmo = True
            outs.multiplexing_is_not_oh = True
            outs.output_per_sample_raw_matrix = True
        elif multiplexing_type == multi_config.CellMultiplexingType.OH:
            outs.multiplexing_is_not_rtl = True
            outs.multiplexing_is_not_cmo = True
            outs.multiplexing_is_not_oh = False
            outs.output_per_sample_raw_matrix = True
        else:
            raise NotImplementedError(
                "Need to pass none or an enum value from multiplexing method."
            )
