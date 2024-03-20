#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine method used for multiplexed data."""

from cellranger.fast_utils import MultiGraph

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
    config = MultiGraph.from_path(args.multi_graph)

    if config.is_cmo_multiplexed():
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo = False
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = False
    elif config.is_rtl_multiplexed():
        outs.multiplexing_is_not_rtl = False
        outs.multiplexing_is_not_cmo = True
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = True
    elif config.is_oh_multiplexed():
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo = True
        outs.multiplexing_is_not_oh = False
        outs.output_per_sample_raw_matrix = True
    else:
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo = True
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = False
