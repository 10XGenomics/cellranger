#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine method used for multiplexed data."""

import cellranger.rna.library as rna_library
from cellranger.fast_utils import MultiGraph

__MRO__ = """
stage MULTIPLEXING_METHOD(
    in  json   multi_graph,
    out bool   multiplexing_is_not_rtl,
    out bool   multiplexing_is_not_cmo_or_hashtag,
    out bool   multiplexing_is_not_oh,
    out bool   output_per_sample_raw_matrix,
    out string multiplexing_method,
    src py     "../rna/stages/multi/multiplexing_method",
) using (
    mem_gb   = 1,
    threads  = 1,
    volatile = strict,
)
"""


def main(args, outs):
    if args.multi_graph is None:
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.multiplexing_is_not_oh = True
        outs.multiplexing_method = None
        return

    # Read in the multi graph
    config = MultiGraph.from_path(args.multi_graph)

    if config.is_cmo_multiplexed():
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo_or_hashtag = False
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = False
        outs.multiplexing_method = rna_library.CellLevel.CMO.value
    elif config.is_hashtag_multiplexed():
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo_or_hashtag = False
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = False
        outs.multiplexing_method = rna_library.CellLevel.Hashtag.value
    elif config.is_rtl_multiplexed():
        outs.multiplexing_is_not_rtl = False
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = True
        outs.multiplexing_method = rna_library.ReadLevel.RTL.value
    elif config.is_oh_multiplexed():
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.multiplexing_is_not_oh = False
        outs.output_per_sample_raw_matrix = True
        outs.multiplexing_method = rna_library.ReadLevel.OH.value
    else:
        outs.multiplexing_is_not_rtl = True
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.multiplexing_is_not_oh = True
        outs.output_per_sample_raw_matrix = False
        outs.multiplexing_method = None
