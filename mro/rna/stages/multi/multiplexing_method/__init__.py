#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine method used for multiplexed data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

import cellranger.rna.library as rna_library
from cellranger.fast_utils import MultiGraph

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage MULTIPLEXING_METHOD(
    in  json   multi_graph,
    out bool   multiplexing_is_not_read_level,
    out bool   multiplexing_is_not_cmo_or_hashtag,
    out bool   output_per_sample_raw_matrix,
    out string multiplexing_method,
    src py     "../rna/stages/multi/multiplexing_method",
) using (
    mem_gb   = 1,
    threads  = 1,
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """MULTIPLEXING_METHOD stage inputs."""

    multi_graph: None | mro_filetypes.json


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """MULTIPLEXING_METHOD stage outputs."""

    multiplexing_is_not_read_level: bool
    multiplexing_is_not_cmo_or_hashtag: bool
    output_per_sample_raw_matrix: bool
    multiplexing_method: None | str


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """MULTIPLEXING_METHOD chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """MULTIPLEXING_METHOD chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    if args.multi_graph is None:
        outs.multiplexing_is_not_read_level = True
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.multiplexing_method = None
        return

    # Read in the multi graph
    config = MultiGraph.from_path(args.multi_graph)

    if config.is_cmo_multiplexed():
        outs.multiplexing_is_not_read_level = True
        outs.multiplexing_is_not_cmo_or_hashtag = False
        outs.output_per_sample_raw_matrix = False
        outs.multiplexing_method = rna_library.CellLevel.CMO.value
    elif config.is_hashtag_multiplexed():
        outs.multiplexing_is_not_read_level = True
        outs.multiplexing_is_not_cmo_or_hashtag = False
        outs.output_per_sample_raw_matrix = False
        outs.multiplexing_method = rna_library.CellLevel.Hashtag.value
    elif config.is_rtl_multiplexed():
        outs.multiplexing_is_not_read_level = False
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.output_per_sample_raw_matrix = True
        outs.multiplexing_method = rna_library.ReadLevel.RTL.value
    elif config.is_oh_multiplexed():
        outs.multiplexing_is_not_read_level = False
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.output_per_sample_raw_matrix = True
        outs.multiplexing_method = rna_library.ReadLevel.OH.value
    elif config.is_cmx_multiplexed():
        outs.multiplexing_is_not_read_level = False
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.output_per_sample_raw_matrix = True
        outs.multiplexing_method = rna_library.ReadLevel.CMX.value
    else:
        outs.multiplexing_is_not_read_level = True
        outs.multiplexing_is_not_cmo_or_hashtag = True
        outs.output_per_sample_raw_matrix = False
        outs.multiplexing_method = None
