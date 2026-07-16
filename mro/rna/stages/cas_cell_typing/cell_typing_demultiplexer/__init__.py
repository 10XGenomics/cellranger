#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Demultiplexing cell typing methods."""

from cellranger.cell_typing.common_cell_typing import AZIMUTH_MODELS, BROAD_TENX_MODELS

__MRO__ = """
stage CELL_TYPING_DEMULTIPLEXER(
    in  string model_family,
    out bool   disable_cloud_cell_types,
    out bool   disable_azimuth_cell_types,
    src py     "stages/cas_cell_typing/cell_typing_demultiplexer",
)
"""


def main(args, outs):
    outs.disable_cloud_cell_types = args.model_family != BROAD_TENX_MODELS
    outs.disable_azimuth_cell_types = args.model_family != AZIMUTH_MODELS
