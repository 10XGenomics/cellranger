#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved
#
"""Summarize antibody analysis."""

from __future__ import annotations

import os
import shutil
from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage SUMMARIZE_ANTIBODY_ANALYSIS(
    in  csv  aggregate_barcodes,
    in  bool is_antibody,
    out path antibody_analysis,
    src py   "stages/feature/summarize_antibody_analysis",
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_ANTIBODY_ANALYSIS stage inputs."""

    aggregate_barcodes: None | mro_filetypes.csv
    is_antibody: bool


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_ANTIBODY_ANALYSIS stage outputs."""

    antibody_analysis: mro_filetypes.path


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_ANTIBODY_ANALYSIS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_ANTIBODY_ANALYSIS chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    os.makedirs(outs.antibody_analysis, exist_ok=True)
    if args.aggregate_barcodes is None:
        return

    shutil.copy(
        args.aggregate_barcodes, os.path.join(outs.antibody_analysis, "aggregate_barcodes.csv")
    )
