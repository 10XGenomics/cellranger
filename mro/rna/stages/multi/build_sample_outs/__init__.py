#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to build CS outs tailored to the multiplexing strategy."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import martian

from cellranger.multi.build_per_sample_outs import build_sample_outs

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage BUILD_SAMPLE_OUTS(
    in  SampleSlfeOuts      sample_slfe_outs,
    in  path                rna_analysis,
    in  path                crispr_analysis,
    in  cloupe              cloupe,
    in  html                web_summary,
    in  csv                 metrics_summary_csv,
    in  VdjOutputsCS        vdj_b_outs,
    in  VdjOutputsCS        vdj_t_outs,
    in  VdjOutputsCS        vdj_t_gd_outs,
    in  bool                output_per_sample_raw_matrix,
    in  BeamAnalyzerOutputs beam_analyzer,
    in  CellTypes           cell_types,
    out SampleOutputsCS     sample_outs,
    src py                  "stages/multi/build_sample_outs",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """BUILD_SAMPLE_OUTS stage inputs."""

    sample_slfe_outs: None | mro_structs.SampleSlfeOuts
    rna_analysis: None | mro_filetypes.path
    crispr_analysis: None | mro_filetypes.path
    cloupe: None | mro_filetypes.cloupe
    web_summary: None | mro_filetypes.html
    metrics_summary_csv: None | mro_filetypes.csv
    vdj_b_outs: None | dict[str, Any]
    vdj_t_outs: None | dict[str, Any]
    vdj_t_gd_outs: None | dict[str, Any]
    output_per_sample_raw_matrix: bool
    beam_analyzer: None | dict[str, Any]
    cell_types: None | dict[str, Any]


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """BUILD_SAMPLE_OUTS stage outputs."""

    sample_outs: dict[str, Any]


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """BUILD_SAMPLE_OUTS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """BUILD_SAMPLE_OUTS chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    """Create the per sample outs.

    Args:
        args:
        outs:

    Returns:
        None
    """
    build_sample_outs(args, outs)
