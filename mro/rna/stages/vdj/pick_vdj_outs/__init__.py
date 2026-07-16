# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Pick T or B output."""

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage PICK_VDJ_OUTS(
    in  bool         disable_vdj_t,
    in  bool         disable_vdj_b,
    in  VdjOutputsCS vdj_t_outs,
    in  html         vdj_t_web_summary,
    in  VdjOutputsCS vdj_b_outs,
    in  html         vdj_b_web_summary,
    out VdjOutputsCS vdj_outs,
    out html         web_summary,
    src py           "stages/vdj/pick_vdj_outs",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """PICK_VDJ_OUTS stage inputs."""

    disable_vdj_t: bool
    disable_vdj_b: bool
    vdj_t_outs: None | mro_structs.VdjOutputsCS
    vdj_t_web_summary: None | mro_filetypes.html
    vdj_b_outs: None | mro_structs.VdjOutputsCS
    vdj_b_web_summary: None | mro_filetypes.html


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """PICK_VDJ_OUTS stage outputs."""

    vdj_outs: None | mro_structs.VdjOutputsCS
    web_summary: None | mro_filetypes.html


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """PICK_VDJ_OUTS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """PICK_VDJ_OUTS chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    assert (
        args.disable_vdj_t ^ args.disable_vdj_b
    )  # Need to be exclusive. Could have taken just one as input

    if not args.disable_vdj_t:
        outs.vdj_outs = args.vdj_t_outs
        outs.web_summary = args.vdj_t_web_summary

    if not args.disable_vdj_b:
        outs.vdj_outs = args.vdj_b_outs
        outs.web_summary = args.vdj_b_web_summary
