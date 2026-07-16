# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Create a config that determines which subpipelines will be run."""

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage _MAKE_VDJ_CONFIG(
    in  VdjInputsCS vdj_t_input,
    in  VdjInputsCS vdj_t_gd_input,
    in  VdjInputsCS vdj_b_input,
    in  bool        disable_vdj,
    in  bool        disable_count,
    in  path        vdj_reference_path,
    out bool        count_is_enabled,
    out bool        disable_vdj_b,
    out bool        disable_vdj_t,
    out bool        disable_vdj_t_gd,
    out bool        has_no_vdj_ref,
    src py          "stages/common/make_vdj_config",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """_MAKE_VDJ_CONFIG stage inputs."""

    vdj_t_input: None | mro_structs.VdjInputsCS
    vdj_t_gd_input: None | mro_structs.VdjInputsCS
    vdj_b_input: None | mro_structs.VdjInputsCS
    disable_vdj: bool
    disable_count: bool
    vdj_reference_path: None | mro_filetypes.path


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """_MAKE_VDJ_CONFIG stage outputs."""

    count_is_enabled: bool
    disable_vdj_b: bool
    disable_vdj_t: bool
    disable_vdj_t_gd: bool
    has_no_vdj_ref: bool


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """_MAKE_VDJ_CONFIG chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """_MAKE_VDJ_CONFIG chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    outs.count_is_enabled = not args.disable_count
    outs.disable_vdj_t = args.vdj_t_input is None or args.disable_vdj
    outs.disable_vdj_t_gd = args.vdj_t_gd_input is None or args.disable_vdj
    outs.disable_vdj_b = args.vdj_b_input is None or args.disable_vdj
    outs.has_no_vdj_ref = args.vdj_reference_path is None
