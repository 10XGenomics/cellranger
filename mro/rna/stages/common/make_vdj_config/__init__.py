# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Create a config that determines which subpipelines will be run."""

__MRO__ = """
stage _MAKE_VDJ_CONFIG(
    in  VdjInputsCS vdj_t_input,
    in  VdjInputsCS vdj_t_gd_input,
    in  VdjInputsCS vdj_b_input,
    in  bool        disable_vdj,
    in  path        vdj_reference_path,
    out bool        disable_vdj_b,
    out bool        disable_vdj_t,
    out bool        disable_vdj_t_gd,
    out bool        has_no_vdj_ref,
    src py          "stages/common/make_vdj_config",
)
"""


def main(args, outs):
    outs.disable_vdj_t = args.vdj_t_input is None or args.disable_vdj
    outs.disable_vdj_t_gd = args.vdj_t_gd_input is None or args.disable_vdj
    outs.disable_vdj_b = args.vdj_b_input is None or args.disable_vdj
    outs.has_no_vdj_ref = args.vdj_reference_path is None
