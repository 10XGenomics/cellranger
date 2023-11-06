# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Pick T or B output."""

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


def main(args, outs):
    assert (
        args.disable_vdj_t ^ args.disable_vdj_b
    )  # Need to be exclusive. Could have taken just one as input

    if not args.disable_vdj_t:
        outs.vdj_outs = args.vdj_t_outs
        outs.web_summary = args.vdj_t_web_summary

    if not args.disable_vdj_b:
        outs.vdj_outs = args.vdj_b_outs
        outs.web_summary = args.vdj_b_web_summary
