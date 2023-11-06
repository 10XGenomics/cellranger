#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to build CS outs tailored to the multiplexing strategy."""
from cellranger.multi.build_per_sample_outs import build_sample_outs

__MRO__ = """
stage BUILD_SAMPLE_OUTS(
    in  SampleSlfeOuts  sample_slfe_outs,
    in  path            rna_analysis,
    in  path            crispr_analysis,
    in  cloupe          cloupe,
    in  html            web_summary,
    in  csv             metrics_summary_csv,
    in  VdjOutputsCS    vdj_b_outs,
    in  VdjOutputsCS    vdj_t_outs,
    in  VdjOutputsCS    vdj_t_gd_outs,
    in  bool            output_per_sample_raw_matrix,
    in  BEAM_ANALYZER   beam_analyzer,
    out SampleOutputsCS sample_outs,
    src py              "stages/multi/build_sample_outs",
)
"""


def main(args, outs):
    """Create the per sample outs.

    Args:
        args:
        outs:

    Returns:
        None
    """
    build_sample_outs(args, outs, False)
