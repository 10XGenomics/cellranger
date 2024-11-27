#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Placeholder."""

__MRO__ = """
stage COMPUTE_CORRECTION_FACTOR(
    in  V1PatternFixArgs v1_pattern_fix,
    in  json             barcodes_under_tissue,
    out float            correction_factor,
    out json             affected_barcodes,
    out bool             disable_downsampling,
    src py               "stages/spatial/compute_correction_factor",
)
"""


# pylint: disable=unused-argument
def main(args, outs):
    outs.correction_factor = None
    outs.affected_barcodes = None
    outs.disable_downsampling = True
