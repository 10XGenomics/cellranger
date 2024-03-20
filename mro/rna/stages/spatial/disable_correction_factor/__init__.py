# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
"""Checks if we should run the Chevron fix."""

__MRO__ = """
stage DISABLE_CORRECTION_FACTOR(
    in  h5   v1_filtered_fbm,
    out bool disable_correction_factor,
    src py   "stages/spatial/disable_correction_factor",
)
"""


def main(args, outs):
    outs.disable_correction_factor = not bool(args.v1_filtered_fbm)
