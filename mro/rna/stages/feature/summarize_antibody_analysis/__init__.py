#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved
#
"""Summarize antibody analysis."""
import os
import shutil

__MRO__ = """
stage SUMMARIZE_ANTIBODY_ANALYSIS(
    in  csv  aggregate_barcodes,
    in  bool is_antibody,
    out path antibody_analysis,
    src py   "stages/feature/summarize_antibody_analysis",
) using (
    mem_gb = 4,
)
"""


def main(args, outs):
    os.makedirs(outs.antibody_analysis, exist_ok=True)
    if args.aggregate_barcodes is None:
        return

    shutil.copy(
        args.aggregate_barcodes, os.path.join(outs.antibody_analysis, "aggregate_barcodes.csv")
    )
