#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

__MRO__ = """
stage CHOOSE_DIMENSION_REDUCTION(
    in  bool batch_alignment,
    out bool disable_run_pca,
    out bool disable_align_batch,
    src py   "stages/analyzer/choose_dimension_reduction",
)
"""

def main(args, outs):
    if args.batch_alignment is None or args.batch_alignment is False:
        outs.disable_run_pca = False
        outs.disable_align_batch = True
    else:
        outs.disable_run_pca = True
        outs.disable_align_batch = False