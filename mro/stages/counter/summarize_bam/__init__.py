#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# This stage only exists to force earlier VDR of prior stages.

import cellranger.utils as cr_utils

__MRO__ = """
stage SUMMARIZE_BAM(
    in json  summary,
    in csv   barcodes_detected,
    out json summary,
    out csv  barcodes_detected,
    sry py   "stages/counter/summarize_bam",
)
"""

def main(args, outs):
    if args.summary is not None:
        cr_utils.copy(args.summary, outs.summary)
    if args.barcodes_detected is not None:
        cr_utils.copy(args.barcodes_detected, outs.barcodes_detected)
