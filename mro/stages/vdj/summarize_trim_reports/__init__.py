#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# This stage only exists to force earlier VDR of prior stages.

import cellranger.utils as cr_utils

__MRO__ = """
stage SUMMARIZE_TRIM_REPORTS(
    in  json trim_reads_summary,
    out json summary,
    src py   "stages/vdj/summarize_trim_reports",
)
"""

def main(args, outs):
    cr_utils.copy(args.trim_reads_summary, outs.summary)
