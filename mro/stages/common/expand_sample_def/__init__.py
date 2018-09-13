#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import martian
import cellranger.preflight as cr_preflight

__MRO__ = """
stage EXPAND_SAMPLE_DEF(
    in  map[]  raw_sample_def,
    out map[]  sample_def,
    src py     "stages/common/expand_sample_def",
)
"""

def main(args, outs):
    ''' Convert sample_def = { "libraries_csv": "/path/to/libraries.csv" } into a
        standard sample_def map used by the rest of the pipeline. Only used by the
        CS pipeline to handle the --libraries cmd-line argument.'''

    if len(args.raw_sample_def) == 1:
        if args.raw_sample_def[0].keys() == ["libraries"]:
            # We've got a 'libraries mode' argument coming in -- load & check the CSV, and expand it into a normal sample def
            try:
                outs.sample_def = cr_preflight.expand_libraries_csv(args.raw_sample_def[0]["libraries"])
                return
            except cr_preflight.PreflightException as e:
                martian.exit(e.msg)

    # Default case -- just copy over the sample_def
    outs.sample_def = args.raw_sample_def
