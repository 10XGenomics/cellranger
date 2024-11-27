#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""This stage extracts a projection from a cloupe file."""

__MRO__ = """
stage ETXRACT_LOUPE_PROJECTION(
    in  cloupe sample_cloupe,
    in  string projection_name,
    out csv    projection,
    src py     "stages/cas_cell_typing/extract_loupe_projection",
)
"""

import subprocess

import martian

import tenkit.log_subprocess as tk_subproc


def main(args, outs):
    if not args.sample_cloupe:
        martian.clear(outs)
        return

    outs.projection = martian.make_path(f"{args.projection_name}.csv").decode()

    call = [
        "exportprojectioncsv",
        args.sample_cloupe,
        outs.projection,
        args.projection_name,
    ]

    unicode_call = [arg.encode("utf-8") for arg in call]

    martian.log_info("Running exportprojectioncsv: {}".format(" ".join(call)))
    try:
        results = tk_subproc.check_output(unicode_call, stderr=subprocess.STDOUT)
        martian.log_info(f"exportprojectioncsv output: {results}")
    except subprocess.CalledProcessError as err:
        martian.clear(outs)
        martian.alarm(
            f"Could not run exportprojectioncsv to get a projection from cloupe file. Error: \n{err.output}"
        )
