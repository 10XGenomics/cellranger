#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Append celltypes to an existing cloupe file."""


import subprocess

import martian

import cellranger.cr_io as cr_io
import tenkit.log_subprocess as tk_subproc

__MRO__ = """
stage APPEND_CELL_TYPES_CLOUPE(
    in  cloupe sample_cloupe,
    in  json   cloupe_cas_types,
    in  string cas_track_name,
    out cloupe cell_annotation_sample_cloupe,
    src py     "stages/cas_cell_typing/append_cell_types_cloupe",
)
"""


def main(args, outs):
    # If we have a cloupe file but no cell types reroute the input to outs
    if args.sample_cloupe and not args.cloupe_cas_types:
        cr_io.hardlink_with_fallback(args.sample_cloupe, outs.cell_annotation_sample_cloupe)
        return

    # Stage gets skipped if no cloupe file is provided. We should still run CAS
    if not args.sample_cloupe or not args.cloupe_cas_types:
        martian.clear(outs)
        return

    cas_track_name = "Cell Types" if not args.cas_track_name else args.cas_track_name

    call = [
        "appendtracks",
        args.cloupe_cas_types,
        args.sample_cloupe,
        outs.cell_annotation_sample_cloupe,
        "--trackName",
        cas_track_name,
        "--emptyName=Unknown",
        "--multiName=Multiple Labels",
    ]

    unicode_call = [arg.encode("utf-8") for arg in call]

    martian.log_info("Running appendtracks: {}".format(" ".join(call)))
    try:
        results = tk_subproc.check_output(unicode_call, stderr=subprocess.STDOUT)
        martian.log_info(f"appendtracks output: {results}")
    except subprocess.CalledProcessError as err:
        outs.cell_annotation_sample_cloupe = None
        martian.throw(f"Could not generate .cloupe file: \n{err.output}")
