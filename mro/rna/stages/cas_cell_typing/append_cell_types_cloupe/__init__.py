#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Append celltypes to an existing cloupe file."""


import subprocess
from dataclasses import dataclass

import martian

import cellranger.cr_io as cr_io
import tenkit.log_subprocess as tk_subproc

__MRO__ = """
struct CellAnnotationLoupeInputs(
    string model_name,
    string cell_annotation_cloupe_name,
    json   cloupe_cas_types,
)


stage APPEND_CELL_TYPES_CLOUPE(
    in  cloupe                      sample_cloupe,
    in  json                        cloupe_cas_types,
    in  string                      cas_track_name,
    in  CellAnnotationLoupeInputs[] cell_annotation_loupe_inputs,
    out cloupe                      cell_annotation_sample_cloupe,
    src py                          "stages/cas_cell_typing/append_cell_types_cloupe",
)
"""


@dataclass
class CellAnnotationLoupeInputs:
    """Struct corresponding to martian struct."""

    model_name: str | None = None
    cell_annotation_cloupe_name: str | None = None
    cloupe_cas_types: str | None = None

    def has_values(self) -> bool:
        return self.cloupe_cas_types is not None


def add_track_to_cloupe(
    sample_cloupe, cloupe_cas_types, cell_annotation_sample_cloupe_path, cas_track_name
):
    """Add track to cloupe files."""
    call = [
        "appendtracks",
        cloupe_cas_types,
        sample_cloupe,
        cell_annotation_sample_cloupe_path,
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
        martian.throw(f"Could not generate .cloupe file: \n{err.output}")


def main(args, outs):
    cell_annotation_loupe_inputs = None
    if args.cell_annotation_loupe_inputs:
        cell_annotation_loupe_inputs = [
            CellAnnotationLoupeInputs(**x) for x in args.cell_annotation_loupe_inputs
        ]

    has_celltypes = bool(args.cloupe_cas_types) or (
        bool(cell_annotation_loupe_inputs)
        and any(x.has_values() for x in cell_annotation_loupe_inputs)
    )
    print(f"{has_celltypes=}")
    # If we have a cloupe file but no cell types reroute the input to outs
    if args.sample_cloupe and not has_celltypes:
        print("hardlinking and exiting")
        cr_io.hardlink_with_fallback(args.sample_cloupe, outs.cell_annotation_sample_cloupe)
        return
    # Stage gets skipped if no cloupe file is provided. We should still run CAS
    if not args.sample_cloupe or not has_celltypes:
        martian.clear(outs)
        return

    if args.cloupe_cas_types:
        cas_track_name = "Cell Types" if not args.cas_track_name else args.cas_track_name
        add_track_to_cloupe(
            args.sample_cloupe,
            args.cloupe_cas_types,
            outs.cell_annotation_sample_cloupe,
            cas_track_name,
        )
    elif cell_annotation_loupe_inputs:
        temp_cloupe_path = martian.make_path("temp.cloupe").decode()
        cr_io.hardlink_with_fallback(args.sample_cloupe, temp_cloupe_path)
        for index, cell_annotation_loupe_input in enumerate(cell_annotation_loupe_inputs):
            if cell_annotation_loupe_input.cloupe_cas_types:
                model_name = (
                    cell_annotation_loupe_input.model_name
                    if cell_annotation_loupe_input.model_name
                    else f"Unnamed mode {index}"
                )
                if args.cas_track_name:
                    track_name_to_use = args.cas_track_name
                elif cell_annotation_loupe_input.cell_annotation_cloupe_name:
                    track_name_to_use = cell_annotation_loupe_input.cell_annotation_cloupe_name
                else:
                    track_name_to_use = "Cell Types"
                cas_track_name = f"{model_name} {track_name_to_use}"
                add_track_to_cloupe(
                    temp_cloupe_path,
                    cell_annotation_loupe_input.cloupe_cas_types,
                    outs.cell_annotation_sample_cloupe,
                    cas_track_name,
                )
                cr_io.hardlink_with_fallback(outs.cell_annotation_sample_cloupe, temp_cloupe_path)
