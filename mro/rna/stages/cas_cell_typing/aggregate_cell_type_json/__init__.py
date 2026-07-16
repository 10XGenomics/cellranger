#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Analyze cell types generated from CALL_CELL_TYPES."""

import functools
import gzip
import operator
from typing import TextIO

import martian

import tenkit.safe_json as tk_safe_json

__MRO__ = """
struct CellTypesAllMethods(
    map<CellTypes> cell_types,
)

struct CellTypes(
    csv     cell_types,
    json.gz cell_annotation_results,
    csv     cell_annotation_differential_expression,
)

stage AGGREGATE_CELL_TYPE_JSON(
    in  map<CellTypesAllMethods> cell_types_all_methods,
    out map<json.gz>             all_cell_annotation_results,
    src py                       "stages/cas_cell_typing/aggregate_cell_type_json",
) using (
    volatile = strict,
)
"""


def read_past_opening_bracket(f: TextIO) -> None:
    """Reads from the file stream `f` until it finds the opening bracket '[' of a JSON array. It skips any leading whitespace.

    Args:
        f: input file.

    """
    while True:
        char = f.read(1)
        if not char:
            raise ValueError("Unexpected end of file while searching for '['")
        # Ignore whitespace
        if char.isspace():
            continue
        # Validate first non-whitespace character is a bracket.
        if char == "[":
            return
        else:
            raise ValueError(
                f"File does not start with a JSON array. Found '{char}' instead of '['."
            )


def stream_copy_until_closing_bracket(in_f: TextIO, out_f: TextIO, chunk_size: int = 10000) -> None:
    """Reads from `in_f` and writes to `out_f` in chunks, stopping just before the final closing bracket ']'.

    This avoids loading the entire file into memory.

    Args:
        in_f: The input stream (opened in text mode).
        out_f: The output stream (opened in text mode).
        chunk_size (int): The size of data to read in each step.
    """
    # Initial chunk
    buffer = in_f.read(chunk_size)

    while buffer:
        # Read the next chunk from the input file
        next_chunk = in_f.read(chunk_size)

        if not next_chunk:
            last_bracket_pos = buffer.rfind("]")
            if last_bracket_pos != -1:
                # Write everything *before* the last bracket
                out_f.write(buffer[:last_bracket_pos])
            else:
                raise ValueError("File does not end with a closing bracket ']'.")
            break  # End of the stream
        else:
            # If not the last chunk, write the entire current buffer to the output file.
            out_f.write(buffer)
            buffer = next_chunk


def aggregate_json_gz_files(files_to_process, out_file_path):
    """Aggregate JSON GZ files."""
    # Open the output file in text mode ('wt')
    with gzip.open(out_file_path, "wt", encoding="utf-8") as out_file:
        out_file.write("[")

        is_first_element = True

        for sample_name, results_path in files_to_process:
            if not is_first_element:
                out_file.write(",")
            else:
                is_first_element = False

            # Dump the sample name
            tk_safe_json.dump_numpy({"sample_name": sample_name}, out_file)

            # Open the source file and stream its contents
            with gzip.open(results_path, "rt", encoding="utf-8") as in_f:
                # Make sure there is content to copy by checking the first element in the array
                first_char = in_f.read(1)
                if first_char and not first_char.isspace() and first_char != "]":
                    # Add a comma to separate the sample_name from the array content.
                    out_file.write(",")
                    # Go back to the start and process the file
                    in_f.seek(0)
                    read_past_opening_bracket(in_f)
                    stream_copy_until_closing_bracket(in_f, out_file)

        out_file.write("]")


def main(args, outs):
    if not args.cell_types_all_methods:
        martian.clear(outs)
        return

    sample_keys = sorted(args.cell_types_all_methods)
    methods_used = sorted(
        set(
            functools.reduce(
                operator.iadd,
                (x.get("cell_types", []) for x in args.cell_types_all_methods.values() if x),
                [],
            )
        )
    )

    if not methods_used:
        martian.clear(outs)
        return

    methods_used_to_files_to_aggr_map = {method: [] for method in methods_used}

    for method in methods_used:
        for sample in sample_keys:
            cell_type_found = (
                args.cell_types_all_methods.get(sample, {})
                .get("cell_types", {})
                .get(method, {})
                .get("cell_annotation_results")
            )
            if cell_type_found:
                methods_used_to_files_to_aggr_map[method].append((sample, cell_type_found))

    num_samples_with_annotations = max(len(x) for x in methods_used_to_files_to_aggr_map.values())
    if num_samples_with_annotations < 2:
        martian.clear(outs)
        return

    outs.all_cell_annotation_results = {}
    for method, files_to_process in methods_used_to_files_to_aggr_map.items():
        if files_to_process:
            out_file_path = martian.make_path(f"{method}.json.gz").decode()
            aggregate_json_gz_files(files_to_process, out_file_path)
            outs.all_cell_annotation_results[method] = out_file_path
