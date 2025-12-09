#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Analyze cell types generated from CALL_CELL_TYPES."""

import gzip
from typing import TextIO

import martian

import tenkit.safe_json as tk_safe_json

__MRO__ = """
stage AGGREGATE_CELL_TYPE_JSON(
    in  map<CellTypes> cell_types,
    out json           all_cell_annotation_results,
    src py             "stages/cas_cell_typing/aggregate_cell_type_json",
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


def main(args, outs):
    if not args.cell_types:
        martian.clear(outs)
        return

    # Check if there are more than one sample with valid cell_annotation_results
    valid_cell_types_count = sum(
        1
        for _, sample_data in args.cell_types.items()
        if sample_data
        and isinstance(sample_data, dict)
        and "cell_annotation_results" in sample_data
        and sample_data["cell_annotation_results"] is not None
    )
    if valid_cell_types_count < 2:
        martian.clear(outs)
        return

    files_to_process = [
        (sample_name, sample_data.get("cell_annotation_results"))
        for sample_name, sample_data in args.cell_types.items()
        if sample_data.get("cell_annotation_results") is not None
    ]

    if files_to_process:
        # Open the output file in text mode ('wt')
        with gzip.open(outs.all_cell_annotation_results, "wt", encoding="utf-8") as out_file:
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
    else:
        martian.clear(outs)
