#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Aggregate cell type CSV files from multiple samples into a single CSV file."""

import functools
import operator

import martian

__MRO__ = """
struct CellTypesAllMethods(
    map<CellTypes> cell_types,
)

struct CellTypes(
    csv     cell_types,
    json.gz cell_annotation_results,
    csv     cell_annotation_differential_expression,
)

stage AGGREGATE_CELL_TYPE_CSV(
    in  map<CellTypesAllMethods> cell_types_all_methods,
    out map<csv>                 all_cell_types,
    src py                       "stages/cas_cell_typing/aggregate_cell_type_csv",
) using (
    volatile = strict,
)
"""


def aggregate_csv_files(
    files_to_process: list[tuple[str, str | bytes]], out_file_path: str | bytes
):
    """Aggregate files."""
    aggregated_data = []
    for idx, (sample_name, csv_path) in enumerate(files_to_process):
        if not csv_path:
            continue

        with open(csv_path) as f:
            # Skip the header for all but the first file
            lines = f.readlines()
            if idx == 0:
                aggregated_data.append(
                    lines[0].strip() + ",sample_name\n"
                )  # Add header with sample_name column
            for line in lines[1:]:
                aggregated_data.append(f"{line.strip()},{sample_name}\n")

    # Write the aggregated data to the output file
    with open(out_file_path, "w") as out_file:
        out_file.writelines(aggregated_data)


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
                .get("cell_types")
            )
            if cell_type_found:
                methods_used_to_files_to_aggr_map[method].append((sample, cell_type_found))

    num_samples_with_annotations = max(len(x) for x in methods_used_to_files_to_aggr_map.values())
    if num_samples_with_annotations < 2:
        martian.clear(outs)
        return

    outs.all_cell_types = {}
    for method, files_to_process in methods_used_to_files_to_aggr_map.items():
        if files_to_process:
            out_file_path = martian.make_path(f"{method}.csv").decode()
            aggregate_csv_files(files_to_process, out_file_path)
            outs.all_cell_types[method] = out_file_path
