#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Aggregate cell type CSV files from multiple samples into a single CSV file."""

import martian

__MRO__ = """
stage AGGREGATE_CELL_TYPE_CSV(
    in  map<CellTypes> cell_types,
    out csv            all_cell_types,
    src py             "stages/cas_cell_typing/aggregate_cell_type_csv",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if not args.cell_types:
        martian.clear(outs)
        return

    # Check if there are more than one sample with valid cell_types
    valid_cell_types_count = sum(
        1
        for _, sample_data in args.cell_types.items()
        if sample_data
        and isinstance(sample_data, dict)
        and "cell_types" in sample_data
        and sample_data["cell_types"] is not None
    )
    if valid_cell_types_count < 2:
        martian.clear(outs)
        return

    files_to_process = [
        (sample_name, sample_data.get("cell_types"))
        for sample_name, sample_data in args.cell_types.items()
        if sample_data.get("cell_types") is not None
    ]

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
    with open(outs.all_cell_types, "w") as out_file:
        out_file.writelines(aggregated_data)
