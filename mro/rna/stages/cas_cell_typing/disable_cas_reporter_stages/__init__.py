#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Disable reporter stages based on cell annotation outs."""

import csv

from cellranger.cell_typing.cas_postprocessing import (
    BARCODE_KEY,
    COARSE_CELL_TYPES_KEY,
)

__MRO__ = """
stage DISABLE_CAS_REPORTER_STAGES(
    in  csv  cell_types,
    out bool disable_differential_expression,
    src py   "stages/cas_cell_typing/disable_cas_reporter_stages",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if not args.cell_types:
        outs.disable_differential_expression = True
        return

    with open(args.cell_types) as read_csv_file:
        reader = csv.DictReader(
            read_csv_file,
        )
        hdrs = reader.fieldnames
        if not hdrs or BARCODE_KEY not in hdrs and COARSE_CELL_TYPES_KEY not in hdrs:
            raise ValueError("Invalid cell type CSV file.")

        set_of_cell_types = {row[COARSE_CELL_TYPES_KEY] for row in reader}
        outs.disable_differential_expression = len(set_of_cell_types) < 2
