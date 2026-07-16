#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Remove meta cells from cell types."""

import csv
import json
import os

import martian

import cellranger.cr_io as cr_io
from cellranger.cell_typing.common_cell_typing import BARCODE_KEY

__MRO__ = """
stage REMOVE_META_CELLS_FROM_CELL_TYPES(
    in  csv  cell_types_in,
    in  json metacell_json,
    out csv  cell_types,
    src py   "stages/cas_cell_typing/remove_meta_cells_from_cell_types",
) split (
) using (
    volatile = strict,
)
"""

NEW_BARCODE_KEY = "new_barcode"
OLD_BARCODES_KEY = "old_barcodes"


def split(args):
    if not args.metacell_json:
        return {"chunks": [], "join": {"__mem_gb": 1}}
    else:
        # Mem to read the JSON ~ 3 x on disk JSON size
        mem_gb = 4 + 3 * os.path.getsize(args.metacell_json) / 1024 / 1024 / 1024
        return {"chunks": [], "join": {"__mem_gb": mem_gb}}


def join(args, outs, _chunk_defs, _chunk_outs):
    if not args.cell_types_in:
        martian.clear(outs=outs)
        return
    elif not args.metacell_json:
        cr_io.hardlink_with_fallback(args.cell_types_in, outs.cell_types)
        return

    with open(args.metacell_json) as f:
        metacell_list = json.load(f)
    metacell_dict = {x[NEW_BARCODE_KEY]: x[OLD_BARCODES_KEY] for x in metacell_list}

    with open(outs.cell_types, "w") as f:
        with open(args.cell_types_in) as g:
            reader = csv.DictReader(g)
            column_names = reader.fieldnames
            writer = csv.DictWriter(f, fieldnames=column_names)
            writer.writeheader()
            for row in reader:
                if row[BARCODE_KEY] not in metacell_dict:
                    writer.writerow(row)
                else:
                    new_row = row.copy()
                    for bc in metacell_dict[row[BARCODE_KEY]]:
                        new_row[BARCODE_KEY] = bc
                        writer.writerow(new_row)
