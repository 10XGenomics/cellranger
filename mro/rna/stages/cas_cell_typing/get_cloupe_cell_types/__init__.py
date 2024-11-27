#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Get JSON to pass into loupe file."""

import csv
import json

from cellranger.cell_typing.cas_postprocessing import BARCODE_KEY, COARSE_CELL_TYPES_KEY

__MRO__ = """
stage GET_CLOUPE_CELL_TYPES(
    in  csv  cell_types,
    out json cloupe_cas_types,
    src py   "stages/cas_cell_typing/get_cloupe_cell_types",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if not args.cell_types:
        outs.cloupe_cas_types = None
        return

    dct = {}

    with open(args.cell_types) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            dct.setdefault(row[COARSE_CELL_TYPES_KEY], []).append(row[BARCODE_KEY])

    with open(outs.cloupe_cas_types, "w") as outfile:
        json.dump(dct, outfile)
