#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Write out cell_types.h5 to use in diff expression."""

__MRO__ = """
stage WRITE_CELL_TYPES_H5(
    in  csv    cell_types,
    in  string cell_types_key,
    out h5     cell_types,
    out json   cell_types_map,
    src py     "stages/cas_cell_typing/write_cell_types_h5",
)
"""

import martian
import polars as pl

import cellranger.analysis.graphclust as cr_graphclust
import cellranger.analysis.io as analysis_io
import tenkit.safe_json as tk_safe_json
from cellranger.cell_typing.common_cell_typing import COARSE_CELL_TYPES_KEY


def main(args, outs):
    if args.cell_types_key is None:
        cell_types_key = COARSE_CELL_TYPES_KEY
    else:
        cell_types_key = args.cell_types_key
    cell_types = pl.read_csv(args.cell_types, null_values=["NA"])
    unique_values = cell_types[cell_types_key].unique().sort()
    mapper = {v: i + 1 for i, v in enumerate(unique_values)}
    cell_types = cell_types.with_columns(
        cell_types[cell_types_key].map_elements(lambda x: mapper[x]).alias(f"{cell_types_key}_int")
    )
    labels = cell_types[f"{cell_types_key}_int"].to_numpy()
    with analysis_io.open_h5_for_writing(outs.cell_types) as f:
        cr_graphclust.save_graphclust_h5(f, labels)

    # Convert DataFrame columns to lists
    cell_types_list = cell_types[cell_types_key].to_list()
    cell_type_ints = cell_types[f"{cell_types_key}_int"].to_list()

    # Zip the lists together and create a dictionary
    cell_dict = dict(zip(cell_types_list, cell_type_ints))
    cell_dict = {k: v for k, v in cell_dict.items() if k is not None}
    with open(martian.make_path(outs.cell_types_map), "w") as f:
        tk_safe_json.dump_numpy(
            cell_dict,
            f,
            indent=4,
            sort_keys=True,
            separators=(",", ": "),
        )
