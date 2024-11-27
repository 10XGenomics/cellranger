#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Write out cell_types.h5 to use in diff expression."""

__MRO__ = """
call WRITE_CELL_TYPES_H5(
    cell_types = POSTPROCESS_CELL_ANNOTATIONS.cell_types,
) using (
    disabled = CALL_CELL_TYPES.skip_cell_annotation,
    volatile = true,
)
"""
import martian
import polars as pl

import cellranger.analysis.graphclust as cr_graphclust
import cellranger.analysis.io as analysis_io
import tenkit.safe_json as tk_safe_json
from cellranger.cell_typing.cas_postprocessing import COARSE_CELL_TYPES_KEY


def main(args, outs):
    cell_types = pl.read_csv(args.cell_types)
    unique_values = cell_types[COARSE_CELL_TYPES_KEY].unique().sort()
    mapper = {v: i + 1 for i, v in enumerate(unique_values)}
    cell_types = cell_types.with_columns(
        cell_types[COARSE_CELL_TYPES_KEY]
        .map_elements(lambda x: mapper[x])
        .alias(f"{COARSE_CELL_TYPES_KEY}_int")
    )
    labels = cell_types[f"{COARSE_CELL_TYPES_KEY}_int"].to_numpy()
    with analysis_io.open_h5_for_writing(outs.cell_types) as f:
        cr_graphclust.save_graphclust_h5(f, labels)

    # Convert DataFrame columns to lists
    cell_types_list = cell_types[COARSE_CELL_TYPES_KEY].to_list()
    cell_type_ints = cell_types[f"{COARSE_CELL_TYPES_KEY}_int"].to_list()

    # Zip the lists together and create a dictionary
    cell_dict = dict(zip(cell_types_list, cell_type_ints))

    with open(martian.make_path(outs.cell_types_map), "w") as f:
        tk_safe_json.dump_numpy(
            cell_dict,
            f,
            indent=4,
            sort_keys=True,
            separators=(",", ": "),
        )
