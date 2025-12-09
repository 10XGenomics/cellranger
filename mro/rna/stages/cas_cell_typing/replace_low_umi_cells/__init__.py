#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Replace low UMI cells from cell_types.csv."""

import csv

import martian
import numpy as np

import cellranger.matrix as cr_matrix
from cellranger.cell_typing.broad_tenx.cas_postprocessing import (
    LOW_UMI_BARCODE_KEY,
    MIN_CELL_TYPE_UMI,
)

__MRO__ = """
stage REPLACE_LOW_UMI_CELLS(
    in  h5  filtered_matrix,
    in  csv cell_types,
    out csv cell_types,
    src py  "stages/cas_cell_typing/replace_low_umi_cells",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    mem_gib = 2 + cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrix, scale=1.1)
    return {
        "chunks": [],
        "join": {"__mem_gb": mem_gib, "__vmem_gb": mem_gib * 2},
    }


def join(args, outs, _chunk_defs, _chunk_outs):
    if any(arg is None for arg in args):
        martian.clear(outs)
        return

    # Get counts per barcode (GEX)
    mat = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrix)
    mat = mat.select_features_by_type("Gene Expression")
    counts = mat.get_counts_per_bc()

    bcs = cr_matrix.CountMatrix.load_bcs_from_h5_file_handle(args.filtered_matrix)
    # Decode barcodes and create a structured array for counts per barcode
    bcs = [b.decode("utf-8") for b in bcs]
    counts_per_bcs = np.array(
        list(zip(bcs, counts)), dtype=[("barcode", "U50"), ("umi_counts", "i4")]
    )

    with open(args.cell_types, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        cell_types = list(reader)

    # Create a dictionary for quick lookup of UMI counts by barcode
    umi_dict = {row["barcode"]: row["umi_counts"] for row in counts_per_bcs}

    # Update the coarse_cell_type based on UMI counts and add new columns
    for cell in cell_types:
        umi_count = umi_dict.get(cell["barcode"])
        cell["coarse_cell_type_unfiltered"] = cell["coarse_cell_type"]
        cell["umi_count"] = umi_count  # Add UMI count column
        if umi_count < MIN_CELL_TYPE_UMI:
            cell["coarse_cell_type"] = f"{LOW_UMI_BARCODE_KEY}"

    # Write the updated cell types back to the output CSV
    with open(outs.cell_types, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cell_types[0].keys())
        writer.writeheader()
        writer.writerows(cell_types)
