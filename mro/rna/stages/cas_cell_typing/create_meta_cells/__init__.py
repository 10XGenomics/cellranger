#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Create meta cells."""

import functools
import json
import operator

import martian
import numpy as np
from scipy.sparse import coo_matrix

from cellranger.cr_io import hardlink_with_fallback
from cellranger.matrix import CountMatrix

__MRO__ = """
stage CREATE_META_CELLS(
    in  h5   filtered_matrix,
    in  json metacell_json,
    out h5   filtered_metacell_matrix,
    src py   "stages/cas_cell_typing/create_meta_cells",
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
        mem_gb = 4 + 2 * CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrix)
        return {"chunks": [], "join": {"__mem_gb": mem_gb}}


def join(args, outs, _chunk_defs, _chunk_outs):
    if not args.filtered_matrix:
        martian.clear(outs=outs)
        return

    if not args.metacell_json:
        print("No meta cells to create; copying filtered matrix to output.")
        outs_path = martian.make_path("filtered_metacell_matrix.h5").decode()
        hardlink_with_fallback(args.filtered_matrix, outs_path)
        outs.filtered_metacell_matrix = outs_path
        return

    cmatrix = CountMatrix.load_h5_file(args.filtered_matrix)

    with open(args.metacell_json) as f:
        metacell_list = json.load(f)

    bcs_in_meta_cells = set(
        functools.reduce(operator.iadd, (x[OLD_BARCODES_KEY] for x in metacell_list), [])
    )
    bcs_not_in_meta_cells = [x for x in cmatrix.bcs if x.decode() not in bcs_in_meta_cells]
    bc_indices_to_keep = cmatrix.bcs_to_ints(bcs_not_in_meta_cells)

    row_indices = bc_indices_to_keep
    new_bcs = list(bcs_not_in_meta_cells)
    col_indices = list(range(len(row_indices)))
    cur_col_index = len(row_indices)
    for new_meta_cell_ind, meta_cell in enumerate(metacell_list):

        new_bcs.append(meta_cell[NEW_BARCODE_KEY].encode())

        row_indices_to_keep = cmatrix.bcs_to_ints([x.encode() for x in meta_cell[OLD_BARCODES_KEY]])
        for row_index in row_indices_to_keep:
            row_indices.append(row_index)
            col_indices.append(cur_col_index + new_meta_cell_ind)

    sparse_coo = coo_matrix(
        (np.ones(len(col_indices)), (row_indices, col_indices)),
        shape=(len(cmatrix.bcs), len(new_bcs)),
        dtype=np.int32,
    )

    new_cmatrix = CountMatrix(
        feature_ref=cmatrix.feature_ref, bcs=new_bcs, matrix=cmatrix.m @ sparse_coo
    )

    new_cmatrix.save_h5_file(outs.filtered_metacell_matrix)
