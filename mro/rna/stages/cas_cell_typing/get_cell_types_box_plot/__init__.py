#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Generate cell types interactive barchart."""

import json

import cellranger.matrix as cr_matrix
import cellranger.websummary.violin_plots as cr_vp

__MRO__ = """
stage GET_CELL_TYPES_BOX_PLOT(
    in  csv  cell_types,
    in  h5   filtered_matrix,
    out json cell_types_box_plot,
    src py   "stages/cas_cell_typing/get_cell_types_box_plot",
) split (
) using (
    volatile = strict,
)
"""
COARSE_CELL_TYPE_LOWER_BOUND = 10


def split(args):
    mem_gib = 2 + cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrix, scale=1.1)
    return {
        "chunks": [],
        "join": {"__mem_gb": mem_gib},
    }


def join(args, outs, _chunk_defs, _chunk_outs):
    if args.pipestance_type.startswith("SC_MULTI_CORE_SAMPLE"):
        plot_width = 450
    else:
        plot_width = 900
    box_plot = cr_vp.make_cell_types_boxplot(
        h5_path=args.filtered_matrix,
        cell_types=args.cell_types,
        barcode_lower_bound=COARSE_CELL_TYPE_LOWER_BOUND,
        final_plot_width=plot_width,
        return_plot=False,
    )
    with open(outs.cell_types_box_plot, "w") as f:
        json.dump(box_plot, f)
