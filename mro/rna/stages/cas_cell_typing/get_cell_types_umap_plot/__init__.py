#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Make umap colored by cell type."""

__MRO__ = """
stage GET_CELL_TYPES_UMAP_PLOT(
    in  path   analysis,
    in  csv    cell_types,
    in  string pipestance_type,
    out json   cell_types_umap_plot,
    src py     "stages/cas_cell_typing/get_cell_types_umap_plot",
) using (
    mem_gb   = 2,
    volatile = strict,
)
"""

import martian

import cellranger.feature.utils as feature_utils
from cellranger.cell_typing.cas_metrics import cell_type_umap


def main(args, outs):
    if args.analysis and args.cell_types:
        ## UMAP of cell types
        umap_plot = cell_type_umap(
            cell_types=args.cell_types, analysis=args.analysis, cloupe_projection=None
        )
        feature_utils.write_json_from_dict(umap_plot, outs.cell_types_umap_plot)
        return
    if args.cloupe_projection and args.cell_types:
        umap_plot = cell_type_umap(
            cell_types=args.cell_types,
            analysis=None,
            cloupe_projection=args.cloupe_projection,
        )
        feature_utils.write_json_from_dict(umap_plot, outs.cell_types_umap_plot)
        return
    else:
        martian.clear(outs)
