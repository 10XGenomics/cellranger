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
from cellranger.cell_typing.broad_tenx.cas_metrics import analysis_has_clustering, cell_type_umap
from cellranger.cell_typing.common_cell_typing import (
    DISPLAY_MAP_KEY,
    FINE_CELL_TYPES_KEY,
)


def main(args, outs):
    if not args.cell_types:
        martian.clear(outs)
        return

    analysis = args.analysis if args.analysis else None

    if (analysis and analysis_has_clustering(analysis)) or args.cloupe_projection:
        umap_plot = cell_type_umap(
            cell_types=args.cell_types,
            analysis=analysis,
            cloupe_projection=args.cloupe_projection,
            configure_legend=True,
            coarse_cell_type=DISPLAY_MAP_KEY,
            fine_cell_type=FINE_CELL_TYPES_KEY,
        )
        feature_utils.write_json_from_dict(umap_plot, outs.cell_types_umap_plot)
    else:
        martian.clear(outs)
