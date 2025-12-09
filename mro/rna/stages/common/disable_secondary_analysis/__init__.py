# Copyright (c) 2023 10x Genomics, Inc. All rights reserved.
"""Disables secondary analysis if there are more than a certain number of barcodes."""

import cellranger.matrix as cr_matrix

__MRO__ = """
stage DISABLE_SECONDARY_ANALYSIS(
    in  bool is_spatial,
    in  h5   filtered_matrices_h5,
    in  bool no_secondary_analysis,
    in  bool is_visium_hd_main_run  "Boolean indicating if this is being called from a main (not-binning) Visium HD run",
    out bool no_secondary_analysis,
    out bool sketch_graphclust      "Boolean indicating if sketch based graphclust should be enabled",
    out bool parallel_umap          "Boolean indicating if parallel UMAP should be enabled",
    src py   "stages/common/disable_secondary_analysis",
) using (
    volatile = strict,
)
"""

# max limit on barcodes at which we run secondary analysis
MAX_BARCODES_FOR_SEC_ANALYSIS = 2_000_000

# max limit for spatial data (visium-hd).
MAX_BARCODES_FOR_SEC_ANALYSIS_SPATIAL = 2_000_000

# max limit for computation on full dataset for spatial data (visium-hd).
SPATIAL_HIGH_BARCODE_LIMIT = 1_000_000


def main(args, outs):
    if args.no_secondary_analysis or args.is_visium_hd_main_run:
        outs.no_secondary_analysis = True
        return

    max_bcs = (
        MAX_BARCODES_FOR_SEC_ANALYSIS_SPATIAL if args.is_spatial else MAX_BARCODES_FOR_SEC_ANALYSIS
    )
    _, num_bcs, _ = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrices_h5)
    outs.no_secondary_analysis = num_bcs >= max_bcs
    spatial_high_barcodes = args.is_spatial and num_bcs >= SPATIAL_HIGH_BARCODE_LIMIT
    outs.sketch_graphclust = spatial_high_barcodes
    outs.parallel_umap = spatial_high_barcodes
