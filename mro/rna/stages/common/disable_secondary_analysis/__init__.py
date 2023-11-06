# Copyright (c) 2023 10x Genomics, Inc. All rights reserved.
"""Disables secondary analysis if there are more than a certain number of barcodes."""

import cellranger.matrix as cr_matrix
import cellranger.molecule_counter as cr_mc

__MRO__ = """
stage DISABLE_SECONDARY_ANALYSIS(
    in  h5   filtered_matrices_h5,
    in  h5   molecule_info,
    in  bool no_secondary_analysis_in,
    out bool no_secondary_analysis,
    src py   "stages/common/disable_secondary_analysis",
)
"""
# max limit on barcodes at which we run secondary analysis
MAX_BARCODES_FOR_SEC_ANALYSIS = 2_000_000
# max limit for spatial data (visium-hd).
MAX_BARCODES_FOR_SEC_ANALYSIS_SPATIAL = 1_000_000


def main(args, outs):
    disable_secondary_analysis = False

    is_spatial = False
    if args.molecule_info is not None:
        mol_info = cr_mc.MoleculeCounter.open(args.molecule_info, "r")
        is_spatial = mol_info.is_spatial_data()

    max_bcs = MAX_BARCODES_FOR_SEC_ANALYSIS_SPATIAL if is_spatial else MAX_BARCODES_FOR_SEC_ANALYSIS

    if args.filtered_matrices_h5 is not None:
        _, num_bcs, _ = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrices_h5)
        disable_secondary_analysis = bool(num_bcs >= max_bcs)
    outs.no_secondary_analysis = disable_secondary_analysis or args.no_secondary_analysis_in
