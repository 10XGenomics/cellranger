#!/usr/bin/env python
#
# Copyright (c) 2026 10x Genomics, Inc. All rights reserved.
#
"""Compute segmentation plots for Visium HD data."""
import json
import os

__MRO__ = """
stage VALIDATE_SEGMENTATION_DIRECTORY(
    in  path                     segmented_outputs,
    out ValidateSegmentationOutputs validate_segmentation_outputs,
    out bool                     disable_segmentation_spatial_plot,
    src py                       "stages/cas_cell_typing/validate_segmentation_directory",
)
"""


def main(
    args,
    outs,
):
    if args.segmented_outputs is None:
        outs.disable_segmentation_spatial_plot = True
        return

    cell_segmentations_geojson = os.path.join(args.segmented_outputs, "cell_segmentations.geojson")
    tissue_hires_image = os.path.join(args.segmented_outputs, "spatial/tissue_hires_image.png")
    scale_factors_json = os.path.join(args.segmented_outputs, "spatial/scalefactors_json.json")
    if os.path.exists(scale_factors_json):
        with open(scale_factors_json) as f:
            scale_factors = json.load(f)
            tissue_hires_image_scale = scale_factors["tissue_hires_scalef"]
    else:
        outs.disable_segmentation_spatial_plot = True
        return

    outs.validate_segmentation_outputs = {
        "cell_segmentations": cell_segmentations_geojson,
        "tissue_hires_image": tissue_hires_image,
        "scale_factors_json": scale_factors_json,
        "tissue_image_scale": tissue_hires_image_scale,
    }
    outs.disable_segmentation_spatial_plot = (
        not os.path.exists(cell_segmentations_geojson)
        or not os.path.exists(tissue_hires_image)
        or not os.path.exists(scale_factors_json)
    )
