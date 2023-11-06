#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Commonly used spatial specific constants."""

from __future__ import annotations

from cellranger.constants import AGG_CLOUPE_FIELD, AGG_H5_FIELD, aggr_files

# Spatial AGGR specific
AGG_SPATIAL_FIELD = "spatial_folder"
AGG_TISSUE_POSITION_FIELD = "tissue_position"
AGG_SCALE_FACTORS_FIELD = "scale_factors"
AGG_HIRES_IMAGES_FIELD = "hires_images"
AGG_LOWRES_IMAGES_FIELD = "lowres_images"


# Definition of spatial aggr files
SPATIAL_AGGR_FILES = {
    AGG_H5_FIELD: aggr_files(paths=["molecule_info.h5"], required=True, default_location=""),
    AGG_CLOUPE_FIELD: aggr_files(paths=["cloupe.cloupe"], required=True, default_location=""),
    AGG_TISSUE_POSITION_FIELD: aggr_files(
        paths=["tissue_positions_list.csv", "tissue_positions.csv"],
        required=True,
        default_location="spatial",
    ),
    AGG_SCALE_FACTORS_FIELD: aggr_files(
        paths=["scalefactors_json.json"], required=True, default_location="spatial"
    ),
    AGG_HIRES_IMAGES_FIELD: aggr_files(
        paths=["tissue_hires_image.png"], required=True, default_location="spatial"
    ),
    AGG_LOWRES_IMAGES_FIELD: aggr_files(
        paths=["tissue_lowres_image.png"], required=True, default_location="spatial"
    ),
}
