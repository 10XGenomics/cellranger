#!/usr/bin/env python
#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Constants associated with the segmentation data."""

FILTERED_CELLS = "filtered_cells"
FRACTION_COUNTS_PER_CELL = "fraction_counts_per_cell"
FRACTION_READS_IN_CELLS = "fraction_reads_in_cells"
MEAN_READS_PER_CELL = "mean_reads_per_cell"
MEAN_COUNTS_PER_CELL = "mean_counts_per_cell"
MEDIAN_GENES_PER_CELL = "median_genes_per_cell"
MEDIAN_COUNTS_PER_CELL = "median_counts_per_cell"
MEDIAN_CELL_AREA = "median_cell_area"
MEDIAN_NUCLEUS_AREA = "median_nucleus_area"
FRACTION_NUCLEI_EXPANDED = "fraction_nuclei_expanded"
SEGMENTATION_SUFFIX = "segmentation_suffix"
MAX_NUCLEUS_DIAMETER = "max_nucleus_diameter_px"

USER_PROVIDED_SEGMENTATION_KEYS = [
    "user_provided_segmentations",
    "square_barcode_to_cell_map",
    "instance_mask_tiff",
    "instance_mask_npy",
]

GEOJSON_CELL_ID_KEY = "cell_id"
GEOJSON_NUCLEUS_CENTROID_KEY = "nucleus_centroid"
GEOJSON_CELL_CENTROID_KEY = "cell_centroid"


GEOJSON_CLASSIFICATION_KEY = "classification"
GEOJSON_CLASSIFICATION_NAME_KEY = "name"
GEOJSON_CLASSIFICATION_COLOR_KEY = "color"
