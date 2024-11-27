#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Functions and constants for all post processing of CAS."""

import csv
import json
from pathlib import Path

# Key in cell annotation out
BARCODE_KEY = "barcode"
SCORE_KEY = "score"
CELL_TYPES_ID_KEY_IN = "cell_type_ontology_term_id"
MATCHES_KEY_IN = "matches"

# Keys in Postprocessed cell annotation out
FINE_CELL_TYPES_KEY = "fine_cell_type"
COARSE_CELL_TYPES_KEY = "coarse_cell_type"


# Key in consolidated cell typing matrix created
CELL_TYPES_KEY = "cas_cell_type"
GROUND_TRUTH_CELL_TYPES_KEY = "gt_cell_type"
CELL_TYPES_AFTER_MAPPING_KEY = "cas_mapped_cell_type"
GROUND_TRUTH_CELL_TYPES_AFTER_MAPPING_KEY = "gt_mapped_cell_type"

# Values when no GT and no cell annotation in cell typing matrix created
NO_GROUND_TRUTH_VALUE = "no ground truth"
NO_CAS_VALUE = "no cas annotation"
DEFAULT_CELL_MAPPING = "cell"
DEFAULT_CELL_MAPPING_ID = "CL:0000000"

# file to get the cell annotation cell map from
CAS_CELL_TYPE_MAP_BASE_DIRECTORY = Path(__file__).parent / "ontology_tree"
CAS_CELL_TYPE_MAP_FILENAME = "V5_mappings.csv"
CAS_CELL_TYPE_MAP_PATH = CAS_CELL_TYPE_MAP_BASE_DIRECTORY / CAS_CELL_TYPE_MAP_FILENAME

# Keys of the cell type maps
SOURCE_NODE_MAP_KEY = "source_node_name"
TARGET_NODE_MAP_KEY = "target_node_name"

# file to get the cell annotation display map from
CAS_DISPLAY_MAP_BASE_DIRECTORY = Path(__file__).parent / "display_map"
CAS_DISPLAY_MAP_FILENAME = "celltype_display_mapping_v6.json"
CAS_DISPLAY_MAP_PATH = CAS_DISPLAY_MAP_BASE_DIRECTORY / CAS_DISPLAY_MAP_FILENAME

# Keys of metrics JSON
AFTER_MAPPING_SUCCESS_RATE_KEY = "after_mapping_success_rate"
BEFORE_MAPPING_SUCCESS_RATE_KEY = "before_mapping_success_rate"
AFTER_MAPPING_SUCCESS_RATE_WITHOUT_DEFAULT_MAPPING_GT_KEY = (
    "after_mapping_success_rate_ignoring_gt_cell"
)
BEFORE_MAPPING_SUCCESS_RATE_WITHOUT_DEFAULT_MAPPING_GT__KEY = (
    "before_mapping_success_rate_ignoring_gt_annotations_mapped_to_cell"
)
SAMPLE_ID_KEY = "sample_id"
SAMPLE_DESCRIPTION_KEY = "sample_desc"
CAS_MODEL_KEY = "cell_annotation_model"
NUM_BCS_SAMPLED_KEY = "number_bcs_to_use"


def get_coarse_cell_type_map() -> dict[str, str]:
    """Get map between OBO celltypes and cell types we keep."""
    with open(CAS_CELL_TYPE_MAP_PATH) as f:
        tree_dct = dict(
            (
                row[SOURCE_NODE_MAP_KEY].replace(",", "").lower(),
                row[TARGET_NODE_MAP_KEY].replace(",", "").lower(),
            )
            for row in csv.DictReader(f)
        )

    with open(CAS_DISPLAY_MAP_PATH) as jsonfile:
        display_long_dct = json.load(jsonfile)

    display_dct = {z.replace(",", "").lower(): x for x, y in display_long_dct.items() for z in y}

    dct = {x: display_dct[y] for x, y in tree_dct.items()}
    dct[NO_CAS_VALUE] = NO_CAS_VALUE
    dct[NO_GROUND_TRUTH_VALUE] = NO_GROUND_TRUTH_VALUE
    return dct


def get_fine_cell_type_map() -> dict[str, str]:
    """Get map between OBO celltypes and fine cell types we keep."""
    with open(CAS_CELL_TYPE_MAP_PATH) as f:
        dct = dict(
            (
                row[SOURCE_NODE_MAP_KEY].replace(",", "").lower(),
                row[TARGET_NODE_MAP_KEY],
            )
            for row in csv.DictReader(f)
        )

    dct[NO_CAS_VALUE] = NO_CAS_VALUE
    dct[NO_GROUND_TRUTH_VALUE] = NO_GROUND_TRUTH_VALUE
    return dct


def get_bucket_dict() -> dict[str, int]:
    """Get cell type to bucket map."""
    dct = get_coarse_cell_type_map()
    level_2_keys = set(x.lower().replace(",", "") for x in dct.values()) - {
        NO_CAS_VALUE,
        NO_GROUND_TRUTH_VALUE,
        DEFAULT_CELL_MAPPING,
    }

    bucket_dct = {}
    for key, value in dct.items():
        if key in [NO_CAS_VALUE, NO_GROUND_TRUTH_VALUE, DEFAULT_CELL_MAPPING]:
            bucket_dct[key] = 0
        elif value in [NO_CAS_VALUE, NO_GROUND_TRUTH_VALUE, DEFAULT_CELL_MAPPING]:
            bucket_dct[key] = 1
        elif key in level_2_keys:
            bucket_dct[key] = 2
        else:
            bucket_dct[key] = 3

    return bucket_dct


def get_tree_version() -> str:
    """Get version of tree being used."""
    return CAS_CELL_TYPE_MAP_PATH.stem


def get_display_map_version() -> str:
    """Get version of tree being used."""
    return CAS_DISPLAY_MAP_PATH.stem
