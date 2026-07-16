#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Cell Typing utility functions."""

import csv
from collections import Counter
from dataclasses import asdict, dataclass

# Method names
BROAD_TENX_MODELS = "10x_Cloud"
AZIMUTH_MODELS = "Azimuth"

# Keys in Postprocessed cell annotation out
FINE_CELL_TYPES_KEY = "fine_cell_type"
BROAD_CELL_TYPES_KEY = "broad_cell_type"
COARSE_CELL_TYPES_KEY = "coarse_cell_type"
DISPLAY_MAP_KEY = "coarse_cell_type"
BARCODE_KEY = "barcode"
FULL_CELL_TYPES_KEY = "full_cell_type"
SUMMARY_SCORE_KEY = "summary_score"
AZIMUTH_SUMMARY_SCORE_KEY = "final_level_softmax_prob"


@dataclass
class CellTypeResults:
    """Common structure for results from remote cell annotation stages."""

    cell_types: str | None = None
    results: str | None = None
    metadata: str | None = None
    frac_returned_bcs: float = 0.0
    model_used: str | None = None
    skip_downstream: bool = False

    def to_outs(self):
        """Return the value this form needs to take to be in stage outs."""
        return asdict(self)


@dataclass
class CellAnnotationMetadata:
    """Dataclass corresponding to metadata for cell annotation."""

    tree_version_used: str | None = None
    display_map_version_used: str | None = None
    fraction_non_informative_annotations: float | None = None
    is_beta_model: bool | None = None
    developer: str | None = None


GENOMES_WITH_SUPPORTED_LOCAL_MODELS = ["hg19", "GRCh38", "GRCh39"]


def get_top_cell_types(cell_types_path: str, level: str, num_cells: int):
    """Get the top cell types from a level from the cell_types.csv file.

    Args:
        cell_types_path: Path to the cell_types.csv file
        level: COARSE_CELL_TYPES_KEY or FINE_CELL_TYPES_KEY cell types
        num_cells: Number of top cell types to return

    Returns:
        dict: Dictionary with numeric keys (starting at 1) and cell type names
            as values, with None for positions without enough cell types
    """
    counter = Counter()

    with open(cell_types_path) as f:
        reader = csv.DictReader(f)
        counter.update(row[level] for row in reader if row.get(level))

    result = {}
    for i in range(1, num_cells + 1):
        if i - 1 < len(counter):
            cell_type, _ = counter.most_common(num_cells)[i - 1]
            result[i] = cell_type
        else:
            result[i] = None

    return result


def get_cell_type_counts(cell_types_path: str, level: str, cell_type: str):
    """Get the count of a specific cell type level from the cell_types.csv file.

    Args:
        cell_types_path: Path to the cell_types.csv file
        level: COARSE_CELL_TYPES_KEY or FINE_CELL_TYPES_KEY cell types
        cell_type: Name of the cell type to count

    Returns:
        int: Count of cell type requested
    """
    with open(cell_types_path) as f:
        reader = csv.DictReader(f)
        return sum(1 for row in reader if row.get(level) == cell_type)
