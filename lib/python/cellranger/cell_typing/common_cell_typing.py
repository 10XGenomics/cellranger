#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Cell Typing utility functions."""

from dataclasses import asdict, dataclass

# Method names
BROAD_TENX_MODELS = "broad_tenx"
AZIMUTH_MODELS = "azimuth"

# Keys in Postprocessed cell annotation out
FINE_CELL_TYPES_KEY = "fine_cell_type"
COARSE_CELL_TYPES_KEY = "coarse_cell_type"
BARCODE_KEY = "barcode"
MEDIUM_CELL_TYPES_KEY = "medium_cell_type"
FULL_CELL_TYPES_KEY = "full_cell_type"


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
