#!/usr/bin/env python3
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#


from __future__ import annotations

import json
import re
from collections.abc import Iterable
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cellranger.sample_bag import SequencingLibrary

GENE_EXPRESSION_LIBRARY_TYPE = "Gene Expression"
CRISPR_LIBRARY_TYPE = "CRISPR Guide Capture"
ANTIBODY_LIBRARY_TYPE = "Antibody Capture"
ANTIGEN_LIBRARY_TYPE = "Antigen Capture"
MULTIPLEXING_LIBRARY_TYPE = "Multiplexing Capture"
FEATURETEST_LIBRARY_TYPE = "FEATURETEST"
CUSTOM_LIBRARY_TYPE = "Custom"
VDJ_LIBRARY_TYPE = "VDJ"
DEFAULT_LIBRARY_TYPE = GENE_EXPRESSION_LIBRARY_TYPE

MULTI_REFS_PREFIX = "multi"

CUSTOM_METRIC_PREFIX = CUSTOM_LIBRARY_TYPE
DISPLAY_PREFIX_CUSTOM = "Custom"
CRISPR_METRIC_PREFIX = "CRISPR"
DISPLAY_PREFIX_CRISPR = "CRISPR:"
ANTIBODY_METRIC_PREFIX = "ANTIBODY"
DISPLAY_PREFIX_ANTIBODY = "Antibody:"
ANTIGEN_METRIC_PREFIX = "ANTIGEN"
DISPLAY_PREFIX_ANTIGEN = "Antigen:"
MULTIPLEXING_METRIC_PREFIX = "MULTIPLEXING"
DISPLAY_PREFIX_MULTIPLEXING = "Multiplexing:"

LIBRARY_TYPE = "library_type"

RECOGNIZED_FEATURE_TYPES = [
    GENE_EXPRESSION_LIBRARY_TYPE,
    CRISPR_LIBRARY_TYPE,
    ANTIBODY_LIBRARY_TYPE,
    ANTIGEN_LIBRARY_TYPE,
    MULTIPLEXING_LIBRARY_TYPE,
]
FEATURE_LIBRARY_TYPES = [
    CRISPR_LIBRARY_TYPE,
    ANTIBODY_LIBRARY_TYPE,
    ANTIGEN_LIBRARY_TYPE,
    MULTIPLEXING_LIBRARY_TYPE,
]

metric_prefix_map = {
    GENE_EXPRESSION_LIBRARY_TYPE: "",
    CRISPR_LIBRARY_TYPE: CRISPR_METRIC_PREFIX,
    ANTIBODY_LIBRARY_TYPE: ANTIBODY_METRIC_PREFIX,
    ANTIGEN_LIBRARY_TYPE: ANTIGEN_METRIC_PREFIX,
    CUSTOM_LIBRARY_TYPE: CUSTOM_METRIC_PREFIX,
    MULTIPLEXING_LIBRARY_TYPE: MULTIPLEXING_METRIC_PREFIX,
}

report_prefix_map = {
    GENE_EXPRESSION_LIBRARY_TYPE: "",
    CRISPR_LIBRARY_TYPE: DISPLAY_PREFIX_CRISPR,
    ANTIBODY_LIBRARY_TYPE: DISPLAY_PREFIX_ANTIBODY,
    ANTIGEN_LIBRARY_TYPE: DISPLAY_PREFIX_ANTIGEN,
    CUSTOM_LIBRARY_TYPE: DISPLAY_PREFIX_CUSTOM,
    MULTIPLEXING_LIBRARY_TYPE: DISPLAY_PREFIX_MULTIPLEXING,
}


# BarcodeMultiplexingType enum from lib/rust/cr_types/src/types.rs
class ReadLevel(Enum):
    RTL = "RTL"
    OH = "OH"


class CellLevel(Enum):
    CMO = "CMO"
    Hashtag = "Hashtag"


class BarcodeMultiplexingType:
    """A class representing the multiplexing type, which can either be cell-level or read-level."""

    def __init__(self, value: str):
        if value == CellLevel.CMO.value:
            self.level = "CellLevel"
            self.type = CellLevel.CMO
        elif value == CellLevel.Hashtag.value:
            self.level = "CellLevel"
            self.type = CellLevel.Hashtag
        elif value == ReadLevel.RTL.value:
            self.level = "ReadLevel"
            self.type = ReadLevel.RTL
        elif value == ReadLevel.OH.value:
            self.level = "ReadLevel"
            self.type = ReadLevel.OH
        else:
            raise ValueError("Invalid type for BarcodeMultiplexingType")

    def __str__(self):
        return self.type.value

    def multiplexing_library_type(self):
        """Returns the library associated with a given cell-level multiplexing type."""
        if self.level == "CellLevel":
            if self.type == CellLevel.CMO:
                return MULTIPLEXING_LIBRARY_TYPE
            elif self.type == CellLevel.Hashtag:
                return ANTIBODY_LIBRARY_TYPE
            else:
                raise ValueError("Invalid CellLevel BarcodeMultiplexingType!")
        else:
            raise ValueError("Multiplexing library is undefined for this BarcodeMultiplexingType!")

    def is_cell_multiplexed(self):
        return self.level == "CellLevel"

    def is_read_multiplexed(self):
        return self.level == "ReadLevel"


# 'target_set_name' should be a key in library_info
TARGET_SET_KEY = "target_set_name"
DEFAULT_TARGET_SETS = ("", None)


def _get_prefix(lib_type: str, sep: str, prefix_map: dict[str, str]) -> str:
    if lib_type == GENE_EXPRESSION_LIBRARY_TYPE:
        return ""
    else:
        value = prefix_map.get(lib_type, lib_type)
        return value + sep


def get_library_type_metric_prefix(lib_type: str) -> str:
    """Get the metric prefix for a given library type."""
    return _get_prefix(lib_type, "_", metric_prefix_map)


def get_library_type_report_prefix(lib_type: str) -> str:
    """Gets the prefix to be used in displayed reports."""
    return _get_prefix(lib_type, " ", report_prefix_map)


def add_multi_prefix(metric: str) -> str:
    """Appends the prefix for cumulative metrics onto a metric name."""
    return f"{MULTI_REFS_PREFIX}_{metric}"


def add_species_prefix(species: str, metric: str) -> str:
    """Append the species/genome name to the front of a metric name."""
    return f"{species}_{metric}"


def get_bam_library_info(bam):
    """Get the library info from a BAM's comment lines.

    Args:
      bam (pysam.AlignmentFile): BAM file

    Returns:
      list of dicts
    """
    comments = bam.header["CO"]
    libraries = []
    for comment in comments:
        m = re.match(r"^library_info:(.+)$", comment)
        if m:
            libraries.append(json.loads(m.group(1)))
    return libraries


def has_genomes(library_type: str):
    """Do genomes make sense for a library type."""
    return library_type == GENE_EXPRESSION_LIBRARY_TYPE


def sorted_library_types(library_info: Iterable[SequencingLibrary]):
    """Sorted list of unique library types in library_info."""
    return sorted({lib[LIBRARY_TYPE] for lib in library_info})


def has_target_set(library):
    return TARGET_SET_KEY in library and library[TARGET_SET_KEY] not in DEFAULT_TARGET_SETS
