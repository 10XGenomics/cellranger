#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Constants used for the targeting pipeline."""


from __future__ import annotations

from typing import NamedTuple

# two types of targeted assays
TARGETING_METHOD_TL = "templated_ligation"
TARGETING_METHOD_TL_FILE_NAME = "probe set"
TARGETING_METHOD_TL_FILE_FORMAT = "probe_set_file_format"
TARGETING_METHOD_TL_OLIGO_NAME = "probe"
TARGETING_METHOD_TL_SPLICED = "spliced"
TARGETING_METHOD_TL_UNSPLICED = "unspliced"
TARGETING_METHOD_TL_ALLOWED_REGIONS = [
    TARGETING_METHOD_TL_SPLICED,
    TARGETING_METHOD_TL_UNSPLICED,
    "other",
]

TARGETING_METHOD_HC = "hybrid_capture"
TARGETING_METHOD_HC_FILE_NAME = "target panel"
TARGETING_METHOD_HC_FILE_FORMAT = "target_panel_file_format"
TARGETING_METHOD_HC_OLIGO_NAME = "bait"

TARGETING_METHOD_FILE_FORMAT_REGEX = (
    "^[0-9]+[.][0-9]$"  # single-digit Y value to allow floating point comparisons
)

GDNA_CONTENT_METRIC = "estimated_gdna_content"
GDNA_UNSPLICED_THRESHOLD = "estimated_gdna_unspliced_threshold"
GDNA_PLOT_NAME = "gDNA_BiPlot"

TARGETING_METHOD_FILE_NAMES = {
    TARGETING_METHOD_TL: TARGETING_METHOD_TL_FILE_NAME,
    TARGETING_METHOD_HC: TARGETING_METHOD_HC_FILE_NAME,
}

TARGETING_METHOD_FILE_FORMATS = {
    TARGETING_METHOD_TL: TARGETING_METHOD_TL_FILE_FORMAT,
    TARGETING_METHOD_HC: TARGETING_METHOD_HC_FILE_FORMAT,
}

OLIGO_NAME = {
    TARGETING_METHOD_TL: TARGETING_METHOD_TL_OLIGO_NAME,
    TARGETING_METHOD_HC: TARGETING_METHOD_HC_OLIGO_NAME,
}


class TargetingMethod(NamedTuple):
    """Aggregate type for holding the constants associated with a probe set or target panel, along with file versions."""

    method: str
    descriptive_name: str
    file_format_tag: str
    oligo_name: str
    file_version: str

    @staticmethod
    def get_file_name(method):
        return TARGETING_METHOD_FILE_NAMES.get(
            method, f"{TARGETING_METHOD_TL_FILE_NAME} or {TARGETING_METHOD_HC_FILE_NAME}"
        )

    @staticmethod
    def get_file_format(method):
        return TARGETING_METHOD_FILE_FORMATS.get(
            method, f"{TARGETING_METHOD_TL_FILE_FORMAT} or {TARGETING_METHOD_HC_FILE_FORMAT}"
        )

    @classmethod
    def get_targeting_method_from_metadata(cls, metadata):
        """Factory takes a metadata dictionary and returns either a completed structure or None.

        if neither file format key is available.
        """
        if TARGETING_METHOD_TL_FILE_FORMAT in metadata:
            return cls(
                TARGETING_METHOD_TL,
                TARGETING_METHOD_TL_FILE_NAME,
                TARGETING_METHOD_TL_FILE_FORMAT,
                TARGETING_METHOD_TL_OLIGO_NAME,
                metadata[TARGETING_METHOD_TL_FILE_FORMAT],
            )
        elif TARGETING_METHOD_HC_FILE_FORMAT in metadata:
            return cls(
                TARGETING_METHOD_HC,
                TARGETING_METHOD_HC_FILE_NAME,
                TARGETING_METHOD_HC_FILE_FORMAT,
                TARGETING_METHOD_HC_OLIGO_NAME,
                metadata[TARGETING_METHOD_HC_FILE_FORMAT],
            )
        return None


ALL_TARGETING_METHODS = list(TARGETING_METHOD_FILE_NAMES.keys())

# List of gene/probe ID prefixes that are excluded from the filtered_feature_bc_matrix.
# Ensure that the corresponding Python and Rust constants are identical.
EXCLUDED_PROBE_ID_PREFIXES = (
    b"DEPRECATED",
    b"Hum-",
    b"IGNORE",
    b"NC-",
    b"VAR_",
    b"VDJ_",
)
