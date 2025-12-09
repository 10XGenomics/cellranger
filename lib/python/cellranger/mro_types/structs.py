#!/usr/bin/env python
#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Martian Python types."""
from typing import TypedDict

from cellranger.mro_types.filetypes import json


class WhitelistSpec(TypedDict):
    """Different ways to specify a whitelist of sequences (typically barcode sequence)."""

    name: str
    translation: bool
    strand: str
    translation_whitelist_path: str
    slide: str
    part: str


class UmiWhitelistSpec(TypedDict):
    """UMI whitelist specification."""

    slide: str
    part: str
    translation: str


class UmiReadComponent(TypedDict):
    """UMI read components."""

    read_type: str
    offset: int
    length: int  # The length of the UMI. At most this number of bases will be extracted for use as a UMI.
    min_length: int  # If a shorter UMI can be used, add it here. None indicates that the full
    whitelist: UmiWhitelistSpec


class RnaReadComponent(TypedDict):
    """RNA read components."""

    read_type: str
    offset: int
    length: int
    min_length: int


class BarcodeReadComponent(TypedDict):
    """Barcode read components."""

    read_type: str
    kind: str
    offset: int
    length: int
    whitelist: WhitelistSpec


class ChemistryDef(TypedDict):
    """Chemistry definition."""

    name: str
    description: str
    endedness: str
    strandedness: str
    barcode: list[BarcodeReadComponent]
    umi: list[UmiReadComponent]
    rna: RnaReadComponent
    rna2: RnaReadComponent
    barcode_extraction: dict


class ProbeBCDef(TypedDict):
    """Probe barcode definition."""

    id: str
    sequence: list[str]
    offset: int
    length: int


class CellCallingParam(TypedDict):
    """Cell calling parameter."""

    per_gem_well: int
    per_sample: dict[str, int]


class CellCalling(TypedDict):
    """Carries options to customize the cell calling mode."""

    recovered_cells: CellCallingParam
    force_cells: CellCallingParam
    emptydrops_minimum_umis: CellCallingParam
    global_minimum_umis: CellCallingParam
    max_mito_percent: CellCallingParam
    cell_barcodes: json
    override_mode: str
    override_library_types: list[str]
    disable_ab_aggregate_detection: bool
    disable_high_occupancy_gem_detection: bool
