#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
#
# Determine the locations of the cell barcode, cDNA sequence, UMI, sample index.
#

from __future__ import annotations

import json
import os
from collections import OrderedDict
from typing import TypedDict

from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

# NOTE: The rust code in `cr_types` also loads chemistry definitions from the same file.
CHEMISTRY_DEFS: dict[str, ChemistryDef] = json.load(
    open(os.path.join(os.path.dirname(__file__), "chemistry_defs.json"))
)

# HT chemistries
HT_CHEMISTRIES = [
    CHEMISTRY_DEFS["SC3Pv3HT-polyA"],
    CHEMISTRY_DEFS["SC3Pv3HT-CS1"],
    CHEMISTRY_DEFS["SC5PHT"],
]

# LT v3 Chemistry
CHEMISTRY_SC3P_LT = CHEMISTRY_DEFS["SC3Pv3LT"]

# Single Cell Feature-Barcoding only
CHEMISTRY_SC_FB = CHEMISTRY_DEFS["SC-FB"]

# Spatial Chemistry
CHEMISTRY_SPATIAL3P_V1 = CHEMISTRY_DEFS["SPATIAL3Pv1"]

# ARC chemistry
CHEMISTRY_ARC_V1 = CHEMISTRY_DEFS["ARC-v1"]

# Single Cell 3' chemistries
CHEMISTRY_SC3P_V3_POLYA = CHEMISTRY_DEFS["SC3Pv3-polyA"]
CHEMISTRY_SC3P_V3_CS1 = CHEMISTRY_DEFS["SC3Pv3-CS1"]

SC3P_V3_CHEMISTRIES = [
    CHEMISTRY_SC3P_V3_POLYA,
    CHEMISTRY_SC3P_V3_CS1,
]

SC3P_V4_PREFIX = "SC3Pv4"

SC3P_V4_CHEMISTRIES = [
    CHEMISTRY_DEFS["SC3Pv4-polyA"],
    CHEMISTRY_DEFS["SC3Pv4-CS1"],
    CHEMISTRY_DEFS["SC3Pv4-polyA-OCM"],
    CHEMISTRY_DEFS["SC3Pv4-CS1-OCM"],
]
SC3P_CHEMISTRIES = [
    CHEMISTRY_DEFS["SC3Pv1"],
    CHEMISTRY_DEFS["SC3Pv2"],
    CHEMISTRY_SC3P_V3_POLYA,
    CHEMISTRY_SC3P_V3_CS1,
] + SC3P_V4_CHEMISTRIES


SC5P_V3_CHEMISTRIES = [
    # v3 versions, with a 3.6M whitelist and 12bp UMI
    CHEMISTRY_DEFS["SC5P-PE-v3"],
    CHEMISTRY_DEFS["SC5P-R2-v3"],
    CHEMISTRY_DEFS["SC5P-R2-OCM-v3"],
    CHEMISTRY_DEFS["SC5P-R1-v3"],
]

# Single Cell 5' chemistries
SC5P_CHEMISTRIES = [
    # Paired-end reads, R2-only, or R1-only
    CHEMISTRY_DEFS["SC5P-PE"],
    CHEMISTRY_DEFS["SC5P-R2"],
    CHEMISTRY_DEFS["SC5P-R1"],
] + SC5P_V3_CHEMISTRIES

SC_GEMX_CHEMISTRIES = SC3P_V4_CHEMISTRIES + SC5P_V3_CHEMISTRIES

# Single Cell V(D)J (5-prime) chemistries
SCVDJ_CHEMISTRIES = [
    CHEMISTRY_DEFS["SCVDJ"],
    CHEMISTRY_DEFS["SCVDJ-R2"],
    # v3 versions, with a 3.6M whitelist and 12bp UMI
    CHEMISTRY_DEFS["SCVDJ-v3"],
    CHEMISTRY_DEFS["SCVDJ-R2-v3"],
]

# Aliases for human usage
CHEMISTRY_ALIASES = OrderedDict([("threeprime", "SC3P_auto"), ("fiveprime", "SC5P_auto")])

DEFINED_CHEMISTRIES = list(CHEMISTRY_DEFS.values())

# User-defined chemistry (use the various read_type/offset/length arguments passed to the pipeline)
CUSTOM_CHEMISTRY_NAME = "custom"

# The field in the JSON
CHEMISTRY_DESCRIPTION_FIELD = "description"
CHEMISTRY_NAME_FIELD = "name"


def get_chemistry(name) -> ChemistryDef:
    """Returns a chemistry definition dict for a given name."""
    name = CHEMISTRY_ALIASES.get(name, name)

    chemistries = [n for n in DEFINED_CHEMISTRIES if n["name"] == name]
    if len(chemistries) == 0:
        raise ValueError(f"Could not find chemistry named {name}")
    if len(chemistries) > 1:
        raise ValueError(f"Found multiple chemistries with name {name}")
    return chemistries[0]


def get_chemistry_from_description(description: str) -> ChemistryDef:
    """Given all the loaded chemistries, return the one with a matching description."""
    chemistries = [n for n in DEFINED_CHEMISTRIES if n[CHEMISTRY_DESCRIPTION_FIELD] == description]
    if len(chemistries) == 0:
        raise ValueError(f"Could not find chemistry with description {description}")
    if len(chemistries) > 1:
        raise ValueError(f"Found multiple chemistries with description {description}")
    return chemistries[0]


def get_whitelist_name_from_chemistry_description(description: str) -> None | str:
    """Try to load a chemistry whitelist from a description of a chemistry.

    Args:
        description:

    Returns:
        the whitelist if it's a simple chemistry or None, since we don't enforce that a whitelist exists for all
        chemistries explicitly in code
        (it's just a JSON file).
    """
    chem = get_chemistry_from_description(description)
    # TODO: We need a refactor here so chemistries aren't a dict but an object instead
    barcode = "barcode"
    barcode_section = chem.get(barcode)
    if barcode_section is not None and len(barcode_section) == 1:
        return barcode_section[0].get("whitelist")
    return None


class ChemistryDef(TypedDict):
    """A chemistry definition."""

    name: str
    description: str
    barcode: list[dict[str, int | str | dict[str, str]]]
    umi: list[dict[str, int | str]]
    rna: dict[str, int | str | None]
    rna2: dict[str, int | str | None] | None
    endedness: str
    strandedness: str


ChemistryDefs = dict[str, ChemistryDef]


def get_primary_chemistry_def(chemistry_defs: ChemistryDefs) -> ChemistryDef:
    """Return the only chemistry definition when there is only one.

    Return the GEX chemistry definition when there are multiple.
    Fail when there are multiple chemistry definitions and no GEX.

    This matches the implementation in the Rust ChemistryDefs type.
    """
    if len(chemistry_defs) == 1:
        return next(iter(chemistry_defs.values()))
    try:
        return chemistry_defs[GENE_EXPRESSION_LIBRARY_TYPE]
    except KeyError as exc:
        raise ValueError("no gene expression chemistry found") from exc


def get_whitelists_from_chemistry_defs(defs: ChemistryDefs) -> set[str]:
    """Get all gel bead whitelists from the provided dict of chemistry defs."""
    whitelists = set()
    for chem in defs.values():
        barcode_read_component = chem["barcode"][0]
        assert barcode_read_component["kind"] == "gel_bead"
        whitelist = barcode_read_component["whitelist"]["name"]
        assert isinstance(whitelist, str)
        whitelists.add(whitelist)
    return whitelists
