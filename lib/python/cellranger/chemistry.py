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

# NOTE: The rust code in `cr_types` also loads chemistry definitions from the same file.
CHEMISTRY_DEFS = json.load(open(os.path.join(os.path.dirname(__file__), "chemistry_defs.json")))


# Single Cell 3' chemistries
CHEMISTRY_SC3P_V1 = CHEMISTRY_DEFS["SC3Pv1"]
CHEMISTRY_SC3P_V2 = CHEMISTRY_DEFS["SC3Pv2"]
CHEMISTRY_SC3P_V3 = CHEMISTRY_DEFS["SC3Pv3"]

# Single Cell V(D)J
CHEMISTRY_SCVDJ = CHEMISTRY_DEFS["SCVDJ"]
# Single Cell V(D)J, transcript on R2 only
CHEMISTRY_SCVDJ_R2 = CHEMISTRY_DEFS["SCVDJ-R2"]
# LT v3 Chemistry
CHEMISTRY_SC3P_LT = CHEMISTRY_DEFS["SC3Pv3LT"]

# Single Cell 5' Gene Expression, paired-end
CHEMISTRY_SC5P_PE = CHEMISTRY_DEFS["SC5P-PE"]
# Single Cell 5' Gene Expression, transcript on R2 only
CHEMISTRY_SC5P_R2 = CHEMISTRY_DEFS["SC5P-R2"]
# Single Cell 5' Gene Expression, transcript on R1 only
CHEMISTRY_SC5P_R1 = CHEMISTRY_DEFS["SC5P-R1"]

# Single Cell Feature-Barcoding only
CHEMISTRY_SC_FB = CHEMISTRY_DEFS["SC-FB"]

# Spatial Chemistry
CHEMISTRY_SPATIAL3P_V1 = CHEMISTRY_DEFS["SPATIAL3Pv1"]

# ARC chemistry
CHEMISTRY_ARC_V1 = CHEMISTRY_DEFS["ARC-v1"]

SC3P_CHEMISTRIES = [CHEMISTRY_SC3P_V1, CHEMISTRY_SC3P_V2, CHEMISTRY_SC3P_V3]

SC5P_CHEMISTRIES = [CHEMISTRY_SC5P_PE, CHEMISTRY_SC5P_R2, CHEMISTRY_SC5P_R1]


DEFINED_CHEMISTRIES = list(CHEMISTRY_DEFS.values())


# Aliases for human usage
CHEMISTRY_ALIASES = OrderedDict([("threeprime", "SC3P_auto"), ("fiveprime", "SC5P_auto")])

# User-defined chemistry (use the various read_type/offset/length arguments passed to the pipeline)
CUSTOM_CHEMISTRY_NAME = "custom"

# The field in the JSON
CHEMISTRY_DESCRIPTION_FIELD = "description"


def get_chemistry(name):
    """Returns a chemistry definition dict for a given name."""
    name = CHEMISTRY_ALIASES.get(name, name)

    chemistries = [n for n in DEFINED_CHEMISTRIES if n["name"] == name]
    if len(chemistries) == 0:
        raise ValueError("Could not find chemistry named %s" % name)
    if len(chemistries) > 1:
        raise ValueError("Found multiple chemistries with name %s" % name)
    return chemistries[0]


def get_chemistry_from_description(description: str):
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
