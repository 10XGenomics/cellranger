#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Chain types used by the VDJ pipeline."""
from __future__ import annotations

import cellranger.vdj.constants as vdj_constants

# Chain types
TR_CHAIN_TYPE = b"TR"
TR_GD_CHAIN_TYPE = b"TR_GD"
IG_CHAIN_TYPE = b"IG"
VDJ_CHAIN_TYPES = [TR_CHAIN_TYPE, TR_GD_CHAIN_TYPE, IG_CHAIN_TYPE]
AUTO_CHAIN_TYPE = b"auto"  # Autodetect the chain type filter
CHAIN_TYPE_SPECS = VDJ_CHAIN_TYPES + [AUTO_CHAIN_TYPE]


def load_chains_from_chain_type(chain_type: bytes | str):
    """Get the chains for a chain type.

    Args:
        chain_type: (Optional) One of vdj_constants.CHAIN_TYPE_SPECS

    Returns:
        Chains corresponding to the chain type. For e.g. TRA, TRB for "TR" chain
    """
    if isinstance(chain_type, str):
        chain_type = chain_type.encode()

    assert chain_type in CHAIN_TYPE_SPECS
    if chain_type == TR_CHAIN_TYPE:
        return vdj_constants.CANONICAL_TR_GENES
    elif chain_type == IG_CHAIN_TYPE:
        return vdj_constants.CANONICAL_IG_GENES
    else:
        return vdj_constants.CANONICAL_VDJ_GENES
