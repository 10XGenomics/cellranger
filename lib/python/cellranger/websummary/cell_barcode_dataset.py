# Copyright (c) 2026 10X Genomics, Inc. All rights reserved.
"""Build the cell-barcode dataset for web summary reports."""

from __future__ import annotations

import hashlib
import random
from collections.abc import Iterable

MAX_BARCODES = 10_000


def build_cell_barcode_dataset(
    barcodes: Iterable[bytes | str],
    sample_id: str,
    max_barcodes: int = MAX_BARCODES,
) -> dict:
    """Return a dataset dict containing cell-associated barcodes for the web summary.

    If len(barcodes) <= max_barcodes every barcode is included; otherwise a
    deterministic random subset of max_barcodes is chosen, seeded from sample_id
    so the same sample always produces the same subset.

    Args:
        barcodes: iterable of barcode sequences (bytes or str).
        sample_id: pipeline sample ID used to seed the RNG.
        max_barcodes: maximum number of barcodes to report.

    Returns:
        dict with keys total_cells, subsampled, max_reported, barcodes.
    """
    decoded = sorted(bc.decode() if isinstance(bc, bytes) else bc for bc in barcodes)
    total = len(decoded)
    if total > max_barcodes:
        seed = int(hashlib.blake2b(sample_id.encode(), digest_size=8).hexdigest(), 16)
        decoded = random.Random(seed).sample(decoded, max_barcodes)
        decoded.sort()
        subsampled = True
    else:
        subsampled = False
    return {
        "total_cells": total,
        "subsampled": subsampled,
        "max_reported": max_barcodes,
        "barcodes": decoded,
    }
