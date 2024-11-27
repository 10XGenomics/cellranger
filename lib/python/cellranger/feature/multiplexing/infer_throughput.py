#!/usr/bin/env python3
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
#
"""Infer Gem well throughput from the barcode rankplot."""

from __future__ import annotations

import cellranger.cell_calling_helpers as cch
from cellranger.feature.throughputs import HT_THROUGHPUT, MT_THROUGHPUT

ANCHOR_BARCODE_IDX = 99
THRESHOLD_BARCODE_IDX = 129999
MT_MAX_COUNT_RATIO = 0.0003
FIRST_BC_IDX_ON_RANKPLOT = 70000
SECOND_BC_IDX_ON_RANKPLOT = 220000


def infer_throughput_from_background_counts(counts_per_bc: list[int]) -> str:
    """Infer throughput (HT vs MT) from the background barcode counts.

    Args:
        counts_per_bc (List[int]): Sorted UMI counts per barcode. Must be in descending order.

    Returns:
        str: "HT" or "MT"
    """
    if len(counts_per_bc) < THRESHOLD_BARCODE_IDX + 1:
        return MT_THROUGHPUT
    if (
        counts_per_bc[THRESHOLD_BARCODE_IDX]
        <= counts_per_bc[ANCHOR_BARCODE_IDX] * MT_MAX_COUNT_RATIO
    ):
        return MT_THROUGHPUT
    else:
        return HT_THROUGHPUT


def infer_throughput_from_rankplot_gradient(counts_per_bc: list[int]) -> str:
    """Infer throughput (HT vs MT) by finding the steepest gradient in the region after cell barcodes."""
    if len(counts_per_bc) < FIRST_BC_IDX_ON_RANKPLOT:
        return None, MT_THROUGHPUT
    outs = cch.filter_cellular_barcodes_gradient(
        counts_per_bc[FIRST_BC_IDX_ON_RANKPLOT:SECOND_BC_IDX_ON_RANKPLOT],
        recovered_cells=None,
        infer_throughput=True,
    )
    if not outs[0].size or not outs[0].any():  # usually indicates atypical rank plots (low depth)
        return None, MT_THROUGHPUT

    slope_bc_idx = FIRST_BC_IDX_ON_RANKPLOT + outs[0][-1]
    inferred_throughput = MT_THROUGHPUT if slope_bc_idx <= THRESHOLD_BARCODE_IDX else HT_THROUGHPUT
    return slope_bc_idx, inferred_throughput
