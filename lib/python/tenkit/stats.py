#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

"""Utils for math operations."""

from __future__ import annotations

from collections.abc import Iterable


def NX(lengths: Iterable[int], fraction: float) -> int:
    """Calculates the N50 for a given set of lengths."""
    lengths = sorted(lengths, reverse=True)
    tot_length = sum(lengths)
    cum_length = 0
    for length in lengths:
        cum_length += length
        if cum_length >= fraction * tot_length:
            return length
    raise ArithmeticError("this should be impossible")


def robust_divide(a: float | int | str | bytes, b: float | int | str | bytes) -> float:
    """Handles 0 division and conversion to floats automatically."""
    a = float(a)
    b = float(b)
    if b == 0:
        return float("NaN")
    else:
        return a / b
