#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import base64
import struct


def array_to_float32_base64(vals):
    """Convert an array of doubles to a little endian base64 encoded ASCII string.

    Args:
        vals: a list of float values

    Returns:
        A byte string.
    """
    assert isinstance(vals, list)
    for val in vals:
        assert isinstance(val, float)
    format_str = "<" + str(len(vals)) + "f"
    return b"data:float32;base64," + base64.b64encode(struct.pack(format_str, *vals))


def _round_float(x: float) -> float:
    """Round a float to the nearest value which can be represented with 4 decimal digits.

    In order to avoid representing floats with full precision, we convert them
    to a lower precision string and then convert that back to a float for use by `json.dumps`
    which will typically then output many fewer digits if possible.
    """
    return float(f"{x:.4g}")


def round_floats_in_list(x):
    """Lower the precision for a whole bunch of floats."""
    return [_round_float(x) for x in x]
