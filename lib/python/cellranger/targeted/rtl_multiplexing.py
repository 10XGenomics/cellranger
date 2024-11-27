#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Utils for dealing with RTL multiplexing."""
from __future__ import annotations


def _get_barcode_length(barcode_def: list[dict[str, str | int]], barcode_kind: str) -> int | None:
    for barcode in barcode_def:
        if barcode["kind"] == barcode_kind:
            return barcode["length"]
    return None


def get_barcode_locations(barcode_def: list[dict[str, str | int]]) -> tuple[
    tuple[int | None, int | None],
    tuple[int | None, int | None],
    tuple[int | None, int | None],
]:
    """Given barcode definitions returns tuples indicating barcode positions.

    The return value ((gb_start, gb_end), (lhs_start, lhs_end), (rhs_start_rhs_end))
    indicates the positions of the barcodes within the overall cell barcode.

    Args:
        barcode_def (List[Dict[str, Union[str, int]]]): list of barcode definitions from chemistry def

    Returns:
        Tuple[Optional[int], Optional[int]]: (gb_start, gb_end)
        Tuple[Optional[int], Optional[int]]: (lhs_start, lhs_end)
        Tuple[Optional[int], Optional[int]]]: (rhs_start_rhs_end)
    """
    gb_start = 0
    gb_end = _get_barcode_length(barcode_def, "gel_bead")

    lhs_length = _get_barcode_length(barcode_def, "left_probe")
    if lhs_length:
        lhs_start = gb_end
        lhs_end = gb_end + lhs_length
    else:
        lhs_start = None
        lhs_end = None

    rhs_length = _get_barcode_length(barcode_def, "right_probe")
    if rhs_length:
        if lhs_length:
            offset = lhs_end
        else:
            offset = gb_end

        rhs_start = offset
        rhs_end = offset + rhs_length
    else:
        rhs_start = None
        rhs_end = None

    return ((gb_start, gb_end), (lhs_start, lhs_end), (rhs_start, rhs_end))


def get_probe_bc_defn(barcode_def: list[dict[str, str | int]]) -> tuple[int | None, int | None]:
    """Returns the offset and length of probe bc if one exists within a cell barcode.

    Args:
        barcode_def (List[Dict[str, Union[str, int]]]): (gb_end, probe_end)

    Returns:
        Tuple[Optional[int], Optional[int]]: [description]
    """
    (_, gb_end), (lhs_start, lhs_end), (rhs_start, rhs_end) = get_barcode_locations(barcode_def)
    if lhs_start is None and rhs_start is None:
        return (None, None)
    elif lhs_end is not None and rhs_end is None:
        probe_end = lhs_end
    elif lhs_end is None and rhs_end is not None:
        probe_end = rhs_end
    else:
        probe_end = max(lhs_end, rhs_end)
    return (gb_end, probe_end - gb_end)


def get_probe_bc_whitelist(barcode_def: list[dict[str, dict[str, str]]]) -> dict:
    """Returns the whitelist spec of the probe barcodes.

    Args:
        barcode_def: list of barcode definitions from chemistry def

    Returns:
        whitelist spec of the probe barcodes
    """
    return next(
        barcode["whitelist"]
        for barcode in barcode_def
        if barcode["kind"] in ("left_probe", "right_probe", "overhang")
    )
