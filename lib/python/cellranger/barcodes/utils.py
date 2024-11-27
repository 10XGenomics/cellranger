#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Barcode handling functions."""

from __future__ import annotations

# Do not add new things to this module.
# Instead, either find or create a module with a name that better describes
# the functionality implemented by the methods or classes you want to add.
import os
from typing import Literal, overload

import numpy as np

import cellranger.constants as cr_constants
import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants

SPATIAL_WHITELIST = [
    "odin-5K-v2",
    "visium-v1",
    "visium-v2",
    "visium-v3",
    "visium-v4",
    "visium-v5",
    "omni-v1",
    "omni-hires-v1",
    "omni-pat7-v1",
    "pseudo-xl-v1",
    "thor-v1",
    "thor-v2",
    "thor-v3",
    "thor-XL-v1",
    "thor-XL-v2",
]
NUM_BARCODES_PER_MEM_GB = 500000


def is_whitelist_spatial(whitelist_name: str):
    """Is the barcode a spatial (visium) whitelist?

    Method returns True if the name is present on a predefined list
    of spatial barcode whitelist names.

    Args:
        whitelist_name: string of the whitelist name

    Returns:
        boolean value indicating if the whitelist is from the spatial assay.
    """
    return whitelist_name in SPATIAL_WHITELIST


@overload
def load_barcode_tsv(filename: str | bytes, as_set: Literal[False] = False) -> list[bytes]: ...


@overload
def load_barcode_tsv(filename: str | bytes, as_set: Literal[True] = ...) -> set[bytes]: ...


@overload
def load_barcode_tsv(filename: str | bytes, as_set: bool = False) -> set[bytes] | list[bytes]: ...


def load_barcode_tsv(filename: str | bytes, as_set=False) -> set[bytes] | list[bytes]:
    """Load barcodes from a tsv file."""
    with cr_io.open_maybe_gzip(filename, "rb") as bc_file:
        barcodes = [x.strip() for x in bc_file if b"#" not in x]
    barcode_set = set(barcodes)
    if len(barcodes) != len(barcode_set):
        raise Exception(f"Duplicates found in barcode whitelist: {filename}")
    return barcode_set if as_set else barcodes


@overload
def get_barcode_whitelist_path(filename: None, translation: bool = False) -> None: ...


@overload
def get_barcode_whitelist_path(filename: str, translation: bool = False) -> str: ...


def get_barcode_whitelist_path(filename: str | None, translation: bool = False) -> str | None:
    # Look for exact path, .txt.gz, or .txt
    if filename is None:
        return None
    elif os.path.exists(filename):
        return filename
    else:
        rootdir = cr_constants.BARCODE_WHITELIST_PATH
        if translation:
            rootdir = os.path.join(rootdir, "translation")

        gz = os.path.join(rootdir, filename + ".txt.gz")
        if os.path.exists(gz):
            return gz

        txt = os.path.join(rootdir, filename + ".txt")
        return txt


def get_all_whitelist_filenames() -> list[str]:
    """Function to return the set of all whitelist names included in the product.

    Uses some heuristics to search through folders in this directory.

    Returns:
        A list of strings, one per filename
    """
    basenames = os.listdir(cr_constants.BARCODE_WHITELIST_PATH)
    allowed_suffixes = (".txt", ".txt.gz")
    all_files = []
    for name in basenames:
        full_name = os.path.join(cr_constants.BARCODE_WHITELIST_PATH, name)
        if (
            name.endswith(allowed_suffixes)
            and name.count("coordinates") == 0
            and os.path.isfile(full_name)
        ):
            all_files.append(full_name)
    return all_files


def get_barcode_whitelist_paths(filenames: str):
    paths = filenames.split(",")
    return ",".join([get_barcode_whitelist_path(p) for p in paths])


class WhitelistNotFound(Exception):
    """Raised when a named whitelist is not found."""


@overload
def load_barcode_whitelist(filename: None, as_set=False) -> None: ...


@overload
def load_barcode_whitelist(
    filename: str, as_set: Literal[False] = False, translation: bool = False
) -> list[bytes]: ...


@overload
def load_barcode_whitelist(
    filename: str, as_set: Literal[True], translation: bool = False
) -> set[bytes]: ...


@overload
def load_barcode_whitelist(
    filename: str, as_set: bool, translation: bool
) -> set[bytes] | list[bytes]: ...


def load_barcode_whitelist(filename: str | None, as_set: bool = False, translation: bool = False):
    path = get_barcode_whitelist_path(filename, translation)

    if path is None:
        return None
    if not os.path.isfile(path):
        raise WhitelistNotFound(f"Unable to find barcode whitelist: {path}")

    return load_barcode_tsv(path, as_set)


def load_probe_barcode_map(
    name: str | None = None, path: str | None = None
) -> dict[str, str] | None:
    """Map probe BC ids to their collaped sequences.

    If the barcode whitelist exists return the mapping dictionary,
    else, return None.
    """
    if name is None and path is None:
        return None

    assert not (name is not None and path is not None)

    if name is not None:
        file_path = None
        for extension in [".txt", ".txt.gz"]:
            file_ext = os.path.join(cr_constants.BARCODE_WHITELIST_TRANSLATE_PATH, name + extension)
            if os.path.exists(file_ext):
                file_path = file_ext
                break
    else:
        file_path = path

    if file_path is None:
        return None
    else:
        translate_map = {}
        for line in cr_io.open_maybe_gzip(file_path, "r"):
            if line.startswith("#"):
                continue
            bcs = line.strip().split()
            translate_map[bcs[2]] = bcs[1]
        return translate_map


def get_mem_gb_request_from_barcode_whitelist(
    barcode_whitelist_fn: str | None, gem_groups=None, use_min=True, double=False
):
    barcode_whitelist = load_barcode_whitelist(barcode_whitelist_fn)

    if use_min:
        if barcode_whitelist is None:
            min_mem_gb = h5_constants.MIN_MEM_GB_NOWHITELIST
        else:
            min_mem_gb = h5_constants.MIN_MEM_GB
    else:
        min_mem_gb = 0

    if barcode_whitelist is None:
        return min_mem_gb

    if gem_groups is not None:
        num_bcs = len(barcode_whitelist) * max(gem_groups)
    else:
        num_bcs = len(barcode_whitelist)

    if double:
        return np.ceil(max(min_mem_gb, 2 * num_bcs // NUM_BARCODES_PER_MEM_GB))
    else:
        return np.ceil(max(min_mem_gb, num_bcs // NUM_BARCODES_PER_MEM_GB))
