#!/usr/bin/env python3
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Functions for interacting with h5 tables."""

from __future__ import annotations

import csv
import os.path
from collections.abc import Callable, Generator, Iterable, Mapping
from contextlib import ExitStack
from typing import NamedTuple, TypeVar, overload

import numpy as np

import cellranger.h5_constants as h5_constants
from cellranger.analysis.analysis_types import PCA
from cellranger.wrapped_tables import tables

# Version for HDF5 format
VERSION_KEY = "version"
VERSION = 2


def __encode_h5_arr(f: tables.File, arr, field: str, group: str, key: str, subgroup: tables.Group):
    # XML encode strings so we can store them as HDF5 ASCII
    if isinstance(arr, str):
        arr = np.bytes_(arr.encode("ascii", "xmlcharrefreplace"))
    elif isinstance(arr, bytes):
        arr = np.bytes_(arr)

    if not hasattr(arr, "dtype"):
        raise ValueError(f"{group}/{key} must be a numpy array or scalar")

    atom = tables.Atom.from_dtype(arr.dtype)
    if len(arr.shape) > 0:
        if arr.size > 0:
            dest_array = f.create_carray(subgroup, field, atom, arr.shape)
        else:
            dest_array = f.create_earray(subgroup, field, atom, arr.shape)
        dest_array[:] = arr
    else:
        f.create_array(subgroup, field, arr)


def save_h5(f: tables.File, group: str, key: str, namedtuple: NamedTuple):
    """Save a namedtuple to an h5 file under a group and subgroup."""
    if VERSION_KEY in f.root:
        version = int(getattr(f.root, VERSION_KEY))
        if version != VERSION:
            raise ValueError(
                "Attempted to write analysis HDF5 version %d data to a version %d file"
                % (VERSION, version)
            )
    else:
        f.create_array(f.root, VERSION_KEY, np.int64(VERSION))

    subgroup = f.create_group(group, "_" + key)
    for field in namedtuple._fields:
        __encode_h5_arr(f, getattr(namedtuple, field), field, group, key, subgroup)


def _string_or_num(value: str | bytes | int | float):
    if isinstance(value, str):
        return value
    elif isinstance(value, bytes):
        return value.decode()
    return value


def _string_or_num_list(
    value: str | bytes | int | float | Iterable[str | bytes | int | float],
) -> list[str | int | float]:
    """Convert a string or collection into a list of str or numbers.

    `bytes` are decoded to `str`s.

    `str` or numbers are wrapped in lists.

    For iterables (other than strings), a list is returned where any `bytes`
    elements have been decoded to `str`.

    Args:
        value (str, bytes, int, float, iterable): The value to convert.

    Returns:
        list: A list containing `str` or numeric values.
    """
    if not isinstance(value, bytes | str | int | float) and hasattr(value, "__iter__"):
        return [_string_or_num(v) for v in value]
    return [_string_or_num(value)]


def save_matrix_csv(
    filename: str | bytes,
    arr: Iterable[str | bytes | int | float | Iterable[str | bytes | int | float]],
    header: Iterable[str],
    prefixes: Iterable[str | bytes | int | float | Iterable[str | bytes | int | float]],
):
    """Save a csv file of the matrix.

    Args:
        filename (Union[str, bytes]): The destination filename.
        arr (Iterable[Union[str, bytes, int, float, Iterable[Union[str, bytes, int, float]]]]): _description_
        header (Iterable[str]): _description_
        prefixes (Iterable[Union[str, bytes, int, float, Iterable[Union[str, bytes, int, float]]]]): _description_
    """
    with open(filename, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(header)
        # Iterate over the given arr, default iteration is by-row
        for row_vec, prefix in zip(arr, prefixes):
            row = _string_or_num_list(prefix)
            row.extend(_string_or_num_list(row_vec))
            writer.writerow(row)


_T1 = TypeVar("_T1", bound=NamedTuple)  # pylint: disable=invalid-name


def __get_table_node(group: tables.Group, field: str):
    try:
        field_value = getattr(group, field).read()
        if field_value.shape == ():
            return field_value.item()
        return field_value
    except tables.NoSuchNodeError:
        return getattr(group._v_attrs, field, None)  # pylint: disable=protected-access


def load_h5_namedtuple(group: tables.Group, namedtuple: type[_T1]) -> _T1:
    """Load a single namedtuple from an h5 group."""
    return namedtuple._make([__get_table_node(group, field) for field in namedtuple._fields])


def load_h5_iter(
    group: tables.Group, namedtuple: type[_T1]
) -> Generator[tuple[str, _T1], None, None]:
    """Iterate through the subgroups of a group, converting each to the given type.

    Args:
        group (tables.Group): _description_
        namedtuple (Type[_T1]): _description_

    Yields:
        Generator[Tuple[str, _T1], None, None]: _description_
    """
    for subgroup in group:
        yield subgroup._v_name[1:], load_h5_namedtuple(  # pylint: disable=protected-access
            subgroup, namedtuple
        )


def h5_path(base_path):
    return os.path.join(base_path, "analysis.h5")


def _combine_h5_group(fins: Iterable[tables.File], fout: tables.File, group: str):
    group_out = fout.create_group(fout.root, group)
    for fin in fins:
        # Skip non-existent input groups
        if group not in fin.root:
            continue

        group_in = fin.root._v_groups[group]  # pylint: disable=protected-access

        # NOTE - this throws an exception if a child group already exists
        fin.copy_children(group_in, group_out, recursive=True)


def combine_h5_files(in_files: Iterable[str], out_file: str, groups: Iterable[str]):
    """Merge a set of groups from a set of h5 files.

    Args:
        in_files (Iterable[str]): The source files.
        out_file (str): The destination file.
        groups (Iterable[str]): The names of the groups.
    """
    with ExitStack() as fin_stack:
        fins: Iterable[tables.File] = [
            fin_stack.enter_context(tables.open_file(filename, "r")) for filename in in_files
        ]
        with tables.open_file(out_file, "w") as fout:
            for group in groups:
                _combine_h5_group(fins, fout, group)


def open_h5_for_writing(filename: str) -> tables.File:
    filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
    return tables.open_file(filename, "w", filters=filters)


def _save_dimension_reduction_h5(
    data_map: Mapping[int, _T1],
    fname: str,
    group_name: str,
    key: Callable[[int], str],
):
    filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(fname, "w", filters=filters) as f:
        group = f.create_group(f.root, group_name)
        for n_components, reduction in data_map.items():
            save_h5(f, group, key(n_components), reduction)


def save_dimension_reduction_h5(data_map: Mapping[int, _T1], fname: str, group_name: str):
    _save_dimension_reduction_h5(data_map, fname, group_name, str)


def save_pca2_dimension_reduction_h5(
    data_map: Mapping[int, _T1],
    fname: str,
    group_name: str,
    library_type: str,
):
    """Save PCA dimension reduction info to HDF5."""
    _save_dimension_reduction_h5(
        data_map, fname, group_name, lambda n_components: f"{library_type}_{n_components}"
    )


@overload
def load_dimension_reduction_from_h5(
    filename: str, group_name: str, ntuple: type[PCA]
) -> list[PCA]: ...


@overload
def load_dimension_reduction_from_h5(filename: str, group_name: str, ntuple: type[_T1]) -> _T1: ...


def load_dimension_reduction_from_h5(filename: str, group_name: str, ntuple: type[_T1]):
    """Iterate over the output of dim reduction and return a list of all projections if PCA."""
    with tables.open_file(filename, "r") as f:
        group: tables.Group = f.root._v_groups[group_name]  # pylint: disable=protected-access
        if ntuple != PCA:  # "lsa", "plsa", etc
            # Just take the first object, assuming we never have multiple
            for _, pca in load_h5_iter(group, ntuple):
                return pca

        # otherwise for PCA, we potentially generate both GEX-, Ab-, and Peak- based projections
        gex_pca, ab_pca, atac_pca = None, None, None
        for pca_type, pca in load_h5_iter(group, ntuple):
            # if pca_type is antibody_capture_10
            if pca_type.startswith("antibody_capture"):
                ab_pca = pca
            # if pca_type is gene_expression_10 or 10
            elif pca_type.startswith("gene_expression") or pca_type.isdigit():
                gex_pca = pca
            elif pca_type.startswith("peaks"):
                atac_pca = pca
        if gex_pca and ab_pca:
            pca_components = [gex_pca, ab_pca]
        elif gex_pca:
            pca_components = [gex_pca]
        elif ab_pca:
            pca_components = [ab_pca]
        elif atac_pca:
            pca_components = [atac_pca]
        else:
            raise ValueError("Expected at least one of gene_expression, antibody_capture, or peaks")
        return pca_components
