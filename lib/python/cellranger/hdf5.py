#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#


from __future__ import annotations

import functools
from collections.abc import Callable, Iterable
from html import unescape

import h5py
import numpy as np
from six import ensure_str

import cellranger.h5_constants as h5_constants
from cellranger.wrapped_tables import tables

STR_DTYPE_CHAR = np.dtype(np.bytes_).char
UNICODE_DTYPE_CHAR = np.dtype(np.str_).char
STRING_DTYPE_CHARS = (STR_DTYPE_CHAR, UNICODE_DTYPE_CHAR)


def is_hdf5(filename, throw_exception=False):
    """A wrapper around h5py.is_hdf5, optionally can throw an exception.

    Args:
        filename: The name of the file to test
        throw_exception: Should we raise an error if not?

    Returns:
        bool as to whether the file is valid
    """
    valid = h5py.is_hdf5(filename)
    if not valid and throw_exception:
        raise OSError(f"File: {filename} is not a valid HDF5 file.")
    return valid


def write_h5(filename, data, append=False):
    filemode = "w"
    if append:
        filemode = "a"
    with h5py.File(filename, filemode) as f:
        for key, value in data.items():
            f[key] = value


def get_h5_filetype(filename):
    with tables.open_file(filename, mode="r") as f:
        try:
            filetype = ensure_str(f.get_node_attr("/", h5_constants.H5_FILETYPE_KEY))
        except AttributeError:
            filetype = None  # older files lack this key
    return filetype


def save_array_h5(filename, name, arr):
    """Save an array to the root of an h5 file."""
    with tables.open_file(filename, "w") as f:
        f.create_carray(f.root, name, obj=arr)


def load_array_h5(filename, name):
    """Load an array from the root of an h5 file."""
    with tables.open_file(filename, "r") as f:
        return getattr(f.root, name).read()


def create_hdf5_string_dataset(group, name, data, **kwargs):
    """Create a dataset of strings under an HDF5 (h5py) group.

    Strings are stored as fixed-length 7-bit ASCII with XML-encoding
    for characters outside of 7-bit ASCII. This is inspired by the
    choice made for the Loom spec:
    https://github.com/linnarsson-lab/loompy/blob/master/doc/format/index.rst

    Args:
        group (h5py.Node): Parent group.
        name (str): Dataset name.
        data (list of str): Data to store. Both None and [] are serialized to an empty dataset.
                            Both elements that are empty strings and elements that are None are
                            serialized to empty strings.
        **kwargs: Additional arguments to `create_dataset`.
    """
    if data is None or hasattr(data, "__len__") and len(data) == 0:
        dtype = np.dtype((np.bytes_, 1))
        group.create_dataset(name, dtype=dtype)
        return

    assert (isinstance(data, np.ndarray) and data.dtype.char == STR_DTYPE_CHAR) or (
        isinstance(data, list) and all(x is None or isinstance(x, bytes | str) for x in data)
    )

    # Convert Nones to empty strings and use XML encoding
    data = [encode_ascii_xml(x) if x else b"" for x in data]

    fixed_len = max(len(x) for x in data)

    # h5py doesn't support strings with zero-length-dtype
    if fixed_len == 0:
        fixed_len = 1
    dtype = np.dtype((np.bytes_, fixed_len))

    group.create_dataset(name, data=data, dtype=dtype, **kwargs)


def decode_ascii_xml(x: str | bytes) -> str:
    """Decode a string from 7-bit ASCII + XML into unicode."""
    if isinstance(x, str):
        return x
    elif isinstance(x, bytes):
        return unescape(x.decode())
    else:
        raise ValueError(f"Expected string type, got type {type(x)!s}")


def decode_ascii_xml_array(data):
    """Decode an array-like container of strings from 7-bit ASCII + XML.

    into unicode.
    """
    if isinstance(data, np.ndarray) and data.dtype.char == UNICODE_DTYPE_CHAR:
        return data

    unicode_data = [decode_ascii_xml(x) for x in data]

    fixed_len = max(len(s) for s in unicode_data)
    # 0 length string type is Ok
    dtype = np.dtype((np.str_, fixed_len))

    # note: python3 would require no.fromiter
    return np.array(unicode_data, dtype=dtype)


def encode_ascii_xml(x: str | bytes):
    """Encode a string as fixed-length 7-bit ASCII with XML-encoding.

    for characters outside of 7-bit ASCII.

    Respect python2 and python3, either unicode or binary.
    """
    if isinstance(x, str):
        return x.encode("ascii", "xmlcharrefreplace")
    elif isinstance(x, bytes):
        return x
    else:
        raise ValueError(f"Expected string type, got type {type(x)!s}")


def encode_ascii_xml_array(
    data: np.ndarray | Iterable[str | bytes],
) -> np.ndarray[tuple[int, int], np.dtype[np.bytes_]]:
    """Encode an array-like container of strings as fixed-length 7-bit ASCII.

    with XML-encoding for characters outside of 7-bit ASCII.
    """
    if (
        isinstance(data, np.ndarray)
        and data.dtype.char == STR_DTYPE_CHAR
        and data.dtype.itemsize > 0
    ):
        return data

    def convert(s):
        return encode_ascii_xml(s) if s is not None else b""

    ascii_data = [convert(x) for x in data]

    fixed_len = max(len(s) for s in ascii_data)
    fixed_len = max(1, fixed_len)
    dtype = np.dtype((np.bytes_, fixed_len))

    # note: python3 would require np.fromiter
    return np.array(ascii_data, dtype=dtype)


_cached_decode: Callable[[bytes | str], str] = functools.lru_cache(maxsize=4)(decode_ascii_xml)


def read_hdf5_string_dataset(dataset: h5py.Dataset, memoize: bool = False) -> list[str]:
    """Read a dataset of strings from HDF5 (h5py).

    Args:
        dataset (h5py.Dataset): Data to read.
        memoize (bool): Whether to use memoization to make string conversion more efficient.

    Returns:
        list[unicode]: Strings in the dataset.
    """
    # h5py doesn't support loading an empty dataset
    if dataset.shape is None:
        return []

    data = dataset[:]
    decode = decode_ascii_xml
    if memoize:
        decode = _cached_decode
    return [decode(x) for x in data]


def set_hdf5_attr(dataset, name: str | bytes, value: str | bytes | Iterable[str | bytes]):
    """Set an attribute of an HDF5 dataset/group.

    Strings are stored as fixed-length 7-bit ASCII with XML-encoding
    for characters outside of 7-bit ASCII. This is inspired by the
    choice made for the Loom spec:
    https://github.com/linnarsson-lab/loompy/blob/master/doc/format/index.rst
    """
    name = encode_ascii_xml(name)

    if isinstance(value, str | bytes):
        value = encode_ascii_xml(value)
    elif isinstance(value, np.ndarray) and value.dtype.char in STRING_DTYPE_CHARS:
        value = encode_ascii_xml_array(value)
    elif isinstance(value, list | tuple) and isinstance(value[0], str | bytes):
        value = encode_ascii_xml_array(value)

    dataset.attrs[name] = value


def combine_h5s_into_one(outfile: str | bytes, files_to_combine: Iterable[str | bytes]):
    assert isinstance(files_to_combine, list)
    with tables.open_file(ensure_str(outfile), "a") as out:
        for fname in files_to_combine:
            with tables.open_file(ensure_str(fname), "r") as data:
                data.copy_children(data.root, out.root, recursive=True)
