#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

"""Functions for fast, memory-efficient operations on masked sparse matrices."""


from __future__ import annotations

from typing import overload

import numpy as np


@overload
def sum_masked(
    matrix,
    row_mask: np.ndarray[int, np.dtype[np.bool_]],
    col_mask: np.ndarray[int, np.dtype[np.bool_]],
    axis: None,
) -> np.generic: ...


@overload
def sum_masked(
    matrix,
    row_mask: np.ndarray[int, np.dtype[np.bool_]],
    col_mask: np.ndarray[int, np.dtype[np.bool_]],
    axis: int,
) -> np.ndarray[int, np.dtype[np.generic]]: ...


def sum_masked(
    matrix,
    row_mask: np.ndarray[int, np.dtype[np.bool_]],
    col_mask: np.ndarray[int, np.dtype[np.bool_]],
    axis: int | None,
):
    """Sum a masked sparse matrix along an axis.

    Args:
      matrix (scipy.sparse.csc_matrix): Matrix
      axis (int): Axis to sum away; None for total sum.
      row_mask (np.array): Boolean mask for rows
      col_mask (np.array): Boolean mask for columns

    Returns:
      np.array
    """
    assert matrix.getformat() == "csc"
    assert len(row_mask) == matrix.shape[0]
    assert len(col_mask) == matrix.shape[1]
    assert row_mask.dtype == "bool"
    assert col_mask.dtype == "bool"

    if axis == 0:
        return matrix.transpose(copy=False).dot(row_mask)[col_mask]

    elif axis == 1 or axis is None:
        row_sums = matrix.dot(col_mask)[row_mask]

        if axis == 1:
            return row_sums
        assert axis is None
        return np.sum(row_sums)

    else:
        raise ValueError("axis must be 0, 1, or None")


@overload
def count_ge_masked(
    m,
    row_mask: np.ndarray[int, np.dtype[np.bool_]],
    col_mask: np.ndarray[int, np.dtype[np.bool_]],
    threshold: int,
    axis: None,
    chunk_size_mb: int = 64,
) -> np.uint64: ...


@overload
def count_ge_masked(
    m,
    row_mask: np.ndarray[int, np.dtype[np.bool_]],
    col_mask: np.ndarray[int, np.dtype[np.bool_]],
    threshold: int,
    axis: int,
    chunk_size_mb: int = 64,
) -> np.ndarray[int, np.dtype[np.uint64]]: ...


def count_ge_masked(
    m,
    row_mask: np.ndarray[int, np.dtype[np.bool_]],
    col_mask: np.ndarray[int, np.dtype[np.bool_]],
    threshold: int,
    axis: int | None,
    chunk_size_mb: int = 64,
):
    """Count values greater than or equal to a threshold in a masked sparse matrix.

    Args:
      m (scipy.sparse.csc_matrix): Matrix
      axis (int): Axis to sum away; None for total count.
      threshold (int): Only count elements exceeding this threshold
      row_mask (np.array): Boolean mask for rows
      col_mask (np.array): Boolean mask for columns
      chunk_size_mb (int): Operate on chunks of this size; don't copy the entire matrix at once.

    Returns:
      np.array
    """
    assert chunk_size_mb > 0
    assert m.getformat() == "csc"
    assert len(row_mask) == m.shape[0]
    assert len(col_mask) == m.shape[1]
    assert row_mask.dtype == "bool"
    assert col_mask.dtype == "bool"

    # Operate on chunks of bounded size
    ind_bytes = m.indices.dtype.itemsize
    data_bytes = m.data.dtype.itemsize
    bytes_per_nz_elem = ind_bytes + data_bytes

    nz_elem_per_col = np.diff(m.indptr)
    bytes_per_col = nz_elem_per_col * bytes_per_nz_elem
    chunk_size_bytes = chunk_size_mb * (1024 * 1024)

    # Split the matrix into approxiamtely equal-RAM barcode chunks
    chunk_starts = [0]
    chunk_size = 0
    for i in range(len(bytes_per_col)):
        chunk_size += bytes_per_col[i]
        if chunk_size >= chunk_size_bytes:
            chunk_starts.append(i)
            chunk_size = 0

    results = []

    if axis not in (0, 1, None):
        raise ValueError("Axis must be 0, 1, or None")

    for chunk_idx, chunk_start in enumerate(chunk_starts):
        chunk_end = (
            chunk_starts[1 + chunk_idx] if chunk_idx < (len(chunk_starts) - 1) else m.shape[1]
        )
        chunk = slice(chunk_start, chunk_end)
        submatrix = (m[:, chunk] >= threshold).astype(np.uint64)
        if axis == 0:
            s = submatrix.transpose(copy=False).dot(row_mask)
        else:
            s = submatrix.dot(col_mask[chunk])
        results.append(s)

    if axis == 0:
        return np.concatenate(results)[col_mask]
    elif axis == 1:
        return sum(results)[row_mask]
    elif axis is None:
        return np.sum(sum(results)[row_mask])
    else:
        raise ValueError("Axis must be 0, 1, or None")
