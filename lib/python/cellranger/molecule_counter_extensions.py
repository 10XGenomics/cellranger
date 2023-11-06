#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""High performance compiled methods for the molecule counter.

In a separate file to avoid loading dependencies during preflights.
"""

from __future__ import annotations

import numpy as np

from cellranger import molecule_counter as cr_mc


def get_indices_for_values(
    mc: cr_mc.MoleculeCounter, col_names, values, chunk_size: int = (1 << 20)
) -> np.ndarray[int, np.dtype[np.bool_]]:
    """Get the molecule indices for those with values (x, y, z) from cols [X, Y, Z].

    Args:
        mc: the input object to work with.
        col_names (list): list of column names values will come from, in some order
        values (list): list of tuples of with values that match columns,
                       in the same order
        chunk_size: the size of each chunk to store.

    Returns:
      np.array_1d(dtype=bool): indices that can be used to lookup any molecules
                               matching a tuple in values
    """
    nmol = mc.nrows()
    idxs = np.zeros((nmol,), dtype=bool)
    if not col_names:
        return idxs
    if isinstance(values, set):
        raise TypeError("Argument cannot be a set, must be list or numpy array")

    # Special case filtering in one dimension, to use numba
    # TODO: Can we use numba on structs? (ran into some compilation issues and had to move on
    # but am not certain it isn't possible).
    if len(col_names) == 1:

        def _get_element_in_set(values):
            col_name = col_names[0]
            single_values_ = np.unique(values)
            dtype = single_values_.dtype
            if col_name == cr_mc.LIBRARY_IDX_COL_NAME:
                max_index = len(mc.get_library_info())
            elif col_name == cr_mc.FEATURE_IDX_COL_NAME:
                max_index = len(mc.feature_reference.feature_defs)
            elif col_name == cr_mc.BARCODE_IDX_COL_NAME:
                max_index = mc.get_barcode_list_size()
            else:
                raise NotImplementedError(
                    "Only barcode, feature and library indexes are currently supported for single indexing."
                )
            element_in_set = np.zeros(max_index, dtype=bool)
            element_in_set[single_values_] = True
            del single_values_
            return element_in_set, mc.get_column_lazy(col_name), dtype

        element_in_set, col, dtype = _get_element_in_set(values)

        def _mask(
            num: int, source_sel: slice, dest_sel: slice
        ) -> np.ndarray[int, np.dtype[np.bool_]]:
            col.read_direct(chunk, source_sel, dest_sel)
            return element_in_set[chunk[0:num]]

    else:
        cols = [(name, mc.get_column_lazy(name)) for name in col_names]
        dtype = np.dtype([(name, col.dtype) for name, col in cols])
        values_ = np.empty((len(values),), dtype=dtype)
        values_[:] = values
        values_ = np.unique(values_)
        # we need a contiguous store for h5py.Dataset.read_direct
        store = {name: np.empty((chunk_size,), dtype=col.dtype) for name, col in cols}

        def _mask(
            num: int, source_sel: slice, dest_sel: slice
        ) -> np.ndarray[int, np.dtype[np.bool_]]:
            for name, col in cols:
                # do not allocate, read directly into the chunk array
                col.read_direct(store[name], source_sel, dest_sel)
            for name in col_names:
                chunk[name][0:num] = store[name][0:num]
            # np.isin should work for us up to most reasonable sizes of values_ (<1M)
            return np.isin(chunk[0:num], values_)

    del values
    # and we need them in struct array layout for np.isin
    chunk = np.empty((chunk_size,), dtype=dtype)
    for begin in range(0, nmol, chunk_size):
        end = min(begin + chunk_size, nmol)
        num = end - begin
        source_sel = np.s_[begin:end]
        dest_sel = np.s_[0:num]
        idxs[begin:end] |= _mask(num, source_sel, dest_sel)
    return idxs
