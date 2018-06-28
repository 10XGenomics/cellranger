#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import tables

UMI_INFO_COLS = {
    'barcode_idx': np.uint32,
    'chain_idx': np.uint32,
    'umi_idx': np.uint32,
    'reads': np.uint32,
}

UMI_REF_COLS = {
    'barcodes': np.str_,
    'chains': np.str_,
}

def get_num_rows(filename):
    with tables.open_file(filename, 'r') as h5:
        return getattr(h5.root, UMI_INFO_COLS.keys()[0]).shape[0]

def get_mem_gb(filename, start_row=None, end_row=None):
    assert (start_row is None) == (end_row is None)

    num_rows = get_num_rows(filename) if start_row is None else end_row - start_row

    return sum(np.dtype(nptype).itemsize for nptype in UMI_INFO_COLS.itervalues()) * num_rows / 1e9

def get_dtype(col):
    if col in UMI_INFO_COLS:
        return np.dtype(UMI_INFO_COLS[col])
    elif col in UMI_REF_COLS:
        return np.dtype(UMI_REF_COLS[col])
    else:
        raise ValueError('Unrecognized column name: %s' % col)

def create_arrays(h5):
    d = {}
    for name, nptype in UMI_INFO_COLS.items():
        d[name] = h5.create_earray('/', name, atom=tables.Atom.from_dtype(np.dtype(nptype)), shape=(0,))
    return d

def set_ref_column(h5, name, data):
    return h5.create_array('/', name, data)

def read_umi_info(filename, start_row=None, end_row=None):
    """ Load data into memory and return a dict """
    assert (start_row is None) == (end_row is None)

    d = {}
    with tables.open_file(filename, 'r') as h5:
        for node in h5.walk_nodes('/', 'Array'):
            if node.name in UMI_INFO_COLS.keys():
                d[node.name] = node[slice(start_row, end_row)]
            elif node.name in UMI_REF_COLS.keys():
                d[node.name] = node[:]

    return d

def get_column(filename, col):
    """ Load a column into memory """
    with tables.open_file(filename, 'r') as h5:
        return getattr(h5.root, col)[:]
