#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import csv
import numpy as np
import os
import six
import tables
import cellranger.h5_constants as h5_constants

# Version for HDF5 format
VERSION_KEY = 'version'
VERSION = 2

def save_h5(f, group, key, namedtuple):
    """ Save a namedtuple to an h5 file under a group and subgroup """
    if VERSION_KEY in f.root:
        version = int(getattr(f.root, VERSION_KEY))
        if version != VERSION:
            raise ValueError("Attempted to write analysis HDF5 version %d data to a version %d file" % (VERSION, version))
    else:
        ds = f.create_array(f.root, VERSION_KEY, np.int64(VERSION))

    subgroup = f.create_group(group, '_'+key)
    for field in namedtuple._fields:
        arr = getattr(namedtuple, field)

        # XML encode strings so we can store them as HDF5 ASCII
        if isinstance(arr, six.string_types):
            arr = np.string_(arr.encode('ascii', 'xmlcharrefreplace'))

        if not hasattr(arr, 'dtype'):
            raise ValueError('%s/%s must be a numpy array or scalar' % (group,key))

        atom = tables.Atom.from_dtype(arr.dtype)
        if len(arr.shape) > 0:
            if arr.size > 0:
                ds = f.create_carray(subgroup, field, atom, arr.shape)
            else:
                ds = f.create_earray(subgroup, field, atom, arr.shape)
            ds[:] = arr
        else:
            ds = f.create_array(subgroup, field, arr)

def save_matrix_csv(filename, arr, header, prefixes):
    with open(filename, 'wb') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(header)
        # Iterate over the given arr, default iteration is by-row
        for i, v in enumerate(arr):
            row = list(prefixes[i]) if hasattr(prefixes[i], '__iter__') else [prefixes[i]]
            row += list(v) if hasattr(v, '__iter__') else [v]
            writer.writerow(row)

def load_h5_namedtuple(group, namedtuple):
    """ Load a single namedtuple from an h5 group """
    args = []
    for field in namedtuple._fields:
        try:
            field_value = getattr(group, field).read()
            if field_value.shape == ():
                field_value = np.asscalar(field_value)
        except tables.NoSuchNodeError:
            try:
                field_value = getattr(group._v_attrs, field)
            except AttributeError:
                field_value = None
        args.append(field_value)
    return namedtuple(*args)

def load_h5_iter(group, namedtuple):
    for subgroup in group:
        yield subgroup._v_name[1:], load_h5_namedtuple(subgroup, namedtuple)

def h5_path(base_path):
    return os.path.join(base_path, "analysis.h5")

def _combine_h5_group(fins, fout, group):
    group_out = fout.create_group(fout.root, group)
    for fin in fins:
        # Skip non-existent input groups
        if group not in fin.root:
            continue

        group_in = fin.root._v_groups[group]

        # NOTE - this throws an exception if a child group already exists
        fin.copy_children(group_in, group_out, recursive=True)

def combine_h5_files(in_files, out_file, groups):
    fins = [tables.open_file(filename, 'r') for filename in in_files]
    with tables.open_file(out_file, 'w') as fout:
        for group in groups:
            _combine_h5_group(fins, fout, group)
    for fin in fins:
        fin.close()

def open_h5_for_writing(filename):
    filters = tables.Filters(complevel = h5_constants.H5_COMPRESSION_LEVEL)
    return tables.open_file(filename, 'w', filters = filters)
