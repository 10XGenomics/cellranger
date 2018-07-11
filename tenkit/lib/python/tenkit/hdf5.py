#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import os.path
import h5py
import numpy as np
import tenkit.pandas as p
import itertools

import pprint
pp = pprint.PrettyPrinter()

CHUNK_SIZE=2**16
vlen_string = h5py.special_dtype(vlen=str)
COMPRESSION = 3
LEVEL_GROUP = "_levels"

def make_index_array(uniques, values, dtype=np.uint16):
    idx_dict = {v:idx for (idx, v) in enumerate(uniques)}

    if type(values) == p.Categorical:
        return make_index_array_cat(idx_dict, values, dtype)
    else:
        return make_index_array_dict(idx_dict, values, dtype)

def make_index_array_dict(idx_dict, values, dtype=np.uint16):

    result = np.zeros(shape=values.shape, dtype=dtype)

    for item in range(len(values)):
        result[item] = idx_dict[values[item]]

    return result

def make_index_array_cat(idx_dict, categorical, dtype=np.uint16):
    # Map codes from categorical_series over to the requested codes in idx_dict

    num_cat_codes = len(categorical.categories)
    trans_table = -1 * np.ones(num_cat_codes)
    for (code, cat) in enumerate(categorical.categories):
        trans_table[code] = idx_dict[categorical.categories[code]]

    new_codes = trans_table[categorical.codes]
    return new_codes

def _pandas_categorical_type():
    """Find Pandas' CategoricalDtype.

    In older pandas, this was pandas.core.common.CategoricalDtype.

    In newer pandas, this moved to pandas.core.dtypes.dtypes.CategoricalDtype
    with a forwarding alias for the constructor which still does not work
    correctly for type comparisons.
    """
    try:
        return p.core.dtypes.dtypes.CategoricalDtype
    except AttributeError:
        return p.core.common.CategoricalDtype

_PANDAS_CATEGORICAL_TYPE = _pandas_categorical_type()

def write_data_column(grp, column):

    # If the dataset already exists, delete it.
    # Overwriting datasets appears to cause some problems
    if column.name in grp.keys():
        del grp[column.name]

    if type(column.dtype) == _PANDAS_CATEGORICAL_TYPE:
        # Pre-made pandas categorical
        uniq = column.cat.categories
        index_array = column.cat.codes

        ds = grp.create_dataset(column.name,
                data = index_array,
                shape = column.shape,
                maxshape = (None,),
                dtype = index_array.dtype,
                compression = COMPRESSION,
                shuffle = True,
                chunks = (CHUNK_SIZE,))

        create_levels(ds, uniq)

    elif column.dtype == np.dtype('O'):
        uniq = column.unique()

        # For repeated strings, encode as 2-byte integer 'factor'
        if len(uniq) < 0.1 * len(column):
            dt = np.uint16

            # Be aggressive about going to 32bits
            if len(uniq) > 255:
                dt = np.uint32

            # Need to pass column.values -- columns is a pandas.Series, which doesn't work for indexing
            uniq = column.unique()
            index_array = make_index_array(uniq, column.values, dt)

            ds = grp.create_dataset(column.name,
                    data = index_array,
                    shape = column.shape,
                    maxshape = (None,),
                    dtype = dt,
                    compression = COMPRESSION,
                    shuffle = True,
                    chunks = (CHUNK_SIZE,))

            # If a dataset has a 'levels' attribute, it is interpreted as an encoded factor

            create_levels(ds, uniq)


        # Otherwise, convert this data to a string array for storage
        else:
            dt = h5py.special_dtype(vlen=str)
            grp.create_dataset(column.name,
                               shape = column.shape,
                               maxshape = (None,),
                               dtype = dt,
                               data  = column,
                               compression = COMPRESSION,
                               shuffle = True,
                               chunks = (CHUNK_SIZE,))

    else:
        # Store as native numpy / h5py types
        try:
            grp.create_dataset(column.name,
                               shape=column.shape,
                               maxshape=(None,),
                               dtype=column.dtype,
                               data=column,
                               compression=COMPRESSION,
                               shuffle=True,
                               chunks=(CHUNK_SIZE,))
        except TypeError:
            print column.name
            print column.dtype
            print type(column.dtype)
            raise


def widen_cat_column(old_ds, new_type):
    name = old_ds.name
    tmp_name = "__tmp_" + old_ds.name
    grp = old_ds.parent

    ds = grp.create_dataset(tmp_name,
            data = old_ds[:],
            shape = old_ds.shape,
            maxshape = (None,),
            dtype = new_type,
            compression = COMPRESSION,
            shuffle = True,
            chunks = (CHUNK_SIZE,))

    del grp[name]
    grp.move(tmp_name, name)
    return ds

def pick_cat_dtype(sz):

    if sz > 2**16 - 1:
        return np.uint32
    if sz > 255:
        return np.uint16


    return np.uint8


def append_data_column(ds, column):

    # Extend the dataset to fit the new data
    new_count = column.shape[0]
    existing_count = ds.shape[0]
    ds.resize((existing_count + new_count,))

    levels = get_levels(ds)

    if levels is not None:
        # update levels if we have new unique values
        if type(column.values) == p.Categorical:
            added_levels = set(column.values.categories) - set(levels)
        elif len(column) == 0:
            # Workaround for bug in pandas - get a crash in .unique() for an empty series
            added_levels = set([])
        else:
            added_levels = set(column.unique()) - set(levels)

        new_levels = list(levels)
        new_levels.extend(added_levels)

        # Check if the new categorical column has more levels
        # than the current bit width supports.
        # If so, rewrite the existing column data w/ more bits
        if len(new_levels) > np.iinfo(ds.dtype).max:
            new_dtype = pick_cat_dtype(len(new_levels))
            ds = widen_cat_column(ds, new_dtype)

        new_levels = np.array(new_levels, dtype=np.object)
        new_data = make_index_array(new_levels, column.values, ds.dtype)

        clear_levels(ds)
        create_levels(ds, new_levels)
    else:
        new_data = column

    # Append new data
    ds[existing_count:(existing_count + new_count)] = new_data


def write_data_frame(fn, df):
    ''' Write the pandas dataframe object to an HDF5 file.  Each column is written as a single 1D dataset at the top
    level of the HDF5 file, using the native pandas datatype'''

    # Always write a fresh file -- the 'w' argument to h5py.File is supposed to truncate an existing file, but it doesn't appear to work correctly
    if os.path.exists(fn):
        os.remove(fn)

    f = h5py.File(fn, "w")

    # To preserve column order, write columns to an attribute
    column_names = np.array(list(df.columns))
    f.attrs.create("column_names", column_names)

    for col in df.columns:
        write_data_column(f, df[col])

    f.close()

def append_data_frame(fn, df):
    ''' Write the pandas dataframe object to an HDF5 file.  Each column is written as a single 1D dataset at the top
    level of the HDF5 file, using the native pandas datatype'''

    if not os.path.exists(fn):
        write_data_frame(fn, df)
        return

    f = h5py.File(fn, "a")

    column_names = f.attrs.get("column_names")

    for col_name in column_names:
        ds = f[col_name]
        col = df[col_name]
        append_data_column(ds, col)

    f.close()


def has_levels(ds):
    ''' Determine if a data column is leveled '''
    if 'levels' in ds.attrs.keys():
        return True

    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp:
        ds_name = ds.name.split("/")[-1]
        level_ds = level_grp.get(ds_name)
        if level_ds:
            return True

    return False

def get_levels(ds):
    ''' Get the level index for a dataset '''
    if 'levels' in ds.attrs.keys():
        return ds.attrs['levels'][:]

    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp:
        ds_name = ds.name.split("/")[-1]
        level_ds = level_grp.get(ds_name)
        if level_ds:
            return level_ds[:]

    return None

def get_col_type(ds):
    ''' Get logical type of column '''
    if 'levels' in ds.attrs.keys():
        return ds.attrs['levels'].dtype

    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp:
        ds_name = ds.name.split("/")[-1]
        level_ds = level_grp.get(ds_name)
        if level_ds:
            return level_ds.dtype

    return ds.dtype

def clear_levels(ds):
    ''' Remove any levels data for a dataset '''
    if 'levels' in ds.attrs.keys():
        del ds.attrs['levels']

    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp:
        ds_name = ds.name.split("/")[-1]
        level_ds = level_grp.get(ds_name)
        if level_ds:
            del level_grp[ds_name]


def create_levels(ds, levels):
    # Create a dataset in the LEVEL_GROUP
    # and store as native numpy / h5py types
    level_grp = ds.file.get(LEVEL_GROUP)
    if level_grp is None:
        # Create a LEVEL_GROUP
        level_grp = ds.file.create_group(LEVEL_GROUP)
    ds_name = ds.name.split("/")[-1]
    dt = h5py.special_dtype(vlen=str)

    level_grp.create_dataset(ds_name,
                       shape = [len(levels)],
                       maxshape = (None,),
                       dtype = dt,
                       data  = levels,
                       compression = COMPRESSION,
                       chunks = (CHUNK_SIZE,))

def get_column_intersection(target_cols, query_cols):
    if len(query_cols) > 0:
        target_cols = sorted(list(set(query_cols) & set(target_cols)))
        if len(target_cols) == 0:
            raise Exception('No valid column specifications.')
    return target_cols

def read_data_frame(fn, query_cols=[]):
    ''' Load a pandas DataFrame from an HDF5 file. If a column list is specified, only load the matching columns '''

    with h5py.File(fn, 'r') as f:

        column_names = f.attrs.get("column_names")
        column_names = get_column_intersection(column_names, query_cols)

        df = p.DataFrame()

        # Add the columns progressively to save memory
        for name in column_names:
            ds = f[name]
            if has_levels(ds):
                indices = ds[:]
                uniques = get_levels(ds)
                # This method of constructing of Categorical avoids copying the indices array
                # which saves memory for big datasets
                df[name] = p.Categorical(indices, categories=uniques, ordered=False, fastpath=True)
            else:
                df[name] = p.Series(ds[:])

        return df

def read_data_frame_limited(fn, query_cols=[], max_rows=None):
    ''' Load a pandas DataFrame from an HDF5 file. If a column list is specified, only load the matching columns '''

    with h5py.File(fn, 'r') as f:

        column_names = f.attrs.get("column_names")
        column_names = get_column_intersection(column_names, query_cols)

        sz = f[column_names[0]].shape[0]
        if max_rows:
            sz = min(sz, max_rows)

        df = p.DataFrame()

        # Add the columns progressively to save memory
        for name in column_names:
            ds = f[name]
            if has_levels(ds):
                indices = ds[:sz]
                uniques = get_levels(ds)
                # This method of constructing of Categorical avoids copying the indices array
                # which saves memory for big datasets
                df[name] = p.Categorical(indices, categories=uniques, ordered=False, fastpath=True)
            else:
                df[name] = p.Series(ds[:sz])

        return df



def read_data_frame_filtered(fn, filter_func, query_cols=[], chunk_size=5000000):
    ''' Load a pandas DataFrame from an HDF5 file. If a column list is specified, only load the matching columns.
        filter_func should take a DataFrame and return a boolean vector of the rows to keep.
        Rows are loaded from the file and filtered in chunks to keep peak memory usage small. '''

    f = h5py.File(fn, 'r')

    column_names = f.attrs.get("column_names")
    column_names = get_column_intersection(column_names, query_cols)
    column_index = p.Index(column_names)

    sz = f[column_names[0]].shape[0]
    starts = np.arange(0,sz, chunk_size)
    ends = np.minimum(sz, starts+chunk_size)

    chunks = []

    for (start,end) in zip(starts, ends):
        cols = {}
        for name in column_names:
            ds = f[name]
            if has_levels(ds):
                indices = ds[start:end]
                uniques = get_levels(ds)
                col = uniques[indices]
            else:
                col = ds[start:end]

            cols[name] = col
        df = p.DataFrame(cols, columns=column_index)
        df = df[filter_func(df)]

        if len(df) > 0 or len(chunks) == 0:
            chunks.append(df)

    f.close()

    result = p.concat(chunks, ignore_index=True)
    return result


def read_data_frame_indexed(fn, queries, query_cols = [], coords = True):
    ''' Read rows from the HDF5 data frame that match each tabix query in the
    queries list.  A tabix query is in the form ('chr1', 100, 200). query_cols
    is a list of columns you want to return. If coords is True, then it it will
    return coordinates regardless of query_cols. If coords is False, it will
    only return the columns specified in query_cols. Returns a concatenated
    pandas DataFrame. WARNING: If the queries overlap in coordinates, the same
    region will appear more than once. In such cases, use
    read_data_frame_indexed_no_concat().'''

    dfs = read_data_frame_indexed_no_concat(fn, queries, query_cols, coords)

    if len(dfs) == 1:
        d = dfs[0]
    else:
        # Return the union of the queries
        d = p.concat(dfs)

    d.reset_index(drop=True, inplace=True)
    return d


def read_data_frame_indexed_no_concat(fn, tabix_queries, query_cols = [], coords = True):
    ''' Read rows from the HDF5 data frame that match each tabix query in the
    queries list.  A tabix query is in the form ('chr1', 100, 200). query_cols
    is a list of columns you want to return. If coords is True, then it it will
    return coordinates regardless of query_cols. If coords is False, it will
    only return the columns specified in query_cols. Returns a list of pandas
    DataFrames, one for each query. '''

    f = h5py.File(fn, 'r')

    # read the index
    tabix_index = read_tabix_index(f)

    dfs = []
    for q in tabix_queries:
        r = _read_data_frame_indexed_sub(f, tabix_index, q, query_cols = query_cols, coords = coords)
        dfs.append(r)

    f.close()

    # Return the union of the queries
    return dfs


class DataFrameReader:
    def __init__(self, fn):
        self.f = h5py.File(fn, 'r')
        self.tabix_index = read_tabix_index(self.f)

    def query(self, tabix_query, query_cols=[], coords=True, id_column=None):
        return _read_data_frame_indexed_sub(self.f, self.tabix_index, tabix_query, query_cols=query_cols, coords=coords, id_column=id_column)

    def query_chroms(self, chroms, query_cols=[], coords=True, id_column=None):
        '''
        hacky optimization for getting multiple whole chromosomes.
        should be merged with _read_data_frame_indexed_sub if possible.
        '''
        column_names = self.f.attrs.get("column_names")
        column_names = get_column_intersection(column_names, query_cols)

        # parse the index
        (index, chrom_col, start_col, end_col) = self.tabix_index

        if len(query_cols) > 0:
            column_names = list(set([chrom_col, start_col, end_col] + column_names))

        # query interval
        slices = get_range_for_chroms(index, chroms)

        # get DataFrame
        df = _get_sliced_df(self.f, column_names, slices, id_column=id_column)

        # Do the final filter of rows
        if not coords:
            df = df[ query_cols ]
        df.reset_index(drop=True, inplace=True)
        return df

    def close(self):
        self.f.close()
        self.f = None
        self.tabix_index = None

def _read_data_frame_indexed_sub(f, tabix_index, tabix_query, query_cols = [], coords = True, id_column=None):
    ''' Read rows from the HDF5 data frame that match the tabix query.
    A tabix query is in the form ('chr1', 100, 200). query_cols
    is a list of columns you want to return. If coords is True, then it it will
    return coordinates regardless of query_cols. If coords is False, it will
    only return the columns specified in query_cols. Returns a pandas DataFrame.'''

    column_names = f.attrs.get("column_names")
    column_names = get_column_intersection(column_names, query_cols)

    # parse the index
    (index, chrom_col, start_col, end_col) = tabix_index

    if len(query_cols) > 0:
        column_names = list(set([chrom_col, start_col, end_col] + column_names))

    # query interval
    (q_chrom, q_start, q_end) = tabix_query

    # get row slices corresponding to this interval.
    # tabix index guarantees that these rows are a superset of the required rows
    slices = find_candidate_ranges(index, (q_chrom, q_start, q_end))

    # get DataFrame
    df = _get_sliced_df(f, column_names, slices, id_column=id_column)

    # Do the final filter of rows
    subset_df = df[ (df[chrom_col] == q_chrom) & (df[end_col] >= q_start) & (df[start_col] < q_end) ]
    if not coords:
        subset_df = subset_df[ query_cols ]
    subset_df.reset_index(drop=True, inplace=True)
    return subset_df

def _get_sliced_df(h5file, column_names, row_slices, id_column=None):
    columns = [ (name, h5file[name], get_levels(h5file[name])) for name in column_names ]

    result_cols = {}

    for (name, ds, translate) in columns:
        if len(row_slices) > 0:
            row_slices.sort()
            rows = np.concatenate([ ds[start:end] for (start, end) in row_slices ])
        else:
            # we'll return an empty data frame if there are no slices
            # np.concatenate fail with 0-length input
            rows = np.array([], dtype=ds.dtype)
        if translate is not None:
            rows = translate[rows]

        result_cols[name] = rows

    if len(row_slices) > 0:
        id_column_values = np.concatenate([np.arange(start,end) for (start, end) in row_slices])
    else:
        id_column_values = np.array([], dtype=np.int32)

    if id_column is not None:
        result_cols[id_column] = id_column_values

    df = p.DataFrame(result_cols)
    df.index = id_column_values

    return df

def chunk_iter(ds, chunk_size):
    l = ds.shape[0]

    for start in xrange(0, l, chunk_size):
        yield ds[start:min(l, start+chunk_size)]

def create_tabix_index(fn, chrom_col, start_col, end_col=None):
    f = h5py.File(fn, 'r+')

    column_names = f.attrs.get("column_names")
    #column_index = p.Index(column_names)

    if end_col is None:
        end_col = start_col

    assert(chrom_col in column_names)
    assert(start_col in column_names)
    assert(end_col in column_names)

    index = {}

    chrom_ds = f[chrom_col]
    start_ds = f[start_col]
    end_ds = f[end_col]

    # check for factors in chrom_ds
    chrom_translate = get_levels(chrom_ds)

    # length
    l = chrom_ds.shape[0]

    # Loop over all the data to create the index
    chunk_iters = itertools.izip(xrange(0, l, CHUNK_SIZE), chunk_iter(chrom_ds, CHUNK_SIZE), chunk_iter(start_ds, CHUNK_SIZE), chunk_iter(end_ds, CHUNK_SIZE))
    for (offset, chrom, start, end) in chunk_iters:

        if chrom_translate is not None:
            chrom = chrom_translate[chrom]

        tabix_index_dataset(index, chrom, start, end, offset)


    # Write the index dataset into the file
    write_tabix_index(f, index, chrom_col, start_col, end_col)
    f.close()



def combine_data_frame_files(output_filename, input_filenames):
    in_files = [ h5py.File(f, 'r') for f in input_filenames ]
    column_names = [ tuple(sorted(f.attrs.get("column_names"))) for f in in_files ]

    uniq = set(column_names)

    if len(uniq) > 1:
        raise Exception("you're attempting to combine incompatible data frames")

    if len(uniq) == 0:
        r = "No input files? output: %s, inputs: %s" % (output_filename, str(input_filenames))
        raise Exception(r)

    column_names = uniq.pop()

    if os.path.exists(output_filename):
        os.remove(output_filename)

    out = h5py.File(output_filename)
    out.attrs.create("column_names", column_names)

    # Write successive columns
    for c in column_names:
        datasets = [f[c] for f in in_files if len(f[c]) > 0]
        num_w_levels = np.sum([has_levels(ds) for ds in datasets if len(ds) > 0])
        fract_w_levels = float(num_w_levels) / (len(datasets) + 1)

        if fract_w_levels > 0.25:
            combine_level_column(out, datasets, c)
            continue

        # filter out empty rows from the type promotion, unless they're all empty
        types = [get_col_type(ds) for ds in datasets if len(ds) > 0]
        if len(types) == 0:
            # Fall back to getting column types from empty data frames
            types = [get_col_type(f[c]) for f in in_files]
        common_type = reduce(np.promote_types, types)

        # numpy doesn't understand vlen strings -- so always promote to vlen strings if anything is using them
        if vlen_string in types:
            common_type = vlen_string

        out_ds = out.create_dataset(c, shape=(0,), maxshape=(None,), dtype=common_type, compression=COMPRESSION, shuffle=True, chunks=(CHUNK_SIZE,))

        item_count = 0
        for ds in datasets:
            new_items = ds.shape[0]
            out_ds.resize((item_count + new_items,))
            data = ds[:]

            if has_levels(ds):
                levels = get_levels(ds)
                data = levels[data]

            out_ds[item_count:(item_count + new_items)] = data
            item_count += new_items

    for in_f in in_files:
        in_f.close()

    out.close()

def combine_level_column(out_h5, datasets, column_name):
    ''' Concatenate columns where one or more of the input columns uses levels into
        an output column that uses levels '''

    def always_get_levels(ds):
        levels = get_levels(ds)
        if levels is not None:
            return levels
        else:
            return np.unique(ds[:])

    level_sets = [always_get_levels(ds) for ds in datasets]

    # Determine the new level set, and a mapping from the old index
    #to the new index for each input dataset
    new_levels = np.unique(np.concatenate(level_sets))
    levels_dict = { value:idx for (idx, value) in enumerate(new_levels) }

    dt = np.uint8
    if len(new_levels) > 2**8:
        dt = np.uint16

    if len(new_levels) > 2**16:
        dt = np.uint32

    out_ds = out_h5.create_dataset(column_name, shape=(0,), maxshape=(None,), dtype=dt, compression=COMPRESSION, shuffle=True, chunks=(CHUNK_SIZE,))
    create_levels(out_ds, new_levels)

    item_count = 0
    for ds in datasets:
        old_levels = get_levels(ds)

        if old_levels is None:
            # Convert over to using levels
            new_data = make_index_array_dict(levels_dict, ds[:], dt)
        else:
            # Dictionary to convert from old levels to new levels
            convert_array = np.array([levels_dict[old_levels[old_index]] for old_index in range(len(old_levels))], dtype=dt)
            # convert to new factor codes
            old_data = ds[:]
            new_data = convert_array[old_data]

        new_items = ds.shape[0]
        out_ds.resize((item_count + new_items,))
        out_ds[item_count:(item_count + new_items)] = new_data
        item_count += new_items



# Tabix magic numbers and functions
bit_levels = [ 11, 14, 17, 20, 23, 26 ]
desc_bit_bins = [ (26, 1), (23, 9), (20, 73), (17, 585), (14, 4681) ]
rev_bit_bins = list(desc_bit_bins)
rev_bit_bins.reverse()


def reg2bin(begin, end):
    '''Compute the tabix bin of the interval'''
    for (bits, bins) in rev_bit_bins:
        if begin >> bits == end >> bits:
            return ((1 << (29 - bits)) - 1) / 7 + (begin >> bits)

    return 0


def reg2bin_vector(begin, end):
    '''Vectorized tabix reg2bin -- much faster than reg2bin'''
    result = np.zeros(begin.shape)

    # Entries filled
    done = np.zeros(begin.shape, dtype=np.bool)

    for (bits, bins) in rev_bit_bins:
        begin_shift = begin >> bits
        new_done = (begin >> bits) == (end >> bits)
        mask = np.logical_and(new_done, np.logical_not(done))
        offset = ((1 << (29 - bits)) - 1) / 7
        result[mask] = offset + begin_shift[mask]

        done = new_done

    return result.astype(np.int32)


def reg2bins(rbeg, rend):
    ''' List the tabix bins that may contain elements intersecting the query interval '''
    bins = [0]  # always include the 0 bin
    rbeg, rend = int(rbeg), int(rend)

    for (bit_shift, prev_bins) in desc_bit_bins:
        bins.extend ( range(prev_bins + (rbeg >> bit_shift), 1 + prev_bins + (rend >> bit_shift)) )

    return bins


def insert_index(idx, tid, bin, beg, end):
    ''' Insert the interval of rows [beg,end] with constant tid and bin into tabix index'''
    tidx = idx.setdefault(tid, {})
    bin_list = tidx.setdefault(bin, [])
    bin_list.append((beg, end))



def find_const_intervals(tid_col, bin_col):
    ''' tid_col is a column indicating is the contig, and bin_col is the bin.
        report tuple of (contig, bin, start_idx, end_idx) for row ranges with
        constant contig and bin '''
    change_points = np.logical_or(tid_col[1:] != tid_col[:-1], bin_col[1:] != bin_col[:-1])
    change_points = np.where(change_points)[0]

    intv_starts = np.concatenate([[0], change_points+1])
    intv_ends =   np.concatenate([change_points, [tid_col.shape[0] - 1]]) + 1

    return [ (tid_col[start], bin_col[start], start, end) for (start,end) in zip(intv_starts, intv_ends) ]


def tabix_index_dataset(index, tid_col, start_col, end_col, offset):
    '''Generate a tabix index for the elements represented by tid_col, start_col,
    and end_col.  The index can be built incrementally by increasing offset
    and calling this function multiple times.  The index variable is updated, not returned'''

    bin_col = reg2bin_vector(start_col, end_col)
    index_intervals = find_const_intervals(tid_col, bin_col)

    for (tid, bin, start_idx, end_idx) in index_intervals:
        insert_index(index, tid, bin, start_idx + offset, end_idx + offset)


def write_tabix_index(h5_group, index, tid_col, start_col, end_col):

    if "_index" in h5_group.keys():
        del h5_group["_index"]

    index_group = h5_group.create_group("_index")

    index_rows = []

    # Support arbitrary tid keys
    # map each template key to an integer
    # and store the map
    tid_levels = list(index.keys())
    tid_map = { key:idx for (idx, key) in enumerate(tid_levels) }

    index_group.create_dataset("tid_levels", data=tid_levels)

    index_group.attrs.create("index_format_version", 3)

    for (tid, bin_lists) in index.items():
        for (bin, intervals) in bin_lists.items():
            for intv in intervals:
                index_rows.append([tid_map[tid], bin, intv[0], intv[1]])

    idx_array = np.array(index_rows, dtype=np.int64)

    index_group.create_dataset("index_rows", data=idx_array,
            maxshape=(None, 4), chunks=(1<<16,4), shuffle=True, compression=1)

    index_group.attrs.create("tid_col", tid_col)
    index_group.attrs.create("start_col", start_col)
    index_group.attrs.create("end_col", end_col)


def read_tabix_index(h5_group):
    index_group = h5_group["_index"]
    index_ds = index_group["index_rows"]

    # Old version of formoat stored the tid_levels in the tid_levels attribute
    # Current version has moved them to a dataset.
    # Check for dataset, fall back to attribute

    if "tid_levels" in index_group:
        tid_levels = index_group["tid_levels"]
    else:
        tid_levels = index_group.attrs.get("tid_levels")

    format_version = index_group.attrs.get("index_format_version", 1)

    tid_code = index_ds[:,0]
    bin = index_ds[:,1]
    row_start = index_ds[:,2]
    row_end = index_ds[:,3]

    # In the original format we used (closed, closed) intervals in the in index.
    # convert to (closed, open)
    if format_version == 1:
        row_end = row_end + 1

    # indexing columns
    tid_col = index_group.attrs.get("tid_col")
    start_col = index_group.attrs.get("start_col")
    end_col = index_group.attrs.get("end_col")

    index = {}

    # index = map with the struture {tid: {bin: ([...record_starts...], [...record_ends...])}}
    for ((_tid, _bin), rows) in itertools.groupby(itertools.izip(tid_code, bin, row_start, row_end), key=lambda x: x[:2]):
        rows = list(rows)
        start_chunk = np.array([x[2] for x in rows], dtype=np.int64)
        end_chunk = np.array([x[3] for x in rows], dtype=np.int64)

        t_idx = index.setdefault(tid_levels[_tid], {})
        t_idx[_bin] = (start_chunk, end_chunk)

    return (index, tid_col, start_col, end_col)

def get_range_for_chroms(index, chroms):
    # return the set of ranges that contain all bins for all specified chroms
    start_chunks = []
    end_chunks = []

    for chrom in chroms:
        if chrom in index:
            t_index = index[chrom]
            for bin in t_index.keys():
                start_chunk, end_chunk = t_index[bin]
                start_chunks.append(start_chunk)
                end_chunks.append(end_chunk)

    return _ranges_from_chunks(start_chunks, end_chunks)

def find_candidate_ranges(index, query):
    '''Find the intervals in the tabix index that may contain entries matching the range query'''
    (tid, start, end) = query

    try:
        t_index = index[tid]
    except KeyError:
        return []

    valid_bins = reg2bins(start, end)

    start_chunks = []
    end_chunks = []

    for bin in valid_bins:
        if t_index.has_key(bin):
            start_chunk, end_chunk = t_index[bin]
            start_chunks.append(start_chunk)
            end_chunks.append(end_chunk)

    return _ranges_from_chunks(start_chunks, end_chunks, 512)

def _ranges_from_chunks(start_chunks, end_chunks, gap_jump=0):
    if len(start_chunks) == 0:
        return []

    start = np.concatenate(start_chunks)
    end = np.concatenate(end_chunks)

    ranges = sorted(zip(start,end))
    ranges.sort()
    new_ranges = []

    cstart = None
    cend = None

    for (start, end) in ranges:
        if cstart is None:
            cstart = start
        if cend is None:
            cend = end

        if start <= cend + gap_jump:
            cend = max(end,cend)
        else:
            new_ranges.append((cstart,cend))
            cstart = start
            cend = end

    new_ranges.append((cstart,cend))

    return new_ranges
