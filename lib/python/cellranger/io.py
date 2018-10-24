#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import h5py
import hashlib
from HTMLParser import HTMLParser
import json
import numpy as np
import os
import pandas as pd
import errno
import shutil
import six
import subprocess
import sys
import tables
import gzip
import lz4.frame as lz4
import _io as io  # this is necessary b/c this module is named 'io' ... :(
import tenkit.log_subprocess as tk_subproc
import cellranger.h5_constants as h5_constants
from cellranger.constants import (MEM_GB_PER_THREAD, STRING_DTYPE_CHARS,
                                  STR_DTYPE_CHAR, UNICODE_DTYPE_CHAR)

def get_thread_request_from_mem_gb(mem_gb):
    """ For systems without memory reservations, reserve multiple threads if necessary to avoid running out of memory"""
    est_threads = round(float(mem_gb) / MEM_GB_PER_THREAD)
    # make sure it's 1, 2, or 4
    for threads in [1, 2, 4]:
        if est_threads <= threads: return threads
    return 4

def fixpath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def get_input_path(oldpath, is_dir=False):
    path = fixpath(oldpath)
    if not os.path.exists(path):
        sys.exit("Input file does not exist: %s" % path)
    if is_dir:
        if not os.path.isdir(path):
            sys.exit("Please provide a directory, not a file: %s" % path)
    else:
        if not os.path.isfile(path):
            sys.exit("Please provide a file, not a directory: %s" % path)
    return path

def get_input_paths(oldpaths):
    paths = []
    for oldpath in oldpaths:
        paths.append(get_input_path(oldpath))
    return paths

def get_output_path(oldpath):
    path = fixpath(oldpath)
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        sys.exit("Output directory does not exist: %s" % dirname)
    if not os.path.isdir(dirname):
        sys.exit("Please provide a directory, not a file: %s" % dirname)
    return path

def write_h5(filename, data):
    with h5py.File(filename, 'w') as f:
        for key, value in data.iteritems():
            f[key] = value

class SubprocessStream(object):
    """ Wrap a subprocess that we stream from or stream to """
    def __init__(self, *args, **kwargs):
        mode = kwargs.pop('mode', 'r')
        if mode == 'r':
            kwargs['stdout'] = subprocess.PIPE
        elif mode == 'w':
            kwargs['stdin'] = subprocess.PIPE
        else:
            raise ValueError('mode %s unsupported' % self.mode)

        kwargs['bufsize'] = 65536  # 64K of buffering in the pipe
        kwargs['preexec_fn'] = os.setsid
        print args[0]
        sys.stdout.flush()
        self.proc = tk_subproc.Popen(*args, **kwargs)

        if mode == 'r':
            self.pipe = self.proc.stdout
        elif mode == 'w':
            self.pipe = self.proc.stdin
        self._read = mode == 'r'
    def __enter__(self):
        return self
    def __iter__(self):
        return self
    def next(self):
        return self.pipe.next()
    def fileno(self):
        return self.pipe.fileno()
    def write(self, x):
        self.pipe.write(x)
    def close(self):
        self.pipe.close()
        if self._read:
            # Just in case the job is stuck waiting for its pipe to get flushed,
            # kill it, since there's nothing listening on the other end of the pipe.
            try:
                self.proc.kill()
            except:
                pass
        self.proc.wait()
    def __exit__(self, tp, val, tb):
        self.close()

def open_maybe_gzip_subprocess(filename, mode='r'):
    compressor = None

    if filename.endswith(h5_constants.GZIP_SUFFIX):
        compressor = 'gzip'
    elif filename.endswith(h5_constants.LZ4_SUFFIX):
        compressor = 'lz4'
    else:
        return open(filename, mode)

    assert compressor is not None

    if mode == 'r':
        return SubprocessStream([compressor, '-c', '-d', filename], mode='r')
    elif mode == 'w':
        f = open(filename, 'w')
        return SubprocessStream([compressor, '-c'], stdout=f, mode='w')
    else:
        raise ValueError("Unsupported mode for compression: %s" % mode)

def open_maybe_gzip(filename, mode='r'):
    # this _must_ be a str
    filename = str(filename)
    if filename.endswith(h5_constants.GZIP_SUFFIX):
        raw = gzip.open(filename, mode + 'b', 2)
    elif filename.endswith(h5_constants.LZ4_SUFFIX):
        raw = lz4.open(filename, mode + 'b')
    else:
        return open(filename, mode)

    bufsize = 1024*1024  # 1MB of buffering
    if mode == 'r':
        return io.BufferedReader(raw, buffer_size=bufsize)
    elif mode == 'w':
        return io.BufferedWriter(raw, buffer_size=bufsize)
    else:
        raise ValueError("Unsupported mode for compression: %s" % mode)

class CRCalledProcessError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def run_command_safely(cmd, args):
    p = tk_subproc.Popen([cmd] + args, stderr=subprocess.PIPE)
    _, stderr_data = p.communicate()
    if p.returncode != 0:
        raise Exception("%s returned error code %d: %s" % (p, p.returncode, stderr_data))

def check_completed_process(p, cmd):
    """ p   (Popen object): Subprocess
        cmd (str):          Command that was run
    """
    if p.returncode is None:
        raise CRCalledProcessError("Process did not finish: %s ." % cmd)
    elif p.returncode != 0:
        raise CRCalledProcessError("Process returned error code %d: %s ." % (p.returncode, cmd))

def mkdir(dst, allow_existing=False):
    """ Create a directory. Optionally succeed if already exists.
        Useful because transient NFS server issues may induce double creation attempts. """
    if allow_existing:
        try:
            os.mkdir(dst)
        except OSError as e:
            if e.errno == errno.EEXIST and os.path.isdir(dst):
                pass
            else:
                raise
    else:
        os.mkdir(dst)

def makedirs(dst, allow_existing=False):
    """ Create a directory recursively. Optionally succeed if already exists.
        Useful because transient NFS server issues may induce double creation attempts. """
    if allow_existing:
        try:
            os.makedirs(dst)
        except OSError as e:
            if e.errno == errno.EEXIST and os.path.isdir(dst):
                pass
            else:
                raise
    else:
        os.makedirs(dst)

def remove(f, allow_nonexisting=False):
    """ Delete a file. Allow to optionally succeed if file doesn't exist.
        Useful because transient NFS server issues may induce double deletion attempts. """
    if allow_nonexisting:
        try:
            os.remove(f)
        except OSError as e:
            if e.errno == errno.ENOENT:
                pass
            else:
                raise
    else:
        os.remove(f)

def copy(src, dst):
    """ Safely copy a file. Not platform-independent """
    run_command_safely('cp', [src, dst])

def move(src, dst):
    """ Safely move a file. Not platform-independent """
    run_command_safely('mv', [src, dst])

def copytree(src, dst, allow_existing=False):
    """ Safely recursively copy a directory. Not platform-independent """
    makedirs(dst, allow_existing=allow_existing)

    for name in os.listdir(src):
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)

        if os.path.isdir(srcname):
            copytree(srcname, dstname)
        else:
            copy(srcname, dstname)

def concatenate_files(out_path, in_paths, mode=''):
    with open(out_path, 'w' + mode) as out_file:
        for in_path in in_paths:
            with open(in_path, 'r' + mode) as in_file:
                shutil.copyfileobj(in_file, out_file)

def concatenate_headered_files(out_path, in_paths, mode=''):
    """ Concatenate files, taking the first line of the first file
        and skipping the first line for subsequent files.
        Asserts that all header lines are equal. """
    with open(out_path, 'w' + mode) as out_file:
        if len(in_paths) > 0:
            # Write first file
            with open(in_paths[0], 'r' + mode) as in_file:
                header = in_file.readline()
                out_file.write(header)
                shutil.copyfileobj(in_file, out_file)

        # Write remaining files
        for in_path in in_paths[1:]:
            with open(in_path, 'r' + mode) as in_file:
                this_header = in_file.readline()
                assert this_header == header
                shutil.copyfileobj(in_file, out_file)

def compute_hash_of_file(filename, block_size_bytes=2**20):
    digest = hashlib.sha1()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(block_size_bytes), b''):
            digest.update(chunk)
    return digest.hexdigest()

def write_empty_json(filename):
    with open(filename, 'w') as f:
        json.dump({}, f)

def load_csv_rownames(csv_file):
    rownames = np.atleast_1d(pd.read_csv(csv_file, usecols=[0]).values.squeeze())
    return rownames

def get_h5_filetype(filename):
    with tables.open_file(filename, mode = 'r') as f:
        try:
            filetype = f.get_node_attr('/', h5_constants.H5_FILETYPE_KEY)
        except AttributeError:
            filetype = None # older files lack this key
    return filetype

def save_array_h5(filename, name, arr):
    """ Save an array to the root of an h5 file """
    with tables.open_file(filename, 'w') as f:
        f.create_carray(f.root, name, obj=arr)

def load_array_h5(filename, name):
    """ Load an array from the root of an h5 file """
    with tables.open_file(filename, 'r') as f:
        return getattr(f.root, name).read()

def merge_jsons_single_level(filenames):
    """ Merge a list of toplevel-dict json files.
        Union dicts at the first level.
        Raise exception when duplicate keys encountered.
        Returns the merged dict.
    """

    merged = {}
    for filename in filenames:
        with open(filename) as f:
            d = json.load(f)

            for key, value in d.iteritems():
                if key not in merged:
                    merged[key] = value
                    continue

                # Handle pre-existing keys
                if merged[key] == value:
                    pass
                elif type(value) == dict and type(merged[key]) == dict:
                    merged[key].update(value)
                else:
                    raise ValueError("No merge strategy for key %s, value of type %s, from %s into %s" % (str(key),
                                                                                                          str(type(value)),
                                                                                                          str(value),
                                                                                                          str(merged[key])))
    return merged

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
    """

    if data is None or hasattr(data, '__len__') and  len(data) == 0:
        dtype = '%s1' % STR_DTYPE_CHAR
        group.create_dataset(name, dtype=dtype)
        return

    assert (isinstance(data, np.ndarray) and data.dtype.char == STR_DTYPE_CHAR) or \
            (isinstance(data, list) and \
            all(x is None or isinstance(x, six.string_types) for x in data))

    # Convert Nones to empty strings and use XML encoding
    data = map(lambda x: x.encode('ascii', 'xmlcharrefreplace') if x is not None else '', data)

    fixed_len = max(len(x) for x in data)

    # h5py doesn't support strings with zero-length-dtype
    if fixed_len == 0:
        fixed_len = 1
    dtype = '%s%d' % (STR_DTYPE_CHAR, fixed_len)

    group.create_dataset(name, data=data, dtype=dtype, **kwargs)

def make_utf8(x):
    """Encode a string as UTF8

    Respect both python2 and python3."""
    if isinstance(x, six.text_type):
        return x
    elif isinstance(x, six.binary_type):
        return x.decode('utf8')
    else:
        raise ValueError('Expected string type, got type %s' % str(type(x)))

def decode_ascii_xml(x):
    """Decode a string from 7-bit ASCII + XML into unicode.
    """
    if isinstance(x, six.text_type):
        return x
    elif isinstance(x, six.binary_type):
        return make_utf8(HTMLParser.unescape.__func__(HTMLParser, x))
    else:
        raise ValueError('Expected string type, got type %s' % str(type(x)))

def decode_ascii_xml_array(data):
    """Decode an array-like container of strings from 7-bit ASCII + XML
    into unicode.
    """
    if isinstance(data, np.ndarray) and \
       data.dtype.char == UNICODE_DTYPE_CHAR:
        return data

    unicode_data = map(decode_ascii_xml, data)

    fixed_len = max(len(s) for s in unicode_data)
    # 0 length string type is Ok
    dtype = '%s%d' % (UNICODE_DTYPE_CHAR, fixed_len)

    # note: python3 would require no.fromiter
    return np.array(unicode_data, dtype=dtype)

def encode_ascii_xml(x):
    """Encode a string as fixed-length 7-bit ASCII with XML-encoding
    for characters outside of 7-bit ASCII.

    Respect python2 and python3, either unicode or binary.
    """
    if isinstance(x, six.text_type):
        return x.encode('ascii', 'xmlcharrefreplace')
    elif isinstance(x, six.binary_type):
        return x
    else:
        raise ValueError('Expected string type, got type %s' % str(type(x)))

def encode_ascii_xml_array(data):
    """Encode an array-like container of strings as fixed-length 7-bit ASCII
    with XML-encoding for characters outside of 7-bit ASCII.
    """
    if isinstance(data, np.ndarray) and \
       data.dtype.char == STR_DTYPE_CHAR and \
       data.dtype.itemsize > 0:
        return data

    convert = lambda s: encode_ascii_xml(s) if s is not None else ''
    ascii_data = map(convert, data)

    fixed_len = max(len(s) for s in ascii_data)
    fixed_len = max(1, fixed_len)
    dtype = '%s%d' % (STR_DTYPE_CHAR, fixed_len)

    # note: python3 would require np.fromiter
    return np.array(ascii_data, dtype=dtype)

def read_hdf5_string_dataset(dataset):
    """Read a dataset of strings from HDF5 (h5py).

    Args:
        dataset (h5py.Dataset): Data to read.
    Returns:
        list of unicode - Strings in the dataset (as utf8)
    """

    # h5py doesn't support loading an empty dataset
    if dataset.shape is None:
        return []

    # Test that the dataset is valid ASCII
    data = dataset[:]
    ascii_strings = [x.decode('ascii', 'ignore') for x in data]

    # Unescape any XML encoded characters
    unescaped = [HTMLParser.unescape.__func__(HTMLParser, x) for x in ascii_strings]

    return map(make_utf8, unescaped)

def set_hdf5_attr(dataset, name, value):
    """Set an attribute of an HDF5 dataset/group.

    Strings are stored as fixed-length 7-bit ASCII with XML-encoding
    for characters outside of 7-bit ASCII. This is inspired by the
    choice made for the Loom spec:
    https://github.com/linnarsson-lab/loompy/blob/master/doc/format/index.rst
    """

    name = encode_ascii_xml(name)

    if isinstance(value, (six.text_type, six.binary_type)):
        value = encode_ascii_xml(value)
    elif isinstance(value, np.ndarray) and value.dtype.char in STRING_DTYPE_CHARS:
        value = encode_ascii_xml_array(value)
    elif isinstance(value, (list, tuple)) and isinstance(value[0], (six.text_type, six.binary_type)):
        value = encode_ascii_xml_array(value)

    dataset.attrs[name] = value
