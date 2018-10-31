#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import struct

def serialize_batch_nearest_neighbor(fp, batch_nearest_neighbor):
    for (a, b), s in batch_nearest_neighbor.iteritems():
        fp.write(struct.pack("qqQ", a, b, len(s)))
        for i, j in s:
            fp.write(struct.pack("qq", i, j))

def deserialize_batch_nearest_neighbor(fp):
    """
    >>> from cStringIO import StringIO
    >>> batch1 = dict()
    >>> batch1[(0, 1)] = set([(1, 2), (3, 4), (5, 6)])
    >>> batch1[(1, 2)] = set([(7, 8), (9, 10)])
    >>> batch1[(3, 4)] = set([(11, 12)])
    >>> fp = StringIO()
    >>> serialize_batch_nearest_neighbor(fp, batch1)
    >>> fp.seek(0)
    >>> batch2 = deserialize_batch_nearest_neighbor(fp)
    >>> batch1 == batch2
    True
    """
    batch_nearest_neighbor = {}
    while True:
        fmt = "qqQ"
        sz = struct.calcsize("qqQ")
        buf = fp.read(sz)
        if len(buf) == 0:
            break
        elif len(buf) != sz:
            raise RuntimeError("corrupted batch_nearest_neighbor stream (key)")
        a, b, slen = struct.unpack(fmt, buf)
        fmt = "qq"
        sz = struct.calcsize("qq")
        s = set()
        for _ in xrange(slen):
            buf = fp.read(sz)
            if len(buf) != sz:
                raise RuntimeError("corrupted batch_nearest_neighbor stream (set)")
            i, j = struct.unpack(fmt, buf)
            s.add((i, j))
        batch_nearest_neighbor[(a, b)] = s
    return batch_nearest_neighbor
