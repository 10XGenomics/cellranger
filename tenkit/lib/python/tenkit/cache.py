#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import resource
import collections

class FileHandleCache:
    ''' LRU cache for file handles. '''

    def __init__(self, mode='w', open_func=open):
        self.mode = mode
        self.open_func = open_func
        self.config_max_files()
        self.have_opened = {}
        self.open_files = collections.OrderedDict()

    def config_max_files(self):
        soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
        self.maxfiles = soft - 100

    def get(self, fn):
        if self.open_files.has_key(fn):
            # move to front
            fh = self.open_files.pop(fn)
            self.open_files[fn] = fh
            return fh
        else:
            if 'w' in self.mode and fn in self.have_opened:
                mode = self.mode.replace('w', 'a')
            else:
                mode = self.mode

            # Close an open file if we are about to fill the cache
            if len(self.open_files) == self.maxfiles - 1:
                close_fn, close_fh = self.open_files.popitem(last=False)
                self.have_opened[close_fn] = close_fh.tell()
                close_fh.close()

            # open the file
            fh = self.open_func(fn, mode)

            # seek to previous position if its been opened before for reading
            if 'r' in mode and fn in self.have_opened:
                fh.seek(self.have_opened[fn])

            # put it on the LRU
            self.open_files[fn] = fh
            return fh

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for f in self.open_files.values():
            f.close()
