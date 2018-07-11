#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
# fadvise support
#
import logging
from cffi import FFI

__all__ = ['fadvise', 'FADV_WILLNEED']

ffi = FFI()
ffi.cdef("int posix_fadvise(int fd, long int offset, long int len, int advice);")
C = ffi.dlopen(None)

FADV_WILLNEED = 3

fadv_error = False
def fadvise(fd, offset, length, advice):
    global fadv_error
    try:
        rc = C.posix_fadvise(fd, offset, length, advice)
        if rc != 0 and not fadv_error:
            fadv_error = True
            logging.warning("posix_fadvise system call failed with return code {}".format(rc))
    except AttributeError:
        if not fadv_error:
            fadv_error = True
            logging.warning("posix_fadvise not found in standard library")
