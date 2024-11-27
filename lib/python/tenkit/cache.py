#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import collections
import os
import resource
from collections.abc import Callable
from typing import BinaryIO, Generic, Literal, TextIO, TypeVar, overload

_T_co = TypeVar("_T_co", bound=TextIO | BinaryIO, covariant=True)


class FileHandleCache(Generic[_T_co]):
    """LRU cache for file handles."""

    @overload
    def __init__(
        self,
        mode: Literal["w", "r"] = ...,
        open_func: Callable[
            [str | bytes | os.PathLike[str] | os.PathLike[bytes], Literal["w", "r", "a"]],
            TextIO,
        ] = ...,
    ): ...

    @overload
    def __init__(
        self,
        mode: Literal["wb", "rb"] = ...,
        open_func: Callable[
            [str | bytes | os.PathLike[str] | os.PathLike[bytes], Literal["wb", "rb", "ab"]],
            BinaryIO,
        ] = ...,
    ): ...

    def __init__(
        self,
        mode: Literal["w", "r", "wb", "rb"] = "w",
        open_func: Callable[
            [
                str | bytes | os.PathLike[str] | os.PathLike[bytes],
                str,
            ],
            _T_co,
        ] = open,
    ):
        self.mode = mode
        self.open_func: Callable[
            [str | bytes | os.PathLike[str] | os.PathLike[bytes], str], _T_co
        ] = open_func
        self.config_max_files()
        self.have_opened: dict[str | bytes | os.PathLike[str] | os.PathLike[bytes], int] = {}
        self.open_files: collections.OrderedDict[
            str | bytes | os.PathLike[str] | os.PathLike[bytes], _T_co
        ] = collections.OrderedDict()

    def config_max_files(self) -> None:
        soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
        self.maxfiles = soft - 100

    def get(self, fn: str | bytes) -> _T_co:
        """Get the cached file object for the given file name, or create one.

        Note:
            While function can take a str or bytes path, they will be considered
            different paths for purposes of caching.
        """
        if fn in self.open_files:
            # move to front
            fh = self.open_files.pop(fn)
            self.open_files[fn] = fh
            return fh
        else:
            if "w" in self.mode and fn in self.have_opened:
                mode = self.mode.replace("w", "a")
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
            if "r" in mode and fn in self.have_opened:
                fh.seek(self.have_opened[fn])

            # put it on the LRU
            self.open_files[fn] = fh
            return fh

    def __enter__(self) -> FileHandleCache[_T_co]:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        for f in self.open_files.values():
            f.close()
