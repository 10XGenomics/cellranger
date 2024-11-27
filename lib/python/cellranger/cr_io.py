#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import errno
import gzip
import io
import os
import shutil
import sys
from collections.abc import Iterable, Sequence
from typing import IO, TYPE_CHECKING, Any, AnyStr, Literal, TextIO, overload

import lz4.frame as lz4
import martian

import cellranger.h5_constants as h5_constants

if TYPE_CHECKING:
    import subprocess


def fixpath(path: AnyStr) -> AnyStr:
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))


def get_input_path(oldpath: AnyStr, is_dir: bool = False) -> AnyStr:
    if not isinstance(oldpath, str | bytes):
        sys.exit(f"'{oldpath}' is not a valid string type and is an invalid path")
    path = fixpath(oldpath)
    if not os.path.exists(path):
        sys.exit(f"Input file does not exist: {path}")
    if not os.access(path, os.R_OK):
        sys.exit(f"Input file path {path} does not have read permissions")
    if is_dir:
        if not os.path.isdir(path):
            sys.exit(f"Please provide a directory, not a file: {path}")
    elif not os.path.isfile(path):
        sys.exit(f"Please provide a file, not a directory: {path}")
    return path


def get_input_paths(paths: list[AnyStr]) -> list[AnyStr]:
    return [get_input_path(path) for path in paths]


def get_output_path(oldpath: AnyStr) -> AnyStr:
    path = fixpath(oldpath)
    dirname = os.path.dirname(path)
    if not os.path.isdir(dirname):
        if isinstance(dirname, bytes):
            dirname = dirname.decode(errors="replace")
        if not os.path.exists(dirname):
            sys.exit(f"Output directory does not exist: {dirname}")
        sys.exit(f"Please provide a directory, not a file: {dirname}")
    return path


@overload
def open_maybe_gzip(filename: str | bytes, mode: Literal["r"] = "r") -> TextIO: ...


@overload
def open_maybe_gzip(filename: str | bytes, mode: Literal["w"] = ...) -> TextIO: ...


@overload
def open_maybe_gzip(filename: str | bytes, mode: Literal["rb"] = ...) -> io.BufferedReader: ...


@overload
def open_maybe_gzip(filename: str | bytes, mode: Literal["wb"] = ...) -> io.BufferedReader: ...


def open_maybe_gzip(filename: str | bytes, mode: str = "r") -> IO[Any]:
    # this _must_ be a bytes
    if not isinstance(filename, bytes):
        filename = str(filename).encode()
    if filename.endswith(h5_constants.GZIP_SUFFIX):
        raw = gzip.open(filename, mode, 2)
    elif filename.endswith(h5_constants.LZ4_SUFFIX):
        raw = lz4.open(filename, mode)
    else:
        return open(filename, mode)

    bufsize = 1024 * 1024  # 1MB of buffering
    if mode == "r":
        return io.TextIOWrapper(io.BufferedReader(raw, buffer_size=bufsize))
    elif mode == "w":
        return io.TextIOWrapper(io.BufferedWriter(raw, buffer_size=bufsize))
    elif mode == "rb":
        return io.BufferedReader(raw, buffer_size=bufsize)
    elif mode == "wb":
        return io.BufferedWriter(raw, buffer_size=bufsize)

    else:
        raise ValueError(f"Unsupported mode for compression: {mode}")


class CRCalledProcessError(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__(msg)
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


def check_completed_process(p: subprocess.CompletedProcess, cmd: str) -> None:
    """Raises an exception if the completed process failed.

    Args:
        p:   Subprocess
        cmd: Command that was run

    Raises:
        CRCalledProcessError
    """
    if p.returncode is None:
        raise CRCalledProcessError(f"Process did not finish: {cmd} .")
    elif p.returncode != 0:
        raise CRCalledProcessError("Process returned error code %d: %s ." % (p.returncode, cmd))


def mkdir(path: str | bytes, exist_ok: bool = True):
    """Create a directory.

    By default succeed if it already exists.

    Useful because transient NFS server issues may induce double creation attempts.
    """
    if exist_ok:
        os.makedirs(path, exist_ok=True)
    else:
        os.mkdir(path)


def remove(f: str | bytes, nonexistent_ok: bool = True):
    """Delete a file. By default succeed if it doesn't exist.

    Useful because transient NFS server issues may induce double deletion attempts.
    """
    if nonexistent_ok:
        try:
            os.remove(f)
        except OSError as e:
            if e.errno == errno.ENOENT:
                pass
            else:
                raise
    else:
        os.remove(f)


def hardlink_with_fallback(src: AnyStr, dst: AnyStr):
    """Hard-links src to dst, falling back to copy if it fails.

    If `src` is a directory, it will attempt to recursively hardlink.
    """

    def _safe_copy(src: AnyStr, dst: AnyStr):
        """Copy a file, like shutil.copy2, but catch errors from copystat."""
        try:
            shutil.copyfile(src, dst, follow_symlinks=False)
        except shutil.SameFileError:
            return
        try:
            shutil.copystat(src, dst, follow_symlinks=False)
        except OSError as ex:
            if ex.errno == errno.EPERM:
                # Some filesystems don't allow changing permissions.
                pass
            else:
                raise

    def _hardlink_file_with_fallback(src: AnyStr, dst: AnyStr):
        """Hardlink a file, fallback to copy if fail."""
        try:
            os.link(src, dst)
        except OSError as ex:
            if ex.errno in [errno.EPERM, errno.EOPNOTSUPP, errno.EXDEV, errno.EMLINK, errno.EEXIST]:
                _safe_copy(src, dst)
            else:
                raise

    if os.path.isdir(src):
        # recursively hardlink a path, fallback to copy if fail
        shutil.copytree(src, dst, copy_function=_hardlink_file_with_fallback, dirs_exist_ok=True)
    else:
        if os.path.isdir(dst):
            dst = os.path.join(dst, os.path.basename(src))
        _hardlink_file_with_fallback(src, dst)


def hard_link(f: AnyStr, relative_path: str | bytes | None = None) -> AnyStr | None:
    """Make a new hard link in a stage directory to the file f, defaulting to the basename of f."""
    if not f:
        return None

    if relative_path is None:
        relative_path = os.path.basename(f)
    assert relative_path

    new_path = martian.make_path(relative_path)
    if isinstance(f, str):
        new_path = new_path.decode("utf8")
    hardlink_with_fallback(f, new_path)

    return new_path


def concatenate_files(
    out_path: str | bytes, in_paths: Iterable[str | bytes], mode: str = ""
) -> None:
    with open(out_path, "w" + mode) as out_file:
        for in_path in in_paths:
            with open(in_path, "r" + mode) as in_file:
                shutil.copyfileobj(in_file, out_file)


def concatenate_headered_files(
    out_path: str | bytes, in_paths: Sequence[str | bytes], mode: str = ""
) -> None:
    """Concatenate files, taking the first line of the first file.

    and skipping the first line for subsequent files.
    Asserts that all header lines are equal.
    """
    with open(out_path, "w" + mode) as out_file:
        if len(in_paths) > 0:
            # Write first file
            with open(in_paths[0], "r" + mode) as in_file:
                header = in_file.readline()
                out_file.write(header)
                shutil.copyfileobj(in_file, out_file)

            # Write remaining files
            for in_path in in_paths[1:]:
                with open(in_path, "r" + mode) as in_file:
                    this_header = in_file.readline()
                    assert this_header == header
                    shutil.copyfileobj(in_file, out_file)


def write_empty_json(filename: str | bytes) -> None:
    with open(filename, "wb") as f:
        f.write(b"{}")


def touch(path: str | bytes) -> None:
    """Create an empty file."""
    with open(path, "wb"):
        pass


def recursive_hard_link_dict(in_files, prefixes=None):
    """Hard link files into this stage directory.

    For a dict with type [String,path_or_file_or_dict], where the keys represent sample IDs
    or other levels of nesting, create a dict with the same keys where all the path_or_file
    values are hardlinked into this stage directory from their old paths.
    If path_or_file_or_dict is a directory, recursively hardlink the contents of the directory.
    When nesting the key used to access each upper level will be used as a prefix
    when constructing the final file name.

    For example,

      {
        "sample1": {
          "gene_expression": "/mydir/gex.txt",
          "antibody": "/mydir/ab.txt"
        },
        "sample2": {
          "gene_expression": "/mydir/gex.txt",
          "antibody": "/mydir/ab.txt"
        },
      }

    would become

      {
        "sample1": {
          "gene_expression": "sample1_gene_expression_gex.txt",
          "antibody": "sample1_antibody.ab.txt"
         },
         "sample2": {
           "gene_expression": "sample2_gene_expression_gex.txt",
           "antibody": "sample2_antibody_ab.txt"
         },
      }
    """
    if in_files is None:
        return None

    if prefixes is None:
        prefixes = []

    out_files = {}
    for k, path_or_dict in in_files.items():
        if path_or_dict is None:
            out_files[k] = None
        elif isinstance(path_or_dict, dict):
            out_files[k] = recursive_hard_link_dict(in_files[k], prefixes + [k])
        elif isinstance(path_or_dict, (str)):
            final_prefixes = prefixes + [k]
            old_path = path_or_dict
            new_path = martian.make_path(
                "_".join(final_prefixes) + "_" + os.path.basename(old_path)
            ).decode("utf8")
            hardlink_with_fallback(old_path, new_path)
            out_files[k] = new_path
        else:
            raise ValueError(
                f"Input dictionary may not contain any elements other than dict and string: {path_or_dict}"
            )
    return out_files
