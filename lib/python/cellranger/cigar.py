#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
"""Code refactored out of align.py to avoid pulling in LibSSW.

Contains utilites that operate on CIGAR strings
"""
from __future__ import annotations

import re
from collections.abc import Iterator
from typing import TYPE_CHECKING

from six import ensure_binary

if TYPE_CHECKING:
    from striped_smith_waterman.ssw_wrap import PyAlignRes


def get_cigar_tuples(cigar_string: bytes) -> Iterator[tuple[int, bytes]]:
    """Get a list of length, CIGAR operation tuples for alignment.

    Args:
        cigar_string (bytes): CIGAR string

    Returns:
        tuple of (int, bytes): tuple of (length CIGAR operation, CIGAR operation).
    """
    cigar_string = ensure_binary(cigar_string)
    cigar_numbers = re.split(b"[A-Za-z]+", cigar_string)[:-1]
    cigar_letters = re.split(b"[0-9]+", cigar_string)[1:]
    return zip((int(number) for number in cigar_numbers), cigar_letters)


def get_max_word_length(alignment: PyAlignRes) -> int:
    """Get the longest M (match or mismatch) operation in a list of cigar tuples."""
    assert alignment.cigar_string is not None
    cigar_tuples = get_cigar_tuples(alignment.cigar_string)
    word_lengths = (op_len for op_len, op in cigar_tuples if op == b"M")

    return max(word_lengths, default=0)
