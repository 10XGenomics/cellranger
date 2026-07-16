#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# General utilities for manipulating nucleotide sequences
#

from __future__ import annotations

NUCS: list[bytes] = [b"A", b"C", b"G", b"T"]
NUCS_INVERSE: dict[bytes, int] = {b"A": 0, b"C": 1, b"G": 2, b"T": 3}


DNA_CONVERT_TABLE: bytes = bytes.maketrans(b"ACGTacgtRYMKBDHVrymkbdhv", b"TGCAtgcaYRKMVHDByrkmvhdb")
RNA_CONVERT_TABLE: bytes = bytes.maketrans(b"ACGUacguRYMKBDHVrymkbdhv", b"UGCAugcaYRKMVHDByrkmvhdb")


def get_rev_comp(seq: str | bytes) -> bytes:
    """Reverse complement for DNA.

    Included ambiguous nucleotides and retains case.
    """
    if isinstance(seq, str):
        seq = seq.encode("ascii")
    return seq.translate(DNA_CONVERT_TABLE)[::-1]


def hamming_distance(a: str | bytes, b: str | bytes) -> int:
    """Compute the Hamming distance between two equal-length sequences."""
    if isinstance(a, bytes):
        a = a.decode("ascii")
    if isinstance(b, bytes):
        b = b.decode("ascii")
    if len(a) != len(b):
        raise ValueError(f"Sequences must be equal length, got {len(a)} and {len(b)}")
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def is_near_palindrome(seq: str | bytes, max_distance: int = 1) -> bool:
    """Check if a sequence is within max_distance of its own reverse complement.

    Near-palindrome sequences cannot be used to distinguish between forward and
    reverse-complement orientations of index reads, since they look the same
    (or nearly the same) in either direction.
    """
    rc = get_rev_comp(seq).decode("ascii")
    if isinstance(seq, bytes):
        seq = seq.decode("ascii")
    return hamming_distance(seq, rc) <= max_distance


def mask(seq: bytes, keep_start: int, keep_end: int) -> bytes:
    """Mask the sequence leaving only [keep_start, keep_end) unmasked."""
    return b"N" * keep_start + seq[keep_start:keep_end] + b"N" * (len(seq) - keep_end)
