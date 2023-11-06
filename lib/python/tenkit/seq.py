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


def mask(seq: bytes, keep_start: int, keep_end: int) -> bytes:
    """Mask the sequence leaving only [keep_start, keep_end) unmasked."""
    return b"N" * keep_start + seq[keep_start:keep_end] + b"N" * (len(seq) - keep_end)
