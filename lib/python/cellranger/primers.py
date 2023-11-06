#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Primer-related helper functions for PD stages."""


from __future__ import annotations

from collections.abc import Iterable
from typing import AnyStr, NamedTuple

from six import ensure_binary

import tenkit.seq as tk_seq


class Primer(NamedTuple):
    name: str
    seq: str


def get_comp(seq: AnyStr) -> bytes:
    return ensure_binary(seq).translate(tk_seq.DNA_CONVERT_TABLE)


def get_rev(seq: AnyStr) -> bytes:
    return ensure_binary(seq)[::-1]


def get_primers_from_dicts(dicts: Iterable[dict[str, str]] | None) -> list[Primer]:
    return [Primer(d["name"], d["seq"]) for d in dicts] if dicts is not None else []


def get_primer_orientation_names(primer: Primer) -> list[str]:
    """Get orientation-related primer metric key names.

    Args:
        primer (Primer): Primer object
    """

    def _get_comp_primer_name(primer: Primer) -> str:
        return primer.name + "_comp"

    def _get_rev_primer_name(primer: Primer) -> str:
        return primer.name + "_rev"

    def _get_rev_comp_primer_name(primer: Primer) -> str:
        return primer.name + "_revcomp"

    assert isinstance(primer.name, str)
    if is_homopolymer_seq(primer.seq):
        return [primer.name]
    else:
        return [
            primer.name,
            _get_comp_primer_name(primer),
            _get_rev_primer_name(primer),
            _get_rev_comp_primer_name(primer),
        ]


def get_primer_orientation_seqs(primer):
    """Generate all combinations of primer sequences that we want to align.

    Args:
        primer (Primer): Primer object

    Returns:
        List[bytes]: List of primer sequences to align
    """
    assert isinstance(primer.name, str)
    seq = ensure_binary(primer.seq)
    if is_homopolymer_seq(seq):
        return [seq]
    else:
        return [
            seq,
            get_comp(seq),
            get_rev(seq),
            tk_seq.get_rev_comp(seq),
        ]


def is_homopolymer_seq(seq):
    return len(set(seq)) == 1
