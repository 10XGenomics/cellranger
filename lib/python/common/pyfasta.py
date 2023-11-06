# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

"""Interface for the pyfasta.Fasta API that uses pyfaidx under the hood."""


from __future__ import annotations

from typing import NamedTuple

import pyfaidx


class BedCoord(NamedTuple):
    """Coordinates in BED format."""

    chrom: str
    start: int
    end: int
    name: str
    strand: str


class FastaIndexed(pyfaidx.Fasta):
    """A wrapper around pyfaidx.Fasta that implements a `pyfasta.Fasta` interface.

    This helps upgrade old code that was dependent on `pyfasta.Fasta`. Any new code should directly use pyfaidx.

    .. deprecated::4.0
       Do not use for new code.  Use `pyfaidx` instead.
    """

    def __init__(self, fasta_name, build_index=False) -> None:
        """Create a pyfaidx.Fasta object from a fasta file. Note that the header after the first.

        " " is ignored. And we use 0-based half-open intervals for indexing.
        """
        super().__init__(filename=fasta_name, one_based_attributes=False, build_index=build_index)

    def get_sequence_region(
        self,
        bed_coord: BedCoord,
    ) -> pyfaidx.Sequence:
        """Get sequence from chrom, [start, end) and reverse complement if strand.

        (either `+` or `-` or `.`) is `-`
        """
        assert bed_coord.strand in ("-", "+", "."), f"Invalid strand = {bed_coord.strand}"
        return self.get_seq(
            name=bed_coord.chrom,
            start=bed_coord.start + 1,
            end=bed_coord.end,
            rc=(bed_coord.strand == "-"),
        )
