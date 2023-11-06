# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

from __future__ import annotations

import os.path
from typing import TYPE_CHECKING

from common.pyfasta import FastaIndexed

if TYPE_CHECKING:
    from pyfaidx import FastaRecord

# 10X generated references that should
# have genes, regions, and snps files
KNOWN_GENOMES: list[str] = ["10X_hg19_ucsc", "10X_b37", "10X_GRCh38_no_alt_decoy"]


def open_reference(reference_path: str) -> dict[str, FastaRecord]:
    """Open a reference fasta and rename the contigs to strip any fasta comments."""
    fasta = FastaIndexed(get_fasta(reference_path))
    new_fasta = {}

    for k, value in fasta.items():
        key_prefix = k.split(maxsplit=1)[0]
        new_fasta[key_prefix] = value

    return new_fasta


def get_fasta(reference_path: str) -> str:
    """Convention for location of reference fasta in a reference path."""
    return os.path.join(reference_path, "fasta", "genome.fa")


def is_tenx(reference_path: str | None) -> bool:
    return not reference_path is None and get_genome(reference_path) in KNOWN_GENOMES


def get_genome(reference_path: str) -> str | None:
    """Load the canonical name of a reference.

    By convention this is stored in <reference_path>/genome
    """
    genome = None

    genome_name_file = os.path.join(reference_path, "genome")

    if os.path.exists(genome_name_file):
        with open(genome_name_file) as f:
            genome = f.readline().strip()

    return genome
