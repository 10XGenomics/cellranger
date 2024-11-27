#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

"""Assorted grab-bag of miscellaneous helper methods.

Do not add to this module.  Instead, find or create a module with a name
that indicates to a potential user what sorts of methods they might find
in that module.
"""

from __future__ import annotations

import json
import os
from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np

import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
from cellranger.targeted.simple_utils import load_target_csv_metadata

if TYPE_CHECKING:
    from pysam import AlignmentFile


def _load_reference_metadata_file(
    reference_path: str,
) -> dict[str, dict | list | str | int]:
    reference_metadata_file = os.path.join(reference_path, cr_constants.REFERENCE_METADATA_FILE)
    with open(reference_metadata_file) as f:
        return json.load(f)


def get_reference_star_path(reference_path: str) -> str:
    return os.path.join(reference_path, cr_constants.REFERENCE_STAR_PATH)


def get_reference_genome_fasta(reference_path: str) -> str:
    return os.path.join(reference_path, cr_constants.REFERENCE_FASTA_PATH)


def get_reference_genomes(
    reference_path: str | None, target_set_path: str | None = None
) -> list[str]:
    """Return the genome names from the reference transcriptome, or target set, or ["NONE"]."""
    if reference_path is not None:
        return _load_reference_metadata_file(reference_path)[cr_constants.REFERENCE_GENOMES_KEY]
    elif target_set_path is not None:
        return [load_target_csv_metadata(target_set_path, "probe set")["reference_genome"]]
    else:
        return ["NONE"]


def is_arc_reference(reference_path: str) -> bool:
    data = _load_reference_metadata_file(reference_path)
    return data.get("mkref_version", "").startswith("cellranger-arc")


def get_reference_mem_gb_request(reference_path: str) -> int:
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_MEM_GB_KEY]


def get_mem_gb_request_from_genome_fasta(reference_path: str) -> float:
    in_fasta_fn = get_reference_genome_fasta(reference_path)
    genome_size_gb = float(os.path.getsize(in_fasta_fn)) / 1e9
    return np.ceil(
        max(
            h5_constants.MIN_MEM_GB,
            cr_constants.BAM_CHUNK_SIZE_GB + max(1, 2 * int(genome_size_gb)),
        )
    )


def chunk_reference(
    input_bam: AlignmentFile, nchunks: int, include_unmapped: bool
) -> list[list[tuple[str, int, int]]]:
    """Chunk up a reference into nchunks roughly equally sized chunks.

    Args:
        input_bam (pysam.AlignmentFile): Source for reference contig names and lengths.
        nchunks (int): The number of chunks to create.
        include_unmapped (bool): If `True` then one of the chunks consists
                                 of all the unmapped reads as specified by
                                 `[("*", None, None)]`.

    Returns:
        list: loci, where each loci is a list of contiguous regions
              `(contig, start, end)`.
    """
    # one chunk is for unmapped reads
    nchunks -= int(include_unmapped)
    assert nchunks >= 1

    chroms = input_bam.references
    chrom_lengths = input_bam.lengths
    genome_size = sum(chrom_lengths)
    chunk_size = int(np.ceil(float(genome_size) / nchunks))
    process = [(chrom, 0, end) for chrom, end in zip(chroms, chrom_lengths)]

    chunks: list[list[tuple[str, int, int]]] = []
    if include_unmapped:
        chunks = [[("*", 0, 0)]]
    chunks.append([])

    current_chunk_size = 0
    while len(process):
        piece = process.pop()
        piece_size = piece[2] - piece[1]
        if current_chunk_size + piece_size < chunk_size:
            chunks[-1].append(piece)
            current_chunk_size += piece_size
        else:
            piece_left_over = current_chunk_size + piece_size - chunk_size
            new_piece = (piece[0], piece[1], piece[2] - piece_left_over)
            chunks[-1].append(new_piece)
            chunks.append([])
            current_chunk_size = 0
            if piece_left_over > 0:
                process.append((piece[0], piece[2] - piece_left_over, piece[2]))
    return [c for c in chunks if len(c)]


def get_ref_name_from_genomes(genomes: Iterable[str]):
    return "_and_".join(genomes)
