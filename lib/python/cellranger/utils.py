# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

"""Assorted grab-bag of miscellaneous helper methods.

Do not add to this module.  Instead, find or create a module with a name
that indicates to a potential user what sorts of methods they might find
in that module.
"""

from __future__ import annotations

import json
import os
from collections.abc import Collection, Generator, Iterable, Sequence
from typing import IO, AnyStr, overload

import numpy as np

import cellranger.constants as cr_constants
import tenkit.seq as tk_seq
from cellranger.fast_utils import (  # pylint: disable=no-name-in-module,unused-import
    FilteredBarcodes,
)


def get_gem_group_from_barcode(barcode: bytes | None) -> int | None:
    """Method to get just the gem group from a barcode.

    Args:
        barcode (bytes): a barcode, a la "ACGTACTAGAC-1"

    Returns:
        The Gem group as an int or None if not present.
    """
    if barcode is None:
        return None
    gem_group = barcode.index(b"-")
    if gem_group > 0:
        return int(barcode[gem_group + 1 :])
    else:
        return None


def format_barcode_seq(barcode: bytes, gem_group: int | None = None) -> bytes:
    if gem_group is not None:
        barcode += b"-%d" % gem_group
    return barcode


def format_barcode_seqs(
    barcode_seqs: Collection[bytes], gem_groups: Iterable[int] | None
) -> list[bytes]:
    """Format a sequence of barcodes as seqs.

    Args:
        barcode_seqs (list[bytes]): _description_
        gem_groups (Optional[Iterable[int]]): _description_

    Returns:
        list[bytes]: _description_
    """
    if gem_groups is None:
        return barcode_seqs
    return [format_barcode_seq(bc, gg) for gg in sorted(set(gem_groups)) for bc in barcode_seqs]


@overload
def split_barcode_seq(
    barcode: None,
) -> tuple[None, None]: ...


@overload
def split_barcode_seq(
    barcode: bytes,
) -> tuple[bytes, int | None]: ...


def split_barcode_seq(barcode: bytes | None):
    """Split a barcode-gem_group.

    Args:
        barcode (bytes): A barcode with an optional gem group suffix.

    Returns:
        _type_: _description_
    """
    if barcode is None:
        return None, None

    assert isinstance(barcode, bytes)
    barcode_parts = barcode.split(b"-", maxsplit=2)

    barcode = barcode_parts[0]
    if len(barcode_parts) > 1:
        gem_group = int(barcode_parts[1])
    else:
        gem_group = None

    return barcode, gem_group


def bcs_suffices_to_names(
    bcs: Iterable[bytes], gg_id_map: dict[int, str | bytes]
) -> list[str | bytes]:
    """Turn list of barcodes into corresponding list of mapped values.

    Args:
        bcs (list): aggr barcodes with suffix corresponding to gem group id
        gg_id_map (dict): mapping gem-group-ids/barcode-suffixes (int) to
            a desired named value (typically `aggr_id` for the library id of that
            one sample, or `batch_id` for the name of a batch it is a part of.)
    """
    mapped_ids = [gg_id_map[split_barcode_seq(bc)[1]] for bc in bcs]
    return mapped_ids


def get_genome_from_str(search: AnyStr | None, genomes: Sequence[AnyStr]) -> AnyStr | None:
    """Return the first genome which is a prefix for the given search string.

    Args:
        search: _description_
        genomes: _description_

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    assert len(genomes) > 0

    if search is None:
        return None

    if len(genomes) == 1:
        return genomes[0]

    for genome in genomes:
        if search.startswith(genome):
            return genome

    raise ValueError(f"{search} does not have valid associated genome")


def remove_genome_from_str(
    seq: AnyStr, genomes: Sequence[AnyStr], prefixes_as_genomes: bool = False
) -> AnyStr | None:
    """_summary_.

    Args:
        seq (_type_): _description_
        genomes (_type_): _description_
        prefixes_as_genomes (bool, optional): _description_. Defaults to False.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    assert len(genomes) > 0
    if seq is None:
        return None
    if len(genomes) == 1:
        return seq

    # Gene names/ids/chroms are padded with N underscores to achieve the same prefix length
    #   for all genomes, e.g., GRCh38_* and mm10___*
    max_len = max(len(g) for g in genomes)

    for genome in genomes:
        if seq.startswith(genome):
            # Strip genome and subsequent underscores
            if prefixes_as_genomes:
                return seq[(1 + len(genome)) :]
            else:
                return seq[(1 + max_len) :]

    raise ValueError(f"{seq} does not have valid associated genome")


def get_read_transcripts_iter(
    read: bytes | None,
) -> Generator[
    (tuple[bytes | None, bytes, int, bytes | None] | tuple[bytes | None, None, None, bytes | None]),
    None,
    None,
]:
    """Iterate over all transcripts compatible with the given read.

    We do this by iterating over the TX tag entries that are of the form
    `TX:<transcript id>,<strand><position>,<cigar>`.

    Note:
        When intronic alignment is turned on, the `TX` tag is set to
        `TX:<gene id>,<strand>` for intronic  reads or exonic reads that are
        not compatible with annotated splice junctions. We ignore these.
    """
    # Method used to work directly with BAM records, since moved to rust/pyo3
    # s = _get_read_tag(read, cr_constants.TRANSCRIPTS_TAG)
    if not read:
        return

    for x in read.split(b";"):
        if len(x) == 0:
            continue

        parts = x.split(b",")
        # len(parts) = 3 for transcriptomic alignments, while for intronic alignments when
        # intron counting mode is enabled len(parts)=2
        if len(parts) != 3:
            continue

        chrom = parts[0] if parts[0] else None
        if parts[1]:
            strand = parts[1][0:1]
            pos = int(parts[1][1:])
        else:
            strand = None
            pos = None
        cigarstring = parts[2] if parts[2] else None
        assert strand in cr_constants.STRANDS or strand is None

        yield chrom, strand, pos, cigarstring


def compress_seq(seq: bytes, bits: int = 64) -> int:
    """Pack a DNA sequence (no Ns!) into a 2-bit format, into a python int.

    Args:
        seq (str): A DNA sequence.
        bits: the bit size of the integer to pack into.

    Returns:
        int: The sequence packed into the bits of an integer.
    """
    assert len(seq) <= (bits // 2 - 1)
    result = 0
    for i in range(len(seq)):
        nuc = seq[i : i + 1]
        assert nuc in tk_seq.NUCS_INVERSE
        result = (result << 2) | tk_seq.NUCS_INVERSE[nuc]
    return result


def merge_jsons_as_dict(in_filenames: Iterable[str | bytes | os.PathLike]) -> dict:
    """Merge a list of json files and return the result as a dictionary."""
    result = {}
    for filename in in_filenames:
        if filename is None:
            continue
        try:
            with open(filename) as f:
                data = json.load(f)
                result.update(data)
        except OSError:
            continue
    return result


def format_barcode_summary_h5_key(
    library_prefix: str, genome: str, region: str, read_type: str
) -> str:
    """Formats the barcode summary into a key to access the `barcode_summary.h5` outs.

    Here, we need to accommodate both accessing old and new versions of the h5 file
    compatible with both v2 and v3 of the matrix.

    Args:
        library_prefix: Name of the library, used as a prefix
        genome: Name of the genome used for alignment
        region: Name of the subset of the genome we are looking at (e.g. transcriptome, regulome,
            epigenome, ...). This should be a controlled vocabulary in `cellranger.constants`
        read_type: Name of the read types we are trying to extract (e.g. `conf_mapped_deduped_barcoded`, ...).
            It should be also controlled in `cellranger.constants`.

    Returns:
        output_key: A string constant with the suffix `reads` appended.
    """
    # Clean the input for non-found or undefined keys. Will break badly if they are not controlled, but will let us
    # keep compatible across v2 and v3 developments.
    str_list = [library_prefix, genome, region, read_type]
    # Append reads
    str_list.append("reads")
    # join
    output_key = "_".join(str_list)
    return output_key


def numpy_groupby(values, keys) -> Generator[tuple[tuple, tuple], None, None]:
    """Group a collection of numpy arrays by key arrays.

    Yields `(key_tuple, view_tuple)` where `key_tuple` is the key grouped
    on and `view_tuple` is a tuple of views into the value arrays.

    Args:
        values (tuple of arrays): tuple of arrays to group.
        keys (tuple): tuple of sorted, numeric arrays to group by.

    Returns:
        sequence of tuple: Sequence of (`key_tuple`, `view_tuple`).
    """
    if len(values) == 0:
        return
    if len(values[0]) == 0:
        return

    for key_array in keys:
        assert len(key_array) == len(keys[0])
    for value_array in values:
        assert len(value_array) == len(keys[0])

    # The indices where any of the keys differ from the previous key become group boundaries
    # pylint: disable=no-member
    key_change_indices = np.logical_or.reduce(
        tuple(np.concatenate(([1], np.diff(key))) != 0 for key in keys)
    )
    group_starts = np.flatnonzero(key_change_indices)
    group_ends = np.roll(group_starts, -1)
    group_ends[-1] = len(keys[0])

    for group_start, group_end in zip(group_starts, group_ends):
        yield tuple(key[group_start] for key in keys), tuple(
            value[group_start:group_end] for value in values
        )


def get_fasta_iter(f: IO[bytes]) -> Generator[tuple[bytes, bytes], None, None]:
    """Iterate through sequences in a fasta file.

    Args:
        f (IO[bytes]): The input file object.

    Yields:
        tuple[bytes, bytes]: _description_
    """
    hdr = b""
    seq = b""
    for line in f:
        line = line.strip()
        if line.startswith(b">"):
            if hdr:
                yield hdr, seq
            hdr = line[1:]
            seq = b""
        else:
            seq += line
    if hdr:
        yield hdr, seq


def load_barcode_csv(barcode_csv: str | bytes | os.PathLike) -> dict[bytes, list[bytes]]:
    """Load a csv file of (genome,barcode)."""
    if isinstance(barcode_csv, bytes):
        barcode_csv = barcode_csv.decode()
    return {
        genome.encode(): barcodes
        for genome, barcodes in FilteredBarcodes(barcode_csv).per_genome_barcodes().items()
    }


def get_cell_associated_barcode_set(
    barcode_csv_filename: str | bytes | os.PathLike, genome: bytes | None = None
) -> set[bytes]:
    """Get set of cell-associated barcode strings.

    Args:
      barcode_csv_filename: TODO
      genome (bytes): Only get cell-assoc barcodes for this genome. If None, disregard genome.

    Returns:
      set of bytes: Cell-associated barcode strings (seq and gem-group).
    """
    if isinstance(genome, str):
        genome = genome.encode()
    cell_bcs_per_genome: dict[bytes, list[bytes]] = load_barcode_csv(barcode_csv_filename)
    cell_bcs: set[bytes] = set()
    for g, bcs in cell_bcs_per_genome.items():
        if genome is None or g == genome:
            cell_bcs.update(bcs)
    return cell_bcs


def string_is_ascii(input_string: bytes | str) -> bool:
    """Returns true if the string can be encoded as ascii.

    Input strings are often stored as ascii in numpy arrays, and we need
    to check that this conversion works.
    """
    if isinstance(input_string, bytes):
        try:
            input_string = input_string.decode(errors="strict")
        except UnicodeDecodeError:
            return False
    if isinstance(input_string, str):
        try:
            input_string.encode("ascii", errors="strict")
        except UnicodeEncodeError:
            return False
        return True
    # Not str or bytes
    return False
