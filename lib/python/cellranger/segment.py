# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

"""Utility functions for working with aligned segments."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import tenkit.constants as tk_constants

if TYPE_CHECKING:
    from pysam import AlignedSegment


def pos_sort_key(read: AlignedSegment) -> tuple:
    return read.tid, read.pos


def get_read_barcode(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.PROCESSED_BARCODE_TAG)


def get_read_raw_barcode(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.RAW_BARCODE_TAG)


def get_read_barcode_qual(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.RAW_BARCODE_QUAL_TAG)


def get_read_umi_qual(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.UMI_QUAL_TAG)


def get_read_umi(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.PROCESSED_UMI_TAG)


def get_read_raw_umi(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.RAW_UMI_TAG)


def get_read_gene_ids(read: AlignedSegment):
    """_summary_.

    Args:
        read (AlignedSegment): _description_

    Returns:
        _type_: _description_
    """
    id_str: str | None = _get_read_tag(read, cr_constants.FEATURE_IDS_TAG)
    if id_str is None:
        return None
    assert isinstance(id_str, str)
    return tuple(id_str.split(";"))


def set_read_tags(read: AlignedSegment, tags) -> None:
    read_tags = read.tags
    read_tags.extend(tags)
    read.tags = read_tags


def _get_read_tag(read: AlignedSegment, tag):
    try:
        r = read.opt(tag)
        if not r:
            r = None
        return r
    except KeyError:
        return None


def get_read_extra_flags(read: AlignedSegment):
    return _get_read_tag(read, cr_constants.EXTRA_FLAGS_TAG) or 0


def get_genome_from_read(read: AlignedSegment, chroms, genomes):
    """_summary_.

    Args:
        read (_type_): _description_
        chroms (_type_): _description_
        genomes (_type_): _description_

    Returns:
        _type_: _description_
    """
    assert len(genomes) > 0

    if read.is_unmapped:
        return None

    if len(genomes) == 1:
        return genomes[0]

    return cr_utils.get_genome_from_str(chroms[read.tid], genomes)


# EXTRA_FLAGS_CONF_MAPPED_TXOME = 1
EXTRA_FLAGS_LOW_SUPPORT_UMI = 2
# EXTRA_FLAGS_GENE_DISCORDANT = 4
EXTRA_FLAGS_UMI_COUNT = 8
EXTRA_FLAGS_CONF_MAPPED_FEATURE = 16
EXTRA_FLAGS_FILTERED_TARGET_UMI = 32


def is_read_low_support_umi(read: AlignedSegment) -> bool:
    return (get_read_extra_flags(read) & EXTRA_FLAGS_LOW_SUPPORT_UMI) > 0


def is_read_filtered_target_umi(read: AlignedSegment) -> bool:
    return (get_read_extra_flags(read) & EXTRA_FLAGS_FILTERED_TARGET_UMI) > 0


def is_read_umi_count(read: AlignedSegment) -> bool:
    return (get_read_extra_flags(read) & EXTRA_FLAGS_UMI_COUNT) > 0


def is_read_conf_mapped_to_feature(read: AlignedSegment) -> bool:
    return (get_read_extra_flags(read) & EXTRA_FLAGS_CONF_MAPPED_FEATURE) > 0


def is_read_dupe_candidate(
    read, high_conf_mapq: int, use_corrected_umi: bool = True, use_umis: bool = True
) -> bool:
    """_summary_.

    Args:
        read: _description_
        high_conf_mapq: _description_
        use_corrected_umi: _description_. Defaults to True.
        use_umis: _description_. Defaults to True.

    Returns:
        bool: _description_
    """
    if use_corrected_umi:
        umi = get_read_umi(read)
    else:
        umi = get_read_raw_umi(read)

    return (
        not read.is_secondary
        and (umi is not None or not use_umis)
        and (get_read_barcode(read) is not None)
        and not is_read_low_support_umi(read)
        and not is_read_filtered_target_umi(read)
        and (
            is_read_conf_mapped_to_transcriptome(read, high_conf_mapq)
            or is_read_conf_mapped_to_feature(read)
        )
    )


def is_read_conf_mapped(read: AlignedSegment, high_conf_mapq: int) -> bool:
    """_summary_.

    Args:
        read (AlignedSegment): _description_
        high_conf_mapq (int): _description_

    Returns:
        bool: _description_
    """
    if read.is_unmapped:
        return False
    elif read.mapq < high_conf_mapq:
        return False
    return True


def is_read_conf_mapped_to_transcriptome(read: AlignedSegment, high_conf_mapq: int) -> bool:
    """_summary_.

    Args:
        read (AlignedSegment): _description_
        high_conf_mapq (int): _description_

    Returns:
        bool: _description_
    """
    if is_read_conf_mapped(read, high_conf_mapq):
        gene_ids = get_read_gene_ids(read)
        return (gene_ids is not None) and (len(gene_ids) == 1)
    return False


__CIGAR_CATEGORIES = b"MIDNSHP=X"


def _cigar_numeric_to_category_map(i: int) -> str:
    """Pysam numeric codes to meaningful categories."""
    return chr(__CIGAR_CATEGORIES[i])


def get_cigar_summary_stats(read: AlignedSegment, strand: bytes) -> dict[str, int]:
    """Get number of mismatches, insertions, deletions, ref skip, soft clip, hard clip bases from a read.

    Returns a dictionary by the element's CIGAR designation. Adds additional
    fields to distinguish between three and five prime soft-clipping for `R1`
    and `R2`: `R1_S_three_prime` and `R1_S_five_prime`, etc. to account for
    soft-clipped local alignments.

    Args:
        read (pysam.AlignedRead): aligned read object
        strand (string): + or - to indicate library orientation (MRO argument strand, for example)

    Returns:
        dict of str,int: Key of base type to base counts for metrics. Adds
                         additional fields to distinguish between three and
                         five prime soft-clipping: `S_three_prime` and
                         `S_five_prime`.
    """
    statistics: dict[str, int] = {}
    cigar_tuples = read.cigar

    for i, (category, count) in enumerate(cigar_tuples):
        # Convert numeric code to category
        category = _cigar_numeric_to_category_map(category)
        count = int(count)

        # Detect 5 prime soft-clipping
        if i == 0:
            if strand == cr_constants.REVERSE_STRAND:
                metric = "R1_S_five_prime" if read.is_read1 else "R2_S_three_prime"
            else:
                metric = "R2_S_five_prime" if read.is_read1 else "R1_S_three_prime"

            if category == "S":
                statistics[metric] = count
            else:
                statistics[metric] = 0

        # Tally up all standard categories from BAM
        statistics[category] = statistics.get(category, 0) + count

        # Detect 3 prime soft-clipping
        if i == len(cigar_tuples):
            if strand == cr_constants.REVERSE_STRAND:
                metric = "R2_S_five_prime" if read.is_read2 else "R1_S_three_prime"
            else:
                metric = "R1_S_five_prime" if read.is_read2 else "R2_S_three_prime"

            if category == "S":
                statistics[metric] = count
            else:
                statistics[metric] = 0

    return statistics


def get_full_alignment_base_quality_scores(
    read: AlignedSegment,
) -> tuple[int, np.ndarray[int, np.dtype[np.byte]]]:
    """Returns base quality scores for the full read alignment.

    Inserts zeroes for deletions and removing inserted and soft-clipped bases.
    Therefore, only returns quality for truly aligned sequenced bases.

    Args:
        read (pysam.AlignedSegment): read to get quality scores for

    Returns:
        np.array: numpy array of quality scores
    """
    # pylint: disable=no-member
    quality_scores = np.fromstring(read.qual, dtype=np.byte) - tk_constants.ILLUMINA_QUAL_OFFSET

    start_pos = 0

    for operation, length in read.cigar:
        operation = _cigar_numeric_to_category_map(operation)

        if operation == "D":
            quality_scores = np.insert(quality_scores, start_pos, [0] * length)
        elif operation in ("I", "S"):
            quality_scores = np.delete(quality_scores, np.s_[start_pos : start_pos + length])

        if not operation in ("I", "S"):
            start_pos += length

    return start_pos, quality_scores
