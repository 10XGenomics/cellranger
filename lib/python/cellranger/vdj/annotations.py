#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
"""Utilities for annotated contigs with gene/chain information, defining clonotypes and more."""

from __future__ import annotations

import itertools
import json
import re
from collections import defaultdict
from collections.abc import Callable, Iterable, Sequence
from typing import TYPE_CHECKING, Any, TextIO, TypedDict, overload

import numpy as np
from six import ensure_str

import cellranger.cigar as cr_cigar
import cellranger.vdj.reference as vdj_reference
import cellranger.vdj.utils as vdj_utils
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
from cellranger.library_constants import MULTI_REFS_PREFIX
from cellranger.vdj.constants import (
    AMBIGUOUS_AA_CODE,
    CODON_TO_AA,
    START_CODONS,
    STOP_CODONS,
    VDJ_5U_FEATURE_TYPES,
    VDJ_ANNOTATION_MATCH_SCORE,
    VDJ_ANNOTATION_MIN_V_OVERLAP_FRAC,
    VDJ_C_FEATURE_TYPES,
    VDJ_CDR3_ALL_END_MOTIFS,
    VDJ_CDR3_COMMON_END_MOTIFS,
    VDJ_CLONOTYPE_TYPES,
    VDJ_D_FEATURE_TYPES,
    VDJ_GENE_PAIRS,
    VDJ_J_FEATURE_TYPES,
    VDJ_MAX_CDR3_LEN,
    VDJ_MIN_CDR3_LEN,
    VDJ_ORDERED_REGIONS,
    VDJ_QUAL_OFFSET,
    VDJ_V_FEATURE_TYPES,
)

# Allow the start codon to shift up or down by this many codons
from tenkit import seq as tk_seq

if TYPE_CHECKING:
    from cellranger.align import SSWAlignmentResult, SSWMultiAligner

START_CODON_SLOP = 1


def codon_to_aa(codon: str | bytes) -> int:
    """Return amino acid corresponding to a codon.

    If the codon is not in the translation table, a default
    AA is returned (see AMBIGUOUS_AA_CODE in vdj.constants)
    """
    assert len(codon) == 3
    if isinstance(codon, bytes):
        codon = codon.decode("ascii")
    assert all(c in "NACGT" for c in codon)
    code = CODON_TO_AA.get(codon, AMBIGUOUS_AA_CODE)
    assert len(code) == 1
    return ord(code)


def filter_alignment(
    alignment_result: SSWAlignmentResult,
    score_ratio: float,
    word_size: int,
    match_score: float = VDJ_ANNOTATION_MATCH_SCORE,
) -> bool:
    """Returns True for a passing alignment and False otherwise.

    Args:
        alignment_result (SSWAlignmentResult): alignment result to filter
        score_ratio: minimum (score / max_possible_score)
        word_size: minimum stretch of contiguous matches/mismatches in the alignment
        match_score: match score used in the alignment

    Returns:
        True if alignment passed filters
    """
    alignment = alignment_result.alignment
    assert alignment is not None

    # alignment is a PyAlignRes object. Ends are inclusive.
    alignment_length = float(alignment.query_end - alignment.query_begin + 1)
    max_score = alignment_length * match_score
    if tk_stats.robust_divide(alignment.score, max_score) < score_ratio:
        return False

    if cr_cigar.get_max_word_length(alignment) < word_size:
        return False

    return True


def filter_v_alignment(
    alignment_result: SSWAlignmentResult,
    score_ratio,
    match_score=VDJ_ANNOTATION_MATCH_SCORE,
    v_overlap_frac=VDJ_ANNOTATION_MIN_V_OVERLAP_FRAC,
):
    """Returns True for a passing alignment and False otherwise.

    Args:
        alignment_result (SSWAlignmentResult): alignment result to filter
        score_ratio: minimum (score / max_possible_score)
        match_score: match score used in the alignment
        v_overlap_frac: Minimum required V gene overlap (alignment length/V gene length)

    Returns:
        True if alignment passed filters
    """
    v_gene_length = len(alignment_result.reference.metadata["feature"].sequence)

    alignment = alignment_result.alignment
    assert alignment is not None

    # alignment is a PyAlignRes object. Ends are inclusive.
    alignment_length = float(alignment.query_end - alignment.query_begin + 1)

    max_score = alignment_length * match_score
    if tk_stats.robust_divide(alignment.score, max_score) < score_ratio:
        return False

    ref_alignment_length = float(alignment.ref_end - alignment.ref_begin + 1)
    min_alignment_length = v_overlap_frac * v_gene_length
    if ref_alignment_length < min_alignment_length:
        return False

    return True


def collect_annotations(
    aligner: SSWMultiAligner,
    ref_seq: bytes,
    seq: bytes,
    filter_func: Callable[[SSWAlignmentResult], bool],
) -> list[Annotation]:
    """Align a sequence against an SSWMultiAligner.

    Return a list of Annotation objects
    for the alignments that passed filter_func.

    Args:
        seq (bytes): the sequence to align.
        ref_seq (bytes): the sequence used to create the returned Annotation objects.
            ref_seq and seq will often be the same but seq could be a partially
            masked version of ref_seq. The two must have the same length.

    Returns:
        list: Annotation objects.
    """
    assert len(ref_seq) == len(seq)
    assert isinstance(ref_seq, bytes)
    assert isinstance(seq, bytes)
    alignment_results = aligner.align(seq)

    annotations = []

    for alignment_result in alignment_results:
        if alignment_result.alignment is None or not filter_func(alignment_result):
            continue

        # seq might be a masked version of ref_seq. The two have the same length, so the cigar strings
        # and match coordinates against seq will make sense against ref_seq as well.
        annotation = Annotation.from_alignment_result(alignment_result, ref_seq)
        annotations.append(annotation)
    return annotations


def find_cdr3_end_motif(j_amino_acids: bytes, allowed_end_motifs: list[str]) -> int | None:
    """Search for the CDR3 end motif in J amino acid sequence.

    Args:
        j_amino_acids: the AA of the J region
        allowed_end_motifs: prioritized list of allowed end motifs. See VDJ_CDR3_ALL_END_MOTIFS
            for the exhaustive list

    Returns:
        end_motif_pos: Index of the end motif in j_amino_acids
    """
    for motif in allowed_end_motifs:
        assert motif in VDJ_CDR3_ALL_END_MOTIFS
        for idx in range(len(j_amino_acids) - len(motif) + 1):
            if all(j_amino_acids[idx + i] == ord(aa) or aa == "X" for i, aa in enumerate(motif)):
                return idx
    return None


def search_cdr3_signature(
    v_region: Annotation, j_region: Annotation, v_frame: int
) -> tuple[tuple[int, int] | None, str]:
    """Search for the CDR3 signature in a sequence.

    Args:
      v_region: Annotation object for the V-Region
      j_region: Annotation object for the J-Region
      v_frame: reading frame (0, 1, 2)

    Returns:
       (a,b): the starting and ending (base) positions of the CDR3
       status_flag: a string indicating what happened within seq or None if no
                    CDR3 could be found.
    """
    flag = ""
    pos = None

    assert v_frame in [0, 1, 2]
    assert j_region.sequence is not None
    assert v_region.sequence is not None

    v_seq = v_region.sequence

    # Start of V within seq
    v_start = v_region.contig_match_start

    # Get the position of the last Cysteine within the V-Region
    v_amino_acids = bytes(
        codon_to_aa(v_seq[i : (i + 3)]) for i in range(v_frame, len(v_seq) - 2, 3)
    )
    last_c_idx = v_amino_acids.rfind(b"C")
    if last_c_idx < 0:
        flag = "GUIDED_NO_C_IN_V"
        return (pos, flag)

    # Translate to a base position within seq. The last Cysteine is
    # included in the CDR3.
    last_c_pos = v_start + v_frame + last_c_idx * 3

    j_seq = j_region.sequence
    j_start = j_region.contig_match_start

    # j_start is not necessarily in frame with the v_start (because of ambiguities
    # in annotations)
    # j_start - v_start - v_frame are the number of bases between the first codon of V
    # and the start of J.
    # Eg V starts at position 0 of seq, frame is 1 and J starts at position 10
    # Amino acids are (1-3), (4-6), (7-9), (10-12), ...
    # j_start - v_start - v_frame % 3 = 0, so we start the first J amino-acid
    # at position 0 of J.
    # If j_start was 11, then:
    # j_start - v_start - v_frame % 3 = 1. I.e. there's a "loose" base between the start
    # of J and the previous amino-acid. We need to skip 3-1 = 2 bases and start the
    # first amino acid of J at position 2 of J.

    j_frame: int = 3 - (j_start - v_start - v_frame) % 3
    if j_frame == 3:
        j_frame = 0
    j_amino_acids = bytes(
        bytearray(
            codon_to_aa(j_seq[i : (i + 3)]) for i in range(j_frame, len(j_region.sequence) - 2, 3)
        )
    )

    # Look for FG(X)G signature or WG(X)G signature
    end_motif_pos = None
    end_motif_idx = find_cdr3_end_motif(j_amino_acids, VDJ_CDR3_ALL_END_MOTIFS)
    if end_motif_idx is not None:
        # The CDR3 includes the first AA of the signature
        end_motif_pos = j_start + j_frame + end_motif_idx * 3 + 3

    if not end_motif_pos:
        flag = "GUIDED_NO_FGXG"

    if end_motif_pos and last_c_pos < end_motif_pos:
        return ((last_c_pos, end_motif_pos), flag)

    return (None, flag)


def search_cdr3_signature_no_vj(seq, v_region=None, j_region=None):
    """Search the CDR3 signature in a sequence without guides from annotations.

    This could lead to more false positive signature hits than the guided version.

    Return value:
    A tuple (CDR3 DNA seq, CDR3 amino-acid seq, start position in seq, end position)
    """
    min_cdr3_aas = VDJ_MIN_CDR3_LEN // 3

    valid_end = False
    cys_pos = 0
    for frame in range(3):
        amino_acids = bytes(
            bytearray(codon_to_aa(seq[i : (i + 3)]) for i in range(frame, len(seq) - 2, 3))
        )

        fgxg_idx = None

        for idx in range(min_cdr3_aas, len(amino_acids) - 3):
            # First try to find the end motif
            for motif in VDJ_CDR3_COMMON_END_MOTIFS:
                valid_end = all(
                    amino_acids[idx + i] == ord(aa) or aa == "X" for i, aa in enumerate(motif)
                )

            if valid_end:
                # The CDR3 includes the first F of the signature
                fgxg_idx = idx
                fgxg_pos = frame + fgxg_idx * 3
                if j_region and (
                    fgxg_pos >= j_region.contig_match_end or fgxg_pos < j_region.contig_match_start
                ):
                    continue

                cys_idx = None
                # Find the Cysteine closer to the end but past a minimum
                # number of amino acids
                for _ in range(fgxg_idx - min_cdr3_aas, 0, -1):
                    if amino_acids[idx] == ord(b"C"):
                        cys_idx = idx
                        cys_pos = frame + cys_idx * 3
                        # Make sure the C of the CDR3 is inside the V
                        if not v_region or (
                            cys_pos >= v_region.contig_match_start
                            and cys_pos < v_region.contig_match_end
                        ):
                            break
                    elif chr(amino_acids[idx]) in STOP_CODONS:
                        # TODO: I don't think this can actually happen.
                        break

                if cys_idx:
                    if fgxg_pos - cys_pos < VDJ_MAX_CDR3_LEN:
                        # include the start of the FGXG
                        cdr3_seq = seq[cys_pos : (fgxg_pos + 3)]
                        cdr3_aas = amino_acids[cys_idx : (fgxg_idx + 1)]
                        return (cdr3_seq, cdr3_aas, cys_pos, fgxg_pos)
    return None


@overload
def _bytes_or_none(maybe_str: None) -> None: ...


@overload
def _bytes_or_none(maybe_str: bytes | str) -> bytes: ...


def _bytes_or_none(maybe_str: bytes | str | None) -> bytes | None:
    if maybe_str is None:
        return None
    if isinstance(maybe_str, bytes):
        return maybe_str
    return maybe_str.encode()


class AnnotatedContigDict(TypedDict):
    barcode: bytes | str | None
    contig_name: Any
    sequence: str
    annotations: Iterable[dict[str, Any]]
    clonotype: Any
    is_cell: bool | None
    read_count: int | None
    umi_count: int | None
    info: dict


class AnnotatedContig:
    """A named sequence with a list of Annotations."""

    # Note: 'frame' is unused
    __slots__ = [
        "contig_name",
        "sequence",
        "quals",
        "annotations",
        "primer_annotations",
        "barcode",
        "clonotype",
        "read_count",
        "umi_count",
        "unannotated_intervals",
        "filtered",
        "cdr3_start",
        "cdr3_stop",
        "cdr3_seq",
        "cdr3",
        "cdr3_flag",
        "start_codon_pos",
        "stop_codon_pos",
        "frame",
        "aa_sequence",
        "productive",
        "info_dict",
        "is_cell",
        "high_confidence",
    ]

    def __init__(
        self,
        name,
        sequence: bytes,
        annotations: list[Annotation] = [],
        primer_annotations: list[Annotation] = [],
        barcode: bytes | None = None,
        is_cell: bool | None = None,
        clonotype=None,
        read_count: int | None = None,
        umi_count: int | None = None,
        unannotated_intervals=[],
        filtered: bool = True,
        info_dict={},
        quals=None,
        high_confidence=False,
    ):
        """The info_dict can contain additional features that are not included in the objects features."""
        self.contig_name = name
        self.sequence: bytes = sequence.upper()
        self.quals = quals
        self.annotations: list[Annotation] = annotations
        self.primer_annotations: list[Annotation] = primer_annotations
        self.barcode: bytes | None = barcode
        self.clonotype = clonotype
        self.read_count: int | None = read_count
        self.umi_count: int | None = umi_count
        self.unannotated_intervals = unannotated_intervals
        # True if contig was selected for downstream analysis
        self.filtered: bool = filtered
        self.cdr3_start: int | None = None
        self.cdr3_stop: int | None = None
        # sequence of CDR3 (self.sequence[cdr_pos[0]:cdr_pos[1]])
        self.cdr3_seq: bytes | None = None
        # CDR3 amino-acid sequence
        self.cdr3: bytes | None = None
        self.cdr3_flag = None
        self.start_codon_pos: int | None = None
        self.stop_codon_pos: int | None = None
        self.frame = None
        self.aa_sequence: bytes | None = None
        self.productive = None
        self.info_dict = info_dict
        # The cell-barcode was called as a cell
        self.is_cell: bool | None = is_cell
        self.high_confidence = high_confidence

    @classmethod
    def from_dict(
        cls,
        annotation_dict: AnnotatedContigDict,
        reference: vdj_reference.VdjReference,
    ):
        """Create an AnnotatedContig object from a dict.

        Does not perform any checks on the self-consistency of the dict features.
        (Eg. do the annotation sequences agree with the contig's sequence,
        do the CDR3 positions make sense etc).
        """
        new_obj = cls(
            name=annotation_dict["contig_name"],
            quals=annotation_dict.get("quals", None),
            sequence=_bytes_or_none(annotation_dict["sequence"]),
            annotations=[
                Annotation.from_dict(hit, reference, contig_seq=annotation_dict["sequence"])
                for hit in annotation_dict["annotations"]
            ],
            primer_annotations=[
                Annotation.from_dict(hit, reference, contig_seq=annotation_dict["sequence"])
                for hit in annotation_dict.get("primer_annotations", [])
            ],
            barcode=_bytes_or_none(annotation_dict["barcode"]),
            clonotype=annotation_dict["clonotype"],
            is_cell=annotation_dict["is_cell"],
            read_count=annotation_dict["read_count"],
            umi_count=annotation_dict["umi_count"],
            info_dict=annotation_dict["info"],
        )
        new_obj.cdr3_start = annotation_dict.get("cdr3_start", None)
        new_obj.cdr3_stop = annotation_dict.get("cdr3_stop", None)
        new_obj.productive = annotation_dict.get("productive", None)
        new_obj.start_codon_pos = annotation_dict.get("start_codon_pos", None)
        new_obj.stop_codon_pos = annotation_dict.get("stop_codon_pos", None)
        new_obj.frame = annotation_dict.get("frame", None)
        new_obj.aa_sequence = _bytes_or_none(annotation_dict.get("aa_sequence", None))
        # Doesn't check that these are consistent with cdr_start, cdr_stop, and sequence.
        new_obj.cdr3_seq = _bytes_or_none(annotation_dict.get("cdr3_seq", None))
        new_obj.cdr3 = annotation_dict.get("cdr3", None)
        new_obj.high_confidence = annotation_dict.get("high_confidence", None)

        return new_obj

    def to_dict(self) -> AnnotatedContigDict:
        return AnnotatedContigDict(
            contig_name=self.contig_name,
            sequence=self.sequence,
            quals=self.quals,
            annotations=[a.to_dict() for a in self.annotations],
            primer_annotations=[a.to_dict() for a in self.primer_annotations],
            barcode=self.barcode,
            clonotype=self.clonotype,
            is_cell=self.is_cell,
            filtered=self.filtered,
            read_count=self.read_count,
            umi_count=self.umi_count,
            cdr3_start=self.cdr3_start,
            cdr3_stop=self.cdr3_stop,
            cdr3_seq=self.cdr3_seq,
            cdr3=self.cdr3,
            start_codon_pos=self.start_codon_pos,
            stop_codon_pos=self.stop_codon_pos,
            frame=self.frame,
            aa_sequence=self.aa_sequence,
            productive=self.productive,
            info=self.info_dict,
            high_confidence=self.high_confidence,
        )

    def hits(self):
        """Iterator over the annotations of the contig."""
        return iter(self.annotations)

    def get_region_hits(self, region_types: Sequence[str]) -> list[Annotation]:
        for annotation in self.annotations:
            assert isinstance(annotation.feature.region_type, bytes)
        for r in region_types:
            assert isinstance(r, str)
        return [
            annotation
            for annotation in self.annotations
            if annotation.feature.region_type.decode() in region_types
        ]

    def get_single_chain(self):
        """Returns 'Multi' if ambiguous."""
        chains = self.contig_chains()
        if len(chains) == 0:
            return None
        elif self.is_single_chain():
            return next(iter(chains))
        else:
            return b"Multi"

    def get_single_gene_display_name(self, region_types):
        """Returns first if ambiguous."""
        for r in region_types:
            assert isinstance(r, str)
        hits = self.get_region_hits(region_types)
        if len(hits) == 0:
            return None
        else:
            return ensure_str(hits[0].feature.display_name)

    def is_single_chain(self):
        return len(self.contig_chains()) == 1

    def contig_chains(self) -> set[bytes]:
        return {annotation.feature.chain for annotation in self.annotations}

    def has_cdr(self):
        return not self.cdr3_seq is None and len(self.cdr3_seq) > 0

    def get_cdr_seq(self):
        """DNA sequence of CDR3."""
        if self.has_cdr():
            return self.cdr3_seq
        return None

    def get_cdr_aa_seq(self):
        """Amino acid sequence of CDR3."""
        if self.has_cdr():
            return self.cdr3
        return None

    def clonotype_seq(self) -> tuple[None, None] | tuple[bytes, bytes]:
        """Clonotype-defining sequence.

        Returns a tuple (chain, <chain_name>_<CDR3_sequence>). If this contig
        does not match a single chain, or does not have a CDR3, returns (None, None).
        """
        if not self.is_single_chain():
            return None, None
        if not self.has_cdr():
            return None, None
        chain = next(iter(self.contig_chains()))
        assert isinstance(chain, bytes), type(chain)
        seq = self.get_cdr_seq()
        assert isinstance(seq, bytes), type(seq)
        cdr_seq = chain + b"_" + seq
        return (chain, cdr_seq)

    def has_full_length_vj_hit(self):
        has_full_len_v_hit = any(
            ensure_str(annotation.feature.region_type) in VDJ_V_FEATURE_TYPES
            for annotation in self.annotations
        )

        # The -3 allows for slop on the 3' end of J in the face of alignment clipping
        has_full_len_j_hit = any(
            ensure_str(annotation.feature.region_type) in VDJ_J_FEATURE_TYPES
            and annotation.annotation_match_end >= annotation.annotation_length - 3
            for annotation in self.annotations
        )
        return has_full_len_v_hit and has_full_len_j_hit

    def spans_v_start(self):
        return any(
            ensure_str(annotation.feature.region_type) in VDJ_V_FEATURE_TYPES
            and annotation.annotation_match_start == 0
            for annotation in self.annotations
        )

    def get_vj_quals(self):
        if self.quals is None:
            return None
        v_regions = [
            annotation
            for annotation in self.annotations
            if ensure_str(annotation.feature.region_type) in VDJ_V_FEATURE_TYPES
            and annotation.annotation_match_start == 0
        ]
        if not v_regions:
            return None
        v_region = v_regions[0]

        j_regions = [
            annotation
            for annotation in self.annotations
            if ensure_str(annotation.feature.region_type) in VDJ_J_FEATURE_TYPES
            and annotation.annotation_match_end == annotation.annotation_length
        ]
        if not j_regions:
            return None
        j_region = j_regions[0]

        quals = np.array(
            [
                ord(q) - VDJ_QUAL_OFFSET
                for q in self.quals[v_region.contig_match_start : j_region.contig_match_end]
            ]
        )
        return quals

    def get_quals(self):
        if not self.quals is None:
            return np.array([ord(q) - VDJ_QUAL_OFFSET for q in self.quals])
        return None

    def get_concat_reference_sequence(self):
        """Return a concatenated reference sequence.

        Return value:
        - (None,None) if this contig isn't annotated with a V and a J segment.
        Otherwise a tuple (seqs, annos) where annos is a list of Annotation objects
        in the order that they should appear in a VDJ sequence and seqs is a list
        of corresponding sequences from the input fasta.
        """
        v_region = self.get_region_hits(VDJ_V_FEATURE_TYPES)
        j_region = self.get_region_hits(VDJ_J_FEATURE_TYPES)

        if not v_region or not j_region:
            return (None, None)

        seqs = []
        ordered_annos = []
        for region_defs in VDJ_ORDERED_REGIONS:
            regions = self.get_region_hits(region_defs)
            if regions:
                seqs.append(regions[0].feature.sequence)
                ordered_annos.append(regions[0])

        return (seqs, ordered_annos)

    def _get_masked_annotations(
        self,
        feature_types: Iterable[str],
        ssw_multi_aligners: Iterable[SSWMultiAligner],
        ssw_filter_funcs: Iterable[Callable[[SSWAlignmentResult], bool]],
    ):
        masked_sequence = self.sequence
        masked_annotations: list[Annotation] = []
        vj_hits: list[Annotation] = []  # Keep track of the V and J annotations

        # V, J, C will be annotated in that order. Skip UTRs and Ds, these will be handled
        # conditionally on the V, J, C annotations.
        for region_type, aligner, align_filter in zip(
            feature_types, ssw_multi_aligners, ssw_filter_funcs
        ):
            assert isinstance(region_type, str)
            if not region_type in VDJ_V_FEATURE_TYPES + VDJ_J_FEATURE_TYPES + VDJ_C_FEATURE_TYPES:
                continue

            # Align against the masked sequence, but create Annotation objects that track the
            # original sequence.
            annotations = collect_annotations(aligner, self.sequence, masked_sequence, align_filter)

            if annotations:
                # Give preference to V and C annotations that are intact at the 5' end. This
                # will effectively give preference to full length Vs/Cs that are
                # more consistent with the Vs/Js annotated so far.
                if region_type in VDJ_V_FEATURE_TYPES + VDJ_C_FEATURE_TYPES:
                    best = max(annotations, key=lambda x: (x.annotation_match_start == 0, x.score))
                else:
                    best = max(annotations, key=lambda x: x.score)

                masked_annotations.append(best)

                if region_type in VDJ_V_FEATURE_TYPES or region_type in VDJ_J_FEATURE_TYPES:
                    vj_hits.append(best)

                # Mask everything before the end of this annotation.
                # contig_match_end is inclusive.
                masked_sequence = tk_seq.mask(
                    self.sequence, best.contig_match_end, len(self.sequence)
                )
        return masked_annotations, vj_hits

    def annotate_features(
        self,
        feature_types: Sequence[str],
        ssw_multi_aligners: Sequence[SSWMultiAligner],
        ssw_filter_funcs: Sequence[Callable[[SSWAlignmentResult], bool]],
    ):
        """Add sequence annotations by comparing against a list of reference segments.

        Args:
            feature_types, ssw_multi_aligners, ssw_filter_funcs: Matched aligners and filters
                for a list of feature types (see setup_feature_aligners)

        Returns:
            A list of Annotation objects.
        """
        masked_annotations, vj_hits = self._get_masked_annotations(
            feature_types, ssw_multi_aligners, ssw_filter_funcs
        )

        # Now try to find a UTR alignment before the beginning of the first thing you've annotated
        utr_aligner = ssw_multi_aligners[feature_types.index(VDJ_5U_FEATURE_TYPES[0])]
        utr_filter = ssw_filter_funcs[feature_types.index(VDJ_5U_FEATURE_TYPES[0])]
        min_start = (
            np.min([anno.contig_match_start for anno in masked_annotations])
            if masked_annotations
            else len(self.sequence)
        )
        annotations = collect_annotations(
            utr_aligner, self.sequence, tk_seq.mask(self.sequence, 0, min_start), utr_filter
        )
        if annotations:
            # Keep best scoring one (doesn't make sense to priorite full length UTRs).
            best_utr = max(annotations, key=lambda x: x.score)
            # If you did find a V region, make sure the gene name of the UTR matches.
            # Otherwise, just pick the best scoring UTR.
            if vj_hits and ensure_str(vj_hits[0].feature.region_type) in VDJ_V_FEATURE_TYPES:
                v_gene = vj_hits[0].feature.display_name
                equally_good = [
                    anno
                    for anno in annotations
                    if anno.score == best_utr.score and anno.feature.display_name == v_gene
                ]
                if equally_good:
                    masked_annotations.append(equally_good[0])
            else:
                masked_annotations.append(best_utr)

        # Only search for D if V and J were found, in order to avoid false hits.
        if len(vj_hits) == 2:
            assert ensure_str(vj_hits[0].feature.region_type) in VDJ_V_FEATURE_TYPES
            assert ensure_str(vj_hits[1].feature.region_type) in VDJ_J_FEATURE_TYPES

            v_end = vj_hits[0].contig_match_end
            j_start = vj_hits[1].contig_match_start
            # Also require a D of the same chain-type. This will also make sure that we
            # don't add a D on a chain that's not supposed to have a D (like TRA).
            if vj_hits[0].feature.chain == vj_hits[1].feature.chain:
                d_aligner = ssw_multi_aligners[feature_types.index(VDJ_D_FEATURE_TYPES[0])]
                d_filter = ssw_filter_funcs[feature_types.index(VDJ_D_FEATURE_TYPES[0])]
                annotations = collect_annotations(
                    d_aligner, self.sequence, tk_seq.mask(self.sequence, v_end, j_start), d_filter
                )
                annotations = [
                    anno for anno in annotations if anno.feature.chain == vj_hits[0].feature.chain
                ]
                if annotations:
                    masked_annotations.append(max(annotations, key=lambda x: x.score))

        return masked_annotations

    def _annotation_from_alignment(
        self,
        alignment_result: SSWAlignmentResult,
        alignment_filter: Callable[[SSWAlignmentResult], bool],
    ):
        alignment = alignment_result.alignment

        if alignment is None:
            return None
        if not alignment_filter(alignment_result):
            return None

        feature = alignment_result.reference.metadata["feature"]

        return Annotation(
            feature=feature,
            cigar=alignment.cigar_string,
            score=alignment.score,
            annotation_length=len(feature.sequence),
            annotation_match_start=alignment.ref_begin,
            annotation_match_end=alignment.ref_end + 1,
            contig_match_start=alignment.query_begin,
            contig_match_end=alignment.query_end + 1,
            mismatches=[],
            contig_seq=self.sequence,
        )

    def annotate_features_by_group(
        self,
        ssw_multi_aligner: SSWMultiAligner,
        alignment_filter: Callable[[SSWAlignmentResult], bool],
    ) -> list[Annotation]:
        """Similar to annotate_features but uses a joint SSWMultiAligner object.

        ssw_multi_aligner can contain references of different types. The sequence
        will be aligned against all of these and then the best alignment by type
        will be returned. The same filter function will be used for all feature types.
        """
        alignment_results = ssw_multi_aligner.align(self.sequence)

        annotations: list[Annotation] = []

        for alignment_result in alignment_results:
            if (
                annotation := self._annotation_from_alignment(alignment_result, alignment_filter)
            ) is not None:
                annotations.append(annotation)

        # Take top scoring hit for each region type
        top_annotations: list[Annotation] = []

        if top_annotations:

            def _by_region(x: Annotation):
                return x.feature.region_type

            for _, group in itertools.groupby(sorted(annotations, key=_by_region), _by_region):
                top_annotations.append(max(group, key=lambda x: x.score))

        return top_annotations

    def get_unannotated_intervals(self):
        """Return a list of Annotation objects corresponding to unannotated regions on the contig."""
        unannotated_intervals = []
        annotation_indicator = np.zeros(len(self.sequence))

        for annotation in self.annotations:
            annotation_indicator[annotation.contig_match_start : annotation.contig_match_end] = 1

        interval_start = 0

        for annotated, region_iter in itertools.groupby(annotation_indicator, lambda x: x == 1):
            region = list(region_iter)

            if not annotated:
                feature = vdj_reference.create_dummy_feature(
                    display_name=b"UNANNOTATED", region_type=b"UNANNOTATED", sequence=None
                )
                unannotated_intervals.append(
                    Annotation(
                        feature=feature,
                        cigar=None,
                        score=0,
                        annotation_length=len(region),
                        annotation_match_start=0,
                        annotation_match_end=len(region),
                        contig_match_start=interval_start,
                        contig_match_end=interval_start + len(region),
                        mismatches=[],
                    )
                )

            interval_start += len(region)

        return unannotated_intervals

    def get_annotations_bed(self):
        """Convert this contig's annotations to a list of BED entries.

        Yields:
            str: each element is a BED file line (no newline)
        """
        for anno_list in (self.annotations, self.unannotated_intervals):
            for anno in anno_list:
                yield "\t".join(
                    [
                        ensure_str(self.contig_name),
                        str(anno.contig_match_start),
                        str(anno.contig_match_end),
                        f"{ensure_str(anno.feature.display_name)}_{ensure_str(anno.feature.region_type)}",
                    ]
                )

    def annotate_cdr3(self):
        v_regions = self.get_region_hits(VDJ_V_FEATURE_TYPES)
        j_regions = self.get_region_hits(VDJ_J_FEATURE_TYPES)

        self.cdr3_start = None
        self.cdr3_stop = None
        self.cdr3_seq = None
        self.cdr3 = None
        self.cdr3_flag = None
        self.frame = None
        self.start_codon_pos = None
        self.stop_codon_pos = None
        self.productive = None
        self.aa_sequence = None

        flags = []

        if len(v_regions) == 1 and len(j_regions) == 1:
            # First try to search in an annotation-guided way
            v_region = v_regions[0]
            j_region = j_regions[0]

            has_v_start = v_region.annotation_match_start == 0
            v_start = v_region.contig_match_start
            v_seq = v_region.sequence
            seq = self.sequence

            cdr_pos = None

            if has_v_start and ensure_str(v_seq[0:3]) in START_CODONS:
                # Full V and V begins with a start codon (because it's an L-REGION+V-REGION).
                # Force the frame to 0.
                # Search for a CDR3 in frame 0.
                (cdr_pos, flag) = search_cdr3_signature(v_region, j_region, 0)
                flags.append(flag)
                flags.append("FULL_V_HAS_START")
                self.start_codon_pos = v_start

            elif has_v_start:
                # Full V but the first annotated codon is not a start.
                # Look for an alternative start nearby that preserves the V frame.
                for i in range(
                    v_start - 3 * START_CODON_SLOP, 1 + v_start + 3 * START_CODON_SLOP, 3
                ):
                    if (
                        i > 0
                        and (i + 3) <= len(seq)
                        and ensure_str(seq[i : (i + 3)]) in START_CODONS
                    ):
                        (cdr_pos, flag) = search_cdr3_signature(v_region, j_region, 0)
                        self.start_codon_pos = i
                        flags.append(flag)
                        flags.append("FULL_V_ALT_START")
                if self.start_codon_pos is None:
                    flags.append("FULL_V_NO_START")

            if not has_v_start or self.start_codon_pos is None:
                # Either we don't contain the start of V or we didn't find a valid start codon above
                # Look for a CDR3 sequence in all frames

                for frame in [0, 1, 2]:
                    ((cdr_pos), flag) = search_cdr3_signature(v_region, j_region, frame)

                    if cdr_pos:
                        flags.append(flag)
                        break
                if cdr_pos is None:
                    flags.append("FAILED_UNGUIDED_SEARCH")

            if cdr_pos and cdr_pos[1] - cdr_pos[0] < VDJ_MAX_CDR3_LEN:
                self.cdr3_start = cdr_pos[0]
                self.cdr3_stop = cdr_pos[1]
                self.cdr3_seq = seq[cdr_pos[0] : cdr_pos[1]]
                self.cdr3 = bytes(
                    bytearray(
                        codon_to_aa(self.cdr3_seq[i : (i + 3)])
                        for i in range(0, len(self.cdr3_seq) - 2, 3)
                    )
                )

                assert (cdr_pos[1] - cdr_pos[0]) % 3 == 0
            elif cdr_pos is None:
                flags.append("NO_CDR3")
            else:
                flags.append("CDR3_TOO_LONG:%d" % (cdr_pos[1] - cdr_pos[0]))

        if not self.cdr3:
            # Either this didn't have both a V and a J, or the annotation-guided search failed to give a valid CDR3.
            # Run the unguided CDR3 search.
            v_region = v_regions[0] if v_regions else None
            j_region = j_regions[0] if j_regions else None

            if v_region is not None:
                res = search_cdr3_signature_no_vj(self.sequence, v_region, j_region)
                if not res is None:
                    (cdr3_seq, cdr3_aas, cys_pos, fgxg_pos) = res
                    self.cdr3 = cdr3_aas
                    self.cdr3_seq = cdr3_seq
                    self.cdr3_start = cys_pos
                    self.cdr3_stop = (
                        fgxg_pos + 3
                    )  # End position is the start position of last amino acid + 3
                    flags.append("FOUND_CDR3_UNGUIDED")

        if self.cdr3:
            # cdr3_frame = self.cdr3_start % 3 # frame wrt start of sequence

            if self.start_codon_pos is None:
                # We don't have a start codon.
                # De novo search for the start codon in the CDR3 frame;
                # Get the leftmost possible start codon match before a stop codon
                for i in range(self.cdr3_start, 0, -3):
                    if ensure_str(self.codon(i)) in START_CODONS:
                        self.start_codon_pos = i
                    if self.start_codon_pos and ensure_str(self.codon(i)) in STOP_CODONS:
                        break

                if self.start_codon_pos:
                    flags.append("FOUND_DE_NOVO_V_START")

            if self.start_codon_pos is not None:
                # If we have a start codon, try to find a stop codon in-frame
                for i in range(self.start_codon_pos, len(self.sequence), 3):
                    if ensure_str(self.codon(i)) in STOP_CODONS:
                        self.stop_codon_pos = i
                        break

            # Determine productivity
            if self.has_full_length_vj_hit():
                has_start = self.start_codon_pos is not None
                cdr3_in_frame = (
                    has_start
                    and self.cdr3_start is not None
                    and (self.start_codon_pos % 3) == (self.cdr3_start % 3)
                )
                vj_nostop = has_start and all(
                    ensure_str(self.codon(i)) not in STOP_CODONS
                    for i in range(self.start_codon_pos, j_region.contig_match_end, 3)
                )
                self.productive = has_start and cdr3_in_frame and vj_nostop

                if not has_start:
                    flags.append("NO_START")
                if has_start and self.cdr3_start is not None and not cdr3_in_frame:
                    flags.append("CDR3_OUT_OF_FRAME")
                if not vj_nostop:
                    flags.append("VJ_STOP")

            elif any(
                ensure_str(self.codon(i)) in STOP_CODONS
                for i in range(self.cdr3_start, self.cdr3_stop, 3)
            ):
                # We don't know where the transcript starts and ends so we'll
                # be cautious about calling productivity.
                flags.append("CDR3_STOP")
                self.productive = False

        # Translate the entire reading frame if possible
        if self.start_codon_pos is not None:
            self.aa_sequence = bytes(
                bytearray(
                    codon_to_aa(self.codon(i))
                    for i in range(self.start_codon_pos, len(self.sequence) - 2, 3)
                )
            )

        self.cdr3_flag = "|".join([f for f in flags if f])

    def codon(self, i):
        return self.sequence[i : (i + 3)] if i + 3 <= len(self.sequence) else b""

    def contains_annotations(self, other_contig):
        """True if this contig contains all the annotations of the other_contig.

        We only check the gene_name to test if the two annotations are the same.
        """
        for other_annotation in other_contig.hits():
            if not any(
                hit.feature.gene_name == other_annotation.feature.gene_name
                for hit in self.annotations
            ):
                return False
        return True

    def is_exact_vj_hit(self, other_contig: AnnotatedContig) -> bool:
        """True if this contig's VJ region is exactly contained in the other_contig."""
        v_regions = self.get_region_hits(VDJ_V_FEATURE_TYPES)
        j_regions = self.get_region_hits(VDJ_J_FEATURE_TYPES)
        if not v_regions or not j_regions:
            return False

        # Get sequence from beginning of V to end of J (this will contain D if present)
        start = np.min([v_region.contig_match_start for v_region in v_regions])
        stop = np.max([j_region.contig_match_end for j_region in j_regions])
        if re.search(self.sequence[start:stop], other_contig.sequence):
            return True
        return False

    def __len__(self):
        return len(self.sequence)

    def annotation_str(self):
        features_to_output = (
            VDJ_V_FEATURE_TYPES + VDJ_J_FEATURE_TYPES + VDJ_C_FEATURE_TYPES + VDJ_D_FEATURE_TYPES
        )

        annotation_strs = []
        for anno in self.annotations:
            if ensure_str(anno.feature.region_type) in features_to_output:
                annotation_strs.append(
                    (
                        anno.contig_match_start,
                        anno.annotation_match_start,
                        anno.contig_match_end - anno.contig_match_start,
                        anno.annotation_length,
                        anno.feature.display_name,
                    )
                )
        annotation_strs = [f"{x[4]}" for x in sorted(annotation_strs, key=lambda x: x[0])]
        return ";".join(annotation_strs)

    def __str__(self):
        annotation_str = ";".join([str(anno) for anno in self.annotations])
        return f"Contig {self.contig_name}: {annotation_str}"

    def __repr__(self):
        return self.to_dict().__repr__()


class AnnotationDict(TypedDict):
    feature: vdj_reference.VdjAnnotationFeatureDict
    cigar: bytes | None
    annotation_length: int
    annotation_match_start: int
    annotation_match_end: int
    contig_match_start: int
    contig_match_end: int


class Annotation:
    """Annotation of a sequence against a reference of segment sequences."""

    __slots__ = [
        "feature",
        "cigar",
        "score",
        "annotation_length",
        "annotation_match_start",
        "annotation_match_end",
        "contig_match_start",
        "contig_match_end",
        "mismatches",
        "sequence",
    ]

    def __init__(
        self,
        feature: vdj_reference.VdjAnnotationFeature,
        cigar: bytes | None,
        score: float,
        annotation_length: int,
        annotation_match_start: int,
        annotation_match_end: int,
        contig_match_start: int,
        contig_match_end: int,
        mismatches,
        contig_seq: bytes | None = None,
    ):
        """Constructor.

        Args:
            feature (vdj.reference.VdjAnnotationFeature): object
        """
        self.feature: vdj_reference.VdjAnnotationFeature = feature
        self.cigar = cigar
        self.score: float = score
        self.annotation_length: int = annotation_length
        self.annotation_match_start: int = annotation_match_start
        self.annotation_match_end: int = annotation_match_end
        self.contig_match_start: int = contig_match_start
        self.contig_match_end: int = contig_match_end
        self.mismatches = mismatches

        if contig_seq is None:
            self.sequence: bytes | None = None
        else:
            self.sequence = contig_seq[contig_match_start:contig_match_end].upper()

    @classmethod
    def from_alignment_result(cls, alignment_result: SSWAlignmentResult, sequence: bytes | None):
        alignment = alignment_result.alignment
        assert alignment is not None

        feature = alignment_result.reference.metadata["feature"]

        # Annotation object points to the input sequence (pre-masking)
        annotation = cls(
            feature=feature,
            cigar=alignment.cigar_string,
            score=alignment.score,
            annotation_length=len(feature.sequence),
            annotation_match_start=alignment.ref_begin,
            # SSW returns inclusive ends - convert to open
            annotation_match_end=alignment.ref_end + 1,
            contig_match_start=alignment.query_begin,
            contig_match_end=alignment.query_end + 1,
            mismatches=[],
            contig_seq=sequence,
        )

        return annotation

    @staticmethod
    def _inferred_values():
        return ["chain_type", "sequence"]

    def annotate_mismatches(self, contig_sequence, feature_sequence):
        """Mark mismatches between the contig and the reference.

        Args:
            contig_sequence (str): contig sequence
            feature_sequence (str): the feature that the annotation derives from

        Returns:
            list of dict: list of annotations for mismatches (each as a dict)
        """
        mismatches = []

        cigar_tuples = cr_cigar.get_cigar_tuples(self.cigar)

        # Generate annotations for Indels
        ## Track position and all operations
        start_position = 0
        all_operations = []

        mismatch_annotation_json_keys = ["region_type", "contig_match_start", "contig_match_end"]

        for length, category in cigar_tuples:
            length = int(length)
            all_operations.extend(length * [category])

            if category == b"D":
                # By convention, make deletion length always equal to 1. Really the length of the deletion
                # in terms of contig coordinates is 0, but this makes the deletions invinsible in the visualization.
                info = dict(
                    zip(
                        mismatch_annotation_json_keys,
                        [category, start_position, start_position + 1],
                    )
                )
                # This is the only case where the length cannot be inferred from the start and end positions
                info["deletion_length"] = length
                mismatches.append(info)
            elif category == b"I":
                mismatches.append(
                    dict(
                        zip(
                            mismatch_annotation_json_keys,
                            [category, start_position, start_position + length],
                        )
                    )
                )

            if category != b"D":
                start_position += length

        # Generate annotations for one base mismatches
        contig_position = 0
        annotation_position = self.annotation_match_start

        for category in all_operations:
            if category == b"S":
                contig_position += 1

            elif category == b"M":
                if (
                    contig_sequence[contig_position : contig_position + 1].upper()
                    != feature_sequence[annotation_position : annotation_position + 1].upper()
                ):
                    mismatches.append(
                        dict(
                            zip(
                                mismatch_annotation_json_keys,
                                ["MISMATCH", contig_position, contig_position + 1],
                            )
                        )
                    )
                annotation_position += 1
                contig_position += 1

            elif category == b"I":
                contig_position += 1

            elif category == b"D":
                annotation_position += 1

        return mismatches

    @classmethod
    def from_dict(
        cls,
        annotation_dict: AnnotationDict,
        reference: vdj_reference.VdjReference,
        contig_seq: str | bytes | None = None,
    ) -> Annotation:
        return cls(
            feature=vdj_reference.convert_dict_to_vdj_feature(
                annotation_dict["feature"], reference
            ),
            cigar=annotation_dict["cigar"],
            score=annotation_dict.get("score", float("nan")),
            annotation_length=annotation_dict["annotation_length"],
            annotation_match_start=annotation_dict["annotation_match_start"],
            annotation_match_end=annotation_dict["annotation_match_end"],
            contig_match_start=annotation_dict["contig_match_start"],
            contig_match_end=annotation_dict["contig_match_end"],
            mismatches=annotation_dict.get("mismatches"),
            contig_seq=_bytes_or_none(contig_seq),
        )

    def to_dict(self) -> AnnotationDict:
        out_dict = {slot: getattr(self, slot) for slot in self.__slots__}
        out_dict["feature"] = vdj_reference.convert_vdj_feature_to_dict(self.feature)
        for val in self._inferred_values():
            out_dict.pop(val, None)
        return out_dict

    def copy(self):
        return self.from_dict(self.to_dict(), self.sequence)

    def __eq__(self, other):
        """Two annotations are the same if they have the same sequence.

        If the sequence is missing, rely on the annotated gene name.
        """
        if not isinstance(other, self.__class__):
            return False
        if not self.sequence is None and self.sequence == other.sequence:
            return True
        return (
            self.feature.gene_name == other.feature.gene_name
            and self.feature.allele_name == other.feature.allele_name
        )

    def __hash__(self):
        """Hash is based on either the sequence or gene name."""
        if not self.sequence is None:
            return hash(self.sequence)
        return hash((self.feature.gene_name, self.feature.allele_name))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return f"{self.feature.display_name}(cigar:{self.cigar})"

    def __repr__(self):
        return self.to_dict().__repr__()


class CellContigs:
    """A list of AnnotatedContigs.

    (eg. representing a clonotype or all the contigs in a cell)
    """

    def __init__(self, name, contigs: Iterable[AnnotatedContig], info_dict={}):
        self.name = name
        self.chains = list(contigs)
        self.info_dict = info_dict

    @classmethod
    def empty(cls, name):
        return cls(name, [], {})

    def contigs(self):
        return iter(self.chains)

    def is_paired(self, require_full_len=True, require_productive=True, require_high_conf=True):
        """True if this cell has contigs matching some VDJ gene pair."""
        good_pairs: list[set[str]] = [set(p.split("_")) for p in VDJ_GENE_PAIRS]
        chains: set[str] = set()
        cl = self.clonotype_tuple(
            require_full_len=require_full_len,
            require_productive=require_productive,
            require_high_conf=require_high_conf,
        )
        for cdr in cl:
            chain = ensure_str(cdr).split("_")[0]
            chains.add(chain)
        return chains in good_pairs

    def clonotype_tuple(
        self,
        require_full_len: bool = True,
        require_productive: bool = True,
        require_high_conf: bool = True,
    ) -> tuple[bytes, ...]:
        """Tuple of unique CDR3s across all productive contigs.

        Args:
            require_productive: Contigs (aka 'chains') need to be productive (productive=True)
            require_full_len: Contigs (aka 'chains') need to be full length (productive or not)
            require_high_conf: Contigs (aka 'chains') need to be high confidence
        """
        cdrs = {
            c.clonotype_seq()[1]
            for c in self.chains
            if (c.productive or not require_productive)
            and (c.has_full_length_vj_hit() or not require_full_len)
            and (c.high_confidence or not require_high_conf)
            and c.has_cdr()
            and c.is_single_chain()
        }

        # Sort for easy comparisons
        return tuple(sorted(cdrs))

    def contains_contig(self, other_contig: AnnotatedContig) -> bool:
        """True if it has a chain that contains all the annotations of the given other contig.

        Args:
            contig (AnnotatedContig):
        """
        return any(contig.contains_annotations(other_contig) for contig in self.chains)

    def has_exact_vj_hit(self, other_contig):
        """Returns True if this cell has an exact hit in the other contig.

        True if this cell has some contig whose V(D)J region is perfectly
        contained in the sequence of the other_contig.
        """
        if other_contig is None:
            return False

        for contig in self.chains:
            if contig.is_exact_vj_hit(other_contig):
                return True

        return False

    def to_dict_list(self):
        """Return as a list of AnnotatedContig dicts.

        Information from this object's info_dict is passed to the contigs.
        """
        out_dicts = []
        for contig in self.chains:
            new_contig = contig.copy()
            new_contig.info_dict.update(dict(self.info_dict))
            contig_out = new_contig.to_dict()
            contig_out["barcode"] = self.name
            out_dicts.append(contig_out)
        return out_dicts

    def copy(self):
        return CellContigs(
            self.name, [contig.copy() for contig in self.chains], dict(self.info_dict)
        )

    def __str__(self):
        return "Cell {}:\n{}".format(self.name, "\n".join([str(chain) for chain in self.chains]))

    def __repr__(self):
        return self.to_dict_list().__repr__()


def load_contig_list_from_json(json_file, reference_path):
    """Returns a list of AnnotatedContig objects from an open json file."""
    reference = vdj_reference.VdjReference(reference_path)
    return [AnnotatedContig.from_dict(x, reference) for x in json.load(json_file)]


def load_cell_contigs_from_json(json_file, reference_path, group_key: str, require_high_conf=True):
    """Returns a list of CellContig objects based on annotations in a json.

    The json is assumed to contain a list of AnnotatedContigs (in dict form).
    The contigs are sorted and grouped by group_key and each such group is put
    into a CellContig object.

    Args:
        group_key: must be 'barcode' or 'clonotype'
    """
    assert group_key in ("barcode", "clonotype")
    with open(json_file) as json_file_obj:
        annotations = load_contig_list_from_json(json_file_obj, reference_path)

    cell_contigs = []

    def key_func(x):
        return getattr(x, group_key)

    anno_iter = itertools.groupby(sorted(annotations, key=key_func), key=key_func)
    for clonotype_name, contig_annotations in anno_iter:  # type: tuple[bytes, Any]
        contigs = []
        for new_contig in contig_annotations:  # type: AnnotatedContig
            # Note, for consensus contigs is_cell=None
            if new_contig.is_cell is not False and (
                new_contig.high_confidence or not require_high_conf
            ):
                contigs.append(new_contig)

        if len(contigs) > 0:
            cell_contigs.append(CellContigs(clonotype_name, contigs))

    return cell_contigs


def save_annotation_list_json(out_file, contigs: Iterable[AnnotatedContig]):
    tk_safe_json.dump_numpy((c.to_dict() for c in contigs), out_file, pretty=True)


def save_contig_list_csv(csv, contigs: Iterable[AnnotatedContig], write_inferred=True):
    """Write contigs to an open csv file."""
    columns = [
        "barcode",
        "is_cell",
        "contig_id",
        "high_confidence",
        "length",
        "chain",
        "v_gene",
        "d_gene",
        "j_gene",
        "c_gene",
        "full_length",
        "productive",
        "cdr3",
        "cdr3_nt",
        "reads",
        "umis",
        "raw_clonotype_id",
        "raw_consensus_id",
    ]

    if write_inferred:
        columns.extend(["inferred_clonotype_id", "inferred_consensus_id"])

    save_annotation_list_csv(csv, contigs, columns)


def save_consensus_list_csv(csv, contigs):
    """Write consensus contigs to an open csv file."""
    columns = [
        "clonotype_id",
        "consensus_id",
        "length",
        "chain",
        "v_gene",
        "d_gene",
        "j_gene",
        "c_gene",
        "full_length",
        "productive",
        "cdr3",
        "cdr3_nt",
        "reads",
        "umis",
    ]
    save_annotation_list_csv(csv, contigs, columns)


def ensure_str_or_none(x):
    if x is None:
        return x
    return ensure_str(x)


def save_annotation_list_csv(csv, contigs, columns):
    """Write AnnotatedContigs to an open csv file."""
    col_set = set(columns)

    vdj_utils.write_csv_row(columns, csv)

    for contig in sorted(contigs, key=lambda x: (not x.is_cell, x.barcode, x.contig_name)):
        if contig.filtered:
            row = {
                "barcode": ensure_str(contig.barcode),
                "is_cell": contig.is_cell,
                "contig_id": ensure_str(contig.contig_name),
                "high_confidence": contig.high_confidence,
                "length": len(contig),
                "chain": ensure_str_or_none(contig.get_single_chain()),
                "v_gene": contig.get_single_gene_display_name(VDJ_V_FEATURE_TYPES),
                "d_gene": contig.get_single_gene_display_name(VDJ_D_FEATURE_TYPES),
                "j_gene": contig.get_single_gene_display_name(VDJ_J_FEATURE_TYPES),
                "c_gene": contig.get_single_gene_display_name(VDJ_C_FEATURE_TYPES),
                "full_length": contig.has_full_length_vj_hit(),
                "productive": contig.productive,
                "cdr3": ensure_str_or_none(contig.cdr3),
                "cdr3_nt": ensure_str_or_none(contig.cdr3_seq),
                "reads": contig.read_count,
                "umis": contig.umi_count,
                "raw_clonotype_id": contig.info_dict.get("raw_clonotype_id"),
                "raw_consensus_id": contig.info_dict.get("raw_consensus_id"),
                "inferred_clonotype_id": contig.info_dict.get("inferred_clonotype_id"),
                "inferred_consensus_id": contig.info_dict.get("inferred_consensus_id"),
                "clonotype_id": contig.clonotype,
                "consensus_id": contig.contig_name,
                "clonotype_frequency": contig.info_dict.get("clonotype_freq"),
                "clonotype_proportion": contig.info_dict.get("clonotype_prop"),
            }
            assert col_set.issubset(row.keys())
            vdj_utils.write_csv_row([row[k] for k in columns], csv)


def save_clonotype_info_csv(csv: TextIO, consensus_contigs: Iterable[AnnotatedContig]):
    """Write a CSV containing clonotype info to an open csv file.

    Takes a list of AnnotatedContigs corresponding to consensuses.
    """
    clonotypes = defaultdict(dict)

    for contig in consensus_contigs:
        chain = contig.get_single_chain()
        if chain is None:
            continue

        clonotype_id = contig.clonotype
        if clonotype_id not in clonotypes:
            clonotypes[clonotype_id]["clonotype_id"] = ensure_str(clonotype_id)
            clonotypes[clonotype_id]["members"] = set(contig.info_dict["cells"])
            clonotypes[clonotype_id]["frequency"] = contig.info_dict["clonotype_freq"]
            clonotypes[clonotype_id]["proportion"] = contig.info_dict["clonotype_prop"]
            clonotypes[clonotype_id]["cdr3s"] = set()
        else:
            assert clonotypes[clonotype_id]["members"] == set(contig.info_dict["cells"])
            assert clonotypes[clonotype_id]["frequency"] == contig.info_dict["clonotype_freq"]
            assert clonotypes[clonotype_id]["proportion"] == contig.info_dict["clonotype_prop"]

        clonotypes[clonotype_id]["cdr3s"].add((chain, contig.cdr3_seq, contig.cdr3))

    # Generate cdr3 annotation strings (for nt and aa), chain:cdr3;chain:cdr3
    # sorted by (chain, cdr3_seq, cdr3_aa)
    def get_cdr3_list_string(chain_cdr3s):
        cdr3s_nt = []
        cdr3s_aa = []
        for chain, nt, aa in sorted(chain_cdr3s):
            cdr3s_nt.append(f"{chain}:{nt}")
            cdr3s_aa.append(f"{chain}:{aa}")
        return ";".join(cdr3s_nt), ";".join(cdr3s_aa)

    for clonotype in clonotypes.values():
        cdr3s_nt, cdr3_aa = get_cdr3_list_string(clonotype["cdr3s"])
        clonotype["cdr3s_nt"] = cdr3s_nt
        clonotype["cdr3s_aa"] = cdr3_aa

    # Sort by frequency, descending
    clonotypes = sorted(clonotypes.values(), key=lambda c: c["frequency"], reverse=True)

    columns = ["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt"]
    col_set = set(columns)

    vdj_utils.write_csv_row(columns, csv)

    for row in clonotypes:
        assert col_set.issubset(set(row))
        vdj_utils.write_csv_row([row[k] for k in columns], csv)


def get_clonotype_vdj_pair(sequence_ids, clonotype_tuple, vdj_gene_pairs):
    """Get the VDJ gene pair associated with the clonotype.

    If the gene pair is not in vdj_gene_pairs, it will return None.
    """
    genes = {sequence_ids[seq].split("_")[0] for seq in clonotype_tuple}
    genes = tuple(sorted(genes))

    return vdj_gene_pairs.get(genes, None)


def report_clonotypes(
    reporter,
    prefix: str,
    cell_barcodes,
    clonotype_ids,
    sequence_ids,
    barcode_contigs,
    bc_clonotype_assignments,
):
    """Group barcodes into clonotypes.

    Args:
        reporter (VdjReporter): object for reporting metrics
        prefix: metrics prefix (must be in VDJ_CLONOTYPE_TYPES)
        cell_barcodes: set of cell barcodes
        clonotype_ids (dict): Mapping from clonotype id to tuple of CDR ids
        sequence_ids (dict): Mapping from CDR id to sequence
        barcode_contigs (list): CellContigs objects
        bc_clonotype_assignments (dict): Mapping from barcode to clonotype id

    Returns:
        dict: clonotype_id -> clonotype_info.
          Clonotype info is itself a dict with the following information::
          - clonotype_id
          - barcodes: barcodes belonging to the clonotype
          - freq: number of cells in the clonotype
          - prop: relative frequency of barcodes in the clonotype
          - consensuses: A dict seq_id -> seq_info with info about the sequences of the
              consensus. seq_info has the following information:
                  - cdr3: AA sequence of CDR3
                  - cdr3_seq: nucleotide sequence of CDR3
                  - chain: chain
                  - cell_contigs: contigs (of the clonotype barcodes) that correspond to this
                  consensus sequence
    """
    clonotype_counts = np.bincount(bc_clonotype_assignments.values())

    freq_metric = reporter._get_metric_attr("vdj_clonotype_freq", prefix)
    prop_metric = reporter._get_metric_attr("vdj_clonotype_prop", prefix)

    # Cells grouped by clonotype assignment. All cells left unassigned will be
    # in the "None" group.
    grouped_barcodes = itertools.groupby(
        sorted(barcode_contigs, key=lambda x: bc_clonotype_assignments.get(x.name, None)),
        key=lambda x: bc_clonotype_assignments.get(x.name, None),
    )

    grouped_barcodes = sorted(
        ((cl, list(b)) for (cl, b) in grouped_barcodes), key=lambda x: len(x[1]), reverse=True
    )
    unassigned_metric = reporter._get_metric_attr("vdj_unassigned_clonotype_bc_frac", prefix)

    observed_barcodes = []

    vdj_gene_pairs = {tuple(sorted(pair.split("_"))): pair for pair in VDJ_GENE_PAIRS}

    out_clonotypes = {}
    clonotype_count = 0

    for cl, bc_contig_lists in grouped_barcodes:
        # Previous operations in the clonotype assignments might have "emptied" out a clonotype.
        if len(bc_contig_lists) == 0:
            continue

        barcode_names = [bc.name for bc in bc_contig_lists]
        # Keep track of barcodes in clonotypes
        observed_barcodes.extend(barcode_names)

        if cl is None:
            vdj_gene_pair = None
            is_cl_paired = False
        else:
            clonotype_count = len(bc_contig_lists)
            assert clonotype_count == clonotype_counts[cl]

            vdj_gene_pair = get_clonotype_vdj_pair(sequence_ids, clonotype_ids[cl], vdj_gene_pairs)

            is_cl_paired = not vdj_gene_pair is None

            assert prefix in VDJ_CLONOTYPE_TYPES
            clonotype_name = vdj_utils.format_clonotype_id(
                len(out_clonotypes), inferred=prefix == "inferred"
            )

            # Compute frequency and proportion
            freq = clonotype_count
            prop = tk_stats.robust_divide(freq, len(cell_barcodes))

            freq_metric.add(clonotype_name, freq)
            prop_metric.add(clonotype_name, prop)

            clonotype_chains = [sequence_ids[cl_seq] for cl_seq in clonotype_ids[cl]]
            clonotype_chain_cdr3_translations = {}

            out_contigs = defaultdict(list)

            for bc in bc_contig_lists:
                for contig in bc.contigs():
                    _, cdr_seq = contig.clonotype_seq()

                    if not cdr_seq in clonotype_chains:
                        continue

                    out_contigs[cdr_seq].append(contig.contig_name)
                    clonotype_chain_cdr3_translations[cdr_seq] = contig.cdr3

            out_clonotype = {
                "clonotype_id": clonotype_name,
                "freq": freq,
                "prop": prop,
                "barcodes": barcode_names,
            }

            consensuses = {}
            for idx, (cdr_seq, members) in enumerate(out_contigs.items()):
                assert len(members) > 0
                consensus_name = f"{clonotype_name}_consensus_{idx + 1}"
                consensus_info = {
                    "chain": cdr_seq.split("_")[0],
                    "cdr3": clonotype_chain_cdr3_translations[cdr_seq],
                    "cdr3_seq": cdr_seq.split("_")[1],
                    "cell_contigs": members,
                }
                consensuses[consensus_name] = consensus_info

            out_clonotype["consensuses"] = consensuses

            out_clonotypes[clonotype_name] = out_clonotype

        unassigned_metric.add(len(bc_contig_lists), filter=(cl is None))

        for gp in reporter.canonical_vdj_gene_pairs:
            # Due to the way we define clonotypes, paired here implies productive.
            # And due to the way we define productive, this also means full length.
            paired_cls_metric = reporter._get_metric_attr("vdj_paired_clonotype_frac", gp, prefix)
            bcs_in_paired_cls_metric = reporter._get_metric_attr(
                "vdj_paired_clonotype_bc_frac", gp, prefix
            )

            diversity_metric = reporter._get_metric_attr("vdj_clonotype_diversity", gp, prefix)
            paired_diversity_metric = reporter._get_metric_attr(
                "vdj_paired_clonotype_diversity", gp, prefix
            )

            if cl is not None:
                paired_cls_metric.add(
                    1, filter=is_cl_paired and (gp in (vdj_gene_pair, MULTI_REFS_PREFIX))
                )
                bcs_in_paired_cls_metric.add(
                    len(bc_contig_lists),
                    filter=is_cl_paired and (gp in (vdj_gene_pair, MULTI_REFS_PREFIX)),
                )

                if gp in (vdj_gene_pair, MULTI_REFS_PREFIX):
                    num_clonotypes_metric = reporter._get_metric_attr(
                        "vdj_clonotype_count", gp, prefix
                    )
                    num_clonotypes_metric.add(1)

                    diversity_metric.add(cl, clonotype_count)
                    if is_cl_paired:
                        paired_diversity_metric.add(cl, clonotype_count)

    # These barcodes have no annotations at all and no clonotype assignments
    remaining_barcodes = cell_barcodes.difference(set(observed_barcodes))
    for gp in reporter.canonical_vdj_gene_pairs:
        bcs_in_paired_cls_metric = reporter._get_metric_attr(
            "vdj_paired_clonotype_bc_frac", gp, prefix
        )
        bcs_in_paired_cls_metric.add(len(remaining_barcodes), filter=False)
    unassigned_metric.add(len(remaining_barcodes), filter=True)

    metric = reporter._get_metric_attr("cdrs_per_bc_histogram", prefix)
    metric.add_many(np.array([len(clonotype_ids[c]) for c in bc_clonotype_assignments.values()]))
    metric.add_many(np.zeros((len(remaining_barcodes),), dtype=int))

    metric = reporter._get_metric_attr("major_clonotype_bc_frac", prefix)
    if len(clonotype_counts) > 0:
        metric.set_value(np.max(clonotype_counts), len(cell_barcodes))
    else:
        metric.set_value(0, len(cell_barcodes))

    return out_clonotypes


def label_contigs_with_consensus(clonotypes, contigs, prefix):
    """Adds clonotype and consensus info to a list of AnnotatedContig objects.

    Args:
        clonotypes: dict like the one output by report_clonotypes.
        contigs: list of AnnotatedContig objects
        prefix: prefix for clonotype ids (must be in VDJ_CLONOTYPE_TYPES)
    """
    assert prefix in VDJ_CLONOTYPE_TYPES
    contig_dict = {contig.contig_name: contig for contig in contigs}

    bc_assignments = {}
    # iterate over clonotypes
    for _, clonotype_info in clonotypes.items():
        # keep track of barcode to clonotype assignments.
        for bc in clonotype_info["barcodes"]:
            bc_assignments[bc] = clonotype_info["clonotype_id"]
        # iterate over consensus sequences of the clonotype and annotate
        # each member-contig of each consensus
        for consensus_name, consensus_info in clonotype_info["consensuses"].items():
            for mem in consensus_info["cell_contigs"]:
                contig_dict[mem].info_dict[f"{prefix}_consensus_id"] = consensus_name

    # Now annotate contigs with their clonotype assignments. We don't do
    # this in the previous pass because we want to annotate contigs whose
    # barcode was assigned to a clonotype but which did not participate in
    # any consensuses (eg. some non-full-length contigs)
    for contig in contigs:
        if contig.barcode in bc_assignments:
            contig.info_dict[f"{prefix}_clonotype_id"] = bc_assignments[contig.barcode]
