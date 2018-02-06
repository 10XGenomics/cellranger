#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Utilities for annotated contigs with gene/chain information,
# defining clonotypes and more.

import itertools
import json
import numpy as np
import re
import cellranger.align as cr_align
from cellranger.vdj.constants import (VDJ_5U_FEATURE_TYPES, VDJ_D_FEATURE_TYPES,
                                      VDJ_V_FEATURE_TYPES, VDJ_J_FEATURE_TYPES,
                                      VDJ_C_FEATURE_TYPES, VDJ_ORDERED_REGIONS,
                                      START_CODONS, STOP_CODONS, CODON_TO_AA,
                                      VDJ_MAX_CDR3_LEN, VDJ_MIN_CDR3_LEN,
                                      VDJ_ANNOTATION_MIN_SCORE_RATIO,
                                      VDJ_ANNOTATION_MIN_WORD_SIZE,
                                      VDJ_ANNOTATION_MATCH_SCORE,
                                      VDJ_ANNOTATION_MISMATCH_PENALTY,
                                      VDJ_ANNOTATION_GAP_OPEN_PENALTY,
                                      VDJ_ANNOTATION_EXTEND_PENALTY,
                                      VDJ_PRIMER_ANNOTATION_MIN_FRACTION_MATCHED,
                                      VDJ_QUAL_OFFSET, VDJ_CLONOTYPE_TYPES,
                                      VDJ_GENE_PAIRS)
from cellranger.constants import MULTI_REFS_PREFIX
import cellranger.vdj.reference as vdj_reference
import cellranger.vdj.utils as vdj_utils
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
import tenkit.seq as tk_seq
from collections import defaultdict

# Allow the start codon to shift up or down by this many codons
START_CODON_SLOP = 1

def filter_alignment(alignment_result, score_ratio, word_size,
                     match_score=VDJ_ANNOTATION_MATCH_SCORE):
    """Returns True for a passing alignment and False otherwise.

    Args:
    - alignment_result (SSWAlignmentResult): alignment result to filter
    - score_ratio: minimum (score / max_possible_score)
    - word_size: minimum stretch of contiguous matches/mismatches in the alignment
    - match_score: match score used in the alignment

    Returns:
        True if alignment passed filters
    """

    alignment = alignment_result.alignment

    # alignment is a PyAlignRes object. Ends are inclusive.
    alignment_length = float(alignment.query_end - alignment.query_begin + 1)
    max_score = alignment_length * match_score
    if tk_stats.robust_divide(alignment.score, max_score) < score_ratio:
        return False

    if cr_align.get_max_word_length(alignment) < word_size:
        return False

    return True

def setup_feature_aligners(reference_path, score_ratios, word_sizes, use_features=None):
    """Returns a list of aligners, one for each type of feature.

    Args:
    - reference_path: Path to VDJ reference.
    - score_ratios: A dict from a feature name to the min score_ratio to use
    for filtering that feature.
    - word_sizes: A dict from a feature name to the min word_size to use
    for filtering that feature.
    - use_features: set(int) of feature IDs to restrict the annotation to

    Return value:
    A tuple (feature_types, aligners, filters). Each element is a list.
    aligners[i] is an aligner (SSWMultiAligner) against elements of type
    feature_types[i]. filters[i] is a filter function to use for alignments
    of type feature_types.
    """
    feature_iter = vdj_reference.get_vdj_feature_iter(reference_path)

    feature_refs = defaultdict(list)

    # Group references by region type. First, match each feature type to
    # a representative name (so, for example, ref sequences with different types
    # which all belong in VDJ_V_FEATURE_TYPES end up in the same group).
    ordered_region_dict = {}
    for names in VDJ_ORDERED_REGIONS:
        for n in names:
            ordered_region_dict[n] = names[0]

    for feature in feature_iter:
        if feature.region_type not in ordered_region_dict:
            continue
        if use_features is not None and feature.feature_id not in use_features:
            continue

        new_ref = cr_align.SSWReference(name=feature.gene_name,
                                        sequence=feature.sequence,
                                        min_score=0,
                                        min_match_length=0,
                                        metadata={'feature': feature})
        feature_refs[ordered_region_dict[feature.region_type]].append(new_ref)

    aligners = []
    filters = []

    for region_types in VDJ_ORDERED_REGIONS:
        feature_aligner = cr_align.SSWMultiAligner(feature_refs[region_types[0]],
                                                   match_score=VDJ_ANNOTATION_MATCH_SCORE,
                                                   mismatch_penalty=VDJ_ANNOTATION_MISMATCH_PENALTY,
                                                   gap_open_penalty=VDJ_ANNOTATION_GAP_OPEN_PENALTY,
                                                   gap_extend_penalty=VDJ_ANNOTATION_EXTEND_PENALTY,
                                                   report_cigar=True)

        score = score_ratios.get(region_types[0], VDJ_ANNOTATION_MIN_SCORE_RATIO)
        word = word_sizes.get(region_types[0], VDJ_ANNOTATION_MIN_WORD_SIZE)
        match = VDJ_ANNOTATION_MATCH_SCORE

        # lambda's variables are resolved when lambda is called.
        # Need to use default arguments, because the default arguments will be evaluated
        # when the lambda is created not when it's called.
        feature_filter_params = lambda x, score=score, word=word, match=match: filter_alignment(x, score, word, match)

        aligners.append(feature_aligner)
        filters.append(feature_filter_params)

    return [v[0] for v in VDJ_ORDERED_REGIONS], aligners, filters


def setup_primer_aligner(primers, score_ratio):
    """ Returns (SSWMultiAligner aligner, lambda aligner_filter) for
    aligning against a set of sequences.

    Sequences are provided as a list of dicts. Each dict must have a "name"
    and a "seq".
    """

    primer_refs = []

    for primer in primers:
        for direction in ['FORWARD', 'REVERSE']:
            if direction == 'REVERSE':
                seq = tk_seq.get_rev_comp(primer['seq'])
            else:
                seq = primer['seq']

            primer_refs.append(cr_align.SSWReference(name=primer['name'],
                                                     sequence=seq,
                                                     min_score=0,
                                                     min_match_length=0,
                                                     metadata={
                                                         'primer': primer['name'] + '_' + direction,
                                                         'feature': vdj_reference.create_dummy_feature(
                                                             display_name=primer['name'] + '_' + direction,
                                                             region_type='PRIMER',
                                                             sequence=seq,
                                                         )
                                                     }))

    primer_aligner = cr_align.SSWMultiAligner(primer_refs,
                                              match_score=VDJ_ANNOTATION_MATCH_SCORE,
                                              mismatch_penalty=VDJ_ANNOTATION_MISMATCH_PENALTY,
                                              gap_open_penalty=VDJ_ANNOTATION_GAP_OPEN_PENALTY,
                                              gap_extend_penalty=VDJ_ANNOTATION_EXTEND_PENALTY,
                                              report_cigar=True)

    primer_filter = lambda x, \
                    score=score_ratio, \
                    match=VDJ_ANNOTATION_MATCH_SCORE: \
                           filter_alignment(x, score,  int(len(x.reference.sequence) * \
                                                           VDJ_PRIMER_ANNOTATION_MIN_FRACTION_MATCHED), match)

    return primer_aligner, primer_filter


def collect_annotations(aligner, ref_seq, seq, filter_func):
    """Align a sequence against an SSWMultiAligner and return a list of Annotation objects
    for the alignments that passed filter_func.

    seq is the sequence to align. ref_seq is the sequence used to create the returned
    Annotation objects. ref_seq and seq will often be the same but seq could be a partially
    masked version of ref_seq. The two must have the same length.
    """
    assert(len(ref_seq) == len(seq))

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


def search_cdr3_signature(seq, v_region, j_region, v_frame):
    """Search for the CDR3 signature in a sequence.

    Args:
    - seq: the sequence of a contig
    - v_region: Annotation object for the V-Region
    - j_region: Annotation object for the J-Region
    - v_frame: reading frame (0, 1, 2)

    Returns:
    A tuple ((a,b), status_flag) where
       - (a,b) are the starting and ending (base) positions of the CDR3
       - status_flag is a string indicating what happened
    within seq or None if no CDR3 could be found.
    """

    flag = ''
    pos = None

    assert v_frame in [0, 1, 2]

    v_seq = v_region.sequence

    # Start of V within seq
    v_start = v_region.contig_match_start

    # Get the position of the last Cysteine within the V-Region
    v_amino_acids = [CODON_TO_AA[v_seq[i:(i+3)]] for i in range(v_frame, len(v_seq) - 3, 3)]
    if not 'C' in v_amino_acids:
        flag = 'GUIDED_NO_C_IN_V'
        return (pos, flag)

    last_c_idx = np.where(np.array(v_amino_acids) == 'C')[0][-1]
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

    j_frame = 3 - (j_start - v_start - v_frame) % 3
    if j_frame == 3:
        j_frame = 0
    j_amino_acids = [CODON_TO_AA[j_seq[i:(i+3)]] for i in range(j_frame, len(j_region.sequence) - 3, 3)]

    # Look for FG(X)G signature or WG(X)G signature
    fgxg_pos = None
    for idx in range(len(j_amino_acids) - 3):
        if (j_amino_acids[idx] == 'F' or j_amino_acids[idx] == 'W') and \
            j_amino_acids[idx + 1] == 'G' and j_amino_acids[idx + 3] == 'G':
            # The CDR3 includes the first F of the signature
            fgxg_pos = j_start + j_frame + idx * 3 + 3
            break

    if not fgxg_pos:
        flag = 'GUIDED_NO_FGXG'

    if fgxg_pos and last_c_pos < fgxg_pos:
        return ((last_c_pos, fgxg_pos), flag)

    return (None, flag)


def search_cdr3_signature_no_vj(seq, v_region=None, j_region=None):
    """Search the CDR3 signature in a sequence without guides from annotations.

    This could lead to more false positive signature hits than the guided version.

    Return value:
    A tuple (CDR3 DNA seq, CDR3 amino-acid seq, start position in seq, end position)
    """
    min_cdr3_aas = VDJ_MIN_CDR3_LEN / 3

    for frame in range(3):
        amino_acids = [CODON_TO_AA.get(seq[i:(i+3)],'') for i in range(frame, len(seq), 3)]

        fgxg_idx = None

        for idx in range(min_cdr3_aas, len(amino_acids) - 3):
            # First try to find the end motif
            if (amino_acids[idx] == 'F' or amino_acids[idx] == 'W') and \
                amino_acids[idx + 1] == 'G' and amino_acids[idx + 3] == 'G':
                # The CDR3 includes the first F of the signature
                fgxg_idx = idx
                fgxg_pos = frame + fgxg_idx * 3
                if j_region and (fgxg_pos >= j_region.contig_match_end or fgxg_pos < j_region.contig_match_start):
                    continue

                cys_idx = None
                # Find the Cysteine closer to the end but past a minimum
                # number of amino acids
                for idx in range(fgxg_idx - min_cdr3_aas, 0, -1):
                    if amino_acids[idx] == 'C':
                        cys_idx = idx
                        cys_pos = frame + cys_idx * 3
                        # Make sure the C of the CDR3 is inside the V
                        if not v_region or (cys_pos >= v_region.contig_match_start and cys_pos < v_region.contig_match_end):
                            break
                    elif amino_acids[idx] in STOP_CODONS:
                        break

                if cys_idx:
                    if fgxg_pos - cys_pos < VDJ_MAX_CDR3_LEN:
                        # include the start of the FGXG
                        cdr3_seq = seq[cys_pos:(fgxg_pos + 3)]
                        cdr3_aas = ''.join(amino_acids[cys_idx:(fgxg_idx + 1)])
                        return (cdr3_seq, cdr3_aas, cys_pos, fgxg_pos)
    return None

class AnnotatedContig(object):
    """A named sequence with a list of Annotations.
    """

    # Note: 'frame' is unused
    __slots__ = ['contig_name', 'sequence', 'quals', 'annotations', 'primer_annotations', 'barcode', 'clonotype', 'read_count', 'umi_count', 'unannotated_intervals', \
                 'filtered', 'cdr3_start', 'cdr3_stop', 'cdr3_seq', 'cdr3', 'cdr3_flag', 'start_codon_pos', 'stop_codon_pos', 'frame', 'aa_sequence', 'productive', \
                 'info_dict', 'is_cell', 'high_confidence' ]

    def __init__(self, name, sequence, annotations=[], primer_annotations=[],
                 barcode=None, is_cell=None,
                 clonotype=None,
                 read_count=None, umi_count=None,
                 unannotated_intervals=[], filtered=True,
                 info_dict={}, quals=None, high_confidence=False):

        """The info_dict can contain additional features that are not included in the
        objects features.
        """

        self.contig_name = name
        self.sequence = sequence.upper()
        self.quals = quals
        self.annotations = annotations
        self.primer_annotations = primer_annotations
        self.barcode = barcode
        self.clonotype = clonotype
        self.read_count = read_count
        self.umi_count = umi_count
        self.unannotated_intervals = unannotated_intervals
        # True if contig was selected for downstream analysis
        self.filtered = filtered
        self.cdr3_start = None
        self.cdr3_stop = None
        # sequence of CDR3 (self.sequence[cdr_pos[0]:cdr_pos[1]])
        self.cdr3_seq = None
        # CDR3 amino-acid sequence
        self.cdr3 = None
        self.start_codon_pos = None
        self.stop_codon_pos = None
        self.frame = None
        self.aa_sequence = None
        self.productive = None
        self.info_dict = info_dict
        # The cell-barcode was called as a cell
        self.is_cell = is_cell
        self.high_confidence = high_confidence


    @classmethod
    def from_dict(cls, annotation_dict, reference):
        """Create an AnnotatedContig object from a dict.

        Does not perform any checks on the self-consistency of the dict features.
        (Eg. do the annotation sequences agree with the contig's sequence,
        do the CDR3 positions make sense etc).
        """
        new_obj = cls(
            name=annotation_dict['contig_name'],
            quals=annotation_dict.get('quals', None),
            sequence=annotation_dict['sequence'],
            annotations=[Annotation.from_dict(hit, reference, contig_seq=annotation_dict['sequence']) for hit in annotation_dict['annotations']],
            primer_annotations=[Annotation.from_dict(hit, reference, contig_seq=annotation_dict['sequence']) for hit in annotation_dict['primer_annotations']],
            barcode=annotation_dict['barcode'],
            clonotype=annotation_dict['clonotype'],
            is_cell=annotation_dict['is_cell'],
            read_count=annotation_dict['read_count'],
            umi_count=annotation_dict['umi_count'],
            info_dict=annotation_dict['info'],
        )
        new_obj.cdr3_start = annotation_dict.get('cdr3_start', None)
        new_obj.cdr3_stop = annotation_dict.get('cdr3_stop', None)
        new_obj.productive = annotation_dict.get('productive', None)
        new_obj.start_codon_pos = annotation_dict.get('start_codon_pos', None)
        new_obj.stop_codon_pos = annotation_dict.get('stop_codon_pos', None)
        new_obj.frame = annotation_dict.get('frame', None)
        new_obj.aa_sequence = annotation_dict.get('aa_sequence', None)
        # Doesn't check that these are consistent with cdr_start, cdr_stop, and sequence.
        new_obj.cdr3_seq = annotation_dict.get('cdr3_seq', None)
        new_obj.cdr3 = annotation_dict.get('cdr3', None)
        new_obj.high_confidence = annotation_dict.get('high_confidence', None)

        return new_obj

    def to_dict(self):
        out = {}
        out['contig_name'] = self.contig_name
        out['sequence'] = self.sequence
        out['quals'] = self.quals
        out['annotations'] = [a.to_dict() for a in self.annotations]
        out['primer_annotations'] = [a.to_dict() for a in self.primer_annotations]
        out['barcode'] = self.barcode
        out['clonotype'] = self.clonotype
        out['is_cell'] = self.is_cell
        out['filtered'] = self.filtered
        out['read_count'] = self.read_count
        out['umi_count'] = self.umi_count
        out['cdr3_start'] = self.cdr3_start
        out['cdr3_stop'] = self.cdr3_stop
        out['cdr3_seq'] = self.cdr3_seq
        out['cdr3'] = self.cdr3
        out['start_codon_pos'] = self.start_codon_pos
        out['stop_codon_pos'] = self.stop_codon_pos
        out['frame'] = self.frame
        out['aa_sequence'] = self.aa_sequence
        out['productive'] = self.productive
        out['info'] = self.info_dict
        out['high_confidence'] = self.high_confidence
        return out

    def copy(self):
        """Returns a deep copy of this object"""
        return self.from_dict(self.to_dict())

    def hits(self):
        """Iterator over the annotations of the contig"""
        return iter(self.annotations)

    def get_region_hits(self, region_types):
        return [annotation for annotation in self.annotations if annotation.feature.region_type in region_types]

    def get_single_chain(self):
        """ Returns 'Multi' if ambiguous """
        chains = self.contig_chains()
        if len(chains) == 0:
            return None
        elif self.is_single_chain():
            return list(chains)[0]
        else:
            return 'Multi'

    def get_single_gene_display_name(self, region_types):
        """ Returns first if ambiguous """
        hits = self.get_region_hits(region_types)
        if len(hits) == 0:
            return None
        else:
            return hits[0].feature.display_name

    def is_single_chain(self):
        return len(self.contig_chains()) == 1

    def contig_chains(self):
        return set([annotation.feature.chain for annotation in self.annotations])

    def has_cdr(self):
        return not self.cdr3_seq is None and len(self.cdr3_seq) > 0

    def get_cdr_seq(self):
        """DNA sequence of CDR3.
        """
        if self.has_cdr():
            return self.cdr3_seq
        return None

    def get_cdr_aa_seq(self):
        """Amino acid sequence of CDR3.
        """
        if self.has_cdr():
            return self.cdr3
        return None

    def clonotype_seq(self):
        """Clonotype-defining sequence.

        Returns a tuple (chain, <chain_name>_<CDR3_sequence>). If this contig
        does not match a single chain, or does not have a CDR3, returns (None, None).
        """
        if not self.is_single_chain():
            return None, None
        if not self.has_cdr():
            return None, None
        chain = list(self.contig_chains())[0]
        cdr_seq = chain + '_' + self.get_cdr_seq()
        return (chain, cdr_seq)

    def has_full_length_vj_hit(self):
        has_full_len_v_hit = any([annotation.feature.region_type in VDJ_V_FEATURE_TYPES and \
                                  annotation.annotation_match_start == 0 for annotation in self.annotations])

        # The -2 allows for slop on the 3' end of J in the face of alignment clipping
        has_full_len_j_hit = any([annotation.feature.region_type in VDJ_J_FEATURE_TYPES  and
                                  annotation.annotation_match_end >= annotation.annotation_length - 2 \
                                  for annotation in self.annotations])
        return has_full_len_v_hit and has_full_len_j_hit

    def get_vj_quals(self):
        if self.quals is None:
            return None
        v_regions = [annotation for annotation in self.annotations if \
                     annotation.feature.region_type in VDJ_V_FEATURE_TYPES and \
                     annotation.annotation_match_start == 0]
        if not v_regions:
            return None
        v_region = v_regions[0]

        j_regions = [annotation for annotation in self.annotations if \
                     annotation.feature.region_type in VDJ_J_FEATURE_TYPES and
                     annotation.annotation_match_end == annotation.annotation_length]
        if not j_regions:
            return None
        j_region = j_regions[0]

        quals = np.array([ord(q) - VDJ_QUAL_OFFSET for q in self.quals[v_region.contig_match_start:j_region.contig_match_end]])
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


    def annotate_features(self, feature_types, ssw_multi_aligners, ssw_filter_funcs):
        """Add sequence annotations by comparing against a list of reference segments.

        Args:
        - feature_types, ssw_multi_aligners, ssw_filter_funcs: Matched aligners and filters
        for a list of feature types (see setup_feature_aligners)

        Return value: A list of Annotation objects.
        """

        masked_sequence = self.sequence
        masked_annotations = []
        vj_hits = [] # Keep track of the V and J annotations

        # V, J, C will be annotated in that order. Skip UTRs and Ds, these will be handled
        # conditionally on the V, J, C annotations.
        for region_type, aligner, align_filter in zip(feature_types, ssw_multi_aligners, ssw_filter_funcs):
            if not region_type in VDJ_V_FEATURE_TYPES + VDJ_J_FEATURE_TYPES + VDJ_C_FEATURE_TYPES:
                continue

            # Align against the masked sequence, but create Annotation objects that track the
            # original sequence.
            annotations = collect_annotations(aligner, self.sequence,
                                              masked_sequence, align_filter)

            if annotations:
                # Give preference to V and C annotations that are intact at the 5' end. This
                # will effectively give preference to full length Vs/Cs that are
                # more consistent with the Vs/Js annotated so far.
                if region_type in VDJ_V_FEATURE_TYPES + VDJ_C_FEATURE_TYPES:
                    sort_key = lambda x:(x.annotation_match_start == 0, x.score)
                else:
                    sort_key = lambda x: x.score

                best = max(annotations, key=sort_key)

                masked_annotations.append(best)

                if region_type in VDJ_V_FEATURE_TYPES or region_type in VDJ_J_FEATURE_TYPES:
                    vj_hits.append(best)

                # Mask everything before the end of this annotation.
                # contig_match_end is inclusive.
                masked_sequence = tk_seq.mask(self.sequence, best.contig_match_end, len(self.sequence))

        # Now try to find a UTR alignment before the beginning of the first thing you've annotated
        utr_aligner = ssw_multi_aligners[feature_types.index(VDJ_5U_FEATURE_TYPES[0])]
        utr_filter = ssw_filter_funcs[feature_types.index(VDJ_5U_FEATURE_TYPES[0])]
        min_start = np.min([anno.contig_match_start for anno in masked_annotations]) if masked_annotations else len(self.sequence)
        annotations = collect_annotations(utr_aligner, self.sequence,
                                          tk_seq.mask(self.sequence, 0, min_start), utr_filter)
        if annotations:
            # Keep best scoring one (doesn't make sense to priorite full length UTRs).
            best_utr = max(annotations, key=lambda x:x.score)
            # If you did find a V region, make sure the gene name of the UTR matches.
            # Otherwise, just pick the best scoring UTR.
            if vj_hits and vj_hits[0].feature.region_type in VDJ_V_FEATURE_TYPES:
                v_gene = vj_hits[0].feature.display_name
                equally_good = [anno for anno in annotations if anno.score == best_utr.score and anno.feature.display_name == v_gene]
                if equally_good:
                    masked_annotations.append(equally_good[0])
            else:
                masked_annotations.append(best_utr)

        # Only search for D if V and J were found, in order to avoid false hits.
        if len(vj_hits) == 2:
            assert(vj_hits[0].feature.region_type in VDJ_V_FEATURE_TYPES)
            assert(vj_hits[1].feature.region_type in VDJ_J_FEATURE_TYPES)

            v_end = vj_hits[0].contig_match_end
            j_start = vj_hits[1].contig_match_start
            # Also require a D of the same chain-type. This will also make sure that we
            # don't add a D on a chain that's not supposed to have a D (like TRA).
            if vj_hits[0].feature.chain == vj_hits[1].feature.chain:
                d_aligner = ssw_multi_aligners[feature_types.index(VDJ_D_FEATURE_TYPES[0])]
                d_filter = ssw_filter_funcs[feature_types.index(VDJ_D_FEATURE_TYPES[0])]
                annotations = collect_annotations(d_aligner, self.sequence,
                                                  tk_seq.mask(self.sequence, v_end, j_start), d_filter)
                annotations = [anno for anno in annotations if anno.feature.chain == vj_hits[0].feature.chain]
                if annotations:
                    masked_annotations.append(max(annotations, key=lambda x:x.score))

        return masked_annotations


    def annotate_features_by_group(self, ssw_multi_aligner, alignment_filter):
        """Similar to annotate_features but uses a joint SSWMultiAligner object.

        ssw_multi_aligner can contain references of different types. The sequence
        will be aligned against all of these and then the best alignment by type
        will be returned. The same filter function will be used for all feature types.
        """

        alignment_results = ssw_multi_aligner.align(self.sequence)

        annotations = []

        for alignment_result in alignment_results:
            alignment = alignment_result.alignment

            if alignment is None:
                continue
            if not alignment_filter(alignment_result):
                continue

            feature = alignment_result.reference.metadata['feature']

            annotation = Annotation(feature=feature,
                                    cigar=alignment.cigar_string,
                                    score=alignment.score,
                                    annotation_length=len(feature.sequence),
                                    annotation_match_start=alignment.ref_begin,
                                    annotation_match_end=alignment.ref_end + 1,
                                    contig_match_start=alignment.query_begin,
                                    contig_match_end=alignment.query_end + 1,
                                    mismatches=[],
                                    contig_seq=self.sequence)

            annotations.append(annotation)

        # Take top scoring hit for each region type
        top_annotations = []

        if top_annotations:
            by_region = lambda x: x.feature.region_type
            for region_type, group in itertools.groupby(sorted(annotations, key=by_region), by_region):
                top_annotations.append(max(group, key=lambda x: x.score))

        return top_annotations


    def get_unannotated_intervals(self):
        """ Return a list of Annotation objects corresponding to unannotated regions on the contig """

        unannotated_intervals = []
        annotation_indicator = np.zeros(len(self.sequence))

        for annotation in self.annotations:
            annotation_indicator[annotation.contig_match_start:annotation.contig_match_end] = 1

        interval_start = 0

        for annotated, region_iter in itertools.groupby(annotation_indicator, lambda x: x == 1):
            region = list(region_iter)

            if not annotated:
                feature = vdj_reference.create_dummy_feature(display_name='UNANNOTATED',
                                                             region_type='UNANNOTATED',
                                                             sequence=None)
                unannotated_intervals.append(Annotation(feature=feature,
                                                        cigar=None,
                                                        score=0,
                                                        annotation_length=len(region),
                                                        annotation_match_start=0,
                                                        annotation_match_end=len(region),
                                                        contig_match_start=interval_start,
                                                        contig_match_end=interval_start + len(region),
                                                        mismatches=[],
                                                    ))

            interval_start += len(region)

        return unannotated_intervals

    def get_annotations_bed(self):
        """
        Convert this contig's annotations to a list of BED entries

        Returns:
            list of str: each element is a BED file line (no newline)
        """

        bed_entries = []

        for anno in self.annotations + self.unannotated_intervals:
            bed_entries.append('\t'.join([self.contig_name,
                                          str(anno.contig_match_start),
                                          str(anno.contig_match_end),
                                          '%s_%s' % (anno.feature.display_name, anno.feature.region_type),
                                      ]))

        return bed_entries

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

            if has_v_start and v_seq[0:3] in START_CODONS:
                # Full V and V begins with a start codon (because it's an L-REGION+V-REGION).
                # Force the frame to 0.
                # Search for a CDR3 in frame 0.
                (cdr_pos, flag) = search_cdr3_signature(seq, v_region, j_region, 0)
                flags.append(flag)
                flags.append('FULL_V_HAS_START')
                self.start_codon_pos = v_start

            elif has_v_start:
                # Full V but the first annotated codon is not a start.
                # Look for an alternative start nearby that preserves the V frame.
                for i in xrange(v_start - 3*START_CODON_SLOP, 1+v_start+3*START_CODON_SLOP, 3):
                    if i > 0 and (i+3) <= len(seq) and seq[i:(i+3)] in START_CODONS:
                        (cdr_pos, flag) = search_cdr3_signature(seq, v_region, j_region, 0)
                        self.start_codon_pos = i
                        flags.append(flag)
                        flags.append('FULL_V_ALT_START')
                if self.start_codon_pos is None:
                    flags.append('FULL_V_NO_START')

            if not has_v_start or self.start_codon_pos is None:
                # Either we don't contain the start of V or we didn't find a valid start codon above
                # Look for a CDR3 sequence in all frames
                cdr_pos = None

                for frame in [0, 1, 2]:
                    ((cdr_pos), flag) = search_cdr3_signature(seq, v_region, j_region, frame)

                    if cdr_pos:
                        flags.append(flag)
                        break
                if cdr_pos is None:
                    flags.append('FAILED_UNGUIDED_SEARCH')

            if cdr_pos and cdr_pos[1] - cdr_pos[0] < VDJ_MAX_CDR3_LEN:
                self.cdr3_start = cdr_pos[0]
                self.cdr3_stop = cdr_pos[1]
                self.cdr3_seq = seq[cdr_pos[0]:cdr_pos[1]]
                self.cdr3 = ''.join([CODON_TO_AA[self.cdr3_seq[i:(i+3)]] for i in xrange(0, len(self.cdr3_seq), 3)])

                assert (cdr_pos[1] - cdr_pos[0]) % 3 == 0
            else:
                if cdr_pos is None:
                    flags.append('NO_CDR3')
                else:
                    flags.append('CDR3_TOO_LONG:%d' % (cdr_pos[1] - cdr_pos[0]))

        if not self.cdr3:
            # Either this didn't have both a V and a J, or the annotation-guided search failed to give a valid CDR3.
            # Run the unguided CDR3 search.
            v_region = v_regions[0] if v_regions else None
            j_region = j_regions[0] if j_regions else None

            res = search_cdr3_signature_no_vj(self.sequence, v_region, j_region)
            if not res is None:
                (cdr3_seq, cdr3_aas, cys_pos, fgxg_pos) = res
                self.cdr3 = cdr3_aas
                self.cdr3_seq = cdr3_seq
                self.cdr3_start = cys_pos
                self.cdr3_stop = fgxg_pos + 3 # End position is the start position of last amino acid + 3
                flags.append('FOUND_CDR3_UNGUIDED')

        if self.cdr3:
            cdr3_frame = self.cdr3_start % 3 # frame wrt start of sequence

            if self.start_codon_pos is None and len(v_regions) == 0:
                # We don't have a V or a start codon.
                # De novo search for the start codon in the CDR3 frame;
                # Get the leftmost possible start codon match
                for i in xrange(cdr3_frame, len(self.sequence), 3):
                    if self.codon(i) in START_CODONS:
                        self.start_codon_pos = i
                        flags.append('NO_V_BUT_FOUND_START')
                        break

            if self.start_codon_pos is not None:
                # If we have a start codon, try to find a stop codon in-frame
                for i in xrange(self.start_codon_pos, len(self.sequence), 3):
                    if self.codon(i) in STOP_CODONS:
                        self.stop_codon_pos = i
                        break

            # Determine productivity
            if self.has_full_length_vj_hit():
                has_start = self.start_codon_pos is not None
                cdr3_in_frame = has_start and self.cdr3_start is not None and \
                                (self.start_codon_pos % 3) == (self.cdr3_start % 3)
                vj_nostop = has_start and \
                            all([self.codon(i) not in STOP_CODONS for i in xrange(self.start_codon_pos,
                                                                                  j_region.contig_match_end, 3)])
                self.productive = has_start and cdr3_in_frame and vj_nostop

                if not has_start:
                    flags.append('NO_START')
                if has_start and self.cdr3_start is not None and not cdr3_in_frame:
                    flags.append('CDR3_OUT_OF_FRAME')
                if not vj_nostop:
                    flags.append('VJ_STOP')

            elif any([self.codon(i) in STOP_CODONS for i in xrange(self.cdr3_start, self.cdr3_stop, 3)]):
                # We don't know where the transcript starts and ends so we'll
                # be cautious about calling productivity.
                flags.append('CDR3_STOP')
                self.productive = False

        # Translate the entire reading frame if possible
        if self.start_codon_pos is not None:
            self.aa_sequence = ''.join([CODON_TO_AA[self.codon(i)] for i in xrange(self.start_codon_pos, len(self.sequence) - 2, 3)])

        self.cdr3_flag = '|'.join([f for f in flags if f is not None and len(f) > 0])



    def codon(self, i):
        return self.sequence[i:(i+3)] if i + 3 <= len(self.sequence) else ''

    def contains_annotations(self, other_contig):
        """True if this contig contains all the annotations of the other_contig.
        We only check the gene_name to test if the two annotations are the same.
        """

        for other_annotation in other_contig.hits():
            if not any([hit.feature.gene_name == other_annotation.feature.gene_name for hit in self.annotations]):
                return False
        return True

    def is_exact_vj_hit(self, other_contig):
        """True if this contig's VJ region is exactly contained in the other_contig.
        """

        v_regions = list(self.get_region_hits(VDJ_V_FEATURE_TYPES))
        j_regions = list(self.get_region_hits(VDJ_J_FEATURE_TYPES))
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
        features_to_output = VDJ_V_FEATURE_TYPES + VDJ_J_FEATURE_TYPES + VDJ_C_FEATURE_TYPES + VDJ_D_FEATURE_TYPES

        annotation_strs = []
        for anno in self.annotations:
            if anno.feature.region_type in features_to_output:
                annotation_strs.append((anno.contig_match_start,
                                        anno.annotation_match_start,
                                        anno.contig_match_end - anno.contig_match_start,
                                        anno.annotation_length,
                                        anno.feature.display_name))
        annotation_strs = ['{}'.format(x[4]) for x in \
                           sorted(annotation_strs, key=lambda x:x[0])]
        return ';'.join(annotation_strs)

    def __str__(self):
        annotation_str = ';'.join([str(anno) for anno in self.annotations])
        return 'Contig {}: {}'.format(self.contig_name, annotation_str)

    def __repr__(self):
        return self.to_dict().__repr__()


class Annotation(object):
    """Annotation of a sequence against a reference of segment sequences."""

    __slots__ = ['feature', 'cigar', 'score', 'annotation_length', 'annotation_match_start', 'annotation_match_end', 'contig_match_start', 'contig_match_end', 'mismatches', 'sequence']


    def __init__(self, feature, cigar, score,
                 annotation_length, annotation_match_start, annotation_match_end,
                 contig_match_start, contig_match_end, mismatches, contig_seq=None):
        """ Args: feature - vdj.reference.VdjAnnotationFeature object """

        self.feature = feature
        self.cigar = cigar
        self.score = score
        self.annotation_length = annotation_length
        self.annotation_match_start = annotation_match_start
        self.annotation_match_end = annotation_match_end
        self.contig_match_start = contig_match_start
        self.contig_match_end = contig_match_end
        self.mismatches = mismatches

        if contig_seq is None:
            self.sequence = None
        else:
            self.sequence = contig_seq[contig_match_start:contig_match_end].upper()

    @classmethod
    def from_alignment_result(cls, alignment_result, sequence):
        alignment = alignment_result.alignment

        feature = alignment_result.reference.metadata['feature']

        # Annotation object points to the input sequence (pre-masking)
        annotation = cls(feature=feature,
                         cigar=alignment.cigar_string,
                         score=alignment.score,
                         annotation_length=len(feature.sequence),
                         annotation_match_start=alignment.ref_begin,
                         # SSW returns inclusive ends - convert to open
                         annotation_match_end=alignment.ref_end + 1,
                         contig_match_start=alignment.query_begin,
                         contig_match_end=alignment.query_end + 1,
                         mismatches=[],
                         contig_seq=sequence)

        return annotation


    @staticmethod
    def _inferred_values():
        return ['chain_type', 'sequence']

    def annotate_mismatches(self, contig_sequence, feature_sequence):
        """
        Mark mismatches between the contig and the reference.
        Args:
            contig_sequence (str): contig sequence
            feature_sequence (str): the feature that the annotation derives from
        Returns:
            list of dict: list of annotations for mismatches (each as a dict)
        """

        mismatches = []

        cigar_tuples = cr_align.get_cigar_tuples(self.cigar)

        # Generate annotations for Indels
        ## Track position and all operations
        start_position = 0
        all_operations = []

        mismatch_annotation_json_keys = ['region_type', 'contig_match_start', 'contig_match_end']

        for length, category in cigar_tuples:
            length = int(length)
            all_operations.extend(length * [category])

            if category == 'D':
                # By convention, make deletion length always equal to 1. Really the length of the deletion
                # in terms of contig coordinates is 0, but this makes the deletions invinsible in the visualization.
                info = dict(zip(mismatch_annotation_json_keys, [category, start_position, start_position + 1]))
                # This is the only case where the length cannot be inferred from the start and end positions
                info['deletion_length'] = length
                mismatches.append(info)
            elif category == 'I':
                mismatches.append(dict(zip(mismatch_annotation_json_keys, [category, start_position, start_position + length])))

            if category != 'D':
                start_position += length

        # Generate annotations for one base mismatches
        contig_position = 0
        annotation_position = self.annotation_match_start

        for category in all_operations:
            if category == 'S':
                contig_position += 1

            elif category == 'M':
                if contig_sequence[contig_position].upper() != feature_sequence[annotation_position].upper():
                    mismatches.append(dict(zip(mismatch_annotation_json_keys, ['MISMATCH', contig_position, contig_position + 1])))
                annotation_position += 1
                contig_position += 1

            elif category == 'I':
                contig_position += 1

            elif category == 'D':
                annotation_position += 1

        return mismatches

    @classmethod
    def from_dict(cls, annotation_dict, reference, contig_seq=None):
        return cls(feature=vdj_reference.convert_dict_to_vdj_feature(annotation_dict['feature'], reference),
                   cigar=annotation_dict['cigar'],
                   score=annotation_dict.get('score', float('nan')),
                   annotation_length=annotation_dict['annotation_length'],
                   annotation_match_start=annotation_dict['annotation_match_start'],
                   annotation_match_end=annotation_dict['annotation_match_end'],
                   contig_match_start=annotation_dict['contig_match_start'],
                   contig_match_end=annotation_dict['contig_match_end'],
                   mismatches=annotation_dict.get('mismatches'),
                   contig_seq=contig_seq)

    def to_dict(self):
        out_dict = dict((slot, self.__getattribute__(slot)) for slot in self.__slots__)
        out_dict['feature'] = vdj_reference.convert_vdj_feature_to_dict(self.feature)
        for val in self._inferred_values():
            out_dict.pop(val, None)
        return out_dict

    def copy(self):
        return self.from_dict(self.to_dict(), self.sequence)

    def __eq__(self, other):
        """Two annotations are the same if they have the same sequence.
        If the sequence is missing, rely on the annotated gene name."""
        if not isinstance(other, self.__class__):
            return False
        if not self.sequence is None and self.sequence == other.sequence:
            return True
        return self.feature.gene_name == other.feature.gene_name and self.feature.allele_name == other.feature.allele_name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return '{}(cigar:{})'.format(self.feature.display_name, self.cigar)

    def __repr__(self):
        return self.to_dict().__repr__()


class CellContigs(object):
    """A  list of AnnotatedContigs (eg. representing a clonotype or all the
    contigs in a cell)
    """
    def __init__(self, name, contigs, info_dict={}):
        self.name = name
        self.chains = []
        for contig in contigs:
            self.chains.append(contig)
        self.info_dict = info_dict

    @classmethod
    def empty(cls, name):
        return cls(name, [], {})

    def contigs(self):
        return iter(self.chains)

    def is_paired(self, require_full_len=True, require_productive=True, require_high_conf=True):
        """True if this cell has contigs matching some VDJ gene pair"""
        good_pairs = [set(p.split('_')) for p in VDJ_GENE_PAIRS]
        chains=set()
        cl = self.clonotype_tuple(require_full_len=require_full_len,
                                  require_productive=require_productive,
                                  require_high_conf=require_high_conf)
        for cdr in cl:
            chain = cdr.split('_')[0]
            chains.add(chain)
        return chains in good_pairs

    def clonotype_tuple(self, require_full_len=True,
                        require_productive=True,
                        require_high_conf=True):
        """Tuple of unique CDR3s across all productive contigs

        - require_productive: Contigs (aka 'chains') need to be productive (productive=True)
        - require_full_len: Contigs (aka 'chains') need to be full length (productive or not)
        - require_high_conf: Contigs (aka 'chains') need to be high confidence
        """
        cdrs = [c.clonotype_seq()[1] for c in self.chains if \
                (c.productive or not require_productive) and \
                (c.has_full_length_vj_hit() or not require_full_len) and \
                (c.high_confidence or not require_high_conf) and \
                c.has_cdr() and \
                c.is_single_chain()]

        # Sort for easy comparisons
        return tuple(sorted(list(set(cdrs))))

    def contains_contig(self, other_contig):
        """True if it has a chain that contains all the annotations
        of the given other contig.

        Args:
        - contig: an AnnotatedContig
        """

        return any([contig.contains_annotations(other_contig) for contig in self.chains])

    def has_exact_hit(self, other_contig):
        """True if this cell has some contig whose sequence is fully contained in the sequence
        of the other_contig. This is stricter than has_exact_vj_hit."""

        for contig in self.chains:
            if re.search(contig.sequence, other_contig.sequence):
                return True
        return False

    def has_exact_vj_hit(self, other_contig):
        """True if this cell has some contig whose V(D)J region is perfectly
        contained in the sequence of the other_contig."""

        if other_contig is None:
            return False

        for contig in self.chains:
            if contig.is_exact_vj_hit(other_contig):
                return True

        return False

    def has_exact_cdr3_hit(self, other_contig):
        """"True if this cell has some contig whose CDR3 matches perfectly
        the CDR3 of the other_contig.
        """

        if other_contig is None:
            return False

        if not other_contig.has_cdr():
            return False

        other_cdr = other_contig.get_cdr_seq()
        return any([contig.get_cdr_seq() == other_cdr for contig in self.chains])

    def to_dict_list(self):
        """Return as a list of AnnotatedContig dicts.

        Information from this object's info_dict is passed to the contigs.
        """
        out_dicts = []
        for contig in self.chains:
            new_contig = contig.copy()
            new_contig.info_dict.update(dict(self.info_dict))
            contig_out = new_contig.to_dict()
            contig_out['barcode'] = self.name
            out_dicts.append(contig_out)
        return out_dicts

    def copy(self):
        return CellContigs(self.name, [contig.copy() for contig in self.chains], dict(self.info_dict))

    def __str__(self):
        return 'Cell {}:\n{}'.format(self.name, '\n'.join([str(chain) for chain in self.chains]))

    def __repr__(self):
        return self.to_dict_list().__repr__()

def load_contig_list_from_json(json_file, reference_path):
    """ Returns a list of AnnotatedContig objects from an open json file """
    reference = vdj_reference.VdjReference(reference_path)
    return [AnnotatedContig.from_dict(x, reference) for x in json.load(json_file)]

def load_cell_contigs_from_json(json_file, reference_path, group_key, require_high_conf=True):
    """Returns a list of CellContig objects based on annotations in a json.

    The json is assumed to contain a list of AnnotatedContigs (in dict form).
    The contigs are sorted and grouped by group_key and each such group is put
    into a CellContig object.

    group_key must be 'barcode' or 'clonotype'
    """

    assert group_key in set(['barcode', 'clonotype'])
    annotations = load_contig_list_from_json(open(json_file), reference_path)

    cell_contigs = []

    key_func = lambda x: x.__getattribute__(group_key)
    anno_iter = itertools.groupby(sorted(annotations, key=key_func), key=key_func)
    for clonotype_name, contig_annotations in anno_iter:

        contigs = []
        for new_contig in contig_annotations:
            # Note, for consensus contigs is_cell=None
            if new_contig.is_cell is not False \
               and (new_contig.high_confidence or not require_high_conf):
                contigs.append(new_contig)

        if len(contigs) > 0:
            cell_contigs.append(CellContigs(clonotype_name, contigs))

    return cell_contigs

def save_annotation_list_json(out_file, contigs):
    tk_safe_json.dump_numpy(tk_safe_json.json_sanitize([c.to_dict() for c in contigs]),
                            out_file, pretty=True)

def save_contig_list_csv(csv, contigs, write_inferred=True):
    """ Write contigs to an open csv file """
    columns = [
        'barcode', 'is_cell', 'contig_id', 'high_confidence', 'length',
        'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene',
        'full_length', 'productive', 'cdr3', 'cdr3_nt',
        'reads', 'umis',
        'raw_clonotype_id', 'raw_consensus_id']

    if write_inferred:
        columns.extend(['inferred_clonotype_id', 'inferred_consensus_id'])

    save_annotation_list_csv(csv, contigs, columns)

def save_consensus_list_csv(csv, contigs):
    """ Write consensus contigs to an open csv file """
    columns = [
        'clonotype_id', 'consensus_id', 'length',
        'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene',
        'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis',
    ]
    save_annotation_list_csv(csv, contigs, columns)

def save_annotation_list_csv(csv, contigs, columns):
    """ Write AnnotatedContigs to an open csv file """
    col_set = set(columns)

    vdj_utils.write_csv_row(columns, csv)

    for contig in sorted(contigs, key=lambda x: (not x.is_cell, x.barcode, x.contig_name)):
        if contig.filtered:
            row = {
                'barcode': contig.barcode,
                'is_cell': contig.is_cell,
                'contig_id': contig.contig_name,
                'high_confidence': contig.high_confidence,
                'length': len(contig),
                'chain': contig.get_single_chain(),
                'v_gene': contig.get_single_gene_display_name(VDJ_V_FEATURE_TYPES),
                'd_gene': contig.get_single_gene_display_name(VDJ_D_FEATURE_TYPES),
                'j_gene': contig.get_single_gene_display_name(VDJ_J_FEATURE_TYPES),
                'c_gene': contig.get_single_gene_display_name(VDJ_C_FEATURE_TYPES),
                'full_length': contig.has_full_length_vj_hit(),
                'productive': contig.productive,
                'cdr3': contig.cdr3,
                'cdr3_nt': contig.cdr3_seq,
                'reads': contig.read_count,
                'umis': contig.umi_count,
                'raw_clonotype_id': contig.info_dict.get('raw_clonotype_id'),
                'raw_consensus_id': contig.info_dict.get('raw_consensus_id'),
                'inferred_clonotype_id': contig.info_dict.get('inferred_clonotype_id'),
                'inferred_consensus_id': contig.info_dict.get('inferred_consensus_id'),
                'clonotype_id': contig.clonotype,
                'consensus_id': contig.contig_name,
                'clonotype_frequency': contig.info_dict.get('clonotype_freq'),
                'clonotype_proportion': contig.info_dict.get('clonotype_prop'),
            }
            assert col_set.issubset(set(row.iterkeys()))
            vdj_utils.write_csv_row([row[k] for k in columns], csv)

def save_clonotype_info_csv(csv, consensus_contigs):
    """ Write a CSV containing clonotype info to an open csv file.
        Takes a list of AnnotatedContigs corresponding to consensuses"""
    clonotypes = defaultdict(dict)

    for contig in consensus_contigs:
        chain = contig.get_single_chain()
        if chain is None:
            continue

        clonotype_id = contig.clonotype
        if clonotype_id not in clonotypes:
            clonotypes[clonotype_id]['clonotype_id'] = clonotype_id
            clonotypes[clonotype_id]['members'] = set(contig.info_dict['cells'])
            clonotypes[clonotype_id]['frequency'] = contig.info_dict['clonotype_freq']
            clonotypes[clonotype_id]['proportion'] = contig.info_dict['clonotype_prop']
            clonotypes[clonotype_id]['cdr3s_aa'] = set()
            clonotypes[clonotype_id]['cdr3s_nt'] = set()
        else:
            assert clonotypes[clonotype_id]['members'] == set(contig.info_dict['cells'])
            assert clonotypes[clonotype_id]['frequency'] == contig.info_dict['clonotype_freq']
            assert clonotypes[clonotype_id]['proportion'] == contig.info_dict['clonotype_prop']

        clonotypes[clonotype_id]['cdr3s_aa'].add((chain, contig.cdr3))
        clonotypes[clonotype_id]['cdr3s_nt'].add((chain, contig.cdr3_seq))

    # Generate cdr3 annotation strings, chain:cdr3;chain:cdr3
    # sorted by (chain, cdr3)
    def get_cdr3_list_string(chain_cdr3s):
        return ';'.join(['%s:%s' % chain_cdr3 for chain_cdr3 in sorted(list(chain_cdr3s))])
    for clonotype in clonotypes.itervalues():
        clonotype['cdr3s_aa'] = get_cdr3_list_string(clonotype['cdr3s_aa'])
        clonotype['cdr3s_nt'] = get_cdr3_list_string(clonotype['cdr3s_nt'])

    # Sort by frequency, descending
    clonotypes = sorted(clonotypes.values(), key=lambda c: c['frequency'], reverse=True)

    columns = ['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt']
    col_set = set(columns)

    vdj_utils.write_csv_row(columns, csv)

    for row in clonotypes:
        assert col_set.issubset(set(row.iterkeys()))
        vdj_utils.write_csv_row([row[k] for k in columns], csv)


def get_clonotype_vdj_pair(sequence_ids, clonotype_tuple, vdj_gene_pairs):
    """Get the VDJ gene pair associated with the clonotype.
    If the gene pair is not in vdj_gene_pairs, it will return None.
    """
    genes = set([sequence_ids[seq].split('_')[0] for seq in list(clonotype_tuple)])
    genes = tuple(sorted(list(genes)))

    return vdj_gene_pairs.get(genes, None)


def report_clonotypes(reporter, prefix, cell_barcodes, clonotype_ids, sequence_ids,
                      barcode_contigs, bc_clonotype_assignments):
    """Group barcodes into clonotypes.

    Args:
    - reporter: VdjReporter object for reporting metrics
    - prefix: metrics prefix (must be in VDJ_CLONOTYPE_TYPES)
    - cell_barcodes: set of cell barcodes
    - clonotype_ids: dict from clonotype id to tuple of CDR ids
    - sequence_ids: dict from CDR id to sequence
    - barcode_contigs: list of CellContigs objects
    - bc_clonotype_assignments: dict from barcode to clonotype id

    Return value:
    A dict clonotype_id -> clonotype_info.
    Clonotype info is itself a dict with the following information:
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

    freq_metric = reporter._get_metric_attr('vdj_clonotype_freq', prefix)
    prop_metric = reporter._get_metric_attr('vdj_clonotype_prop', prefix)

    # Cells grouped by clonotype assignment. All cells left unassigned will be
    # in the "None" group.
    grouped_barcodes = itertools.groupby(sorted(barcode_contigs,
                                                key=lambda x:bc_clonotype_assignments.get(x.name, None)),
                                         key=lambda x: bc_clonotype_assignments.get(x.name, None))

    grouped_barcodes = sorted([(cl, list(b)) for (cl, b) in grouped_barcodes],
                              key=lambda x:len(x[1]), reverse=True)
    unassigned_metric = reporter._get_metric_attr('vdj_unassigned_clonotype_bc_frac', prefix)

    observed_barcodes = []

    vdj_gene_pairs = {tuple(sorted(pair.split('_'))):pair for pair in VDJ_GENE_PAIRS}

    out_clonotypes = {}

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
            clonotype_name = vdj_utils.format_clonotype_id(len(out_clonotypes), inferred=prefix == 'inferred')

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
                'clonotype_id': clonotype_name,
                'freq': freq,
                'prop': prop,
                'barcodes': barcode_names,
            }

            consensuses = {}
            for idx, (cdr_seq, members) in enumerate(out_contigs.iteritems()):
                assert len(members) > 0
                consensus_name = '{}_consensus_{}'.format(clonotype_name, idx + 1)
                consensus_info = {
                    'chain': cdr_seq.split('_')[0],
                    'cdr3': clonotype_chain_cdr3_translations[cdr_seq],
                    'cdr3_seq': cdr_seq.split('_')[1],
                    'cell_contigs': members,
                }
                consensuses[consensus_name] = consensus_info

            out_clonotype['consensuses'] = consensuses

            out_clonotypes[clonotype_name] = out_clonotype

        unassigned_metric.add(len(bc_contig_lists), filter=(cl is None))

        for gp in reporter.canonical_vdj_gene_pairs:
            # Due to the way we define clonotypes, paired here implies productive.
            # And due to the way we define productive, this also means full length.
            paired_cls_metric = reporter._get_metric_attr('vdj_paired_clonotype_frac', gp, prefix)
            bcs_in_paired_cls_metric = reporter._get_metric_attr('vdj_paired_clonotype_bc_frac', gp, prefix)

            diversity_metric = reporter._get_metric_attr('vdj_clonotype_diversity', gp, prefix)
            paired_diversity_metric = reporter._get_metric_attr('vdj_paired_clonotype_diversity', gp, prefix)

            if not cl is None:

                paired_cls_metric.add(1, filter=is_cl_paired and (gp == vdj_gene_pair or gp == MULTI_REFS_PREFIX))
                bcs_in_paired_cls_metric.add(len(bc_contig_lists), filter=is_cl_paired and (gp == vdj_gene_pair or gp == MULTI_REFS_PREFIX))

                if gp == vdj_gene_pair or gp == MULTI_REFS_PREFIX:
                    num_clonotypes_metric = reporter._get_metric_attr('vdj_clonotype_count', gp, prefix)
                    num_clonotypes_metric.add(1)

                    diversity_metric.add(cl, clonotype_count)
                    if is_cl_paired:
                        paired_diversity_metric.add(cl, clonotype_count)

    # These barcodes have no annotations at all and no clonotype assignments
    remaining_barcodes = cell_barcodes.difference(set(observed_barcodes))
    for gp in reporter.canonical_vdj_gene_pairs:
        bcs_in_paired_cls_metric = reporter._get_metric_attr('vdj_paired_clonotype_bc_frac', gp, prefix)
        bcs_in_paired_cls_metric.add(len(remaining_barcodes), filter=False)
    unassigned_metric.add(len(remaining_barcodes), filter=True)

    metric = reporter._get_metric_attr('cdrs_per_bc_histogram', prefix)
    metric.add_many(np.array([len(clonotype_ids[c]) for c in bc_clonotype_assignments.values()]))
    metric.add_many(np.zeros((len(remaining_barcodes),), dtype=np.int))

    metric = reporter._get_metric_attr('major_clonotype_bc_frac', prefix)
    if len(clonotype_counts) > 0:
        metric.set_value(np.max(clonotype_counts), len(cell_barcodes))
    else:
        metric.set_value(0, len(cell_barcodes))

    return out_clonotypes


def label_contigs_with_consensus(clonotypes, contigs, prefix):
    """Adds clonotype and consensus info to a list of AnnotatedContig objects.

    Args:
    - clonotypes: dict like the one output by report_clonotypes.
    - contigs: list of AnnotatedContig objects
    - prefix: prefix for clonotype ids (must be in VDJ_CLONOTYPE_TYPES)
    """
    assert prefix in VDJ_CLONOTYPE_TYPES
    contig_dict = {contig.contig_name:contig for contig in contigs}

    bc_assignments = {}
    # iterate over clonotypes
    for _, clonotype_info in clonotypes.iteritems():
        # keep track of barcode to clonotype assignments.
        for bc in clonotype_info['barcodes']:
            bc_assignments[bc] = clonotype_info['clonotype_id']
        # iterate over consensus sequences of the clonotype and annotate
        # each member-contig of each consensus
        for consensus_name, consensus_info in clonotype_info['consensuses'].iteritems():
            for mem in consensus_info['cell_contigs']:
                contig_dict[mem].info_dict['{}_consensus_id'.format(prefix)] = consensus_name

    # Now annotate contigs with their clonotype assignments. We don't do
    # this in the previous pass because we want to annotate contigs whose
    # barcode was assigned to a clonotype but which did not participate in
    # any consensuses (eg. some non-full-length contigs)
    for contig in contigs:
        if contig.barcode in bc_assignments:
            contig.info_dict['{}_clonotype_id'.format(prefix)] = bc_assignments[contig.barcode]
