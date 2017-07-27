#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import collections
import ctypes
import itertools
import re
import striped_smith_waterman.ssw_wrap as ssw

# A single reference sequence and associated metadata
SSWReference = collections.namedtuple('SSWReference', ['name', 
                                                       'sequence',
                                                       'min_score', 'min_match_length',
                                                       'metadata'])

SSWAlignmentResult = collections.namedtuple('SSWAlignmentResult', ['reference', 'alignment', 'query_len'])

class SSWMultiAligner(object):
    """ Store multiple references for speed (single setup/teardown per-query) """

    def __init__(self, ssw_refs, match_score=2, mismatch_penalty=2,
                 gap_open_penalty=3, gap_extend_penalty=1, report_cigar=True):
        self.references = []
        self.aligners = []
        for ssw_ref in ssw_refs:
            self.references.append(ssw_ref)
            self.aligners.append(ssw.Aligner(ref_seq=ssw_ref.sequence,
                                             report_cigar=report_cigar,
                                             mismatch=mismatch_penalty,
                                             match=match_score,
                                             gap_open=gap_open_penalty,
                                             gap_extend=gap_extend_penalty))

    def align(self, query_seq):
        """ Re-wrap ssw aligner for speedup w/ multiple refs;
        only initialize and dealloc the query sequence once-per-query, not once-per-ref.
        NOTE: this code was largely kept intact from striped_smith_waterman/src/ssw_wrap.py

        Args: ssw_references: list of SSWReference objects
        """
        if len(self.aligners) == 0:
            return []

        first_aligner = self.aligners[0]

        _query_seq = first_aligner._DNA_to_int_mat(query_seq, len(query_seq))
        profile = first_aligner.ssw_init(_query_seq, # Query seq in c type integers
                                         ctypes.c_int32(len(query_seq)), # Length of Queryseq in bytes
                                         first_aligner.mat, # Score matrix
                                         5, # Square root of the number of elements in mat
                                         2) # flag = no estimation of the best alignment score
        if len(query_seq) > 30:
            mask_len = len(query_seq) / 2
        else:
            mask_len = 15

        results = []
        for ssw_ref, aligner in itertools.izip(self.references, self.aligners):
            c_result = aligner.ssw_align(profile, # Query profile
                                         aligner.ref_seq, # Ref seq in c type integers
                                         ctypes.c_int32(len(aligner.ref_seq)), # Length of Refseq in bytes
                                         aligner.gap_open, # Absolute value of gap open penalty
                                         aligner.gap_extend, # absolute value of gap extend penalty
                                         1, # Bitwise FLAG for output values = return all
                                         0, # Score filter = return all
                                         0, # Distance filter = return all
                                         mask_len) # Distance between the optimal and suboptimal alignment

            # Transform the Cstructure into a python object if score and length match the requirements
            score = c_result.contents.score
            match_len = c_result.contents.query_end - c_result.contents.query_begin + 1

            if score >= ssw_ref.min_score and match_len >= ssw_ref.min_match_length:
                alignment = ssw.PyAlignRes(Res=c_result,
                                           query_len=len(query_seq),
                                           report_secondary=False,
                                           report_cigar=aligner.report_cigar)
            else:
                alignment = None

            results.append(SSWAlignmentResult(ssw_ref, alignment, len(query_seq)))

            # Free reserved space by ssw.init and ssw_init methods.
            aligner._align_destroy(c_result)

        first_aligner._init_destroy(profile)

        return results

def get_cigar_tuples(cigar_string):
    """
    Get a list of length, CIGAR operation tuples for alignment
    Args:
        cigar_string str): CIGAR string
    Returns:
        tuple of (int, str): tuple of (length CIGAR operation, CIGAR operation). Returns none if no cigar
    """
    cigar_numbers = re.split('[A-Za-z]+', cigar_string)[:-1]
    cigar_letters = re.split('[0-9]+', cigar_string)[1:]
    return zip([int(number) for number in cigar_numbers], cigar_letters)

def get_max_word_length(alignment):
    """ Get the longest M (match or mismatch) operation in a list of cigar tuples """
    cigar_tuples = get_cigar_tuples(alignment.cigar_string)
    word_lengths = [op_len for op_len, op in cigar_tuples if op == 'M']

    if len(word_lengths) == 0:
        return 0
    else:
        return max(word_lengths)
