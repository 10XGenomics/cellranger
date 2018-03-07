#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Unit tests for cellranger.vdj.annotations.py.
#

from collections import namedtuple
import cellranger.align as cr_align
import cellranger.vdj.annotations as cr_annotations
import tenkit.test as tk_test
from cellranger.vdj.test import in_path
import cellranger.vdj.reference as vdj_reference

class TestAnnotations(tk_test.UnitTestBase):
    def setUp(self):
        pass


    def test_filter(self):
        SSWAlignment = namedtuple('SSWAlignment', ['query_begin', 'query_end', 'score', 'cigar_string'])
        alignment = SSWAlignment(100, 109, 20, '100S10M')
        res = cr_align.SSWAlignmentResult(None, alignment, 110)
        assert(cr_annotations.filter_alignment(res, 1.0, 5, 2.0))
        assert(not cr_annotations.filter_alignment(res, 1.0, 100))

        alignment = SSWAlignment(0, 0, 0, '100S')
        res = cr_align.SSWAlignmentResult(None, alignment, 100)
        assert(cr_annotations.filter_alignment(res, 0, 0))
        assert(not cr_annotations.filter_alignment(res, 0, 1))


    def test_setup_feature_aligner(self):
        score_ratios = {"5U":0.8,"C":0.8,"D":0.5,"J":0.8,"V":0.8}
        word_sizes = {"5U":2,"C":3,"D":3,"J":5,"V":6}
        test_ref = in_path('annotation_setup_test_ref')
        features, aligners, filters = cr_annotations.setup_feature_aligners(test_ref, score_ratios, word_sizes)

        # This tests that the lambda constants are not dynamically determined.
        #score = 2.0
        #word = 100

        seq = 'TTAAAAAAAATTTTCCCC'
        for t, al, f in zip(features, aligners, filters):
            alignments = cr_annotations.collect_annotations(al, seq, seq, f)
            if t == '5U' or t == 'V' or t == 'D':
                assert(alignments)
            else:
                # J doesn't have a match because we require at least 5 bases
                # C doesn't have a good enough hit
                assert(not alignments)


    def test_coordinates(self):
        """Test that coordinates for matches are 0-based, half-open.
        """
        score_ratios = {"5U":0.8,"C":0.8,"D":0.5,"J":0.8,"V":0.8}
        word_sizes = {"5U":2,"C":3,"D":3,"J":5,"V":6}
        test_ref = in_path('annotation_setup_test_ref')
        features, aligners, filters = cr_annotations.setup_feature_aligners(test_ref, score_ratios, word_sizes)

        aligner = [(al, f) for (t, al, f) in zip(features, aligners, filters) if t == 'V'][0]

        seq = 'AAAAAAAA'
        alignments = cr_annotations.collect_annotations(aligner[0], seq, seq, aligner[1])
        self.assertEqual(len(alignments), 1)
        anno = alignments[0]
        self.assertEqual(len(anno.mismatches), 0)
        self.assertEqual(anno.annotation_match_end, len(seq))
        self.assertEqual(anno.contig_match_end, len(seq))

        score_ratios = {"5U":0.8,"C":0.8,"D":0.5,"J":0.8,"V":0.5}
        word_sizes = {"5U":2,"C":3,"D":3,"J":5,"V":3}
        seq = 'AAATAAAA'
        features, aligners, filters = cr_annotations.setup_feature_aligners(test_ref, score_ratios, word_sizes)
        aligner = [(al, f) for (t, al, f) in zip(features, aligners, filters) if t == 'V'][0]
        alignments = cr_annotations.collect_annotations(aligner[0], seq, seq, aligner[1])
        anno = alignments[0]
        mismatches = anno.annotate_mismatches(seq, anno.feature.sequence)
        mismatch = {'region_type':'MISMATCH',
                    'contig_match_start':3,
                    'contig_match_end':4}
        self.assertEqual(mismatches[0], mismatch)

        test_ref = in_path('annotation_setup_test_ref')
        features, aligners, filters = cr_annotations.setup_feature_aligners(test_ref, score_ratios, word_sizes)

        aligner = [(al, f) for (t, al, f) in zip(features, aligners, filters) if t == 'V'][0]

        test_ref = in_path('annotation_setup_test_ref2')
        fasta = vdj_reference.get_vdj_reference_fasta(test_ref)
        features, aligners, filters = cr_annotations.setup_feature_aligners(test_ref, score_ratios, word_sizes)
        aligner = [(al, f) for (t, al, f) in zip(features, aligners, filters) if t == 'V'][0]

        with open(fasta, 'r') as f:
            true_seq = f.readlines()[1].strip()

        # Test deletion coordinates
        middle = len(true_seq) / 2
        seq = true_seq[0:middle] + true_seq[middle+3:]

        alignments = cr_annotations.collect_annotations(aligner[0], seq, seq, aligner[1])
        anno = alignments[0]
        mismatches = anno.annotate_mismatches(seq, anno.feature.sequence)
        mismatch = {'region_type':'D',
                    'contig_match_start':middle,
                    'contig_match_end':middle + 1, # convention
                    'deletion_length':3}
        self.assertEqual(mismatches[0], mismatch)

        # Test insertion coordinates
        seq = true_seq[0:middle] + 'AAA' + true_seq[middle:]

        alignments = cr_annotations.collect_annotations(aligner[0], seq, seq, aligner[1])
        anno = alignments[0]
        mismatches = anno.annotate_mismatches(seq, anno.feature.sequence)
        mismatch = {'region_type':'I',
                    'contig_match_start':middle,
                    'contig_match_end':middle + 3}
        self.assertEqual(mismatches[0], mismatch)
