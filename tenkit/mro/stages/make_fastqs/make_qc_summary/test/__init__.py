#!/usr/bin/env python
#
# Copyright (c) 2016 10x Genomics, Inc.  All rights reserved.
#
# Code for testing various make_qc_summary functionality

import tenkit.test as tk_test
import os
import json
from .. import shim_standard_qc_synonyms, ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP


class TestShimStandardArguments(tk_test.UnitTestBase):
    def setUp(self):
        super(TestShimStandardArguments, self).setUp()
        self.base_file = os.path.join(os.path.dirname(__file__), "analyze_run_pd.json")

    def test_shim_standard_qc_synonyms(self):
        with open(self.base_file, 'r') as base_file:
            base_dict = json.load(base_file)

        replaced_dict = shim_standard_qc_synonyms(base_dict, replace=True)
        for old_attr in ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP.keys():
            self.assertNotIn(old_attr, replaced_dict)
        for new_attr in ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP.values():
            if 'fract' in new_attr:
                self.assertIn(new_attr, replaced_dict)

        supplemented_dict = shim_standard_qc_synonyms(base_dict, replace=False)
        for old_attr in ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP.keys():
            if 'fraction' in old_attr:
                self.assertIn(old_attr, supplemented_dict)
        for new_attr in ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP.values():
            if 'fract' in new_attr:
                self.assertIn(new_attr, replaced_dict)

        # test a sample_qc
        replaced_sample_qc = replaced_dict['sample_qc']['20486']['all']
        self.assertNotIn('total_reads', replaced_sample_qc)
        self.assertNotIn('gems', replaced_sample_qc)
        self.assertIn('number_reads', replaced_sample_qc)
        self.assertIn('gems_detected', replaced_sample_qc)

        supplemented_sample_qc = supplemented_dict['sample_qc']['20486']['all']
        self.assertIn('total_reads', supplemented_sample_qc)
        self.assertIn('gems', supplemented_sample_qc)
        self.assertIn('number_reads', supplemented_sample_qc)
        self.assertIn('gems_detected', supplemented_sample_qc)
