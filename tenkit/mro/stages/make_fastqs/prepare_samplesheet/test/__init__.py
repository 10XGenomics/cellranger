#!/usr/bin/env python
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Code for testing prepare_samplesheet stage
#

import os
from functools import partial

import tenkit.test as tk_test
import csv
from .. import make_csv_from_specs


class TestLaneCount(tk_test.UnitTestBase):

    def setUp(self):
        super(TestLaneCount, self).setUp()
        self.input_dir = tk_test.in_path("lane")
        self.output_dir = tk_test.out_path("prepare_samplesheet")
        os.makedirs(self.output_dir)

    def get_layout_rows(self):
        with open(os.path.join(self.output_dir, "layout.csv"), 'r') as infile:
            reader = csv.reader(infile)
            return [line for line in reader]

    def test_make_csv_from_specs(self):
        f = make_csv_from_specs
        single_lane_spec = {
            "lanes": [1],
            "sample": "Sample",
            "indices": ["SI-GA-A03"]
        }
        miseq_xml = os.path.join(self.input_dir, "miseq.xml")
        f([single_lane_spec], miseq_xml, self.output_dir)
        rows = self.get_layout_rows()
        self.assertEqual(len(rows), 2, msg="Single row count check")
        self.assertEqual(rows[1][0], "1", msg="Lane check")
        self.assertEqual(rows[1][1], "Sample", msg="Sample check")
        self.assertEqual(rows[1][2], "SI-GA-A03", msg="index check")

        empty_lane_spec = {
            "sample": "Sample",
            "indices": ["SI-GA-A03"]
        }
        xten_xml = os.path.join(self.input_dir, "xten.xml")
        f([empty_lane_spec], xten_xml, self.output_dir)
        rows = self.get_layout_rows()
        self.assertEqual(len(rows), 9, msg="Eight row count check")
        self.assertEqual(rows[1][0], "1", msg="1 Lane check")
        self.assertEqual(rows[1][1], "Sample", msg="1 Sample check")
        self.assertEqual(rows[1][2], "SI-GA-A03", msg="1 index check")
        self.assertEqual(rows[5][0], "5", msg="5 Lane check")
        self.assertEqual(rows[5][1], "Sample", msg="5 Sample check")
        self.assertEqual(rows[5][2], "SI-GA-A03", msg="5 index check")
        self.assertEqual(rows[8][0], "8", msg="8 Lane check")
        self.assertEqual(rows[8][1], "Sample", msg="8 Sample check")
        self.assertEqual(rows[8][2], "SI-GA-A03", msg="8 index check")

        multi_specs = [
            {
                "lanes": [1],
                "sample": "Sample1",
                "indices": ["SI-GA-A01","SI-GA-A02"]
            },
            {
                "lanes": [2,3],
                "sample": "Sample2",
                "indices": ["SI-GA-A03"]
            }
        ]

        f(multi_specs, xten_xml, self.output_dir)
        rows = self.get_layout_rows()
        self.assertEqual(len(rows), 5, msg="complicated row count check")
        self.assertEqual(rows[1][0], "1")
        self.assertEqual(rows[1][1], "Sample1")
        self.assertEqual(rows[1][2], "SI-GA-A01")
        self.assertEqual(rows[2][0], "1")
        self.assertEqual(rows[2][1], "Sample1")
        self.assertEqual(rows[2][2], "SI-GA-A02")
        self.assertEqual(rows[3][0], "2")
        self.assertEqual(rows[3][1], "Sample2")
        self.assertEqual(rows[3][2], "SI-GA-A03")
        self.assertEqual(rows[4][0], "3")
        self.assertEqual(rows[4][1], "Sample2")
        self.assertEqual(rows[4][2], "SI-GA-A03")



if __name__ == '__main__':
    tk_test.run_tests()
