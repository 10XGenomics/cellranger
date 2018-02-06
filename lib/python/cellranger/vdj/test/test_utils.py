#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Unit tests for cellranger.vdj.utils.py
#

import json
import StringIO
import tenkit.test as tk_test
import cellranger.vdj.utils as vdj_utils

class TestVdjUtils(tk_test.UnitTestBase):
    def setUp(self):
        pass


    def test_json_stream(self):
        s = '''[{"x": "abc"},
                {"x": "abc\\""},
                {"x": "abc\\\\"},
                {"x": "abc\\\\\\""},
                {"y": 123}
               ]'''

        self.assertEqual(json.loads(s), list(vdj_utils.get_json_obj_iter(StringIO.StringIO(s))))
