#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
from cellranger.vdj.constants import TEST_FILE_IN_DIR, TEST_FILE_OUT_DIR
import os.path

def in_path(filename):
    return os.path.join(TEST_FILE_IN_DIR, filename)

def out_path(filename):
    if not os.path.isdir(TEST_FILE_OUT_DIR):
        os.mkdir(TEST_FILE_OUT_DIR)
    return os.path.join(TEST_FILE_OUT_DIR, filename)
