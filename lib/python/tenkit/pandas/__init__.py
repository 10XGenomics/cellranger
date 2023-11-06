# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""Prevent pandas frm using eval when subsetting data frames.

We ran into a problem with a bug in glibc and haswell processors that causes
a lock in some situations. We first found it in bwa and fixed it by linking to
jmalloc. But now it is turning up in pandas in numexpr eval.
"""


# import only what we need, both for efficiency and to keep pylint relatively
# happy.
from pandas import DataFrame, set_option  # noqa: F401

set_option("compute.use_numexpr", False)
