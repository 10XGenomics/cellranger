#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Determine if the cloupe was not generated because this is a single sample run of MULTI.

In this case we hard-link the single sample cloupe as the library-level cloupe.
"""

from shutil import copyfile

import cellranger.cr_io as cr_io

__MRO__ = """
stage CHOOSE_CLOUPE(
    in  cloupe library_cloupe,
    in  map<cloupe> sample_cloupe,
    out cloupe cloupe,
    src py     "stages/multi/choose_cloupe",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if args.library_cloupe is None:
        if args.sample_cloupe is None:
            outs.cloupe = None
            return

        sample_cloupe = next(iter(args.sample_cloupe.values()))
        if sample_cloupe is None:
            outs.cloupe = None
            return
        # changed this from a hardlink to a copy to avoid issues with moving/deleting directories etc
        copyfile(sample_cloupe, outs.cloupe)

    else:
        cr_io.hardlink_with_fallback(args.library_cloupe, outs.cloupe)
