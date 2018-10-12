
#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

__MRO__ = """
stage CHOOSE_DIMENSION_REDUCTION(
    in  bool chemistry_batch_correction,
    out bool disable_run_pca,
    out bool disable_correct_chemistry_batch,
    src py   "stages/analyzer/choose_dimension_reduction",
)
"""

def main(args, outs):
    if args.chemistry_batch_correction is None or args.chemistry_batch_correction is False:
        outs.disable_run_pca = False
        outs.disable_correct_chemistry_batch = True
    else:
        outs.disable_run_pca = True
        outs.disable_correct_chemistry_batch = False
