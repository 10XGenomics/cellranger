#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to assign a GEM well to a sample def.

This stage will be unnecessary once it is enforced that
GEM wells are properly assigned at every entry point to the pipeline.
"""


__MRO__ = """
stage _FORCE_SAMPLE_DEF_GEM_WELL(
    in  map[]  sample_def,
    in  int    gem_group,
    out map[]  sample_def,
    src py     "stages/multi/_force_sample_def_gem_well",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if args.sample_def is None:
        outs.sample_def = None
    else:
        outs.sample_def = []
        for sample_def in args.sample_def:
            sample_def["gem_group"] = args.gem_group
            outs.sample_def.append(sample_def)
