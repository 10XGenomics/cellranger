#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""When merging metrics and files from multiple gem wells, the aggr pipeline is used.

to do a lot of the heavy lifting of generating combined molecule info, matrices etc.
This stage generates the aggr sample defs and aggr CSV that are passed to aggr to do this.
Each sample in the sample def specifies a different gem well.
"""
import cellranger.constants as cr_constants

__MRO__ = """
struct CountAggrSampleDef(
    string library_id,
    h5     molecule_h5,
)

stage MAKE_MULTI_GEM_RNA_AGGR_SAMPLE_DEFS(
    in  int[]                gem_groups,
    in  h5[]                 molecule_info,
    out CountAggrSampleDef[] sample_defs,
    out csv                  aggr_csv,
    src py                   "stages/multi/make_multi_gem_rna_aggr_sample_defs",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    assert len(args.gem_groups) == len(args.molecule_info)

    outs.sample_defs = []
    with open(outs.aggr_csv, "w") as outfile:
        print(
            f"{cr_constants.AGG_ID_FIELD},{cr_constants.AGG_H5_FIELD}",
            file=outfile,
        )

        for gem_group, molecule_info in zip(args.gem_groups, args.molecule_info):
            gem_group_label = f"gem_group_{gem_group}"
            # write CSV line
            print(f"{gem_group_label},{molecule_info}", file=outfile)

            # add sample def
            sample_def = {
                cr_constants.AGG_ID_FIELD: gem_group_label,
                cr_constants.AGG_H5_FIELD: molecule_info,
            }
            outs.sample_defs.append(sample_def)
