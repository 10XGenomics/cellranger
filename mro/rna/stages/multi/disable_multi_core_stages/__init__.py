#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to disable stages in SC_MULTI_CORE."""

import martian

__MRO__ = """
stage DISABLE_MULTI_CORE_STAGES(
    in  bool is_pd,
    in  bool disable_gex,
    in  bool disable_multi_count,
    in  bool disable_count,
    in  bool no_secondary_analysis,
    in  bool skip_cell_annotation,
    out bool disable_sample_cas_celltyping,
    out bool disable_library_cas_celltyping,
    src py   "stages/multi/disable_multi_core_stages",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    disable_cell_typing_totally = (
        bool(args.disable_gex)
        or bool(args.no_secondary_analysis)
        or bool(args.skip_cell_annotation)
    )
    if args.disable_gex:
        martian.alarm(
            "Not running cell annotation as there is no gene expression library in sample!"
        )
    if args.no_secondary_analysis:
        martian.alarm("Not running cell annotation as secondary analysis has been disabled!")

    outs.disable_sample_cas_celltyping = disable_cell_typing_totally or bool(
        args.disable_multi_count
    )
    outs.disable_library_cas_celltyping = disable_cell_typing_totally or bool(args.disable_count)

    # dont run cell typing at library level if it is run at sample level
    outs.disable_library_cas_celltyping = outs.disable_library_cas_celltyping or (
        not outs.disable_sample_cas_celltyping
    )
