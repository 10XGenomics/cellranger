#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Structify aggregated cell typing outputs."""

from cellranger.cr_io import recursive_hard_link_dict

__MRO__ = """
struct AggregatedCellTypes(
    csv     all_cell_types,
    json.gz all_cell_annotation_results,
)

stage STRUCTIFY_AGGREGATED_CELLTYPES(
    in  map<csv>                 all_cell_types,
    in  map<json.gz>             all_cell_annotation_results,
    out map<AggregatedCellTypes> library_cell_types,
    src py                       "stages/cas_cell_typing/structify_aggregated_celltypes",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    all_methods = (
        args.all_cell_types
        if args.all_cell_types
        else {} | args.all_cell_annotation_results if args.all_cell_annotation_results else {}
    )

    library_cell_types = {}
    for method in all_methods:
        library_cell_types[method] = {
            "all_cell_types": args.all_cell_types.get(method) if args.all_cell_types else None,
            "all_cell_annotation_results": (
                args.all_cell_annotation_results.get(method)
                if args.all_cell_annotation_results
                else None
            ),
        }

    outs.library_cell_types = recursive_hard_link_dict(library_cell_types)
