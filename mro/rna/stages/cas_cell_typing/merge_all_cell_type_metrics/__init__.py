#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Get legacy cell type outs."""
import json

import martian

__MRO__ = """
stage MERGE_ALL_CELL_TYPE_METRICS(
    in  json cell_annotation_primary_metrics,
    in  json cell_annotation_analysis_metrics,
    out json merged_cell_annotation_metrics,
    src py   "stages/cas_cell_typing/merge_all_cell_type_metrics",
) using (
    volatile = strict,
) retain (
    merged_cell_annotation_metrics,
)
"""


def main(args, outs):
    if not args.cell_annotation_primary_metrics and not args.cell_annotation_analysis_metrics:
        martian.clear(outs)
        return

    merged_metrics = {}
    if args.cell_annotation_primary_metrics:
        with open(args.cell_annotation_primary_metrics) as f:
            merged_metrics.update(json.loads(f.read()))
    if args.cell_annotation_analysis_metrics:
        with open(args.cell_annotation_analysis_metrics) as f:
            merged_metrics.update(json.loads(f.read()))

    with open(outs.merged_cell_annotation_metrics, "w") as f:
        json.dump(merged_metrics, f)
