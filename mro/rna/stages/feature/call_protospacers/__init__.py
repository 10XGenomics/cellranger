#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#
"""Assign protospacers to cells.

Identify which cells express which protospacers above background.
"""

import csv
import json
import math
import operator
from collections import defaultdict
from collections.abc import Iterable
from functools import reduce
from statistics import median, pstdev

import cellranger.rna.library as rna_library
from cellranger.feature.feature_assigner import GuideAssigner
from cellranger.feature.utils import write_json_from_dict
from cellranger.matrix import CountMatrix

__MRO__ = """
stage CALL_PROTOSPACERS(
    in  h5   filtered_feature_counts_matrix,
    in  int  min_crispr_umi_threshold,
    out csv  protospacer_calls_summary,
    out csv  protospacer_calls_per_cell,
    out json protospacer_call_metrics_json,
    out json cells_per_protospacer,
    out json protospacer_umi_thresholds_json,
    out csv  protospacer_umi_thresholds_csv,
    src py   "stages/feature/call_protospacers",
) split (
    in  int  chunk_start,
    in  int  chunk_end,
    out json chunk_cells_per_protospacer,
    out json chunk_protospacer_umi_thresholds,
) using (
    volatile = strict,
)
"""


def split(args):
    num_features = CountMatrix.load_feature_ref_from_h5_file(
        args.filtered_feature_counts_matrix
    ).get_count_of_feature_type(rna_library.CRISPR_LIBRARY_TYPE)
    _num_features, num_barcodes, nnz = CountMatrix.load_dims_from_h5(
        args.filtered_feature_counts_matrix
    )
    mem_gib = CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=1)
    print(f"{num_features=},{num_barcodes=},{nnz=},{mem_gib=}")

    if num_barcodes == 0:
        return {"chunks": []}

    num_chunks = min(100, num_features, math.ceil(max(num_features / 100, num_barcodes / 10_000)))
    num_features_per_chunk = math.ceil(num_features / num_chunks)
    return {
        "chunks": [
            {
                "chunk_start": start,
                "chunk_end": min(num_features, start + num_features_per_chunk),
                "__mem_gb": mem_gib,
            }
            for start in range(0, num_features, num_features_per_chunk)
        ],
        "join": {"__mem_gb": mem_gib},
    }


def main(args, outs):
    filtered_feature_counts_matrix = (
        CountMatrix.load_h5_file(args.filtered_feature_counts_matrix)
        .select_features_by_type(rna_library.CRISPR_LIBRARY_TYPE)
        .select_features(range(args.chunk_start, args.chunk_end))
    )

    guide_assigner = GuideAssigner(
        matrix=filtered_feature_counts_matrix,
        min_crispr_umi_threshold=args.min_crispr_umi_threshold,
    )

    guide_assigner.assignments = guide_assigner.get_feature_assignments()
    write_json_from_dict(guide_assigner.get_cells_per_feature(), outs.chunk_cells_per_protospacer)
    write_json_from_dict(
        guide_assigner.compute_assignment_metadata().umi_thresholds,
        outs.chunk_protospacer_umi_thresholds,
    )


def dump_json(obj, filename, *, indent=4):
    """Dump a JSON object to a file with a trailing newline."""
    with open(filename, "w") as f:
        json.dump(obj, f, indent=indent)
        f.write("\n")


def get_counts_for_features(
    matrix: CountMatrix, barcode: bytes, features: list[str]
) -> Iterable[int]:
    """Return the counts for the given features in the given barcode."""
    barcode_index = matrix.bc_to_int(barcode)
    return (matrix.m[matrix.feature_id_to_int(x.encode()), barcode_index] for x in features)


def join(args, outs, _chunk_defs, chunk_outs):
    matrix = CountMatrix.load_h5_file(args.filtered_feature_counts_matrix).select_features_by_type(
        rna_library.CRISPR_LIBRARY_TYPE
    )
    num_cells = matrix.bcs_dim
    if num_cells == 0:
        for out in outs.items():
            setattr(outs, out, None)
        return

    # cells_per_protospacer.json
    cells_per_protospacer: dict[str, list[str]] = reduce(
        operator.ior,
        (json.load(open(chunk.chunk_cells_per_protospacer)) for chunk in chunk_outs),
    )
    dump_json(cells_per_protospacer, outs.cells_per_protospacer, indent=2)

    # protospacer_calls_per_cell.csv
    protospacer_calls_per_cell: dict[bytes, list[str]] = defaultdict(list)
    for protospacer, cells in cells_per_protospacer.items():
        for cell in cells:
            protospacer_calls_per_cell[cell.encode()].append(protospacer)
    del cells_per_protospacer
    with open(outs.protospacer_calls_per_cell, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("cell_barcode", "num_features", "feature_call", "num_umis"))
        for barcode, protospacers in protospacer_calls_per_cell.items():
            counts = "|".join(map(str, get_counts_for_features(matrix, barcode, protospacers)))
            writer.writerow((barcode.decode(), len(protospacers), "|".join(protospacers), counts))

    # protospacer_call_metrics_json.json
    num_cells_with_single_protospacer: float = sum(
        len(x) == 1 for x in protospacer_calls_per_cell.values()
    )
    num_cells_with_multiple_protospacer = sum(
        len(x) >= 2 for x in protospacer_calls_per_cell.values()
    )
    num_cells_with_protospacer = (
        num_cells_with_single_protospacer + num_cells_with_multiple_protospacer
    )
    num_cells_without_protospacer = num_cells - num_cells_with_protospacer
    num_cells_with_no_guide_molecules = sum(x == 0 for x in matrix.get_counts_per_bc())
    frac_cells_with_single_protospacer = num_cells_with_single_protospacer / num_cells
    frac_cells_with_multiple_protospacer = num_cells_with_multiple_protospacer / num_cells
    frac_cells_with_protospacer = num_cells_with_protospacer / num_cells
    frac_cells_without_protospacer = num_cells_without_protospacer / num_cells
    frac_cells_with_no_guide_molecules = num_cells_with_no_guide_molecules / num_cells
    dump_json(
        {
            "CRISPR_frac_cells_with_single_protospacer": frac_cells_with_single_protospacer,
            "CRISPR_frac_cells_with_multiple_protospacer": frac_cells_with_multiple_protospacer,
            "CRISPR_frac_cells_with_protospacer": frac_cells_with_protospacer,
        },
        outs.protospacer_call_metrics_json,
    )

    # protospacer_calls_summary.csv
    protospacer_calls_summary: dict[str, list[int]] = defaultdict(list)
    for barcode, protospacers in protospacer_calls_per_cell.items():
        protospacer_calls_summary["|".join(protospacers)].append(
            int(sum(get_counts_for_features(matrix, barcode, protospacers)))
        )
    del matrix, protospacer_calls_per_cell

    with open(outs.protospacer_calls_summary, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("feature_call", "num_cells", "pct_cells", "median_umis", "stddev_umis"))
        for name, num, frac in (
            (
                "No guide molecules",
                num_cells_with_no_guide_molecules,
                frac_cells_with_no_guide_molecules,
            ),
            ("No confident call", num_cells_without_protospacer, frac_cells_without_protospacer),
            (
                "1 protospacer assigned",
                num_cells_with_single_protospacer,
                frac_cells_with_single_protospacer,
            ),
            (
                "More than 1 protospacer assigned",
                num_cells_with_multiple_protospacer,
                frac_cells_with_multiple_protospacer,
            ),
        ):
            writer.writerow((name, num, frac, "None", "None"))

        for feature_call, umis in protospacer_calls_summary.items():
            writer.writerow(
                (
                    feature_call,
                    len(umis),
                    len(umis) / num_cells,
                    median(umis),
                    pstdev(umis) if len(umis) >= 2 else 0,
                )
            )
        del protospacer_calls_summary

    # protospacer_umi_thresholds_json.json
    protospacer_umi_thresholds = reduce(
        operator.ior,
        (json.load(open(chunk.chunk_protospacer_umi_thresholds)) for chunk in chunk_outs),
    )
    dump_json(protospacer_umi_thresholds, outs.protospacer_umi_thresholds_json)

    # protospacer_umi_thresholds_csv.csv
    with open(outs.protospacer_umi_thresholds_csv, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("Protospacer", "UMI threshold"))
        for protospacer, umi_threshold in protospacer_umi_thresholds.items():
            writer.writerow((protospacer, umi_threshold))
