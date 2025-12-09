#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#

import csv
import itertools
import math
import os
from collections import Counter

import martian

from cellranger import cr_io
from cellranger.analysis.diffexp import save_differential_expression_csv
from cellranger.constants import FILTER_LIST
from cellranger.feature.crispr.measure_perturbations import (
    MIN_NUMBER_CELLS_PER_PERTURBATION,
    NON_TARGETING,
    NUM_THREADS,
    TargetInfo,
    get_perturbation_efficiency,
    read_and_validate_feature_ref,
    save_perturbation_efficiency_summary,
)
from cellranger.matrix import CountMatrix
from cellranger.rna.library import CRISPR_LIBRARY_TYPE, GENE_EXPRESSION_LIBRARY_TYPE

__MRO__ = """
stage MEASURE_PERTURBATIONS(
    in  csv  protospacer_calls_per_cell,
    in  h5   filtered_feature_counts_matrix,
    in  bool by_feature,
    out csv  perturbation_efficiencies,
    out path perturbation_effects_path,
    src py   "stages/feature/measure_perturbations",
) split (
    in  int  chunk_start,
    in  int  chunk_end,
    out csv  chunk_perturbation_efficiencies,
    out csv  chunk_perturbation_effects,
) using (
    volatile = strict,
)
"""


def get_target(perturbation: str, target_info: dict[str, TargetInfo]) -> str:
    """Return the target gene names of the specified pipe-separated perturbation.

    Deduplicate by gene ID and return the gene names.
    Match the behaviour of `_get_ps_clusters_by_target`.
    """
    gene_id_names = sorted(
        set((target_info[x].gene_id, target_info[x].gene_name) for x in perturbation.split("|"))
    )
    if all(gene_id == NON_TARGETING for gene_id, _ in gene_id_names):
        return NON_TARGETING

    target = "|".join(
        gene_name for gene_id, gene_name in gene_id_names if gene_id not in FILTER_LIST
    )
    return target if target else "Ignore"


def get_perturbations(
    protospacer_calls_per_cell: str,
    target_info: dict[str, TargetInfo],
    by_feature: bool,
) -> list[str]:
    """Return a sorted list of perturbations from the protospacer calls per cell file.

    Translate protospacer IDs to target gene names when by_feature is false.
    """
    with open(protospacer_calls_per_cell) as f:
        feature_calls = (x["feature_call"] for x in csv.DictReader(f))
        perturbation_counts = (
            Counter(feature_calls)
            if by_feature
            else Counter(get_target(x, target_info) for x in feature_calls)
        )
        return sorted(
            perturbation
            for perturbation, cells in perturbation_counts.items()
            if cells >= MIN_NUMBER_CELLS_PER_PERTURBATION
        )


def split(args):
    feature_ref = CountMatrix.load_feature_ref_from_h5_file(args.filtered_feature_counts_matrix)
    target_info = read_and_validate_feature_ref(feature_ref)
    if args.protospacer_calls_per_cell is None or target_info is None:
        return {"chunks": []}

    num_perturbations = len(
        get_perturbations(args.protospacer_calls_per_cell, target_info, args.by_feature)
    )
    if num_perturbations == 0:
        return {"chunks": []}

    num_gex_features = feature_ref.get_count_of_feature_type(GENE_EXPRESSION_LIBRARY_TYPE)
    num_crispr_features = feature_ref.get_count_of_feature_type(CRISPR_LIBRARY_TYPE)
    num_features, num_barcodes, nnz = CountMatrix.load_dims_from_h5(
        args.filtered_feature_counts_matrix
    )
    mem_gib = 1 + CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=1.5)

    print(
        f"{num_perturbations=},{num_gex_features=},{num_crispr_features=},"
        f"{num_features=},{num_barcodes=},{nnz=},{mem_gib=}"
    )
    num_chunks = min(100, num_perturbations, math.ceil(num_perturbations / 50))
    num_perturbations_per_chunk = math.ceil(num_perturbations / num_chunks)

    return {
        "chunks": [
            {
                "chunk_start": start,
                "chunk_end": min(num_perturbations, start + num_perturbations_per_chunk),
                "__mem_gb": mem_gib,
                "__threads": NUM_THREADS,
            }
            for start in range(0, num_perturbations, num_perturbations_per_chunk)
        ]
    }


def main(args, outs):
    target_info = read_and_validate_feature_ref(
        CountMatrix.load_feature_ref_from_h5_file(args.filtered_feature_counts_matrix)
    )
    assert target_info is not None

    gex_count_matrix = CountMatrix.load_data_for_library_type_from_h5(
        args.filtered_feature_counts_matrix, GENE_EXPRESSION_LIBRARY_TYPE
    )

    which_perturbations = set(
        get_perturbations(args.protospacer_calls_per_cell, target_info, args.by_feature)[
            args.chunk_start : args.chunk_end
        ]
    )

    perturbation_result = get_perturbation_efficiency(
        which_perturbations,
        target_info,
        args.protospacer_calls_per_cell,
        gex_count_matrix,
        by_feature=args.by_feature,
    )
    if perturbation_result is None:
        outs.chunk_perturbation_efficiencies = None
        outs.chunk_perturbation_effects = None
        return

    (
        perturbation_names,
        results_all_perturbations,
        fold_change_per_perturbation,
    ) = perturbation_result

    save_perturbation_efficiency_summary(
        outs.chunk_perturbation_efficiencies,
        fold_change_per_perturbation,
        args.by_feature,
    )
    save_differential_expression_csv(
        None,
        results_all_perturbations,
        gex_count_matrix,
        cluster_names=perturbation_names,
        base_dir=os.path.dirname(outs.chunk_perturbation_effects),
        file_name=os.path.basename(outs.chunk_perturbation_effects).removesuffix(".csv"),
    )


def join(_args, outs, _chunk_defs, chunk_outs):
    if chunk_outs == [] or all(
        chunk.chunk_perturbation_efficiencies is None for chunk in chunk_outs
    ):
        martian.clear(outs)
        return

    cr_io.concatenate_headered_files(
        outs.perturbation_efficiencies,
        [
            chunk.chunk_perturbation_efficiencies
            for chunk in chunk_outs
            if chunk.chunk_perturbation_efficiencies is not None
        ],
    )

    # Sort by column `Log2 Fold Change`
    with open(outs.perturbation_efficiencies) as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header[2] == "Log2 Fold Change"
        perturbation_efficiencies = [header] + sorted(reader, key=lambda x: float(x[2]))
    with open(outs.perturbation_efficiencies, "w") as f:
        csv.writer(f, lineterminator="\n").writerows(perturbation_efficiencies)

    os.makedirs(outs.perturbation_effects_path, exist_ok=True)
    with open(
        os.path.join(outs.perturbation_effects_path, "transcriptome_analysis.csv"), "w"
    ) as outfile:
        writer = csv.writer(outfile, lineterminator="\n")
        for first_reader, *other_readers in zip(
            *(
                csv.reader(open(chunk.chunk_perturbation_effects))
                for chunk in chunk_outs
                if chunk.chunk_perturbation_effects is not None
            )
        ):
            # The first two columns `Feature ID` and `Feature Name` are identical in all chunks.
            for line in other_readers:
                assert line[0:2] == first_reader[0:2]
            writer.writerow(itertools.chain(first_reader, *(line[2:] for line in other_readers)))
