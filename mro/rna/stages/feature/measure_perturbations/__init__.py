#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#

from __future__ import annotations

import csv
import itertools
import json
import os
from typing import TYPE_CHECKING

import martian
import numpy as np
import scipy.sparse as sp_sparse

from cellranger import cr_io
from cellranger.analysis.diffexp import save_differential_expression_csv
from cellranger.feature.crispr.measure_perturbations import (
    NUM_THREADS,
    get_perturbation_efficiency,
    read_and_validate_feature_ref,
    read_bootstrap_count_total,
    read_perturbation_stats,
    save_perturbation_efficiency_summary,
)
from cellranger.matrix import CountMatrix
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage MEASURE_PERTURBATIONS(
    in  csv   protospacer_calls_per_cell,
    in  h5    control_stats,
    in  h5    control_bootstrap,
    in  h5[]  chunk_stats,
    in  h5[]  chunk_bootstrap,
    in  json  chunk_perturbations,
    in  bool  by_feature,
    out csv   perturbation_efficiencies,
    out path  perturbation_effects_path,
    src py    "stages/feature/measure_perturbations",
) split (
    in  h5       chunk_stats_slice,
    in  h5       chunk_bootstrap_slice,
    in  string[] which_perturbations,
    out csv      chunk_perturbation_efficiencies,
    out csv      chunk_perturbation_effects,
) using (
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """MEASURE_PERTURBATIONS stage inputs.

    SPLIT_PERTURBATION_MATRIX has reduced the matrix upstream: `control_stats` holds the shared
    non-targeting control's sSeq sufficient statistics (+ the full feature reference), `chunk_stats[k]`
    holds chunk k's perturbations' statistics, and `control_bootstrap` / `chunk_bootstrap[k]` hold the
    small target-gene matrices the empirical-CI bootstrap resamples. `chunk_perturbations` (from
    COMPUTE_PRESPLIT_ASSIGNMENT) lists the perturbation names for each chunk, in slice order -- the
    single source of chunk membership, so it is never recomputed here.
    """

    protospacer_calls_per_cell: None | mro_filetypes.csv
    control_stats: mro_filetypes.h5
    control_bootstrap: mro_filetypes.h5
    chunk_stats: list[mro_filetypes.h5]
    chunk_bootstrap: list[mro_filetypes.h5]
    chunk_perturbations: mro_filetypes.json
    by_feature: bool


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """MEASURE_PERTURBATIONS stage outputs."""

    perturbation_efficiencies: mro_filetypes.csv
    perturbation_effects_path: mro_filetypes.path


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """MEASURE_PERTURBATIONS chunk inputs."""

    chunk_stats_slice: mro_filetypes.h5
    chunk_bootstrap_slice: mro_filetypes.h5
    which_perturbations: list[str]


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """MEASURE_PERTURBATIONS chunk outputs."""

    chunk_perturbation_efficiencies: None | mro_filetypes.csv
    chunk_perturbation_effects: None | mro_filetypes.csv


def split(args: StageInputs):
    with open(args.chunk_perturbations) as f:
        chunk_perturbations = json.load(f)

    if not chunk_perturbations:
        return {"chunks": []}

    assert (
        len(chunk_perturbations) == len(args.chunk_stats) == len(args.chunk_bootstrap)
    ), "chunk_perturbations, chunk_stats and chunk_bootstrap must align position-for-position"

    # Each chunk loads only the shared control bootstrap + its bootstrap slice (target-gene rows
    # only) and the compact statistics, so size memory from the bootstrap dims (the matrix is gone).
    _, control_barcodes, control_nnz = CountMatrix.load_dims_from_h5(args.control_bootstrap)

    chunks = []
    for chunk_stats_slice, chunk_bootstrap_slice, which_perturbations in zip(
        args.chunk_stats, args.chunk_bootstrap, chunk_perturbations
    ):
        _, slice_barcodes, slice_nnz = CountMatrix.load_dims_from_h5(chunk_bootstrap_slice)
        # The scale=1.5 matrix term already exceeds observed peak (~1.8x peak on the by-target chunks),
        # so it carries the data scaling; a +0.5 cushion + 1 GiB floor covers
        # the DE result arrays and clustering on [control ++ slice] for the small base-dominated chunks.
        mem_gib = max(
            1.0,
            0.5
            + CountMatrix.get_mem_gb_from_matrix_dim(
                control_barcodes + slice_barcodes, control_nnz + slice_nnz, scale=1.5
            ),
        )
        chunks.append(
            {
                "chunk_stats_slice": chunk_stats_slice,
                "chunk_bootstrap_slice": chunk_bootstrap_slice,
                "which_perturbations": which_perturbations,
                "__mem_gb": mem_gib,
                "__threads": NUM_THREADS,
            }
        )
    return {"chunks": chunks}


def main(args: ChunkInputs, outs: ChunkOutputs):
    full_feature_ref = CountMatrix.load_feature_ref_from_h5_file(args.control_stats)
    target_info = read_and_validate_feature_ref(full_feature_ref)
    assert target_info is not None
    # Full GEX feature reference: gene_id -> result-array row, and the gene axis of the output CSVs.
    gex_feature_ref = full_feature_ref.select_features_by_type(GENE_EXPRESSION_LIBRARY_TYPE)

    control_stats = read_perturbation_stats(args.control_stats)
    chunk_stats = read_perturbation_stats(args.chunk_stats_slice)

    # The bootstrap matrix = target-gene rows x [control cells ++ this chunk's perturbed cells], in
    # storage order -- the same cell order the pre-split path got from hstack([control, slice]), so
    # the bootstrap draws stay bit-identical. combined_col_sums is each cell's all-GEX count_total in
    # that same order (the size-factor numerator), read alongside each bootstrap matrix.
    control_boot = CountMatrix.load_h5_file(args.control_bootstrap)
    chunk_boot = CountMatrix.load_h5_file(args.chunk_bootstrap_slice)
    combined_matrix = sp_sparse.hstack([control_boot.m, chunk_boot.m], format="csc")
    combined_bcs = list(control_boot.bcs) + list(chunk_boot.bcs)
    bootstrap_matrix = CountMatrix(control_boot.feature_ref, combined_bcs, combined_matrix)
    combined_col_sums = np.concatenate(
        [
            read_bootstrap_count_total(args.control_bootstrap),
            read_bootstrap_count_total(args.chunk_bootstrap_slice),
        ]
    ).astype(np.float64)

    perturbation_result = get_perturbation_efficiency(
        set(args.which_perturbations),
        target_info,
        args.protospacer_calls_per_cell,
        combined_bcs,
        gex_feature_ref,
        control_stats,
        chunk_stats,
        bootstrap_matrix,
        combined_col_sums,
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
    # save_differential_expression_csv only needs the gene axis (feature_ref.feature_defs); the DE
    # result rows are the full GEX genes, so pass an empty matrix carrying the GEX feature reference.
    save_differential_expression_csv(
        None,
        results_all_perturbations,
        CountMatrix.empty(gex_feature_ref, []),
        cluster_names=perturbation_names,
        base_dir=os.path.dirname(outs.chunk_perturbation_effects),
        file_name=os.path.basename(outs.chunk_perturbation_effects).removesuffix(".csv"),
    )


def join(
    _args: StageInputs,
    outs: StageOutputs,
    _chunk_defs: list[ChunkInputs],
    chunk_outs: list[ChunkOutputs],
):
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
