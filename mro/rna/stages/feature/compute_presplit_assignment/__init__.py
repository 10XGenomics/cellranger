#!/usr/bin/env python
#
# Copyright (c) 2026 10X Genomics, Inc. All rights reserved
#

"""Compute the per-barcode reduction plan for the CRISPR perturbation sufficient-statistics split.

`SPLIT_PERTURBATION_MATRIX` (Rust) streams the feature-barcode matrix once and emits, per
MEASURE chunk, the sSeq differential-expression *sufficient statistics* (per gene x cell-group
aggregates) plus a small target-gene "bootstrap" matrix — never the full per-chunk matrix the
pre-split MEASURE used to hstack. The *biology* that decides which barcode belongs to which
perturbation stays here in Python, reusing the exact perturbation/target clustering
(`_get_ps_clusters`, `get_perturbations`) so the control set, the per-cell perturbation identity,
and the per-chunk perturbation membership remain a single source of truth shared with
`MEASURE_PERTURBATIONS`.

We emit:
  * a single `reduce_plan` JSON consumed by the Rust splitter, holding the global target-gene
    GEX rows (the bootstrap matrix's feature rows) and, per mode (by-feature and by-target):
      - `group_assignment`: one int per barcode in matrix storage order — -1 (drop), 0
        (control / non-targeting) or `p` in `1..=num_perturbations` (the p-th perturbation, in
        `get_perturbations` order). This is the per-cell *perturbation* identity the streaming
        reducer accumulates per (gene x group);
      - `group_to_chunk`: for each perturbation group `p` (1-based), the MEASURE chunk (1-based)
        that owns it, so the splitter buckets each perturbation's statistics into a chunk;
      - `num_perturbations`, `num_chunks`.
  * per-mode `*_chunks` JSONs `[[perturbation_name, ...], ...]` (length = num_chunks) giving the
    perturbation names handled by each chunk, in slice order; consumed by
    `MEASURE_PERTURBATIONS.split` to set each chunk's `which_perturbations` (unchanged contract).
"""

from __future__ import annotations

import json
import math
from typing import TYPE_CHECKING

import martian

from cellranger.feature.crispr.measure_perturbations import (
    NON_TARGETING,
    TargetInfo,
    _get_ps_clusters,
    get_perturbations,
    read_and_validate_feature_ref,
)
from cellranger.matrix import CountMatrix
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage COMPUTE_PRESPLIT_ASSIGNMENT(
    in  csv  protospacer_calls_per_cell,
    in  h5   filtered_feature_counts_matrix,
    out json reduce_plan,
    out json by_feature_chunks,
    out json by_target_chunks,
    src py   "stages/feature/compute_presplit_assignment",
) using (
    mem_gb   = 3,
    volatile = strict,
)
"""
# mem_gb note: this stage loads the barcode list + protospacer calls (not the matrix counts), so its
# working set is linear in cell count -- ~0.67 GiB at ~1M cells. 3 GiB stays
# >=2x peak for any realistic CRISPR cell count (<~1.5M); it is a fixed request, not data-sized,
# because adding a split() purely to size memory is not worth it for a sub-GiB-peak stage.


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """COMPUTE_PRESPLIT_ASSIGNMENT stage inputs."""

    protospacer_calls_per_cell: None | mro_filetypes.csv
    filtered_feature_counts_matrix: mro_filetypes.h5


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """COMPUTE_PRESPLIT_ASSIGNMENT stage outputs."""

    reduce_plan: mro_filetypes.json
    by_feature_chunks: mro_filetypes.json
    by_target_chunks: mro_filetypes.json


def _compute_mode(
    target_info: dict[str, TargetInfo],
    calls: str,
    barcodes: list[str],
    by_feature: bool,
) -> tuple[dict, list[list[str]]]:
    """Return (mode plan dict, chunks list[list[str]]) for one mode.

    The plan dict has `group_assignment` (-1 drop / 0 control / 1..P perturbation, per barcode),
    `group_to_chunk` (chunk 1..N for each perturbation group 1..P), `num_perturbations`,
    `num_chunks`. `chunks[k - 1]` is the perturbation-name list handled by chunk k, preserving the
    exact membership of `get_perturbations(...)[start:end]` so the chunk boundaries — and therefore
    each `MEASURE_PERTURBATIONS` chunk's `group_a`/`group_b` — cannot drift from the matrix path.
    """
    target_calls, perturbation_keys = _get_ps_clusters(target_info, calls, barcodes, by_feature)
    perturbations = get_perturbations(calls, target_info, by_feature)

    num_perturbations = len(perturbations)
    if num_perturbations == 0:
        empty = {
            "group_assignment": [-1] * len(barcodes),
            "group_to_chunk": [],
            "perturbation_names": [],
            "num_perturbations": 0,
            "num_chunks": 0,
        }
        return empty, []

    num_chunks = min(100, num_perturbations, math.ceil(num_perturbations / 50))
    per_chunk = math.ceil(num_perturbations / num_chunks)
    chunks = [perturbations[s : s + per_chunk] for s in range(0, num_perturbations, per_chunk)]

    # Global perturbation group ids (1..P) in get_perturbations order — the order the streaming
    # reducer's per-group accumulators are indexed by, and the order chunk slices partition.
    name_to_group = {name: i + 1 for i, name in enumerate(perturbations)}

    # group -> chunk (1-based). Non-Targeting never owns a perturbed slice even if it appears in a
    # chunk's name list (it is filtered when the DE is computed, exactly as in the full-matrix path),
    # so its group has no chunk and no cells.
    group_to_chunk = [0] * num_perturbations
    for k, chunk_perturbations in enumerate(chunks, start=1):
        for name in chunk_perturbations:
            if name != NON_TARGETING:
                group_to_chunk[name_to_group[name] - 1] = k

    # cluster -> per-cell group id: 0 control (Non-Targeting), p perturbation, -1 drop.
    cluster_to_group: dict[int, int] = {}
    for cluster, name in perturbation_keys.items():
        if name == NON_TARGETING:
            cluster_to_group[cluster] = 0
        else:
            cluster_to_group[cluster] = name_to_group.get(name, -1)

    plan = {
        "group_assignment": [cluster_to_group.get(c, -1) for c in target_calls.tolist()],
        "group_to_chunk": group_to_chunk,
        "perturbation_names": perturbations,
        "num_perturbations": num_perturbations,
        "num_chunks": num_chunks,
    }
    return plan, chunks


def _target_gene_rows(feature_ref, target_info: dict[str, TargetInfo]) -> list[int]:
    """GEX rows (0-based among GEX features, in storage order) of every targeted gene present.

    These are the feature rows of the small per-chunk "bootstrap" matrix the splitter writes; the
    MEASURE chunk selects a single target gene out of it for the empirical-CI bootstrap via
    `select_features_by_ids`, so every gene any perturbation could target must be a row. The GEX
    row indexing matches the splitter's (both enumerate GENE_EXPRESSION features in storage order)
    and the MEASURE chunk's `gene_id_to_row` (built from the same GEX feature defs).
    """
    gene_id_to_gex_row: dict[str, int] = {}
    gex_row = 0
    for feature_def in feature_ref.feature_defs:
        if feature_def.feature_type == GENE_EXPRESSION_LIBRARY_TYPE:
            gene_id_to_gex_row.setdefault(feature_def.id.decode(), gex_row)
            gex_row += 1
    rows = {
        gene_id_to_gex_row[info.gene_id]
        for info in target_info.values()
        if info.gene_id in gene_id_to_gex_row
    }
    return sorted(rows)


def main(args: StageInputs, outs: StageOutputs):
    feature_ref = CountMatrix.load_feature_ref_from_h5_file(args.filtered_feature_counts_matrix)
    target_info = read_and_validate_feature_ref(feature_ref)
    barcodes = [
        bc.decode() for bc in CountMatrix.load_bcs_from_h5(args.filtered_feature_counts_matrix)
    ]

    if args.protospacer_calls_per_cell is None or target_info is None:
        empty_mode = {
            "group_assignment": [-1] * len(barcodes),
            "group_to_chunk": [],
            "perturbation_names": [],
            "num_perturbations": 0,
            "num_chunks": 0,
        }
        _write_json(
            outs.reduce_plan,
            {"target_gene_rows": [], "by_feature": empty_mode, "by_target": empty_mode},
        )
        _write_json(outs.by_feature_chunks, [])
        _write_json(outs.by_target_chunks, [])
        return

    by_feature_plan, by_feature_chunks = _compute_mode(
        target_info, args.protospacer_calls_per_cell, barcodes, by_feature=True
    )
    by_target_plan, by_target_chunks = _compute_mode(
        target_info, args.protospacer_calls_per_cell, barcodes, by_feature=False
    )

    _write_json(
        outs.reduce_plan,
        {
            "target_gene_rows": _target_gene_rows(feature_ref, target_info),
            "by_feature": by_feature_plan,
            "by_target": by_target_plan,
        },
    )
    _write_json(outs.by_feature_chunks, by_feature_chunks)
    _write_json(outs.by_target_chunks, by_target_chunks)


def _write_json(path: str, obj: dict | list) -> None:
    with open(path, "w") as f:
        json.dump(obj, f)
