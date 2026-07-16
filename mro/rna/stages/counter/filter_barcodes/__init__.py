#!/usr/bin/env python3
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Identify partitions that contain cells."""

from __future__ import annotations

import dataclasses
import json
import math
import os
from collections import defaultdict
from typing import TYPE_CHECKING, Any, NamedTuple

import martian
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

from cellranger.cell_barcodes import InvalidCellBarcode
from cellranger.cell_calling_helpers import (
    CellCalling,
    FilterMethod,
    HighOccupancyGemSummary,
    MetricGroups,
    apply_global_minimum_umis_threshold,
    apply_mitochondrial_threshold,
    call_additional_cells,
    call_initial_cells,
    combine_initial_metrics,
    get_emptydrops_minimum_umis,
    get_filter_method_from_string,
    get_filter_method_name,
    get_force_cells,
    get_global_minimum_umis,
    get_max_mito_percent,
    get_recovered_cells,
    remove_antibody_antigen_aggregates,
    remove_bcs_from_high_occupancy_gems,
    remove_cells_with_zero_targeted_counts,
)
from cellranger.chemistry import CHEMISTRY_DESCRIPTION_FIELD, get_primary_chemistry_def
from cellranger.constants import DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP
from cellranger.csv_utils import (
    combine_csv,
    write_filtered_barcodes,
    write_isotype_normalization_csv,
)
from cellranger.fast_utils import MultiGraph, load_filtered_bcs_groups, save_filtered_bcs_groups
from cellranger.feature.antibody.analysis import (
    FRACTION_CORRECTED_READS,
    FRACTION_SAMPLE_READS,
    FRACTION_TOTAL_READS,
    augment_correction_table_with_corrected_reads_fraction,
    augment_correction_table_with_read_fraction,
)
from cellranger.matrix import CountMatrix, make_matrix_attrs_count
from cellranger.metrics import BarcodeFilterResults
from cellranger.read_level_multiplexing import get_overhang_bc_defn, get_sample_tag_barcodes
from cellranger.rna.library import (
    ANTIBODY_LIBRARY_TYPE,
    ANTIGEN_LIBRARY_TYPE,
    GENE_EXPRESSION_LIBRARY_TYPE,
    MULTIPLEXING_LIBRARY_TYPE,
)
from cellranger.rna.matrix import save_mex
from cellranger.targeted.rtl_multiplexing import get_probe_bc_defn
from tenkit.safe_json import dump_numpy

if TYPE_CHECKING:
    from collections.abc import Mapping

    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

PROBE_BC_SAMPLE_ID = "id"
PROBE_BC_SEQS = "sequence"
PROBE_BC_OFFSET = "offset"
PROBE_BC_LEN = "length"

# FBc control tag name
FBC_CONTROL_TAG_NAME = "isotype_control"

# FBc Scale factors
FBC_SCALE_FACTOR = 10_000

# Number of FBc s needed in a count matrix to normalise
MINIMUM_FBCS_TO_NORMALIZE = 10

np.random.seed(0)


# TODO: We can get rid of all the gem-group logic here,
# because we run this only on data from a single gem-group

__MRO__ = """
struct ProbeBCDef(
    string   id,
    string[] sequence,
    int      offset,
    int      length,
)

struct CellCallingParam(
    float      per_gem_well,
    map<float> per_sample,
)


struct CellCalling(
    CellCallingParam recovered_cells,
    CellCallingParam force_cells,
    CellCallingParam emptydrops_minimum_umis,
    CellCallingParam global_minimum_umis,
    json             cell_barcodes,
    string           override_mode,
    string[]         override_library_types,
    bool             disable_ab_aggregate_detection,
    bool             disable_high_occupancy_gem_detection,
)

stage FILTER_BARCODES(
    in  map<ChemistryDef> chemistry_defs,
    in  string            sample_id,
    in  h5                matrices_h5,
    in  csv               barcode_correction_csv,
    in  bool              is_antibody_only,
    in  int[]             gem_groups,
    in  CellCalling       cell_calling_config,
    in  json              multi_graph,
    in  csv               per_barcode_metrics,
    in  bool              is_spatial,
    out json              summary,
    out csv               filtered_barcodes,
    out csv               aggregate_barcodes,
    out h5                filtered_matrices_h5,
    out path              filtered_matrices_mex,
    out csv               nonambient_calls,
    out csv               mitochondrial_summary,
    out csv               isotype_normalization_factors,
    src py                "stages/counter/filter_barcodes",
) split (
    in  ProbeBCDef        probe_bc_def,
    out json              filtered_metrics_groups,
    out bincode           filtered_bcs_groups,
    out csv               co_mitochondrial_summary,
) using (
    mem_gb   = 8,
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """FILTER_BARCODES stage inputs."""

    chemistry_defs: dict[str, mro_structs.ChemistryDef]
    sample_id: str
    matrices_h5: mro_filetypes.h5
    barcode_correction_csv: None | mro_filetypes.csv
    is_antibody_only: bool
    gem_groups: list[int]
    cell_calling_config: mro_structs.CellCalling
    multi_graph: None | mro_filetypes.json
    per_barcode_metrics: mro_filetypes.csv
    is_spatial: bool


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """FILTER_BARCODES stage outputs."""

    summary: None | mro_filetypes.json
    filtered_barcodes: mro_filetypes.csv
    aggregate_barcodes: None | mro_filetypes.csv
    filtered_matrices_h5: mro_filetypes.h5
    filtered_matrices_mex: mro_filetypes.path
    nonambient_calls: None | mro_filetypes.csv
    mitochondrial_summary: None | mro_filetypes.csv
    isotype_normalization_factors: None | mro_filetypes.csv


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """FILTER_BARCODES chunk inputs."""

    probe_bc_def: None | mro_structs.ProbeBCDef


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """FILTER_BARCODES chunk outputs."""

    filtered_metrics_groups: None | mro_filetypes.json
    filtered_bcs_groups: None | mro_filetypes.bincode
    co_mitochondrial_summary: None | mro_filetypes.csv


class CellCallingArgs(NamedTuple):
    """Arguments for cell calling initialization."""

    chemistry_description: str
    config: CellCalling
    feature_types: list[str]
    force_cells: int | None
    genomes: list[str]
    gg_recovered_cells: int | None
    has_cmo_data: bool
    matrix: CountMatrix
    method: FilterMethod
    num_probe_barcodes: int | None
    probe_barcode_sample_id: str | None
    summary: dict[Any, Any]
    target_features: list[int] | None
    unique_gem_groups: list[int]


# Heuristic constants - calibrated against 29-case baselines + orthogonal
# OOM set (see deck/cellranger/fb_smoke/heuristic_eval.py and
# FILTER_BARCODES_IDEAS.md). Aims for 10-25% headroom over measured peak;
# never under-predicts on the calibration set.
_HEUR_B_PER_NNZ = 12
_HEUR_B_PER_BC = 208
_HEUR_B_PER_FEAT = 1500
_HEUR_GIB = 1024**3
_HEUR_MIN_MEM_GIB = 4
_HEUR_MAIN_SAFETY_FACTOR = 1.29
_HEUR_JOIN_SAFETY_FACTOR = 1.00


def _heuristic_mem_gib(
    num_features: int,
    num_barcodes: int,
    nnz: int,
    *,
    has_ab: bool,
    is_spatial: bool,
    n_samples: int,
) -> tuple[int, int]:
    """Return (main_chunk_gib, join_gib) for FB using cheap split-time inputs."""
    load = _HEUR_B_PER_NNZ * nnz + _HEUR_B_PER_BC * num_barcodes + _HEUR_B_PER_FEAT * num_features
    main_scaled = (1.85 if n_samples > 1 else 4.10) * load
    # Visium HD carries extra fixed overhead (spot grid + barcode<->coordinate
    # maps + spatial harness state) that the load*scale term doesn't capture.
    # Calibrated against 1893385_1984732 (30M-spot HD): bare main ratio was
    # 3.75x load with 4.70 GiB residual; this 4 GiB bump restores 1.20x+ on
    # the Python-overhead piece.
    main_fixed_gib = 5.0 if is_spatial else 1.0
    main_bytes = main_scaled + (_HEUR_GIB if has_ab else 0) + main_fixed_gib * _HEUR_GIB
    join_scaled = 2.25 * load
    if is_spatial:
        join_scaled += load
    if has_ab:
        join_scaled += 0.5 * _HEUR_GIB
    join_bytes = join_scaled + 1.5 * _HEUR_GIB
    return (
        max(_HEUR_MIN_MEM_GIB, math.ceil(main_bytes * _HEUR_MAIN_SAFETY_FACTOR / _HEUR_GIB)),
        max(_HEUR_MIN_MEM_GIB, math.ceil(join_bytes * _HEUR_JOIN_SAFETY_FACTOR / _HEUR_GIB)),
    )


def split(args: StageInputs):
    chemistry_def = get_primary_chemistry_def(args.chemistry_defs)
    chemistry_barcode = chemistry_def["barcode"]

    num_features, num_barcodes, nnz = CountMatrix.load_dims_from_h5(args.matrices_h5)
    lib_types = (
        MultiGraph.from_path(args.multi_graph).library_types()
        if args.multi_graph
        else CountMatrix.load_library_types_from_h5_file(args.matrices_h5)
    )
    has_ab = ANTIBODY_LIBRARY_TYPE in lib_types

    # Plan chunks first so n_samples is known for the heuristic.
    chunks_planned: list[dict] = []
    if args.multi_graph is None:
        chunks_planned.append({"probe_bc_def": None})
    else:
        multi_graph = MultiGraph.from_path(args.multi_graph)
        if multi_graph.is_rtl_multiplexed():
            offset, length = get_probe_bc_defn(chemistry_barcode)
        elif multi_graph.is_oh_multiplexed():
            offset, length = get_overhang_bc_defn(chemistry_barcode)
        else:
            offset, length = None, None

        if length is not None:
            probe_bcs = get_sample_tag_barcodes(multi_graph, chemistry_barcode)
            # Sort by sample_id so the chunk order (and therefore the row order of
            # the join-combined per-chunk CSVs: nonambient_calls, aggregate_barcodes,
            # mitochondrial_summary) is deterministic. `get_sample_tag_barcodes`
            # returns a dict whose key order derives from a Rust HashMap
            # (MultiGraph.sample_tag_ids), which iterates in a per-process-random
            # order; without this sort the combined CSVs differ run-to-run only in
            # row order. Cell calling is independent per chunk, so reordering chunks
            # cannot change any computed value (filtered barcodes are merged through
            # a sorted BTreeSet; metrics are content-keyed; summary is summed).
            for sample_id, probe_bc_seqs in sorted(probe_bcs.items()):
                chunks_planned.append(
                    {
                        "probe_bc_def": {
                            "id": sample_id,
                            "sequence": probe_bc_seqs,
                            "offset": offset,
                            "length": length,
                        },
                    }
                )
        else:
            chunks_planned.append({"probe_bc_def": None})

    n_samples = len(chunks_planned)
    chunk_mem_gib, join_mem_gib = _heuristic_mem_gib(
        num_features,
        num_barcodes,
        nnz,
        has_ab=has_ab,
        is_spatial=args.is_spatial,
        n_samples=n_samples,
    )
    print(
        f"FILTER_BARCODES.split: num_features={num_features} num_barcodes={num_barcodes} "
        f"nnz={nnz} n_samples={n_samples} has_ab={has_ab} "
        f"is_spatial={args.is_spatial} chunk_mem_gib={chunk_mem_gib} join_mem_gib={join_mem_gib}"
    )

    chunks = [{**c, "__mem_gb": chunk_mem_gib} for c in chunks_planned]
    return {"chunks": chunks, "join": {"__mem_gb": join_mem_gib}}


def main(args: ChunkInputs, outs: ChunkOutputs):
    filter_barcodes(args, outs)


def join(
    args: StageInputs,
    outs: StageOutputs,
    _chunk_defs: list[ChunkInputs],
    chunk_outs: list[ChunkOutputs],
):
    chemistry_def = get_primary_chemistry_def(args.chemistry_defs)

    filtered_metrics_groups = parse_chunked_metrics(chunk_outs)
    filtered_bcs, genome_filtered_bcs = parse_chunked_filtered_bcs(chunk_outs)
    summary = parse_chunked_summary(chunk_outs)

    # high occupancy GEM filtering to mitigate any GEMs that appear overloaded
    # in initial cell calls. Only applies if input has probe bcs and not disabled.
    config = CellCalling(args.cell_calling_config)
    chemistry_barcode = chemistry_def["barcode"]
    probe_bc_offset, _ = get_probe_bc_defn(chemistry_barcode)
    if probe_bc_offset and not config.disable_high_occupancy_gem_detection:
        (
            filtered_bcs,
            genome_filtered_bcs,
            high_occupancy_gem_metrics,
        ) = remove_bcs_from_high_occupancy_gems(
            filtered_bcs,
            genome_filtered_bcs,
            args.per_barcode_metrics,
            probe_bc_offset,
        )
        high_occupancy_gem_metrics = dataclasses.asdict(high_occupancy_gem_metrics)
    else:
        high_occupancy_gem_metrics = dataclasses.asdict(HighOccupancyGemSummary())

    # Select cell-associated barcodes
    raw_matrix = CountMatrix.load_h5_file(args.matrices_h5)

    # If non-spatial targeted ensure all filtered barcodes have nonzero targeted UMI counts
    is_targeted = raw_matrix.feature_ref.has_target_features()
    if (not args.is_spatial) and is_targeted:
        filtered_bcs, genome_filtered_bcs = remove_cells_with_zero_targeted_counts(
            raw_matrix, filtered_bcs, genome_filtered_bcs, summary
        )

    # Build the filtered matrix from raw_matrix.
    # On the targeted path (Opt-3) we fuse `select_barcodes_by_seq` and
    # `remove_genes_not_on_list` into a single CSC walk: the intermediate
    # bcs-only matrix is never materialised, saving ~36 GiB of peak co-residency
    # with `remove_genes_not_on_list`'s fancy-index transients at 384-plex Flex
    # scale. Non-targeted path keeps the simple `select_barcodes_by_seq` because
    # no gene-restriction follows.
    if is_targeted:
        target_features = raw_matrix.feature_ref.get_target_feature_indices()
        assert target_features is not None
        bc_indices = raw_matrix.bcs_to_ints(filtered_bcs, return_sorted=False)
        filtered_matrix = raw_matrix.select_barcodes_and_genes_targeted(bc_indices, target_features)
    else:
        filtered_matrix = raw_matrix.select_barcodes_by_seq(filtered_bcs)

    # raw_matrix.m is only needed past this point on the spatial AB-normalization path.
    # On non-spatial paths, free it now so the post-load process phase does not
    # co-reside with raw_matrix's ~38 GiB at 384-plex Flex scale.
    if not args.is_spatial:
        del raw_matrix

    # Write the filtered barcodes file
    write_filtered_barcodes(outs.filtered_barcodes, genome_filtered_bcs)

    # Normalise matrix on ABs in spatial
    libraries_in_h5 = filtered_matrix.get_library_types()
    # Check if there is a spatial FBc run in the matrix and if so normalise by the
    # control FBc counts
    normalisation_metrics: dict[str, int | str] = {"number_of_overflows_in_normalization": 0}
    outs.isotype_normalization_factors = None
    # This accounts for when the --no-libraries flag is passed in but a feature ref is defined.
    # An empty Antibody matrix will be passed in and we don't want to go through the waste of normalizing
    # when it's not needed. Gated on is_spatial so the non-spatial path (where raw_matrix
    # has already been freed) does not scan it.
    antibody_counts = (
        sum(raw_matrix.view().select_features_by_type(ANTIBODY_LIBRARY_TYPE).get_counts_per_bc())
        if args.is_spatial
        else 0
    )
    if ANTIBODY_LIBRARY_TYPE in libraries_in_h5 and args.is_spatial and antibody_counts != 0:
        # Get names of negative controls
        fbc_control_names = [
            x.id
            for x in filtered_matrix.view().feature_ref.feature_defs
            if x.tags.get(FBC_CONTROL_TAG_NAME, "FALSE").upper() == "TRUE"
        ]
        number_of_fbcs_in_matrix = filtered_matrix.get_count_of_feature_type(ANTIBODY_LIBRARY_TYPE)

        if fbc_control_names or number_of_fbcs_in_matrix > MINIMUM_FBCS_TO_NORMALIZE:
            # Normalise the filtered matrix
            normalisation_metrics["number_of_overflows_in_normalization"] = (
                filtered_matrix.normalise_library_type_by_control_features(
                    ANTIBODY_LIBRARY_TYPE, FBC_SCALE_FACTOR, fbc_control_names
                )
            )

            # Compute normalisation factor for all barcodes in raw matrix
            if fbc_control_names:
                normalization_factor = FBC_SCALE_FACTOR / (
                    raw_matrix.view().select_features_by_ids(fbc_control_names).get_counts_per_bc()
                    + 1
                )
                normalisation_metrics["ANTIBODY_isotype_normalized"] = "On"
            else:
                normalization_factor = FBC_SCALE_FACTOR / (
                    raw_matrix.view()
                    .select_features_by_type(ANTIBODY_LIBRARY_TYPE)
                    .get_counts_per_bc()
                    + 1
                )
                normalisation_metrics["ANTIBODY_isotype_normalized"] = "Off"

            outs.isotype_normalization_factors = martian.make_path(
                "isotype_normalization_factors.csv"
            ).decode()
            filtered_bcs_set = set(filtered_bcs)
            in_tissue = [x in filtered_bcs_set for x in raw_matrix.bcs]
            write_isotype_normalization_csv(
                outs.isotype_normalization_factors, raw_matrix.bcs, in_tissue, normalization_factor
            )

    if args.is_spatial:
        del raw_matrix

    matrix_attrs = make_matrix_attrs_count(
        args.sample_id, args.gem_groups, chemistry_def[CHEMISTRY_DESCRIPTION_FIELD]
    )
    filtered_matrix.save_h5_file(
        outs.filtered_matrices_h5,
        extra_attrs=matrix_attrs,
        sw_version=martian.get_pipelines_version(),
    )

    save_mex(filtered_matrix, outs.filtered_matrices_mex, martian.get_pipelines_version())

    parse_filtered_bcs_method(filtered_metrics_groups, summary)

    summary = combine_initial_metrics(filtered_metrics_groups, genome_filtered_bcs, summary)

    # Add keys in from high occupancy GEM filtering in if it was performed
    summary.update(high_occupancy_gem_metrics)

    # Add normalisation metrics
    summary.update(normalisation_metrics)

    # Write metrics json
    assert outs.summary is not None
    with open(ensure_binary(outs.summary), "w") as f:
        dump_numpy(summary, f, indent=4, sort_keys=True)

    chunk_aggregate_barcodes = [
        chunk.aggregate_barcodes
        for chunk in chunk_outs
        if chunk.aggregate_barcodes is not None and os.path.exists(chunk.aggregate_barcodes)
    ]
    if chunk_aggregate_barcodes:
        combine_csv(chunk_aggregate_barcodes, outs.aggregate_barcodes)
    else:
        outs.aggregate_barcodes = None

    chunk_nonambient_calls = [
        chunk.nonambient_calls
        for chunk in chunk_outs
        if chunk.nonambient_calls is not None and os.path.exists(chunk.nonambient_calls)
    ]
    combine_csv(chunk_nonambient_calls, outs.nonambient_calls)

    chunk_mitochondrial_summary = [
        x.co_mitochondrial_summary for x in chunk_outs if x.co_mitochondrial_summary is not None
    ]
    if chunk_mitochondrial_summary:
        combine_csv(chunk_mitochondrial_summary, outs.mitochondrial_summary)
    else:
        outs.mitochondrial_summary = None


def get_cell_calling_args(args: ChunkInputs, outs: ChunkOutputs) -> CellCallingArgs:
    """Initializes cell calling arguments."""
    chemistry_description = get_primary_chemistry_def(args.chemistry_defs)[
        CHEMISTRY_DESCRIPTION_FIELD
    ]
    # Apply the cell calling config in args.cell_calling_config
    config = CellCalling(args.cell_calling_config)
    assert config.force_cells is not None
    force_cells = get_force_cells(
        config.force_cells,
        args.probe_bc_def.get(PROBE_BC_SAMPLE_ID, None) if args.probe_bc_def else None,
    )
    assert config.recovered_cells is not None
    recovered_cells = get_recovered_cells(
        config.recovered_cells,
        args.probe_bc_def.get(PROBE_BC_SAMPLE_ID, None) if args.probe_bc_def else None,
    )

    if args.barcode_correction_csv is None:
        correction_data = None
    else:
        try:
            correction_data = pd.read_csv(
                ensure_str(args.barcode_correction_csv), converters={"barcode": ensure_binary}
            )
            correction_data["barcode"] = correction_data["barcode"].astype("S")
            augment_correction_table_with_corrected_reads_fraction(correction_data)
            augment_correction_table_with_read_fraction(correction_data)
        except pd.errors.EmptyDataError:
            correction_data = None  # will circumvent aggregate detection below

    raw_matrix = CountMatrix.load_h5_file(args.matrices_h5)
    has_cmo_data = raw_matrix.feature_ref.has_feature_type(MULTIPLEXING_LIBRARY_TYPE)

    # If there is a probe barcode def, restrict to data from this specific sample
    if args.probe_bc_def and args.probe_bc_def[PROBE_BC_SEQS]:
        probe_bc_seq = {seq.encode() for seq in args.probe_bc_def[PROBE_BC_SEQS]}
        offset = args.probe_bc_def[PROBE_BC_OFFSET]
        ending_idx = offset + args.probe_bc_def[PROBE_BC_LEN]
        mask = [bc[offset:ending_idx] in probe_bc_seq for bc in raw_matrix.bcs]
        bcs_idx = [i for i, x in enumerate(mask) if x]
        raw_matrix = raw_matrix.select_barcodes(bcs_idx)
        if correction_data is not None:
            correction_data = correction_data[
                correction_data["barcode"].str[offset:ending_idx].isin(probe_bc_seq)
            ]
    # Always do this whether multiplexed or not so the column is always present.
    # For singleplex data, this will ultimately be the same as FRACION_TOTAL_READS.
    if correction_data is not None:
        augment_correction_table_with_read_fraction(correction_data, FRACTION_SAMPLE_READS)
    method: FilterMethod
    if config.cell_barcodes is not None:
        method = FilterMethod.MANUAL
    elif force_cells is not None:
        method = FilterMethod.TOP_N_BARCODES
    elif config.override_mode is not None:
        method = get_filter_method_from_string(config.override_mode)
    elif args.is_antibody_only:
        method = FilterMethod.ORDMAG
    else:
        method = FilterMethod.ORDMAG_NONAMBIENT

    # Get unique gem groups
    unique_gem_groups = sorted(list(set(args.gem_groups)))

    if config.override_library_types is None:
        if args.is_antibody_only:
            feature_types = [ANTIBODY_LIBRARY_TYPE]
        else:
            feature_types = [GENE_EXPRESSION_LIBRARY_TYPE]
    else:
        feature_types = config.override_library_types

    # For targeted GEX, retrieve target gene indices for cell calling
    if raw_matrix.feature_ref.has_target_features():
        target_features = raw_matrix.feature_ref.get_target_feature_indices()
        assert target_features is not None
        assert all(
            x in range(raw_matrix.features_dim) for x in target_features
        ), "Invalid index value found in target_features."
    else:
        target_features = None

    num_probe_barcodes = len(args.probe_bc_def[PROBE_BC_SEQS]) if args.probe_bc_def else None
    if (
        method not in [FilterMethod.MANUAL]
        and correction_data is not None
        and any(
            x in correction_data.library_type.unique()
            for x in [ANTIBODY_LIBRARY_TYPE, ANTIGEN_LIBRARY_TYPE]
        )
    ):
        lib_types = (
            MultiGraph.from_path(args.multi_graph).library_types()
            if args.multi_graph
            else raw_matrix.get_library_types()
        )
        matrix, summary, removed_bcs_df = remove_antibody_antigen_aggregates(
            correction_data,
            raw_matrix,
            lib_types,
            config.disable_ab_aggregate_detection,
            num_probe_barcodes=num_probe_barcodes,
        )
        # report all identified aggregate barcodes, together with their umis,
        # umi corrected reads, fraction of corrected reads, and fraction of total reads
        removed_bcs_df = removed_bcs_df.round(
            {
                FRACTION_CORRECTED_READS: 3,
                FRACTION_TOTAL_READS: 3,
                FRACTION_SAMPLE_READS: 3,
            }
        )
        if len(removed_bcs_df) != 0:
            removed_bcs_df.to_csv(outs.aggregate_barcodes)
    else:
        matrix = raw_matrix
        summary = {}

    if config.cell_barcodes is not None:
        assert method == FilterMethod.MANUAL
    elif force_cells is not None:
        assert method == FilterMethod.TOP_N_BARCODES

    total_diversity_key = (
        args.probe_bc_def[PROBE_BC_SAMPLE_ID] + "_total_diversity"
        if args.probe_bc_def
        else "total_diversity"
    )
    summary[total_diversity_key] = matrix.bcs_dim

    # Get per-gem group cell load
    if recovered_cells is not None:
        gg_recovered_cells = int(float(recovered_cells) / float(len(unique_gem_groups)))
    elif method in (FilterMethod.ORDMAG, FilterMethod.ORDMAG_NONAMBIENT):
        gg_recovered_cells = None
    else:
        gg_recovered_cells = DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    genomes = CountMatrix.get_genomes_from_h5(args.matrices_h5) or ["NONE"]
    probe_barcode_sample_id = args.probe_bc_def[PROBE_BC_SAMPLE_ID] if args.probe_bc_def else None

    return CellCallingArgs(
        chemistry_description=chemistry_description,
        config=config,
        feature_types=feature_types,
        force_cells=force_cells,
        genomes=genomes,
        gg_recovered_cells=gg_recovered_cells,
        has_cmo_data=has_cmo_data,
        matrix=matrix,
        method=method,
        num_probe_barcodes=num_probe_barcodes,
        probe_barcode_sample_id=probe_barcode_sample_id,
        summary=summary,
        target_features=target_features,
        unique_gem_groups=unique_gem_groups,
    )


def initial_cell_calling(cell_calling_args: CellCallingArgs) -> tuple[
    dict[MetricGroups, BarcodeFilterResults],
    dict[tuple[int, str], list[bytes]],
]:
    """Perform initial cell calling."""
    try:
        filtered_metrics_groups, filtered_bcs_groups = call_initial_cells(
            cell_calling_args.matrix,
            cell_calling_args.genomes,
            cell_calling_args.probe_barcode_sample_id,
            cell_calling_args.unique_gem_groups,
            cell_calling_args.method,
            cell_calling_args.gg_recovered_cells,
            cell_calling_args.config.cell_barcodes,
            cell_calling_args.force_cells,
            cell_calling_args.feature_types,
            cell_calling_args.chemistry_description,
            cell_calling_args.target_features,
            cell_calling_args.has_cmo_data,
            num_probe_barcodes=cell_calling_args.num_probe_barcodes,
        )
        return filtered_metrics_groups, filtered_bcs_groups
    except InvalidCellBarcode as ex:
        # This is thrown deeper in the code, but caught here with a martian.exit call
        msg = (
            "Cell Barcodes did not match list of valid barcodes with observed reads.  "
            "If attempting manual cell calling, please be sure all input "
            f"barcodes are expected for that chemistry.\n\nError: {ex}"
        )
        martian.exit(msg)


def run_empty_drops(
    cell_calling_args: CellCallingArgs,
    filtered_bcs_groups: dict[tuple[int, str], list[bytes]],
) -> tuple[dict[tuple[int, str], list[bytes]], pd.DataFrame, dict[tuple[str, int], int]]:
    """Run the empty drops algorithm to call additional cells."""
    assert cell_calling_args.config.emptydrops_minimum_umis is not None
    (
        filtered_bcs_groups,
        nonambient_summary,
        emptydrops_minimum_umis,
    ) = call_additional_cells(
        cell_calling_args.matrix,
        cell_calling_args.unique_gem_groups,
        cell_calling_args.genomes,
        filtered_bcs_groups,
        cell_calling_args.feature_types,
        cell_calling_args.chemistry_description,
        cell_calling_args.probe_barcode_sample_id,
        num_probe_barcodes=cell_calling_args.num_probe_barcodes,
        emptydrops_minimum_umis=get_emptydrops_minimum_umis(
            cell_calling_args.config.emptydrops_minimum_umis,
            cell_calling_args.probe_barcode_sample_id,
        ),
    )
    return filtered_bcs_groups, nonambient_summary, emptydrops_minimum_umis


def filter_barcodes(args: ChunkInputs, outs: ChunkOutputs):
    """Identify cell-associated partitions."""
    cell_calling_args = get_cell_calling_args(args, outs)

    filtered_metrics_groups, filtered_bcs_groups = initial_cell_calling(cell_calling_args)

    if cell_calling_args.method == FilterMethod.ORDMAG_NONAMBIENT:
        (
            filtered_bcs_groups,
            nonambient_summary,
            emptydrops_minimum_umis,
        ) = run_empty_drops(
            cell_calling_args,
            filtered_bcs_groups,
        )
        if nonambient_summary.empty:
            outs.nonambient_calls = None
        else:
            nonambient_summary.to_csv(outs.nonambient_calls, index=False)
        for (genome, gg), value in emptydrops_minimum_umis.items():
            cell_calling_args.summary[
                f"{genome}_gem_group_{gg}_{get_filter_method_name(FilterMethod.ORDMAG_NONAMBIENT)}_threshold"
            ] = value
    else:
        outs.nonambient_calls = None

    ### Apply a minimum UMI threshold on cell calls
    if cell_calling_args.method not in [FilterMethod.MANUAL, FilterMethod.TOP_N_BARCODES]:
        assert cell_calling_args.config.global_minimum_umis is not None
        filtered_bcs_groups = apply_global_minimum_umis_threshold(
            cell_calling_args.matrix,
            cell_calling_args.unique_gem_groups,
            cell_calling_args.genomes,
            filtered_bcs_groups,
            cell_calling_args.feature_types,
            global_minimum_umis=get_global_minimum_umis(
                cell_calling_args.config.global_minimum_umis,
                cell_calling_args.probe_barcode_sample_id,
            ),
        )

        ### Apply a max mitochondrial read percentage threshold on cell calls
        max_mito_percent = get_max_mito_percent(
            cell_calling_args.config.max_mito_percent,
            cell_calling_args.probe_barcode_sample_id,
        )
        filtered_bcs_groups, mitochondrial_summary = apply_mitochondrial_threshold(
            cell_calling_args.matrix,
            cell_calling_args.genomes,
            cell_calling_args.unique_gem_groups,
            filtered_bcs_groups,
            max_mito_percent,
        )
        mitochondrial_summary.to_csv(outs.co_mitochondrial_summary)
    else:
        outs.co_mitochondrial_summary = None

    def remap_keys_metrics(groups):
        return [{"key": k, "value": v.__dict__} for k, v in groups.items()]

    assert outs.filtered_metrics_groups is not None
    with open(outs.filtered_metrics_groups, "w") as f:
        dump_numpy(
            remap_keys_metrics(filtered_metrics_groups),
            f,
            indent=4,
            sort_keys=True,
        )
    assert outs.filtered_bcs_groups is not None
    save_filtered_bcs_groups(filtered_bcs_groups, outs.filtered_bcs_groups)

    assert outs.summary is not None
    with open(outs.summary, "w") as f:
        dump_numpy(cell_calling_args.summary, f, indent=4, sort_keys=True)


def parse_chunked_metrics(
    chunk_outs: list[ChunkOutputs],
) -> dict[MetricGroups, BarcodeFilterResults]:
    filtered_metrics_groups_list = []
    for chunk_out in chunk_outs:
        if chunk_out.filtered_metrics_groups is not None:
            with open(chunk_out.filtered_metrics_groups) as infile:
                filtered_metrics_groups_list.append(json.load(infile))

    filtered_metrics_groups: dict[MetricGroups, BarcodeFilterResults] = {}
    for filtered_metrics in filtered_metrics_groups_list:
        for groups in filtered_metrics:
            key_tuple = MetricGroups(
                groups["key"][0], groups["key"][1], groups["key"][2], groups["key"][3]
            )
            result = BarcodeFilterResults(0)
            result.filtered_bcs = groups["value"]["filtered_bcs"]
            result.filtered_bcs_var = groups["value"]["filtered_bcs_var"]
            result.filtered_bcs_cv = groups["value"]["filtered_bcs_cv"]
            result.filtered_bcs_lb = groups["value"]["filtered_bcs_lb"]
            result.filtered_bcs_ub = groups["value"]["filtered_bcs_ub"]
            result.filtered_bcs_cutoff = groups["value"]["filtered_bcs_cutoff"]
            filtered_metrics_groups[key_tuple] = result
    return filtered_metrics_groups


def parse_chunked_filtered_bcs(
    chunk_outs: list[ChunkOutputs],
) -> tuple[list[bytes], dict[str, list[bytes]]]:
    filtered_bcs_groups_list = [
        chunk_out.filtered_bcs_groups
        for chunk_out in chunk_outs
        if chunk_out.filtered_bcs_groups is not None
    ]

    return load_filtered_bcs_groups(filtered_bcs_groups_list)


def parse_chunked_summary(chunk_outs: list[ChunkOutputs]) -> dict[str, int]:
    summary: dict[str, int] = defaultdict(int)
    for chunk_out in chunk_outs:
        if chunk_out.summary is not None:
            with open(chunk_out.summary) as infile:
                per_chunk_summary = json.load(infile)
                for k, v in per_chunk_summary.items():
                    summary[k] += v
    return dict(summary)


def parse_filtered_bcs_method(
    filtered_metrics_groups: Mapping[MetricGroups, BarcodeFilterResults],
    summary: dict,
) -> None:
    for key_tuple in filtered_metrics_groups:
        if key_tuple.sample is not None:
            summary.update({key_tuple.sample + "_filter_barcodes_method": key_tuple.method})
        else:
            summary.update({"filter_barcodes_method": key_tuple.method})
