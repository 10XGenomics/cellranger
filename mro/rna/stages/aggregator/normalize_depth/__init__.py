#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import array
import copy
import json
from collections import OrderedDict, defaultdict
from typing import Any

import martian
import numpy as np
import scipy.sparse as sp_sparse

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.molecule_counter as cr_mol_counter
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
from cellranger.fast_utils import count_usable_reads
from cellranger.logperf import LogPerf
from cellranger.matrix import CountMatrix
from cellranger.molecule_counter import BarcodeInfo, MoleculeCounter

__MRO__ = """
stage NORMALIZE_DEPTH(
    in  map        gem_group_index,
    in  h5         molecules,
    in  string     normalization_mode,
    in  map<int[]> gem_group_barcode_ranges,
    in  float      targeted_depth_factor,
    out h5[]       raw_matrices_h5,
    out int        raw_nnz,
    out h5[]       filtered_matrices_h5,
    out int        filtered_nnz,
    out json       summary,
    src py         "stages/aggregator/normalize_depth",
) split (
    in  float[]    frac_reads_kept,
    in  int[]      num_cells,
    in  int        chunk_start,
    in  int        chunk_len,
    in  json       reads_per_library,
    out json       chunk_summary,
    out h5         raw_matrix_h5,
    out h5         filtered_matrix_h5,
) using (
    mem_gb = 4,
    volatile = strict,
)
"""

TARGETED_READ_PAIRS_METRIC = "targeted_read_pairs"
NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK_RUST = 400000000


def _get_non_zero_min_rpc_by_lt(library_info: list[dict[str, Any]], usable_rpc: np.ndarray):
    """Determine lowest non-zero depth for each library type."""
    lt_rpcs = {lib[rna_library.LIBRARY_TYPE]: array.array("d") for lib in library_info}
    for lib, rpc in zip(library_info, usable_rpc):
        if rpc > 0:
            lt_rpcs[lib[rna_library.LIBRARY_TYPE]].append(rpc)
    return {lt: min(rpcs) if len(rpcs) else 0 for lt, rpcs in lt_rpcs.items()}


def _adjust_frac_kept(
    targeted_depth_factor: float, library_info: list[dict[str, Any]], frac_reads_kept: np.ndarray
):
    """Try to multiply all targeted sample `frac_reads_kept` by `targeted_depth_factor` if possible.

    this avoids throwing out very large fractions of the targeted samples when they are paired
    with control samples.
    """
    frac_reads_kept_adjusted = [
        (targeted_depth_factor if MoleculeCounter.is_targeted_library(lib) else 1.0)
        * frac_reads_kept[i]
        for i, lib in enumerate(library_info)
    ]
    if all(frac <= 1.0 for frac in frac_reads_kept_adjusted):
        return frac_reads_kept_adjusted
    return frac_reads_kept


def any_incorrect_usable_read_pairs(mc: MoleculeCounter) -> bool:
    """Return whether the metric usable_read_pairs is incorrect for any library.

    Return whether feature_read_pairs exceeds usable_read_pairs for any library.
    CR7.0 and CR7.1 incorrectly set usable_read_pairs equal to feature_read_pairs for MFRP.
    """
    return any(
        feature_read_pairs >= usable_read_pairs
        for feature_read_pairs, usable_read_pairs in zip(
            mc.get_transcriptomic_read_pairs_per_library(), mc.get_usable_read_pairs_per_library()
        )
    )


def split(args):
    # default to downsampling by mapped reads
    downsample = args.normalization_mode != cr_constants.NORM_MODE_NONE
    reads_per_library = None

    # compute downsample rates for each library
    with MoleculeCounter.open(args.molecules, "r") as mc:
        library_info = mc.get_library_info()
        usable_reads = mc.get_usable_read_pairs_per_library()

        assert mc.feature_reference is not None
        is_targeted_aggr = mc.feature_reference.has_target_features()
        if is_targeted_aggr:
            try:
                if any_incorrect_usable_read_pairs(mc):
                    raise ValueError
                # when the only GEX libraries present are targeted, can just fetch the on target usable reads
                gex_indices = mc.get_library_indices_by_type().get(
                    rna_library.GENE_EXPRESSION_LIBRARY_TYPE, []
                )
                all_libraries_usable_reads = mc.get_on_target_usable_read_pairs_per_library()
                for idx in gex_indices:
                    usable_reads[idx] = all_libraries_usable_reads[idx]
            except ValueError:
                # likely an aggr of WGS+Targeted libraries. Cannot fetch usable targeted reads from all of them,
                # instead calculate that separately for each library by iterating over the molecule info file
                print("Re-calculating on_target_usable_read_pairs")
                usable_reads_per_library = count_usable_reads(args.molecules)
                for library, on_target_usable_read_pairs in usable_reads_per_library.items():
                    usable_reads[library] = on_target_usable_read_pairs

        reads_per_library = martian.make_path("reads_per_library.json")
        with open(reads_per_library, "w") as f:
            json.dump(usable_reads, f, indent=4)

        cells = mc.get_num_filtered_barcodes_for_libraries(np.arange(0, len(library_info)))

        print(f"Libraries: {library_info}")
        print(f"Usable reads: {usable_reads}")
        print(f"Cells: {cells}")

        assert len(library_info) == len(usable_reads)
        usable_reads = np.array(usable_reads, dtype=np.float64)
        usable_rpc = np.divide(
            np.array(usable_reads, dtype=np.float64),
            cells.astype(np.float64),
            out=np.zeros(len(library_info), dtype=np.float64),
            where=cells > 0,
        )

    # Determine lowest depth for each library type
    min_rpc_by_lt = _get_non_zero_min_rpc_by_lt(library_info, usable_rpc)

    for lib_idx, lib in enumerate(library_info):
        lib_type = lib[rna_library.LIBRARY_TYPE]
        print(f"{lib_type} Usable read pairs per cell: {usable_rpc[lib_idx]}")
        print("%s Minimum read pairs usable per cell: %d" % (lib_type, min_rpc_by_lt[lib_type]))

    if not downsample:
        frac_reads_kept = np.ones(len(library_info), dtype=float)
    else:
        frac_reads_kept = np.zeros(len(library_info), dtype=float)
        for i, lib in enumerate(library_info):
            lib_type = lib[rna_library.LIBRARY_TYPE]
            min_rpc = min_rpc_by_lt[lib_type]
            if min_rpc != 0:
                frac_reads_kept[i] = tk_stats.robust_divide(min_rpc, usable_rpc[i])

        if is_targeted_aggr:
            # If we're not downsampling, it's all zeroes anyway.
            frac_reads_kept = _adjust_frac_kept(
                args.targeted_depth_factor, library_info, frac_reads_kept
            )

    # Split the molecule info h5 into equi-RAM chunks, preserving (barcode, gem_group) boundaries
    # Assumes the molecule_info is sorted by (gem_group, barcode)
    tgt_chunk_len = NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK_RUST

    chunks = []

    with MoleculeCounter.open(args.molecules, "r") as mc:
        # Number of barcodes in the full matrix
        num_barcodes = mc.get_ref_column_lazy("barcodes").shape[0]

        for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
            # 1024**3 / (1e9 / 28) = 30.07 bytes per molecule
            mol_mem_gib = MoleculeCounter.estimate_mem_gb(chunk_len, scale=1.0, cap=False)

            # Worst case number of nonzero elements in chunk matrix
            # 1024**3 / 50e6 = 21.47 bytes per molecule
            # 1024**3 / 10e6 = 107.37 bytes per barcode
            matrix_mem_gib = CountMatrix.get_mem_gb_from_matrix_dim(
                num_barcodes, chunk_len, scale=1.0
            )

            barcodes_mem_gib = round(140 * num_barcodes / 1024**3, 1)

            # bytes per molecule = 30.07 + 21.47 = 51.54
            # bytes per barcode = 107.37 + 140 = 247.37
            mem_gib = 2 + matrix_mem_gib + mol_mem_gib + barcodes_mem_gib
            print(
                f"chunk={len(chunks)},{chunk_len=},{num_barcodes=},{mol_mem_gib=},{matrix_mem_gib=},{barcodes_mem_gib=},{mem_gib=}"
            )

            chunks.append(
                {
                    "frac_reads_kept": list(frac_reads_kept),
                    "num_cells": [int(x) for x in cells],
                    "chunk_start": chunk_start,
                    "chunk_len": chunk_len,
                    "reads_per_library": reads_per_library,
                    # Request enough for two copies
                    "__mem_gb": mem_gib,
                }
            )

    # Join is not loading the merged matrix, so it doesn't need as much memory as main.
    # Nonetheless, it needs to load the barcodes, which can get large when many samples
    # are being aggregated.
    # WRITE_MATRICES will use the precise nnz counts to make an appropriate mem request.
    join_mem_gib = 1 + round(56 * num_barcodes / 1024**3, 1)
    print(f"{num_barcodes=},{join_mem_gib=}")
    return {"chunks": chunks, "join": {"__mem_gb": join_mem_gib, "__threads": 2}}


def summarize_read_matrix(matrix, library_info, barcode_info, barcode_seqs):
    """Summarize matrix of read-pair counts."""
    lib_types = rna_library.sorted_library_types(library_info)

    view = matrix.view()
    summary = {}

    for lib_type in lib_types:
        if rna_library.has_genomes(lib_type):
            sum_genomes = [str(x) for x in barcode_info.genomes]
        else:
            sum_genomes = [rna_library.MULTI_REFS_PREFIX]

        for genome in sum_genomes:
            m = view.select_features_by_type(lib_type)
            if rna_library.has_genomes(lib_type):
                m = m.select_features_by_genome(genome)
                genome_idx = barcode_info.genomes.index(genome)
            else:
                genome_idx = None

            prefix = f"{rna_library.get_library_type_metric_prefix(lib_type)}{genome}"
            summary[f"{prefix}_raw_mapped_reads"] = m.sum()

            filtered_bcs = MoleculeCounter.get_filtered_barcodes(
                barcode_info,
                library_info,
                barcode_seqs,
                genome_idx=genome_idx,
                library_type=lib_type,
            )
            filtered_m = m.select_barcodes_by_seq(filtered_bcs)
            summary[f"{prefix}_flt_mapped_reads"] = filtered_m.sum()

    return summary


def _make_metrics_out(args, library_info: list[dict[str, Any]], mc: MoleculeCounter):
    metrics_in = mc.get_all_metrics()
    metrics_out = copy.deepcopy(metrics_in)

    # Compute subsampling rate and approximate new total readpair count
    frac_reads_kept = np.array(args.frac_reads_kept, dtype=float)
    assert len(frac_reads_kept) == len(library_info)

    for lib_idx, total_reads in enumerate(
        mc.get_raw_read_pairs_per_library(metrics_in) * frac_reads_kept
    ):
        metrics_out[cr_mol_counter.LIBRARIES_METRIC][str(lib_idx)][
            cr_mol_counter.DOWNSAMPLED_READS_METRIC
        ] = total_reads

    # this is a relatively new metric and might not exist in older versions of molecule info
    try:
        total_feature_reads_in = mc.get_transcriptomic_read_pairs_per_library()
    except ValueError:
        total_feature_reads_in = np.zeros(len(library_info))
    total_feature_reads_out = total_feature_reads_in * frac_reads_kept

    for lib_idx, _ in enumerate(library_info):
        metrics_out[cr_mol_counter.LIBRARIES_METRIC][str(lib_idx)][
            cr_mol_counter.DOWNSAMPLED_FEATURE_READS_METRIC
        ] = total_feature_reads_out[lib_idx]

    assert mc.feature_reference is not None
    if mc.feature_reference.has_target_features():
        with open(args.reads_per_library) as f:
            usable_reads = json.load(f)
        for lib_idx, lib in enumerate(library_info):
            lib_metrics = metrics_out[cr_mol_counter.LIBRARIES_METRIC][str(lib_idx)]
            if (
                lib_metrics[cr_mol_counter.FEATURE_READS_METRIC]
                >= lib_metrics[cr_mol_counter.USABLE_READS_METRIC]
                and lib[rna_library.LIBRARY_TYPE] == rna_library.ANTIBODY_LIBRARY_TYPE
            ):
                # CR7.0 and CR7.1 incorrectly set usable_read_pairs equal to feature_read_pairs.
                lib_metrics[cr_mol_counter.USABLE_READS_METRIC] = usable_reads[lib_idx]

            lib_metrics[TARGETED_READ_PAIRS_METRIC] = (
                usable_reads[lib_idx]
                if lib[rna_library.LIBRARY_TYPE] == rna_library.GENE_EXPRESSION_LIBRARY_TYPE
                else 0
            )

    return metrics_out, frac_reads_kept


def _get_new_read_pairs(
    chunk: slice,
    mc: MoleculeCounter,
    frac_reads_kept: np.ndarray,
) -> np.ndarray:
    return np.random.binomial(
        mc.get_column_lazy("count")[chunk],
        frac_reads_kept[mc.get_column_lazy("library_idx")[chunk]],
    )


def _get_matrix(
    chunk: slice,
    mc: MoleculeCounter,
    new_read_pairs: np.ndarray,
    keep_mol: np.ndarray,
    gg_barcode_idx_start: np.ndarray,
    gg_barcode_idx_len: np.ndarray,
    num_features: int,
    num_bcs: int,
):
    mol_gem_group = mc.get_column_lazy("gem_group")[chunk][keep_mol]
    mol_feature_idx = mc.get_column_lazy("feature_idx")[chunk][keep_mol]
    mol_barcode_idx = mc.get_column_lazy("barcode_idx")[chunk][keep_mol]

    # Convert molecule barcode indices into matrix barcode indices
    # The molecule info barcode_idx is in this space:
    #  [W_0, W_1, ...] where W_i is distinct original whitelist i.
    # The matrix is in, e.g., this space:
    #  [w_0-1, w_1-2, w_0-3, ...] where w_i-j is a copy of whitelist i for gem group j.

    # Return to the original whitelist index
    mol_barcode_idx -= gg_barcode_idx_start.astype(np.uint64)[mol_gem_group]

    # Offset by the cumulative whitelist length up to a barcode's gem group
    gg_barcode_matrix_start = np.cumsum(gg_barcode_idx_len).astype(np.uint64)
    mol_barcode_idx += gg_barcode_matrix_start[mol_gem_group - 1]

    umi_matrix = sp_sparse.coo_matrix(
        (
            np.ones(len(mol_barcode_idx), dtype=cr_matrix.DEFAULT_DATA_DTYPE),
            (mol_feature_idx, mol_barcode_idx),
        ),
        shape=(num_features, num_bcs),
    )
    print("created umi matrix")
    LogPerf.mem()

    # Create a read-count matrix so we can summarize reads per barcode
    read_matrix = sp_sparse.coo_matrix(
        (new_read_pairs, (mol_feature_idx, mol_barcode_idx)), shape=(num_features, num_bcs)
    )
    return read_matrix, umi_matrix


def _get_barcode_idxs(gem_groups, gem_group_barcode_ranges):
    # Get the range of possible barcode indices for each gem group.
    idx_ranges = np.array(
        [(gg, idx[0], idx[1]) for gg, idx in gem_group_barcode_ranges], dtype=np.int64
    )
    # Sort by gem group
    idx_ranges = idx_ranges[np.argsort(idx_ranges[:, 0]), :]
    # Allocate as a single block for improved memory locality.
    gg_barcode_idx: np.ndarray = np.zeros((1 + len(gem_groups), 2), dtype=np.int64)
    gg_barcode_idx[idx_ranges[:, 0]] = idx_ranges[:, 1:3]
    gg_barcode_idx[:, 1] = np.subtract(
        gg_barcode_idx[:, 1], gg_barcode_idx[:, 0], out=gg_barcode_idx[:, 1]
    )
    return gg_barcode_idx[:, 0], gg_barcode_idx[:, 1]


def _update_metrics(
    args,
    frac_reads_kept: np.ndarray,
    mc: MoleculeCounter,
    library_info: list[dict[str, Any]],
    barcode_info: BarcodeInfo,
):
    # downsample molecule info
    chunk = slice(args.chunk_start, args.chunk_start + args.chunk_len)
    new_read_pairs = _get_new_read_pairs(chunk, mc, frac_reads_kept)
    keep_mol = np.flatnonzero(new_read_pairs)
    new_read_pairs = new_read_pairs[keep_mol]

    # Assert that gem groups start at 1 and are contiguous.  If they are,
    # then the sorted set of unique groups will be identically range [1, N].
    gem_groups = np.array(sorted({lib["gem_group"] for lib in library_info}), dtype=np.int32)
    assert gem_groups[0] == 1
    assert gem_groups[-1] == len(gem_groups)

    feature_ref = mc.get_feature_ref()
    num_features = feature_ref.get_num_features()

    # Compute matrix dimensions
    gg_barcode_idx_start, gg_barcode_idx_len = _get_barcode_idxs(
        gem_groups, args.gem_group_barcode_ranges.items()
    )

    num_bcs = gg_barcode_idx_len.sum()

    print("downsampled")
    LogPerf.mem()

    read_matrix, umi_matrix = _get_matrix(
        chunk,
        mc,
        new_read_pairs,
        keep_mol,
        gg_barcode_idx_start,
        gg_barcode_idx_len,
        num_features,
        num_bcs,
    )
    del new_read_pairs

    # Get all barcode strings and all features for the raw matrix
    barcode_seqs = mc.get_barcodes()
    feature_ref = mc.get_feature_ref()

    print(len(barcode_seqs), len(gem_groups))
    print("creating barcode strings")
    LogPerf.mem()

    barcode_list = []
    for gg, idx_start, idx_end in zip(
        gem_groups,
        gg_barcode_idx_start[1:],
        gg_barcode_idx_start[1:] + gg_barcode_idx_len[1:],
    ):
        chnksize = idx_end - idx_start
        if chnksize > 0:
            cur_bcs = barcode_seqs[idx_start:idx_end]  # this can be of size 0
            assert cur_bcs.shape[0] == chnksize, "Expected elements missing from BC array."
            new_bc = cr_utils.format_barcode_seq(cur_bcs[0], gg)
            barcode_dtype = np.dtype("S%d" % len(new_bc))
            barcode_list.append(
                np.fromiter(
                    (cr_utils.format_barcode_seq(bc, gg) for bc in cur_bcs),
                    dtype=barcode_dtype,
                    count=len(cur_bcs),
                )
            )
            del cur_bcs
    barcodes = np.concatenate(barcode_list)
    del barcode_list

    barcodes.flags.writeable = False

    print("created barcode strings")
    LogPerf.mem()

    # Get mapped reads per barcode per library,genome
    read_matrix = CountMatrix(feature_ref, barcodes, read_matrix)
    read_matrix.m = read_matrix.m.tocsc(copy=True)
    read_summary = summarize_read_matrix(read_matrix, library_info, barcode_info, barcode_seqs)

    del read_matrix

    print("created read matrix")
    LogPerf.mem()
    # Construct the raw UMI matrix
    raw_umi_matrix = CountMatrix(feature_ref, barcodes, umi_matrix)
    return raw_umi_matrix, barcode_seqs, read_summary


def main(args, outs):
    np.random.seed(0)

    LogPerf.mem()

    with MoleculeCounter.open(args.molecules, "r") as mc:
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()

        metrics_out, frac_reads_kept = _make_metrics_out(args, library_info, mc)

        raw_umi_matrix, barcode_seqs, read_summary = _update_metrics(
            args, frac_reads_kept, mc, library_info, barcode_info
        )
        raw_umi_matrix.save_h5_file(outs.raw_matrix_h5, sw_version=martian.get_pipelines_version())
        outs.raw_nnz = raw_umi_matrix.m.nnz

        # Construct the filtered UMI matrix
        filtered_bcs = MoleculeCounter.get_filtered_barcodes(
            barcode_info, library_info, barcode_seqs
        )
        filtered_umi_matrix = raw_umi_matrix.select_barcodes_by_seq(filtered_bcs)
        assert mc.feature_reference is not None
        if mc.feature_reference.has_target_features():
            target_feature_indices = filtered_umi_matrix.feature_ref.get_target_feature_indices()
            assert target_feature_indices is not None
            feature_indices_to_select = [
                fd.index
                for fd in filtered_umi_matrix.feature_ref.feature_defs
                if fd.feature_type != rna_library.GENE_EXPRESSION_LIBRARY_TYPE
                or fd.index in target_feature_indices
            ]
            filtered_umi_matrix = filtered_umi_matrix.select_features(feature_indices_to_select)
        filtered_umi_matrix.save_h5_file(
            outs.filtered_matrix_h5, sw_version=martian.get_pipelines_version()
        )
        outs.filtered_nnz = filtered_umi_matrix.m.nnz

        print("created filtered umi matrix")
        LogPerf.mem()

        summary = {
            "read_summary": read_summary,
            "mol_metrics": metrics_out,
        }

        with open(outs.chunk_summary, "w") as f:
            tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)

    # Don't write MEX from chunks.
    outs.raw_matrices_mex = None
    outs.filtered_matrices_mex = None


def join(args, outs, chunk_defs, chunk_outs):
    IS_SPATIAL = MoleculeCounter.open(args.molecules, "r").is_spatial_data()

    # Pass through the matrix chunks and nnz counts
    outs.raw_matrices_h5 = [o.raw_matrix_h5 for o in chunk_outs]
    outs.raw_nnz = sum(o.raw_nnz for o in chunk_outs)
    outs.filtered_matrices_h5 = [o.filtered_matrix_h5 for o in chunk_outs]
    outs.filtered_nnz = sum(o.filtered_nnz for o in chunk_outs)

    if IS_SPATIAL:
        summary = {
            "frac_reads_kept": chunk_defs[0].frac_reads_kept,
            "num_spots_by_library": chunk_defs[0].num_cells,
        }
    else:
        summary = {
            "frac_reads_kept": chunk_defs[0].frac_reads_kept,
            "num_cells_by_library": chunk_defs[0].num_cells,
        }

    with MoleculeCounter.open(args.molecules, "r") as mc:
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()
        barcode_seqs = mc.get_barcodes()

    lib_types = rna_library.sorted_library_types(library_info)
    for lib_type in lib_types:
        if rna_library.has_genomes(lib_type):
            genomes = [str(x) for x in barcode_info.genomes]
            genome_indices = [barcode_info.genomes.index(g) for g in genomes]
        else:
            genomes, genome_indices = [rna_library.MULTI_REFS_PREFIX], [None]

        for genome, genome_idx in zip(genomes, genome_indices):
            prefix = f"{rna_library.get_library_type_metric_prefix(lib_type)}{genome}"
            filtered_bcs = MoleculeCounter.get_filtered_barcodes(
                barcode_info,
                library_info,
                barcode_seqs,
                genome_idx=genome_idx,
                library_type=lib_type,
            )
            summary[f"{prefix}_filtered_bcs"] = len(filtered_bcs)

    # Merge read summary metrics
    read_summary = defaultdict(int)
    for filename in [co.chunk_summary for co in chunk_outs]:
        with open(filename) as f:
            d = json.load(f)
            for k, v in d["read_summary"].items():
                read_summary[k] += v
    summary.update(read_summary)

    # Get summary metrics
    with open(chunk_outs[0].chunk_summary) as f:
        mol_metrics = json.load(f)["mol_metrics"]
    summary.update((k, v) for k, v in mol_metrics.items() if k.startswith(("chemistry", "target")))
    print(json.dumps(mol_metrics, indent=4, sort_keys=True))

    # Report normalization metrics
    all_batches = OrderedDict()

    # These are all per-library-type
    min_frac_reads_kept = np.ones(len(lib_types), dtype="float")
    total_raw_read_pairs = np.zeros(len(lib_types), dtype="uint64")
    total_feature_read_pairs = np.zeros(len(lib_types), dtype="uint64")
    total_ds_raw_read_pairs = np.zeros(len(lib_types), dtype="uint64")
    total_ds_feature_read_pairs = np.zeros(len(lib_types), dtype="uint64")
    total_cells = np.zeros(len(lib_types), dtype="uint64")

    for lib_type_idx, lib_type in enumerate(lib_types):
        lib_inds = [i for i, lib in enumerate(library_info) if lib["library_type"] == lib_type]
        for lib_idx in lib_inds:
            aggr_id = library_info[lib_idx]["aggr_id"]
            old_gg = library_info[lib_idx]["old_gem_group"]
            batch = aggr_id + ("-%d" % old_gg if old_gg > 1 else "")
            all_batches[batch] = None
            if IS_SPATIAL:
                n_cells = summary["num_spots_by_library"][lib_idx]
            else:
                n_cells = summary["num_cells_by_library"][lib_idx]
            total_cells[lib_type_idx] += n_cells

            lib_metrics = mol_metrics[cr_mol_counter.LIBRARIES_METRIC][str(lib_idx)]
            raw_read_pairs = lib_metrics[cr_mol_counter.TOTAL_READS_METRIC]
            mapped_read_pairs = lib_metrics[cr_mol_counter.USABLE_READS_METRIC]
            # Cast string "NaN" to float("NaN")
            ds_read_pairs = float(lib_metrics[cr_mol_counter.DOWNSAMPLED_READS_METRIC])
            # these metrics are relatively new and might not exist in older versions of molecule info
            feature_reads = lib_metrics.get(cr_mol_counter.FEATURE_READS_METRIC, 0)
            ds_feature_reads = float(
                lib_metrics.get(cr_mol_counter.DOWNSAMPLED_FEATURE_READS_METRIC, 0)
            )

            total_raw_read_pairs[lib_type_idx] += raw_read_pairs
            total_ds_raw_read_pairs[lib_type_idx] += (
                ds_read_pairs if not np.isnan(ds_read_pairs) else 0
            )
            total_feature_read_pairs[lib_type_idx] += feature_reads
            total_ds_feature_read_pairs[lib_type_idx] += (
                ds_feature_reads if not np.isnan(ds_feature_reads) else 0
            )

            frac_reads_kept = summary["frac_reads_kept"][lib_idx]
            if not frac_reads_kept is None:
                min_frac_reads_kept[lib_type_idx] = min(
                    min_frac_reads_kept[lib_type_idx], frac_reads_kept
                )

            pre_norm_raw_rppc = tk_stats.robust_divide(raw_read_pairs, n_cells)
            pre_norm_mapped_rppc = tk_stats.robust_divide(mapped_read_pairs, n_cells)
            pre_norm_feature_rppc = tk_stats.robust_divide(feature_reads, n_cells)

            # analysis parameter
            gg = library_info[lib_idx]["gem_group"]
            gg_metrics = mol_metrics[cr_mol_counter.GEM_GROUPS_METRIC][str(gg)]
            intron_mode_param = gg_metrics.get(
                cr_mol_counter.INTRON_MODE_PARAM, cr_mol_counter.INTRON_MODE_HISTORIC_DEFAULT
            )

            # Prefix with batch and library type
            lib_prefix = rna_library.get_library_type_metric_prefix(lib_type)

            p = (batch, lib_prefix)
            summary.update(
                {
                    "{}_{}frac_reads_kept".format(*p): frac_reads_kept,
                    "{}_{}filtered_bcs".format(*p): n_cells,
                    "{}_{}pre_normalization_feature_reads".format(*p): feature_reads,
                    "{}_{}pre_normalization_raw_reads_per_filtered_bc".format(
                        *p
                    ): pre_norm_raw_rppc,
                    "{}_{}pre_normalization_cmb_reads_per_filtered_bc".format(
                        *p
                    ): pre_norm_mapped_rppc,
                    "{}_{}pre_normalization_feature_reads_per_filtered_bc".format(
                        *p
                    ): pre_norm_feature_rppc,
                    "{}_{}introns_included".format(*p): intron_mode_param,
                }
            )

            if TARGETED_READ_PAIRS_METRIC in lib_metrics.keys():
                targeted_read_pairs = lib_metrics[TARGETED_READ_PAIRS_METRIC]
                pre_norm_targeted_rpc = tk_stats.robust_divide(targeted_read_pairs, n_cells)
                summary.update(
                    {
                        "{}_{}pre_normalization_targeted_cmb_reads_per_filtered_bc".format(
                            *p
                        ): pre_norm_targeted_rpc,
                    }
                )

    summary["batches"] = list(all_batches)

    for lib_type_idx, lib_type in enumerate(lib_types):
        mean_rppc = tk_stats.robust_divide(
            total_raw_read_pairs[lib_type_idx], total_cells[lib_type_idx]
        )
        ds_mean_rppc = tk_stats.robust_divide(
            total_ds_raw_read_pairs[lib_type_idx], total_cells[lib_type_idx]
        )
        mean_feature_reads_pc = tk_stats.robust_divide(
            total_feature_read_pairs[lib_type_idx], total_cells[lib_type_idx]
        )
        ds_mean_feature_reads_pc = tk_stats.robust_divide(
            total_ds_feature_read_pairs[lib_type_idx], total_cells[lib_type_idx]
        )

        p = rna_library.get_library_type_metric_prefix(lib_type)
        summary.update(
            {
                f"{p}pre_normalization_total_reads": total_raw_read_pairs[lib_type_idx],
                f"{p}post_normalization_total_reads": total_ds_raw_read_pairs[lib_type_idx],
                f"{p}pre_normalization_total_feature_reads": total_feature_read_pairs[lib_type_idx],
                f"{p}post_normalization_total_feature_reads": total_ds_feature_read_pairs[
                    lib_type_idx
                ],
                f"{p}filtered_bcs_transcriptome_union": total_cells[lib_type_idx],
                f"{p}pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc": mean_rppc,
                f"{p}post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc": ds_mean_rppc,
                f"{p}pre_normalization_multi_transcriptome_total_feature_reads_per_filtered_bc": mean_feature_reads_pc,
                f"{p}post_normalization_multi_transcriptome_total_feature_reads_per_filtered_bc": ds_mean_feature_reads_pc,
                f"{p}lowest_frac_reads_kept": min_frac_reads_kept[lib_type_idx],
            }
        )

    with open(outs.summary, "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)
