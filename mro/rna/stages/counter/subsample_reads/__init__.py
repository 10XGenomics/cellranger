#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import math
import os
import pickle as cPickle

import martian

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mc
import cellranger.subsample as cr_ss
import tenkit.safe_json as tk_safe_json
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

__MRO__ = """
stage SUBSAMPLE_READS(
    in  h5     molecule_info,
    in  csv    filtered_barcodes,
    in  string target_mode,
    out json   summary,
    out pickle merged_metrics,
    src py     "stages/counter/subsample_reads",
) split (
    in  int    chunk_start,
    in  int    chunk_len,
    in  map[]  subsample_info,
    out pickle metrics,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""

MEM_GB_MAIN_PER_FILT_BCS_GB = 60.0
MEM_GB_JOIN_PER_FILT_BCS_GB = 2.0 * MEM_GB_MAIN_PER_FILT_BCS_GB


def split(args):
    # construct a range of subsampling depths and write out metadata
    if args.target_mode not in cr_ss.SUBSAMPLE_TARGET_MODES:
        martian.exit(f"target mode must be one of {cr_ss.SUBSAMPLE_TARGET_MODES}")
    is_targeted = args.target_mode is not None
    subsamplings = cr_ss.construct_all_subsamplings(
        args.molecule_info, args.filtered_barcodes, is_targeted=is_targeted, include_bulk=True
    )

    if len(subsamplings) == 0:
        return {"chunks": []}

    # Dynamically adjust mem_gb based on size of input filtered barcodes file
    filt_bcs_bytes = os.path.getsize(args.filtered_barcodes) if args.filtered_barcodes else 0
    mem_gb_main = math.ceil(MEM_GB_MAIN_PER_FILT_BCS_GB * float(filt_bcs_bytes) / 1e9)
    mem_gb_join = math.ceil(MEM_GB_JOIN_PER_FILT_BCS_GB * float(filt_bcs_bytes) / 1e9)

    # Split the molecule info h5 into equi-RAM chunks
    chunks = []
    # Use smaller chunks because this code is a slower than other users of mol_info
    tgt_chunk_len = cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK // 4

    with cr_mc.MoleculeCounter.open(args.molecule_info, "r") as mc:
        for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
            chunks.append(
                {
                    "chunk_start": int(chunk_start),
                    "chunk_len": int(chunk_len),
                    "subsample_info": subsamplings,
                    # The estimate_mem_gb only count the memory usage for the MoleculeCounter object, which is
                    # under-estimated the actual memory usage.
                    # Based on memory profiling with test case fuzzer_114, actual memory usage is ~4x more
                    # than estimate_mem_gb (without cap), here set scale = 6.
                    "__mem_gb": max(
                        mem_gb_main, cr_mc.MoleculeCounter.estimate_mem_gb(chunk_len, scale=6)
                    ),
                }
            )
    return {"chunks": chunks, "join": {"__mem_gb": max(mem_gb_join, 6)}}


def main(args, outs):
    if args.target_mode is None:
        feature_indices = None
    else:
        with cr_mc.MoleculeCounter.open(args.molecule_info, "r") as mc:
            feature_indices = mc.feature_reference.get_target_feature_indices()
            if args.target_mode == cr_constants.OFF_TARGET_SUBSAMPLE:
                # subsampling to off-target genes -- invert selection
                all_gene_indices = [
                    fd.index
                    for fd in mc.feature_reference.feature_defs
                    if fd.feature_type == GENE_EXPRESSION_LIBRARY_TYPE
                ]
                feature_indices = sorted(set(all_gene_indices) - set(feature_indices))

    data = cr_ss.run_subsampling(
        args.molecule_info,
        args.subsample_info,
        args.filtered_barcodes,
        feature_indices,
        args.chunk_start,
        args.chunk_len,
    )
    with open(outs.metrics, "wb") as f:
        cPickle.dump(data, f, protocol=cPickle.HIGHEST_PROTOCOL)


def join(args, outs, chunk_defs, chunk_outs):
    # Merge tallies
    metrics = [chunk.metrics for chunk in chunk_outs]
    data = cr_ss.join_metrics(metrics)

    subsample_info = chunk_defs[0].subsample_info if len(chunk_defs) > 0 else []
    summary = cr_ss.calculate_subsampling_metrics(
        data, args.molecule_info, args.filtered_barcodes, subsample_info, args.target_mode
    )

    with open(outs.summary, "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)

    with open(outs.merged_metrics, "wb") as f:
        cPickle.dump([data, subsample_info], f, protocol=cPickle.HIGHEST_PROTOCOL)
