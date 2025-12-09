#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import gzip
import pickle

import martian

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mc
import cellranger.subsample as cr_ss
import tenkit.safe_json as tk_safe_json
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

__MRO__ = """
stage SUBSAMPLE_READS(
    in  h5        molecule_info,
    in  csv       filtered_barcodes,
    in  string    target_mode,
    out json      summary,
    src py        "stages/counter/subsample_reads",
) split (
    in  int       chunk_start,
    in  int       chunk_len,
    in  map[]     subsample_info,
    out pickle.gz metrics,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


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

    with cr_mc.MoleculeCounter.open(args.molecule_info, "r") as mc:
        num_filtered_barcodes = mc.get_num_filtered_barcodes()
        num_subsamplings = len(subsamplings)
        mem_gib_chunk = 7 + round(100 * num_filtered_barcodes * num_subsamplings / 1024**3, 1)
        mem_gib_join = 2 + round(112 * num_filtered_barcodes * num_subsamplings / 1024**3, 1)
        print(f"{num_filtered_barcodes=},{num_subsamplings=},{mem_gib_chunk=},{mem_gib_join=}")
        chunks = [
            {
                "chunk_start": chunk_start,
                "chunk_len": chunk_len,
                "subsample_info": subsamplings,
                "__mem_gb": mem_gib_chunk,
                # cr_ss.run_subsampling requires additional VMEM
                "__vmem_gb": 3 + 3 * mem_gib_chunk,
            }
            for chunk_start, chunk_len in mc.get_chunks(
                cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK, preserve_boundaries=True
            )
        ]

    return {
        "chunks": chunks,
        "join": {"__mem_gb": mem_gib_join},
    }


def main(args, outs):
    if args.target_mode is None:
        feature_indices = None
    else:
        with cr_mc.MoleculeCounter.open(args.molecule_info, "r") as mc:
            assert mc.feature_reference is not None
            feature_indices = mc.feature_reference.get_target_feature_indices()
            assert feature_indices is not None
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
    with gzip.open(
        outs.metrics,
        "wb",
        compresslevel=1,  # Data is very sparse, so fast but inefficient compression is fine
    ) as outfile:
        pickle.dump(data, outfile, protocol=pickle.HIGHEST_PROTOCOL)


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
