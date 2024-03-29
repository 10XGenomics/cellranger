#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import math

import cellranger.analysis.multigenome as cr_mg_analysis
import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants
from cellranger.matrix import CountMatrix

__MRO__ = """
stage RUN_MULTIGENOME_ANALYSIS(
    in  h5   filtered_matrices_h5,
    in  bool is_multi_genome,
    out path multi_genome_csv,
    out path multi_genome_json,
    out json multi_genome_summary,
    src py   "stages/analyzer/run_multigenome_analysis",
) split using (
)
"""


def split(args):
    if not args.is_multi_genome:
        return {"chunks": []}

    mem_gb = CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrices_h5)

    chunks = [{"__mem_gb": int(math.ceil(max(h5_constants.MIN_MEM_GB, mem_gb)))}]
    return {"chunks": chunks}


def main(args, outs):
    filtered_matrix = CountMatrix.load_h5_file(args.filtered_matrices_h5)
    analysis = cr_mg_analysis.MultiGenomeAnalysis(filtered_matrix)
    analysis.run_all()
    analysis.save_summary_json(outs.multi_genome_summary)
    analysis.save_gem_class_json(outs.multi_genome_json)
    analysis.save_gem_class_csv(outs.multi_genome_csv)


def join(args, outs, chunk_defs, chunk_outs):
    if not args.is_multi_genome:
        return

    chunk_out = chunk_outs[0]
    cr_io.hardlink_with_fallback(chunk_out.multi_genome_summary, outs.multi_genome_summary)
    cr_io.hardlink_with_fallback(chunk_out.multi_genome_csv, outs.multi_genome_csv)
    cr_io.hardlink_with_fallback(chunk_out.multi_genome_json, outs.multi_genome_json)
