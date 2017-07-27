#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.multigenome as cr_mg_analysis
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils

__MRO__ = """
stage RUN_MULTIGENOME_ANALYSIS(
    in  h5   raw_matrices_h5,
    in  h5   filtered_matrices_h5,
    in  bool is_multi_genome,
    in  bool skip,
    out path multi_genome_csv,
    out path multi_genome_json,
    out json multi_genome_summary,
    src py   "stages/analyzer/run_multigenome_analysis",
) split using (
)
"""

def split(args):
    if args.skip or not args.is_multi_genome:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    chunks = [{
        '__mem_gb': round(1.5 * cr_constants.MIN_MEM_GB)
    }]
    return {'chunks': chunks}

def main(args, outs):
    if args.skip or not args.is_multi_genome:
        return

    raw_matrices = cr_matrix.GeneBCMatrices.load_h5(args.raw_matrices_h5)
    filtered_matrices = cr_matrix.GeneBCMatrices.load_h5(args.filtered_matrices_h5)
    analysis = cr_mg_analysis.MultiGenomeAnalysis(raw_matrices, filtered_matrices)
    analysis.run_all()
    analysis.save_summary_json(outs.multi_genome_summary)
    analysis.save_gem_class_json(outs.multi_genome_json)
    analysis.save_gem_class_csv(outs.multi_genome_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or not args.is_multi_genome:
        return

    chunk_out = chunk_outs[0]
    cr_utils.copy(chunk_out.multi_genome_summary, outs.multi_genome_summary)
    cr_utils.copytree(chunk_out.multi_genome_csv, outs.multi_genome_csv)
    cr_utils.copytree(chunk_out.multi_genome_json, outs.multi_genome_json)
