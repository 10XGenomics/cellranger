#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.multigenome as cr_mg_analysis
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io

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
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    chunks = [{
        '__mem_gb': round(1.5 * h5_constants.MIN_MEM_GB)
    }]
    return {'chunks': chunks}

def main(args, outs):
    if args.skip or not args.is_multi_genome:
        return

    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.raw_matrices_h5)
    filtered_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrices_h5)
    analysis = cr_mg_analysis.MultiGenomeAnalysis(raw_matrix, filtered_matrix)
    analysis.run_all()
    analysis.save_summary_json(outs.multi_genome_summary)
    analysis.save_gem_class_json(outs.multi_genome_json)
    analysis.save_gem_class_csv(outs.multi_genome_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or not args.is_multi_genome:
        return

    chunk_out = chunk_outs[0]
    cr_io.copy(chunk_out.multi_genome_summary, outs.multi_genome_summary)
    cr_io.copytree(chunk_out.multi_genome_csv, outs.multi_genome_csv)
    cr_io.copytree(chunk_out.multi_genome_json, outs.multi_genome_json)
