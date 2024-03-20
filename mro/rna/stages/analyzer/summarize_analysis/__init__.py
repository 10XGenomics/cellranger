#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#


import json
import os

import h5py as h5
from six import ensure_binary

import cellranger.analysis.io as analysis_io
import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants
import cellranger.hdf5 as cr_h5

__MRO__ = """
stage SUMMARIZE_ANALYSIS(
    in  h5   matrix_h5,
    in  h5   pca_h5,
    in  h5   clustering_h5,
    in  h5   diffexp_h5,
    in  h5   tsne_h5,
    in  h5   umap_h5,
    in  path pca_csv,
    in  path clustering_csv,
    in  path diffexp_csv,
    in  path tsne_csv,
    in  path umap_csv,
    in  json multi_genome_summary,
    in  path multi_genome_csv,
    in  path multi_genome_json,
    in  bool is_multi_genome,
    in  bool chemistry_batch_correction,
    in  float batch_score_before_correction,
    in  float batch_score_after_correction,
    out path analysis,
    out path analysis_csv,
    out json summary,
    src py   "stages/analyzer/summarize_analysis",
) split using (
)
"""


def split(args):
    chunks = [{"__mem_gb": h5_constants.MIN_MEM_GB}]
    return {"chunks": chunks}


def main(args, outs):
    # Multi genome spatial jobs do not have a barnyard analysis.
    # the CSV is produced whenever the json is
    if args.is_multi_genome and args.multi_genome_json:
        cr_io.hardlink_with_fallback(args.multi_genome_json, outs.analysis)
        cr_io.hardlink_with_fallback(args.multi_genome_csv, outs.analysis_csv)

    analysis_h5 = analysis_io.h5_path(outs.analysis)
    os.makedirs(os.path.dirname(analysis_h5), exist_ok=True)

    # Pytables doesn't support variable len strings, so use h5py first
    with h5.File(ensure_binary(args.matrix_h5), "r") as matrix, h5.File(analysis_h5, "w") as out:
        # TODO: copy the first group; fixme when we have a key
        name = next(iter(matrix.keys()))
        matrix.copy(matrix[name], out, name="matrix")

    h5s_to_combine = [args.pca_h5, args.clustering_h5, args.diffexp_h5, args.umap_h5]
    if args.tsne_h5:
        h5s_to_combine.append(args.tsne_h5)
    cr_h5.combine_h5s_into_one(analysis_h5, h5s_to_combine)

    pca_dir = os.path.join(outs.analysis_csv, "pca")
    cr_io.hardlink_with_fallback(args.pca_csv, pca_dir)

    clustering_dir = os.path.join(outs.analysis_csv, "clustering")
    cr_io.hardlink_with_fallback(args.clustering_csv, clustering_dir)

    diffexp_dir = os.path.join(outs.analysis_csv, "diffexp")
    cr_io.hardlink_with_fallback(args.diffexp_csv, diffexp_dir)

    if args.tsne_csv:
        tsne_dir = os.path.join(outs.analysis_csv, "tsne")
        cr_io.hardlink_with_fallback(args.tsne_csv, tsne_dir)

    umap_dir = os.path.join(outs.analysis_csv, "umap")
    cr_io.hardlink_with_fallback(args.umap_csv, umap_dir)


def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]
    cr_io.hardlink_with_fallback(chunk_out.analysis, outs.analysis)
    cr_io.hardlink_with_fallback(chunk_out.analysis_csv, outs.analysis_csv)

    summary = {}

    # batch correction summary
    if args.chemistry_batch_correction is True:
        summary["batch_effect_score_before_correction"] = args.batch_score_before_correction
        summary["batch_effect_score_after_correction"] = args.batch_score_after_correction

    # Multigenome analysis is disabled in spatial even if the reference has multiple genomes
    if args.is_multi_genome and args.multi_genome_summary:
        with open(args.multi_genome_summary) as reader:
            multi_genome_summary = json.load(reader)
        summary.update(multi_genome_summary)

    with open(outs.summary, "w") as f:
        json.dump(summary, f, indent=4, sort_keys=True)
