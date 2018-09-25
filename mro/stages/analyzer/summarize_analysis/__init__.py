#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io

import h5py as h5
import os
import json
import tables

__MRO__ = """
stage SUMMARIZE_ANALYSIS(
    in  h5   matrix_h5,
    in  h5   pca_h5,
    in  h5   clustering_h5,
    in  h5   diffexp_h5,
    in  h5   tsne_h5,
    in  path pca_csv,
    in  path clustering_csv,
    in  path diffexp_csv,
    in  path tsne_csv,
    in  json multi_genome_summary,
    in  path multi_genome_csv,
    in  path multi_genome_json,
    in  bool is_multi_genome,
    in  bool batch_alignment,
    in  float batch_score_before_alignment,
    in  float batch_score_after_alignment,
    in  bool skip,
    out path analysis,
    out path analysis_csv,
    out json summary,
    src py   "stages/analyzer/summarize_analysis",
) split using (
)
"""

def split(args):
    chunks = [{
        '__mem_gb': h5_constants.MIN_MEM_GB
    }]
    return {'chunks': chunks}

def main(args, outs):
    if args.skip:
        return

    if args.is_multi_genome:
        cr_io.copytree(args.multi_genome_json, outs.analysis)
        cr_io.copytree(args.multi_genome_csv, outs.analysis_csv)

    analysis_h5 = analysis_io.h5_path(outs.analysis)
    cr_io.makedirs(os.path.dirname(analysis_h5), allow_existing=True)

    # Pytables doesn't support variable len strings, so use h5py first
    with h5.File(args.matrix_h5, 'r') as matrix,\
         h5.File(analysis_h5, 'w') as out:
        # TODO: copy the first group; fixme when we have a key
        name = matrix.keys()[0]
        matrix.copy(matrix[name], out, name='matrix')

    with tables.open_file(args.pca_h5, 'r') as pca,\
         tables.open_file(args.clustering_h5, 'r') as clustering,\
         tables.open_file(args.diffexp_h5, 'r') as diffexp,\
         tables.open_file(args.tsne_h5, 'r') as tsne,\
         tables.open_file(analysis_h5, 'a') as out:

         pca.copy_children(pca.root, out.root, recursive=True)
         clustering.copy_children(clustering.root, out.root, recursive=True)
         diffexp.copy_children(diffexp.root, out.root, recursive=True)
         tsne.copy_children(tsne.root, out.root, recursive=True)

    pca_dir = os.path.join(outs.analysis_csv, 'pca')
    cr_io.copytree(args.pca_csv, pca_dir)

    clustering_dir = os.path.join(outs.analysis_csv, 'clustering')
    cr_io.copytree(args.clustering_csv, clustering_dir)

    diffexp_dir = os.path.join(outs.analysis_csv, 'diffexp')
    cr_io.copytree(args.diffexp_csv, diffexp_dir)

    tsne_dir = os.path.join(outs.analysis_csv, 'tsne')
    cr_io.copytree(args.tsne_csv, tsne_dir)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        outs.analysis = None
        outs.analysis_csv = None
        outs.summary = None
        return

    chunk_out = chunk_outs[0]
    cr_io.copytree(chunk_out.analysis, outs.analysis)
    cr_io.copytree(chunk_out.analysis_csv, outs.analysis_csv)

    # batch alignment summary
    summary = {
        'batch_alignment': args.batch_alignment,
        'batch_effect_score_before_alignment': args.batch_score_before_alignment, 
        'batch_effect_score_after_alignment': args.batch_score_after_alignment,        
    }

    if args.is_multi_genome:
        with open(args.multi_genome_summary) as reader:
            multi_genome_summary = json.load(reader)
        summary.update(multi_genome_summary)
    else:
        summary = summary

    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)
