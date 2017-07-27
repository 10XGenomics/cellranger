#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

import os
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
        '__mem_gb': cr_constants.MIN_MEM_GB
    }]
    return {'chunks': chunks}

def main(args, outs):
    if args.skip:
        return

    if args.is_multi_genome:
        cr_utils.copytree(args.multi_genome_json, outs.analysis)
        cr_utils.copytree(args.multi_genome_csv, outs.analysis_csv)
        return

    analysis_h5 = cr_io.h5_path(outs.analysis)
    cr_utils.makedirs(os.path.dirname(analysis_h5), allow_existing=True)

    with tables.open_file(args.matrix_h5, 'r') as matrix,\
         tables.open_file(args.pca_h5, 'r') as pca,\
         tables.open_file(args.clustering_h5, 'r') as clustering,\
         tables.open_file(args.diffexp_h5, 'r') as diffexp,\
         tables.open_file(args.tsne_h5, 'r') as tsne,\
         tables.open_file(analysis_h5, 'w') as out:

         # NOTE - genome name is replaced with 'matrix'
         mat_groups = [m for m in matrix.root]
         matrix.copy_node(mat_groups[0], out.root, recursive=True, newname='matrix')

         pca.copy_children(pca.root, out.root, recursive=True)
         clustering.copy_children(clustering.root, out.root, recursive=True)
         diffexp.copy_children(diffexp.root, out.root, recursive=True)
         tsne.copy_children(tsne.root, out.root, recursive=True)

    pca_dir = os.path.join(outs.analysis_csv, 'pca')
    cr_utils.copytree(args.pca_csv, pca_dir)

    clustering_dir = os.path.join(outs.analysis_csv, 'clustering')
    cr_utils.copytree(args.clustering_csv, clustering_dir)

    diffexp_dir = os.path.join(outs.analysis_csv, 'diffexp')
    cr_utils.copytree(args.diffexp_csv, diffexp_dir)

    tsne_dir = os.path.join(outs.analysis_csv, 'tsne')
    cr_utils.copytree(args.tsne_csv, tsne_dir)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        outs.analysis = None
        outs.analysis_csv = None
        outs.summary = None
        return

    chunk_out = chunk_outs[0]
    cr_utils.copytree(chunk_out.analysis, outs.analysis)
    cr_utils.copytree(chunk_out.analysis_csv, outs.analysis_csv)

    if args.is_multi_genome:
        cr_utils.copy(args.multi_genome_summary, outs.summary)
    else:
        outs.summary = None
