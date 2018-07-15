#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import os
import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.io as cr_io

__MRO__ = """
stage COMBINE_CLUSTERING(
    in  bool skip,
    in  h5   kmeans_h5,
    in  path kmeans_csv,
    in  h5   graphclust_h5,
    in  path graphclust_csv,
    out h5   clustering_h5,
    out path clustering_csv,
    src py   "stages/analyzer/combine_clustering",
)
"""

def copy_subdirs(src_dir, dest_dir):
    for subdir in os.listdir(src_dir):
        cr_io.copytree(os.path.join(src_dir, subdir), os.path.join(dest_dir, subdir))

def main(args, outs):
    if args.skip:
        return

    analysis_io.combine_h5_files([args.kmeans_h5, args.graphclust_h5], outs.clustering_h5,
                           [analysis_constants.ANALYSIS_H5_KMEANS_GROUP,
                            analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP])

    csv_path = os.path.join(outs.clustering_csv)
    cr_io.makedirs(csv_path, allow_existing=True)
    copy_subdirs(args.kmeans_csv, csv_path)
    copy_subdirs(args.graphclust_csv, csv_path)
