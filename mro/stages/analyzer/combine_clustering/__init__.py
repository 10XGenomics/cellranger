#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import os
import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

__MRO__ = """
stage COMBINE_CLUSTERING(
    in  bool skip,
    in  bool is_multi_genome,
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
        cr_utils.copytree(os.path.join(src_dir, subdir), os.path.join(dest_dir, subdir))

def main(args, outs):
    if args.skip or args.is_multi_genome:
        return

    cr_io.combine_h5_files([args.kmeans_h5, args.graphclust_h5], outs.clustering_h5,
                           [cr_constants.ANALYSIS_H5_KMEANS_GROUP,
                            cr_constants.ANALYSIS_H5_CLUSTERING_GROUP])

    csv_path = os.path.join(outs.clustering_csv)
    cr_utils.makedirs(csv_path, allow_existing=True)
    copy_subdirs(args.kmeans_csv, csv_path)
    copy_subdirs(args.graphclust_csv, csv_path)
