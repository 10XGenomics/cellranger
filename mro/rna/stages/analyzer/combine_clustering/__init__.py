#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#


import os

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
import cellranger.cr_io as cr_io

__MRO__ = """
stage COMBINE_CLUSTERING(
    in  h5   kmeans_h5,
    in  path kmeans_csv,
    in  h5   graphclust_h5,
    in  path graphclust_csv,
    in  h5   hclust_h5,
    in  path hclust_csv,
    out h5   clustering_h5,
    out path clustering_csv,
    src py   "stages/analyzer/combine_clustering",
) using (
    volatile = strict,
)
"""


def copy_subdirs(src_dir, dest_dir):
    for subdir in os.listdir(src_dir):
        cr_io.hardlink_with_fallback(
            os.path.join(src_dir, subdir),
            os.path.join(dest_dir, subdir),
        )


def main(args, outs):
    list_of_hdfs = [args.kmeans_h5, args.graphclust_h5]
    if args.hclust_h5:
        list_of_hdfs.append(args.hclust_h5)
    analysis_io.combine_h5_files(
        list_of_hdfs,
        outs.clustering_h5,
        [
            analysis_constants.ANALYSIS_H5_KMEANS_GROUP,
            analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP,
        ],
    )

    csv_path = os.path.join(outs.clustering_csv)
    os.makedirs(csv_path, exist_ok=True)
    copy_subdirs(args.kmeans_csv, csv_path)
    copy_subdirs(args.graphclust_csv, csv_path)
    if args.hclust_csv:
        copy_subdirs(args.hclust_csv, csv_path)
