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


def hardlink_files_in_dir(src_dir, dest_dir):
    """Hardlink all files in src_dir to dest_dir."""
    for filename in os.listdir(src_dir):
        cr_io.hardlink_with_fallback(
            os.path.join(src_dir, filename),
            os.path.join(dest_dir, filename),
        )


def main(args, outs):
    analysis_io.combine_h5_files(
        [x for x in (args.kmeans_h5, args.graphclust_h5, args.hclust_h5) if x is not None],
        outs.clustering_h5,
        [
            analysis_constants.ANALYSIS_H5_KMEANS_GROUP,
            analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP,
        ],
    )

    os.makedirs(outs.clustering_csv, exist_ok=True)
    for csv_directory in [args.kmeans_csv, args.graphclust_csv, args.hclust_csv]:
        if csv_directory is not None:
            hardlink_files_in_dir(csv_directory, outs.clustering_csv)
