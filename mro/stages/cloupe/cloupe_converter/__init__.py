#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

import martian
import os
import subprocess

__MRO__ = """
stage CLOUPE_CONVERTER(
    in  string pipestance_type,
    in  string sample_id,
    in  path   pipestance_path,
    out cloupe output_for_cloupe,
    src py     "stages/cloupe/cloupe_converter",
)
"""

def main(args, outs):
    if args.pipestance_type != "count" and args.pipestance_type != "aggr":
        martian.exit("The type argument must be one of: count, aggr")

    if args.pipestance_type == "count":
        pname = "SC_RNA_COUNTER_CS"
    if args.pipestance_type == "aggr":
        pname = "SC_RNA_AGGREGATOR_CS"

    pipestance_exists = os.path.exists(args.pipestance_path)
    if not pipestance_exists:
        martian.exit("Invalid pipestance path: %s" % args.pipestance_path)

    # check to see if an analysis file exists.  If it doesn't, then
    # this is likely a barnyard sample, and we cannot generate a
    # .loupe file (CELLRANGER-773);
    analysis_h5_path = os.path.join(args.pipestance_path, "outs/analysis/analysis.h5")

    # 1.2.0 location only
    internal_count_h5_path = os.path.join(
        args.pipestance_path,
        "SC_RNA_COUNTER_CS/SC_RNA_COUNTER/SC_RNA_ANALYZER/SUMMARIZE_ANALYSIS/fork0/files/analysis/analysis.h5"
    )

    internal_aggr_h5_path = os.path.join(
        args.pipestance_path,
        "SC_RNA_AGGREGATOR_CS/SC_RNA_AGGREGATOR/SC_RNA_ANALYZER/SUMMARIZE_ANALYSIS/fork0/files/analysis/analysis.h5"
    )

    if not os.path.exists(analysis_h5_path) \
            and not os.path.exists(internal_count_h5_path) \
            and not os.path.exists(internal_aggr_h5_path):
        martian.exit("Could not find single-species analysis HDF5 file. " +
                     "Loupe Cell Browser files are not generated for multi-species experiments.")

    # has to be 1.2 or higher
    cellranger_pd_before_1_2_path = os.path.join(args.pipestance_path, "CELLRANGER_PD")
    cellranger_cs_before_1_2_path = os.path.join(args.pipestance_path, "CELLRANGER_CS")
    if os.path.exists(cellranger_pd_before_1_2_path) or os.path.exists(cellranger_cs_before_1_2_path):
        martian.exit("mkloupe is only supported for Cell Ranger 1.2 and later.")

    call = ["crconverter",
            args.sample_id,
            pname,
            "--pipestance", args.pipestance_path,
            "--output", outs.output_for_cloupe]

    martian.log_info("Running crconverter: %s" % " ".join(call))
    try:
        results = subprocess.check_output(call)
        martian.log_info("crconverter output: %s" % results)
    except subprocess.CalledProcessError, e:
        outs.output_for_cloupe = None
        martian.throw("Could not generate .cloupe file: \n%s" % e.output)