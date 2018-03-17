#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc.  All rights reserved.
#
# Compute the run-specific demux parameters required to guide
# downstream operation.
#
import os
import martian
import tenkit.bcl as tk_bcl
from tenkit.constants import DEMULTIPLEX_BARCODE_LENGTH, DEMULTIPLEX_DEFAULT_SAMPLE_INDEX_LENGTH

__MRO__ = """
stage COMPUTE_DEMUX_PARAMS(
    in  path   run_path,
    out bool   rc_i2_read,
    out string si_read_type,
    out bool   split_by_tile,
    src py     "stages/bcl_processor/compute_demux_params",
)
"""
def main(args, outs):
    """
    run_path must be the top-level Illumina flowcell directory
    """
    if not os.path.exists(args.run_path):
        martian.throw("Run directory does not exist: %s" % args.run_path)
    
    run_info_xml = os.path.join(args.run_path, "RunInfo.xml")
    read_info, flowcell = tk_bcl.load_run_info(run_info_xml)
    outs.si_read_type = get_si_read_type(read_info)

    (rta_version, rc_i2_read, bcl_params) = tk_bcl.get_rta_version(args.run_path)
    martian.log_info("BCL folder RTA Version: %s" % rta_version)
    martian.log_info("BCL params: %s" % str(bcl_params))
    martian.log_info("RC'ing i2 read: %s" % str(rc_i2_read))
    outs.rc_i2_read = rc_i2_read
    
    split_by_tile = _split_by_tile(args)
    martian.log_info("Splitting by tile: %s" % str(split_by_tile))
    outs.split_by_tile = split_by_tile


def _split_by_tile(args):
    """
    Whether to split bcl2fastq into a number of stages by tile.
    """
    sequencer_type = tk_bcl.get_sequencer_type(args.run_path)
    return sequencer_type in ("NovaSeq",)


def get_si_read_type(read_info):
    si_read_type = None
    # First search for an index read of the default length
    for read in read_info:
        if read["index_read"] and read["read_length"] == DEMULTIPLEX_DEFAULT_SAMPLE_INDEX_LENGTH:
            si_read_type = read["read_name"]
            break
    # Fall back to searching for an index read that is not the old 10x barcode length
    if si_read_type is None:
        for read in read_info:
            if read["index_read"] and read["read_length"] != DEMULTIPLEX_BARCODE_LENGTH:
                si_read_type = read["read_name"]
                break

    if si_read_type:
        martian.log_info("SI read type: %s" % si_read_type)
    else:
        martian.log_info("No SI read detected")

    return si_read_type

