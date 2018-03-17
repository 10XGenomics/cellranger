#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc.  All rights reserved.
#
# If BCL files were processed and demultiplex in parallel by tile, reassemble
# them (in tile suffix order).
#
import os
import glob
import martian
import shutil
import subprocess
from tenkit.fasta import BclProcessorFastqFile
from tenkit.constants import DEMULTIPLEX_INVALID_SAMPLE_INDEX
import tenkit.log_subprocess as tk_proc

MRO = """
stage MERGE_FASTQS_FROM_TILES(
    in  path     demultiplexed_fastq_path,
    in  string[] common_bcs,
    in  bool     split_by_tile,
    out path     demultiplexed_fastq_path
) split using (
    in int      lane,
    in string[] bcs,
)
"""

BARCODE_GROUP_SIZE = 16
CHUNK_MEM_GB = 2

def split(args):
    # if the files have not been split by tile, we're done, just bail
    if not args.split_by_tile:
        return {'chunks': [{'lane': None, 'bcs': []}]}

    # from here forward, assume that we're dealing with FASTQs separated
    # by tile
    file_glob = os.path.join(args.demultiplexed_fastq_path, "Tile*", "*.fastq*")
    files = glob.glob(file_glob)

    if len(files) == 0:
        martian.throw("No FASTQ files were found.")

    # find the unique # of lanes there are in all
    file_info = [ BclProcessorFastqFile(x) for x in files ]
    lanes = sorted(set([fi.lane for fi in file_info]))

    # lexicographically sort barcodes (incoming order is in reverse frequency
    # order) in order to spread work around more evenly
    sorted_bcs = sorted(args.common_bcs)
    bc_groups, bc_remgroup = divmod(len(sorted_bcs), BARCODE_GROUP_SIZE)

    chunks = []
    for group_index in range(bc_groups):
        bcs = sorted_bcs[group_index*BARCODE_GROUP_SIZE:(group_index+1)*BARCODE_GROUP_SIZE]
        for lane in lanes:
            chunks.append({'__mem_gb': CHUNK_MEM_GB, 'lane': lane, 'bcs': bcs})
    if bc_remgroup > 0:
        bcs = sorted_bcs[bc_groups*BARCODE_GROUP_SIZE:]
        for lane in lanes:
            chunks.append({'__mem_gb': CHUNK_MEM_GB, 'lane': lane, 'bcs': bcs})

    # finally, the leftovers (si_X)
    for lane in lanes:
        chunks.append({'__mem_gb': CHUNK_MEM_GB, 'lane': lane, 'bcs': [DEMULTIPLEX_INVALID_SAMPLE_INDEX]})

    return {'chunks': chunks}


def main(args, outs):
    if not args.split_by_tile:
        return

    os.makedirs(outs.demultiplexed_fastq_path)

    demux_read_types = ("RA", "I1", "I2")  # covering the bases
    # like tenkit.fasta.find_input_fastq_files_10x_preprocess but allow Ns
    # from combined barcode list
    for read_type in demux_read_types:
        for barcode in args.bcs:
            file_glob = "read-%s_si-%s_lane-%03d[_\-]*.fastq*" % (read_type, barcode, args.lane)
            dir_glob = os.path.join(args.demultiplexed_fastq_path, "Tile*", file_glob)
            files = glob.glob(dir_glob)
            # assuming here that all files are already gzipped
            out_path = os.path.join(outs.demultiplexed_fastq_path,
                                    "read-%s_si-%s_lane-%03d-chunk-001.fastq.gz" % (read_type, barcode, args.lane))
            if files:
                subprocess_args = ["cat"] + files + [">", out_path]
                tk_proc.check_call(" ".join(subprocess_args), shell=True)


def join(args, outs, chunk_defs, chunk_outs):
    os.makedirs(outs.demultiplexed_fastq_path)
    # if not originally split by tile, just move files from
    # the demultiplex stage into this stage
    if not args.split_by_tile:
        for f in os.listdir(args.demultiplexed_fastq_path):
            in_file = os.path.join(args.demultiplexed_fastq_path, f)
            shutil.move(in_file, outs.demultiplexed_fastq_path)
        return

    for chunk_out in chunk_outs:
        for f in os.listdir(chunk_out.demultiplexed_fastq_path):
            in_file = os.path.join(chunk_out.demultiplexed_fastq_path, f)
            shutil.move(in_file, outs.demultiplexed_fastq_path)
