#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import json
import os
import martian
import tenkit.log_subprocess as tk_subproc

__MRO__ = """
stage CHUNK_READS(
    in  map[] chunks,
    in  int   reads_per_file,
    out map[] out_chunks,
    src py    "stages/common/chunk_reads",
) split using (
    in  map   read_chunks,
)
"""

def split(args):
    # One chunk per input file set
    chunks = []
    for chunk in args.chunks:
        chunk_def = {}
        chunk_def['read_chunk'] = chunk
        chunk_def['__threads'] = 4
        chunk_def['__mem_gb'] = 1
        chunks.append(chunk_def)

    return {'chunks': chunks, 'join': {'__mem_gb' : 1}}


def join(args, outs, chunk_defs, chunk_outs):
    # Combine chunk JSONs
    all_out_chunks = []
    for out in chunk_outs:
        all_out_chunks.extend(out.out_chunks)

    outs.out_chunks = all_out_chunks


def main(args, outs):

    # Write read_chunk for consumption by Rust
    with open("chunk_args.json", "w") as f:
        json.dump(args.read_chunk, f)

    output_path = martian.make_path("")
    prefix = "fastq_chunk"
    chunk_reads_args = ['chunk_reads',  '--reads-per-fastq', str(args.reads_per_file), output_path, prefix, "--martian-args", "chunk_args.json", '--compress', 'lz4']
    print "running chunk reads: [%s]" % str(chunk_reads_args)
    tk_subproc.check_call(chunk_reads_args)

    with open(os.path.join(output_path, "read_chunks.json")) as f:
        chunk_results = json.load(f)

    outs.out_chunks = []

    # Write out a new chunk entry for each resulting chunk
    for chunk in chunk_results:
        print args.read_chunk
        chunk_copy = args.read_chunk.copy()
        print chunk_copy
        chunk_copy['read_chunks'] = chunk
        outs.out_chunks.append(chunk_copy)
