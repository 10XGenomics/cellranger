#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import resource
import tenkit.bam as tk_bam
import cellranger.utils as cr_utils

__MRO__ = """
stage SORT_BY_POS(
    in  bam[]   inputs,
    in  int     num_threads,
    in  int     mem_gb,
    out bam     output,
    out bam.bai index,
    src py      "stages/counter/sort_reads_by_pos",
) split using (
    in  bam     chunk_input,
)
"""

def split(args):
    chunks = []
    for chunk_input in args.inputs:
        chunks.append({
            'chunk_input': chunk_input,
        })
    join = {
        '__threads': args.num_threads,
        '__mem_gb': args.mem_gb,
    }
    return {'chunks': chunks, 'join': join}

def main(args, outs):
    args.coerce_strings()
    bam_prefix, ext = os.path.splitext(outs.output)
    tk_bam.sort(str(args.chunk_input), str(bam_prefix))

def merge(input_bams, output_bam, threads=1):
    ''' Merge the sorted bam chunks hierarchically to conserve open file handles '''
    soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    soft -= 100

    tmp_dir = os.path.dirname(output_bam)
    while len(input_bams) > 1:
        new_bams = []
        for i in range(0, len(input_bams), soft):
            bam_chunk = input_bams[i:i+soft]
            if len(bam_chunk) > 1:
                new_bam = os.path.join(tmp_dir, "%d-%d.bam" % (i, len(input_bams)))
                tk_bam.merge(new_bam, bam_chunk, threads)
                new_bams.append(new_bam)
            else:
                new_bams.append(input_bams[i])
        input_bams = new_bams
    cr_utils.move(input_bams[0], output_bam)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.output) for chunk in chunk_outs]
    merge(input_bams, outs.output, args.__threads)
    outs.index = outs.output + '.bai'
    tk_bam.index(outs.output)
