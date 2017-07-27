#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import tenkit.bam as tk_bam
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils

__MRO__ = '''
stage ALIGN_READS(
    in  fastq[]  reads,
    in  fastq[]  read2s,
    in  string[] read_groups,
    in  path     reference_path,
    in  int      threads,
    in  int      max_hits_per_read,
    out bam[]    genome_output,
    out bool     paired_end,
    src py       "stages/counter/align_reads",
) split using (
    in  fastq    read_chunk,
    in  fastq    read2_chunk,
    in  string[] read_group,
)
'''

def split(args):
    chunks = []

    mem_gb = cr_utils.get_reference_mem_gb_request(args.reference_path)
    assert len(args.read2s) == 0 or (len(args.reads) == len(args.read2s))
    for read_chunk, read2_chunk, read_group in itertools.izip_longest(args.reads, args.read2s, args.read_groups):
        chunks.append({
            'read_chunk': read_chunk,
            'read2_chunk': read2_chunk,
            'read_group': read_group,
            '__mem_gb': mem_gb,
            '__threads': args.threads,
        })
    return {'chunks': chunks}

def join(args, outs, chunk_defs, chunk_outs):
    outs.genome_output = [chunk_out.genome_output for chunk_out in chunk_outs]
    outs.paired_end = len(args.read2s) > 0

def main(args, outs):
    reference_star_path = cr_utils.get_reference_star_path(args.reference_path)
    star = cr_reference.STAR(reference_star_path)

    star.align(args.read_chunk, args.read2_chunk, outs.genome_output,
               max_report_alignments_per_read=args.max_hits_per_read,
               threads=args.threads,
               read_group_tags=tk_bam.make_star_rg_header(args.read_group))
