#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import cStringIO
import os
import pysam
import re
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import cellranger.constants as cr_constants
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils
import cellranger.io as cr_io

__MRO__ = '''
stage ALIGN_READS(
    in  fastq[]  reads,
    in  fastq[]  read2s,
    in  string[] read_groups,
    in  string[] library_types,
    in  path     reference_path,
    in  int      threads,
    in  int      max_hits_per_read,
    out bam[]    genome_output,
    src py       "stages/counter/align_reads",
) split using (
    in  fastq    read_chunk,
    in  fastq    read2_chunk,
    in  string   read_group,
    in  string   library_type,
)
'''

def split(args):
    chunks = []

    star_mem_gb = cr_utils.get_reference_mem_gb_request(args.reference_path)
    assert len(args.read2s) == 0 or (len(args.reads) == len(args.read2s))
    for read_chunk, read2_chunk, read_group, library_type in \
        itertools.izip_longest(args.reads, args.read2s, args.read_groups, args.library_types):

        star_align = library_type in cr_constants.ALIGN_LIBRARY_TYPES

        chunks.append({
            'read_chunk': read_chunk,
            'read2_chunk': read2_chunk,
            'read_group': read_group,
            'library_type': library_type,
            '__mem_gb': star_mem_gb if star_align else 1,
            '__threads': args.threads if star_align else 1,
        })
    return {'chunks': chunks}

def join(args, outs, chunk_defs, chunk_outs):
    outs.genome_output = [chunk_out.genome_output for chunk_out in chunk_outs]

''' Align reads using STAR '''
def align_reads(args, outs):
    star_ref_path = cr_utils.get_reference_star_path(args.reference_path)
    star = cr_reference.STAR(star_ref_path)

    star.align(args.read_chunk, args.read2_chunk, outs.genome_output,
               max_report_alignments_per_read=args.max_hits_per_read,
               threads=args.threads,
               read_group_tags=tk_bam.make_star_rg_header(args.read_group))

''' Write reads as unaligned BAM '''
def create_unaligned_bam(args, outs):
    star_ref_path = cr_utils.get_reference_star_path(args.reference_path)

    header_buf = cStringIO.StringIO()

    header_buf.write('@HD\tVN:1.4\n')

    # SQ header lines
    with open(os.path.join(star_ref_path, 'chrNameLength.txt')) as f:
        for line in f:
            chr_name, chr_len = line.strip().split('\t')
            header_buf.write('@SQ\tSN:{}\tLN:{}\n'.format(chr_name, chr_len))

    # RG header lines
    for packed_rg in args.read_groups:
        header_buf.write(re.sub('\\\\t', '\t', tk_bam.make_rg_header(packed_rg)) + '\n')

    # Get read group ID for this chunk of reads
    read_group = args.read_group

    # pysam doesn't support reading SAM from a StringIO object
    with open('tmphdr', 'w') as f:
        f.write(header_buf.getvalue())
    samfile = pysam.AlignmentFile('tmphdr', 'r', check_sq=False)

    outbam = pysam.AlignmentFile(outs.genome_output, 'wb', template=samfile)

    fastq_file1 = cr_io.open_maybe_gzip(args.read_chunk)
    fastq_file2 = cr_io.open_maybe_gzip(args.read2_chunk) if args.read2_chunk else None
    read1s = tk_fasta.read_generator_fastq(fastq_file1)
    read2s = tk_fasta.read_generator_fastq(fastq_file2) if fastq_file2 else []

    record = pysam.AlignedSegment()
    record.flag = 4

    for read1, read2 in itertools.izip_longest(read1s, read2s):
        name, seq, qual = read1
        record.query_name, record.query_sequence = name.split(' ')[0], seq
        record.query_qualities = tk_fasta.get_qvs(qual)
        record.set_tag('RG', read_group, 'Z')
        outbam.write(record)

        if read2:
            name, seq, qual = read2
            record.query_name, record.query_sequence = name.split(' ')[0], seq
            record.query_qualities = tk_fasta.get_qvs(qual)
            record.set_tag('RG', read_group, 'Z')
            outbam.write(record)

    samfile.close()
    fastq_file1.close()
    if fastq_file2 is not None:
        fastq_file2.close()
    outbam.close()


def main(args, outs):
    if args.library_type in cr_constants.ALIGN_LIBRARY_TYPES:
        align_reads(args, outs)
    else:
        create_unaligned_bam(args, outs)
