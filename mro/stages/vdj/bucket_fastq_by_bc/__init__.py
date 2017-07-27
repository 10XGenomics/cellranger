#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# 1) Bucket augmented FASTQs by gem group and then by barcode.
# 2) Chunk the reads.
#
import itertools
import json
import martian
import math
import tenkit.fasta as tk_fasta
import cellranger.utils as cr_utils
import cellranger.vdj.utils as vdj_utils

__MRO__ = """
stage BUCKET_FASTQ_BY_BC(
    in  int     nbases,
    in  fastq[] read1s,
    in  fastq[] read2s,
    in  int[]   gem_groups,
    in  json    reads_summary,
    in  int     readpairs_per_chunk,
    out map[]   buckets,
    src py      "stages/vdj/bucket_fastq_by_bc",
) split using (
    in  fastq   read1s_chunk,
    in  fastq   read2s_chunk,
    in  map     chunks_per_gem_group,
)
"""

def split(args):
    assert len(args.read1s) == len(args.read2s)

    chunks = []

    # Determine the number of buckets required to achieve
    # the given chunk size.
    chunks_per_gem_group = {}
    with open(args.reads_summary) as f:
        reads_summary = json.load(f)
        for gg in args.gem_groups:
            readpairs = reads_summary['%d_total_reads_per_gem_group' % gg]
            chunks_per_gem_group[str(gg)] = max(2,
                                                int(math.ceil(float(readpairs) / \
                                                              args.readpairs_per_chunk)))

    for fastq1, fastq2 in itertools.izip(args.read1s, args.read2s):
        chunks.append({
            'read1s_chunk': fastq1,
            'read2s_chunk': fastq2,
            'chunks_per_gem_group': chunks_per_gem_group,
        })
    return {'chunks': chunks}

def enumerate_bucket_names(chunks_per_gem_group):
    """ yield (gem_group, bucket_name) """
    for gem_group, num_chunks in chunks_per_gem_group.iteritems():
        for chunk_idx in xrange(num_chunks):
            yield gem_group, str(gem_group) + '-' + str(chunk_idx)

def get_bucket_name(gem_group, barcode, num_buckets):
    """ gem_group - integer
        barcode - barcode sequence
        num_buckets - desired number of buckets for gem group """

    # NOTE: Python modulo returns non-negative numbers here, which we want
    return str(gem_group) + '-' + str(hash(barcode) % num_buckets)

def parse_bucket_name(bucket_name):
    """ Returns (int,str) - (gem_group, hash_str) """
    gem_group_str, hash_str = bucket_name.split('-')
    return int(gem_group_str), hash_str


def main(args, outs):
    # Martian coerces dict keys to string
    # Coerce keys back to int
    args.chunks_per_gem_group = {int(k): v for k, v in args.chunks_per_gem_group.iteritems()}

    with open(args.read1s_chunk) as f1:
        read1s = [read for read in tk_fasta.read_generator_fastq(f1)]

    with open(args.read2s_chunk) as f2:
        read2s = [read for read in tk_fasta.read_generator_fastq(f2)]

    assert len(read1s) == len(read2s)

    fastqs_out = {}
    buckets = {}

    outs.buckets = {}

    for gem_group, bucket_name in enumerate_bucket_names(args.chunks_per_gem_group):
        filename = martian.make_path("%s.fastq" % bucket_name)
        fastqs_out[bucket_name] = open(filename, 'w')
        outs.buckets[bucket_name] = filename
        buckets[bucket_name] = []

    for read1, read2 in itertools.izip(read1s, read2s):
        barcode = vdj_utils.get_fastq_read_barcode(read1)

        # Exclude unbarcoded reads
        if barcode is None:
            continue

        assert barcode == vdj_utils.get_fastq_read_barcode(read2)

        barcode_seq, gem_group = cr_utils.split_barcode_seq(barcode)
        bucket_name = get_bucket_name(gem_group, barcode_seq, args.chunks_per_gem_group[gem_group])

        buckets[bucket_name].append(read1)
        buckets[bucket_name].append(read2)

    # Sort and write each bucket
    for bucket_name, bucket in buckets.iteritems():
        bucket.sort(key=vdj_utils.fastq_barcode_sort_key)

        fastq_out = fastqs_out[bucket_name]
        for read in bucket:
            tk_fasta.write_read_fastq(fastq_out, *read)

        fastq_out.close()


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    buckets = {}

    # Populate buckets from chunks
    for out in chunk_outs:
        for bucket_name, filename in out.buckets.iteritems():
            if bucket_name not in buckets:
                gem_group, _ = parse_bucket_name(bucket_name)

                buckets[bucket_name] = {
                    'bucket_name': bucket_name,
                    'gem_group': gem_group,
                    'fastqs': []
                }

            buckets[bucket_name]['fastqs'].append(filename)

    # Return list of buckets sorted lexicographically by bucket name.
    # This ensures contiguous gem groups.
    outs.buckets = sorted(buckets.itervalues(), key=lambda x: x['bucket_name'])
