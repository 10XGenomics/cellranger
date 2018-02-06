#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# 1) Attach corrected barcodes
# 2) Discard reads that are too short
# 3) Bucket augmented FASTQs by gem group and then by barcode.
# 4) Chunk the reads.
#
import itertools
import json
import martian
import math
import tenkit.fasta as tk_fasta
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.utils as cr_utils
import cellranger.vdj.utils as vdj_utils

__MRO__ = """
stage BUCKET_FASTQ_BY_BC(
    in  fastq[] read1s,
    in  fastq[] read2s,
    in  tsv[]   corrected_bcs,
    in  int[]   gem_groups,
    in  json    reads_summary,
    in  int     readpairs_per_chunk,
    in  map     chemistry_def,
    out map[]   buckets,
    src py      "stages/vdj/bucket_fastq_by_bc",
) split using (
    in  fastq   read1s_chunk,
    in  fastq   read2s_chunk,
    in  fastq   bcs,
    in  map     chunks_per_gem_group,
)
"""

# Discard read pairs where either read is shorter than this. If we don't, bad things will happen in vdj_asm
# when reads fall below the kmer length. Namely reads will go blank.
MIN_READ_LENGTH = 50

def split(args):
    paired_end = cr_chem.is_paired_end(args.chemistry_def)

    if paired_end:
        assert len(args.read1s) == len(args.read2s)

    assert len(args.corrected_bcs) == len(args.read1s)

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

    for fastq1, fastq2, bcs in itertools.izip_longest(args.read1s, args.read2s, args.corrected_bcs):
        chunks.append({
            'read1s_chunk': fastq1,
            'read2s_chunk': fastq2 if paired_end else None,
            'bcs': bcs,
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

    paired_end = args.read2s_chunk is not None

    # Lazy load R1
    r1_file = cr_utils.open_maybe_gzip(args.read1s_chunk)
    read1s = tk_fasta.read_generator_fastq(r1_file)

    # Lazy load R2
    if paired_end:
        r2_file = cr_utils.open_maybe_gzip(args.read2s_chunk)
        read2s = tk_fasta.read_generator_fastq(r2_file)
    else:
        read2s = []

    # Lazy load corrected BCs
    bc_file = cr_utils.open_maybe_gzip(args.bcs)
    bcs = (line.strip() for line in bc_file)

    buckets = {}

    bucket_filenames = {}

    for gem_group, bucket_name in enumerate_bucket_names(args.chunks_per_gem_group):
        filename = martian.make_path("%s.fastq" % bucket_name)
        bucket_filenames[bucket_name] = filename
        buckets[bucket_name] = []

    for read1, read2, barcode in itertools.izip_longest(read1s, read2s, bcs):
        # Exclude unbarcoded reads
        if barcode == '':
            continue

        # Exclude short reads
        if len(read1[1]) < MIN_READ_LENGTH or (read2 is not None and len(read2[1]) < MIN_READ_LENGTH):
            continue

        # Attach processed barcode to reads
        r1_hdr = cr_fastq.AugmentedFastqHeader(read1[0])
        r1_hdr.set_tag(cr_constants.PROCESSED_BARCODE_TAG, barcode)
        r1_new_qname = r1_hdr.to_string()

        if paired_end:
            r2_hdr = cr_fastq.AugmentedFastqHeader(read2[0])
            r2_hdr.set_tag(cr_constants.PROCESSED_BARCODE_TAG, barcode)
            r2_new_qname = r2_hdr.to_string()

        barcode_seq, gem_group = cr_utils.split_barcode_seq(barcode)
        bucket_name = get_bucket_name(gem_group, barcode_seq, args.chunks_per_gem_group[gem_group])

        buckets[bucket_name].append((r1_new_qname, read1[1], read1[2]))
        if paired_end:
            buckets[bucket_name].append((r2_new_qname, read2[1], read2[2]))

    outs.buckets = {}

    # Sort and write each bucket
    for bucket_name, bucket in buckets.iteritems():
        bucket.sort(key=vdj_utils.fastq_barcode_sort_key)

        # Don't create empty bucket files.
        # This is common when the reads are ordered by gem group
        # And a chunk sees only a single gem group.
        if len(bucket) == 0:
            continue

        filename = bucket_filenames[bucket_name]
        with cr_utils.open_maybe_gzip(filename, 'w') as f:
            for read in bucket:
                tk_fasta.write_read_fastq(f, *read)

        outs.buckets[bucket_name] = bucket_filenames[bucket_name]

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    buckets = {}

    # Merge bucket FASTQ filenames across chunks
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
