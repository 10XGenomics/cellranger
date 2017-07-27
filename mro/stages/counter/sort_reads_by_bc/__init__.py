#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file by sorting buckets, then concatenating bucket files
#
import heapq
import pysam
import tenkit.bam as tk_bam
import tenkit.cache as tk_cache
import cellranger.utils as cr_utils

__MRO__ = """
stage SORT_BY_BC(
    in  map    buckets,
    out int    total_reads,
    out bam,
    src py     "stages/counter/sort_reads_by_bc",
) split using (
    in  string prefix,
    in  bam[]  bucket,
)
"""

def split(args):
    chunks = []
    for prefix, bucket in args.buckets.iteritems():
        chunks.append({
            'prefix': prefix,
            'bucket': bucket,
        })
    return {'chunks': chunks}

def main(args, outs):
    outs.coerce_strings()
    bam_in = tk_bam.create_bam_infile(args.bucket[0])
    bam_out, _ = tk_bam.create_bam_outfile(outs.default, None, None, template=bam_in)
    bam_in.close()

    outs.total_reads = merge_by_key(args.bucket, bc_and_qname_sort_key, bam_out)
    bam_out.close()

def bc_and_qname_sort_key(read):
    # Maintain qname ordering within each BC
    bc_key = cr_utils.barcode_sort_key(read)
    return (bc_key, read.qname)

def merge_by_key(bam_filenames, key_func, bam_out):
    file_cache = tk_cache.FileHandleCache(mode='rb', open_func=pysam.Samfile)
    total_reads = 0
    heap  = []

    for bam_filename in bam_filenames:
        try:
            bam = file_cache.get(bam_filename)
            first_read = bam.next()
            heapq.heappush(heap, (key_func(first_read), first_read, bam_filename))
        except StopIteration:
            pass

    while len(heap) > 0:
        # Get the minimum item and write it to the bam.
        key, read, bam_filename = heapq.heappop(heap)
        bam = file_cache.get(bam_filename)
        bam_out.write(read)
        total_reads += 1

        # Get the next read from the source bam we just wrote from
        # If that BAM file is out of reads, then we leave that one out
        try:
            next_read = bam.next()
            heapq.heappush(heap, (key_func(next_read), next_read, bam_filename))
        except StopIteration:
            pass

    return total_reads

def join(args, outs, chunk_defs, chunk_outs):
    chunks = zip(chunk_defs, chunk_outs)
    chunks.sort(key=lambda chunk: chunk[0].prefix)

    buckets = []
    outs.total_reads = 0
    for chunk in chunks:
        buckets.append(chunk[1].default)
        outs.total_reads += chunk[1].total_reads

    tk_bam.concatenate(outs.default, buckets)
