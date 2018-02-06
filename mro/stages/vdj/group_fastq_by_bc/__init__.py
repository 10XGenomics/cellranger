#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Group FASTQ file by concatenating buckets into chunks.
#
import heapq
import json
import martian
import tenkit.fasta as tk_fasta
import tenkit.cache as tk_cache
import tenkit.safe_json as tk_safe_json
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.utils as vdj_utils

__MRO__ = """
stage GROUP_FASTQ_BY_BC(
    in  map[]   buckets,
    in map chemistry_def,
    out fastq[] read1s,
    out fastq[] read2s,
    out int[]   chunk_gem_groups,
    out json[]  chunk_barcodes,
    src py      "stages/vdj/group_fastq_by_bc",
) split using (
    in  string  bucket_name,
    in  int     gem_group,
    in  fastq[] fastqs,
)
"""

def split(args):
    return {'chunks': args.buckets}

def main(args, outs):
    outs.coerce_strings()

    paired_end = cr_chem.is_paired_end(args.chemistry_def)

    outs.read1s = martian.make_path('reads_1.fastq' + cr_constants.LZ4_SUFFIX)
    r1_fq_out = cr_utils.open_maybe_gzip(outs.read1s, 'w')

    if paired_end:
        outs.read2s = martian.make_path('reads_2.fastq' + cr_constants.LZ4_SUFFIX)
        r2_fq_out = cr_utils.open_maybe_gzip(outs.read2s, 'w')
    else:
        outs.read2s = None
        r2_fq_out = None

    barcodes_out = cr_utils.open_maybe_gzip(outs.chunk_barcodes, 'w')

    merge_by_barcode(args.fastqs, r1_fq_out, r2_fq_out,
                     barcodes_out, paired_end)

    r1_fq_out.close()
    if r2_fq_out is not None:
        r2_fq_out.close()
    barcodes_out.close()

def merge_by_barcode(in_filenames, r1_out_file, r2_out_file, bcs_out_file, paired_end):
    barcodes = set()

    # Note: The filehandle cache precludes the use of compressed files
    file_cache = tk_cache.FileHandleCache(mode='r', open_func=open)
    heap = []

    key_func = vdj_utils.fastq_barcode_sort_key

    for filename in in_filenames:
        try:
            fastq = tk_fasta.read_generator_fastq(file_cache.get(filename), paired_end=paired_end)
            first_readpair = fastq.next()

            key = key_func(first_readpair[0:3])
            barcode = key[0]
            barcodes.add(barcode)

            heapq.heappush(heap, (key, first_readpair, filename))

        except StopIteration:
            pass

    while len(heap) > 0:
        # Get the minimum item and write it.
        key, readpair, in_filename = heapq.heappop(heap)

        fastq = tk_fasta.read_generator_fastq(file_cache.get(in_filename), paired_end=paired_end)

        tk_fasta.write_read_fastq(r1_out_file, *readpair[0:3])
        if paired_end:
            tk_fasta.write_read_fastq(r2_out_file, *readpair[3:6])

        # Get the next item from the source file we just wrote from
        # If that file is out of items, then we leave that one out
        try:
            next_readpair = fastq.next()

            key = key_func(next_readpair[0:3])
            barcode = key[0]
            barcodes.add(barcode)

            heapq.heappush(heap, (key, next_readpair, in_filename))

        except StopIteration:
            pass

    json.dump(tk_safe_json.json_sanitize(list(barcodes)), bcs_out_file)


def join(args, outs, chunk_defs, chunk_outs):
    outs.read1s = [out.read1s for out in chunk_outs]
    outs.read2s = [out.read2s for out in chunk_outs]
    outs.chunk_gem_groups = [chunk.gem_group for chunk in chunk_defs]
    outs.chunk_barcodes = [out.chunk_barcodes for out in chunk_outs]
