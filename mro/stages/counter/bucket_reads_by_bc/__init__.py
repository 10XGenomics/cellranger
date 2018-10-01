#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Put BAM file into buckets by (gem group, barcode prefix)
#
from collections import OrderedDict
import itertools
import martian
import math
import os
import tenkit.bam as tk_bam
import tenkit.stats as tk_stats
import cellranger.h5_constants as h5_constants
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils

__MRO__ = """
stage BUCKET_BY_BC(
    in  int   nbases,
    in  bam[] inputs,
    in  int[] num_alignments,
    out map   buckets,
    src py    "stages/counter/bucket_reads_by_bc",
) split using (
    in  bam   chunk_input,
    in  map[] read_groups,
)
"""

# empirically, how much memory does a set of alignments use relative to the BAM size?
# the average number is closer to 7.5, but say 10 to be safe
PYSAM_TO_BAM_COMPRESSION_RATIO = 10.0

# Empirical observation + 1.25x safety factor
BAM_ALIGNMENTS_PER_MEM_GB = 750000

def get_mem_gb_request_from_bam(bam_path):
    bytes_on_disk = os.path.getsize(bam_path)
    bytes_in_ram = round(PYSAM_TO_BAM_COMPRESSION_RATIO * bytes_on_disk)
    return max(bytes_in_ram / 1e9, h5_constants.MIN_MEM_GB)

def get_mem_gb_request_from_num_alignments(num_alignments):
    mem_gb = tk_stats.robust_divide(num_alignments, BAM_ALIGNMENTS_PER_MEM_GB)
    return max(int(math.ceil(mem_gb)), h5_constants.MIN_MEM_GB)

def get_read_group_union(in_filenames):
    """ Get union of read group (@RG) headers from a list of BAM filenames """
    if len(in_filenames) == 0:
        return []

    # Get all read groups (collapsed by ID)
    read_groups = OrderedDict()
    for bam_fn in in_filenames:
        bam = tk_bam.create_bam_infile(bam_fn)
        for rg in bam.header['RG']:
            if rg['ID'] not in read_groups:
                read_groups[rg['ID']] = rg
        bam.close()

    return read_groups.values()

def split(args):
    chunks = []

    assert args.num_alignments is None or len(args.inputs) == len(args.num_alignments)

    read_groups = get_read_group_union(args.inputs)

    for chunk_input, num_alignments in itertools.izip_longest(args.inputs, args.num_alignments or []):
        # Base it on compression ratios if we don't have alignment counts
        if num_alignments is None:
            mem_gb = get_mem_gb_request_from_bam(chunk_input)
        else:
            mem_gb = get_mem_gb_request_from_num_alignments(num_alignments)

        chunks.append({
            'chunk_input': chunk_input,
            'read_groups': read_groups,
            '__mem_gb': mem_gb,
            # HACK: Prevent memory oversubscription on poorly configured clusters
            '__threads': 2,
        })

    return {'chunks': chunks}

def main(args, outs):
    bam_in = tk_bam.create_bam_infile(args.chunk_input)

    # Get gem groups
    library_info = rna_library.get_bam_library_info(bam_in)
    gem_groups = sorted(list(set(lib['gem_group'] for lib in library_info)))

    # Define buckets
    bucket_names = []
    prefixes = cr_utils.get_seqs(args.nbases)
    for gg in gem_groups:
        for prefix in prefixes:
            bucket_names.append('%s-%d' % (prefix, gg))
    bucket_names.append('')

    # Read all records
    reads = [read for read in bam_in]

    # Bucket the records
    bams_out = {}
    outs.buckets = {}
    buckets = {}
    for bucket_name in bucket_names:
        filename = martian.make_path("bc-%s.bam" % bucket_name)
        bam_out, _ = tk_bam.create_bam_outfile(filename, None, None, template=bam_in, rgs=args.read_groups, replace_rg=True)

        bams_out[bucket_name] = bam_out
        outs.buckets[bucket_name] = filename
        buckets[bucket_name] = []

    for r in reads:
        barcode = cr_utils.get_read_barcode(r)
        if barcode is None:
            bucket_name = ''
        else:
            barcode_seq, gem_group = cr_utils.split_barcode_seq(barcode)
            prefix = barcode_seq[:args.nbases]
            bucket_name = '%s-%d' % (prefix, gem_group)
        buckets[bucket_name].append(r)

    for bucket_name, bucket in buckets.iteritems():
        bucket.sort(key=cr_utils.barcode_sort_key)
        bam_out = bams_out[bucket_name]
        for r in bucket:
            bam_out.write(r)
        bam_out.close()

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.buckets = {}
    for out in chunk_outs:
        for bucket_name, filename in out.buckets.iteritems():
            if bucket_name not in outs.buckets:
                outs.buckets[bucket_name] = []
            outs.buckets[bucket_name].append(filename)
