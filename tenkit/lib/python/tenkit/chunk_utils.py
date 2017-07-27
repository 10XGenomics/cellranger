#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import tenkit.bio_io as tk_io
from tenkit.constants import PARALLEL_LOCUS_SIZE

def pack_loci(loci):
    packed_loci = []

    current_group = []
    current_size = 0
    for locus in loci:
        (chrom, start, end) = tk_io.get_locus_info(locus)
        current_group.append(locus)
        current_size += end - start

        if current_size >= 0.25 * PARALLEL_LOCUS_SIZE:
            packed_loci.append(current_group)
            current_group = []
            current_size = 0

    if len(current_group) > 0:
        packed_loci.append(current_group)

    return packed_loci

# Choose chromosome slices so that boundaries don't split target regions
def adjust_start(chrom, start, regions):
    if regions is None or chrom not in regions:
        return start
    r = regions[chrom].get_region_containing_point(start)
    if r is None or r[0] == start:
        return start
    else:
        return r[1]

"""
Bins the chromosome into windows of (up to) the given chunk size.
If target regions are given, then the window starts never split a region.
An overlap can be added between the chunks. This might be useful if you
want to mitigate boundary effects.
"""

def get_parallel_locus_size(targets_file):
    if targets_file is None:
        return PARALLEL_LOCUS_SIZE
    else:
        return PARALLEL_LOCUS_SIZE * 3

def generate_chrom_loci(target_regions, chrom, chrom_length, chunk_size, overlap = 0):
    starts = [adjust_start(chrom, s, target_regions) for s in range(0, chrom_length, chunk_size)]
    ends = [min(chrom_length, s + overlap) for s in starts[1:]] + [chrom_length]
    chunks = [tk_io.create_locus_info(chrom, s, e) for (s,e) in zip(starts, ends)]
    return chunks

def chunk_by_locus(chroms, chrom_lengths, chunk_size, overlap = 0, contig_whitelist = None, target_regions = None, extra_args = {}):
    chunks = []
    for (chrom, length) in zip(chroms, chrom_lengths):
        if contig_whitelist is None or chrom in contig_whitelist:
            chrom_loci = generate_chrom_loci(target_regions, chrom, length, chunk_size, overlap = overlap)
            for locus in chrom_loci:
                chunk = {'locus': locus}
                chunk.update(extra_args)
                chunks.append(chunk)
    return chunks
