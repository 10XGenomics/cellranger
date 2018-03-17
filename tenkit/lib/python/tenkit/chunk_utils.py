#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import os.path
import math
import pysam
from itertools import groupby
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


def get_sized_bam_chunks(bam_fn, gb_per_chunk, contig_whitelist=None, target_regions=None, extra_args = {}):
    '''Divide a BAM file into disjoint loci with a max compressed size of roughly gb_per_chunk. If contig_whitelist is supplied,
       those contigs will not be included. If target_regions is supplied, boundaries will be adjusted to avoid on-target regions'''
    total_size = os.path.getsize(bam_fn)
    bam = pysam.Samfile(bam_fn)

    file_starts = []
    last_pos = 0
    for chrom in bam.references:
        offset = get_voffset(bam, chrom, 0)
        if offset is None:
            offset = last_pos
        else:
            last_pos = offset

        file_starts.append(offset)

    file_sizes = []
    for i in range(len(file_starts) - 1):
        file_sizes.append(file_starts[i+1] - file_starts[i])
    file_sizes.append(total_size - file_starts[-1])

    loci = []
    for (chrom, file_start, file_size, chrom_size) in zip(bam.references, file_starts, file_sizes, bam.lengths):
        if contig_whitelist is None or chrom in contig_whitelist:
        
            n_chunks = max(1, int(math.ceil(float(file_size) / 1e9 / gb_per_chunk)))
            chunk_size = int(file_size / n_chunks)
            chunk_starts = []#

            for i in range(n_chunks):
                if i == 0:
                    pos = 0
                else:
                    voffset = file_start + chunk_size * i
                    _pos = find_pos_of_voffset(bam, chrom, voffset, err=chunk_size/20)
                    pos = adjust_start(chrom, _pos, target_regions)

                # Don't create very small chunks
                if len(chunk_starts) > 0 and pos - chunk_starts[-1] < 100:
                    continue

                chunk_starts.append(pos)

            for i in range(len(chunk_starts)):
                if i < len(chunk_starts) - 1:
                    locus = tk_io.create_locus_info(chrom, chunk_starts[i], chunk_starts[i+1])
                else:
                    locus = tk_io.create_locus_info(chrom, chunk_starts[i], chrom_size)

                loci.append(locus)

    validate_loci(bam, loci, contig_whitelist)
    chunks = []
    for l in loci:
        chunk = {'locus': l}
        chunk.update(extra_args)
        chunks.append(chunk)

    return chunks


def validate_loci(bam, loci, whitelist):

    loci = sorted([tk_io.get_locus_info(l) for l in loci])

    good_chroms = []
    for (chrom, items) in groupby(loci, lambda x: x[0]):
        sorted_items = sorted(items, key=lambda x:x[1])
        last_end = 0
        for (_,s,e) in sorted_items:
            assert(e-s > 0)
            assert(s == last_end)
            last_end = e

        assert(last_end == chrom_size(bam, chrom))
        good_chroms.append(chrom)

    if whitelist == None:
        assert(set(good_chroms) == set(bam.references))
    else:
        assert(set(good_chroms) == set(whitelist))


def chrom_size(bam, chrom):
    sz = [sz for (c, sz) in zip(bam.references, bam.lengths) if c == chrom][0]
    return sz


def get_voffset(bam, chrom, pos):
    itr = bam.fetch(chrom, pos, chrom_size(bam, chrom))

    
    try:
        first_read = itr.next()
    except StopIteration:
        first_read = None

    if first_read is not None:
        return bam.tell() >> 16
    else:
        return None

def find_pos_of_voffset(bam, chrom, voffset, err=1e6):

    start = 0
    end = chrom_size(bam, chrom)
    probe_offset = -1e12

    while probe_offset is None or abs(voffset - probe_offset) > err:
        probe_pos = start / 2 + end / 2
        probe_offset = get_voffset(bam, chrom, probe_pos)

        if probe_offset is None or probe_offset > voffset:
            end = probe_pos
        else:
            start = probe_pos

        if end - start < 1000:
            break

    return probe_pos
