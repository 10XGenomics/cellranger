#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# This stage only exists to force earlier VDR of prior stages.

import itertools
import cellranger.chemistry as cr_chem
import cellranger.utils as cr_utils

__MRO__ = """
stage SUMMARIZE_READ_REPORTS(
    in  json     extract_reads_summary,
    in  json     barcode_counts,
    in  int[]    gem_groups,
    in  string[] read_groups,
    in  map      align,
    in  string[] bam_comments,
    in  fastq[]  read1s,
    in  fastq[]  read2s,
    in  bool     retain_fastqs,
    in  map      chemistry_def,
    out json     summary,
    out json     barcode_counts,
    out h5       barcode_summary,
    out int[]    gem_groups,
    out string[] read_groups,
    out map      align,
    out string[] bam_comments,
    out fastq[]  read1s,
    out fastq[]  read2s,
    src py       "stages/counter/summarize_read_reports",
) split using (
    in  fastq    read1,
    in  fastq    read2,
)
"""

def split(args):
    paired_end = cr_chem.is_paired_end(args.chemistry_def)
    if paired_end:
        assert len(args.read1s) == len(args.read2s)

    chunks = []

    if args.retain_fastqs:
        for read1, read2 in itertools.izip_longest(args.read1s, args.read2s):
            chunks.append({
                'read1': read1,
                'read2': read2 if paired_end else None,
            })

    return {'chunks': chunks}

def main(args, outs):
    if args.read1 is not None:
        # Ensure same extension
        out_path, _ = cr_utils.splitexts(outs.read1s)
        _, in_ext = cr_utils.splitexts(args.read1)
        outs.read1s = out_path + in_ext
        cr_utils.copy(args.read1, outs.read1s)
    if args.read2 is not None:
        out_path, _ = cr_utils.splitexts(outs.read2s)
        _, in_ext = cr_utils.splitexts(args.read2)
        outs.read2s = out_path + in_ext
        cr_utils.copy(args.read2, outs.read2s)

def join(args, outs, chunk_defs, chunk_outs):
    cr_utils.copy(args.extract_reads_summary, outs.summary)
    cr_utils.copy(args.barcode_counts, outs.barcode_counts)

    outs.gem_groups = args.gem_groups
    outs.read_groups = args.read_groups
    outs.align = args.align
    outs.bam_comments = args.bam_comments

    outs.read1s = [co.read1s for co in chunk_outs]
    outs.read2s = [co.read2s for co in chunk_outs]
