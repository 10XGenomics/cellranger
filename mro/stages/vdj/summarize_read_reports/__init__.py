#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# This stage only exists to force earlier VDR of prior stages.

import cellranger.report as cr_report
import cellranger.utils as cr_utils

__MRO__ = """
stage SUMMARIZE_READ_REPORTS(
    in  json     extract_reads_summary,
    in  json     correct_barcodes_summary,
    in  json     raw_barcode_counts,
    in  json     corrected_barcode_counts,
    in  h5       barcode_summary,
    in  int[]    gem_groups,
    in  string[] read_groups,
    in  map      align,
    in  string[] bam_comments,
    in  fastq[]  bc_corrected_read1s,
    in  fastq[]  bc_corrected_read2s,
    in  bool     retain_fastqs,
    out json     summary,
    out json     raw_barcode_counts,
    out json     corrected_barcode_counts,
    out h5       barcode_summary,
    out int[]    gem_groups,
    out string[] read_groups,
    out map      align,
    out string[] bam_comments,
    out fastq[]  bc_corrected_read1s,
    out fastq[]  bc_corrected_read2s,
    src py       "stages/vdj/summarize_read_reports",
) split using (
    in  fastq    read1s,
    in  fastq    read2s,
)
"""

def split(args):
    chunks = []

    if args.retain_fastqs:
        for read1s, read2s in zip(args.bc_corrected_read1s, args.bc_corrected_read2s):
            chunks.append({
                'read1s': read1s,
                'read2s': read2s,
            })
    else:
        chunks = [{'read1s': None, 'read2s': None}]

    return {'chunks': chunks}


def main(args, outs):
    if args.read1s is not None:
        cr_utils.copy(args.read1s, outs.bc_corrected_read1s)
    if args.read2s is not None:
        cr_utils.copy(args.read2s, outs.bc_corrected_read2s)

def join(args, outs, chunk_defs, chunk_outs):
    summary_files = [
        args.extract_reads_summary,
        args.correct_barcodes_summary,
    ]

    summary_files = [sum_file for sum_file in summary_files if not sum_file is None]

    cr_report.merge_jsons(summary_files, outs.summary)

    cr_utils.copy(args.raw_barcode_counts, outs.raw_barcode_counts)
    cr_utils.copy(args.corrected_barcode_counts, outs.corrected_barcode_counts)
    cr_utils.copy(args.barcode_summary, outs.barcode_summary)
    outs.gem_groups = args.gem_groups
    outs.read_groups = args.read_groups
    outs.align = args.align
    outs.bam_comments = args.bam_comments

    outs.bc_corrected_read1s = [out.bc_corrected_read1s for out in chunk_outs]
    outs.bc_corrected_read2s = [out.bc_corrected_read2s for out in chunk_outs]
