#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# This stage only exists to force earlier VDR of prior stages.

import itertools
import cellranger.chemistry as cr_chem
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.io as cr_io

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
    in  fastq[]  read1s,
    in  fastq[]  read2s,
    in  bool     retain_fastqs,
    in  json     trim_reads_summary,
    in  map      chemistry_def,
    out json     summary,
    out json     raw_barcode_counts,
    out json     corrected_barcode_counts,
    out h5       barcode_summary,
    out int[]    gem_groups,
    out string[] read_groups,
    out map      align,
    out string[] bam_comments,
    out fastq[]  read1s,
    out fastq[]  read2s,
    src py       "stages/vdj/summarize_read_reports",
) split using (
    in  fastq    read1,
    in  fastq    read2,
    in  tsv      bcs,
)
"""

def split(args):
    paired_end = cr_chem.is_paired_end(args.chemistry_def)
    if paired_end:
        assert len(args.read1s) == len(args.read2s)
    assert len(args.corrected_bcs) == len(args.read1s)

    chunks = []

    if args.retain_fastqs:
        for read1, read2, bcs in itertools.izip_longest(args.read1s, args.read2s, args.corrected_bcs):
            chunks.append({
                'read1': read1,
                'read2': read2 if paired_end else None,
                'bcs': bcs,
            })

    return {'chunks': chunks}


def main(args, outs):
    if args.read1 is not None:
        # Ensure same extension
        out_path, _ = cr_utils.splitexts(outs.read1s)
        _, in_ext = cr_utils.splitexts(args.read1)
        outs.read1s = out_path + in_ext
        cr_io.copy(args.read1, outs.read1s)

    if args.read2 is not None:
        out_path, _ = cr_utils.splitexts(outs.read2s)
        _, in_ext = cr_utils.splitexts(args.read2)
        outs.read2s = out_path + in_ext
        cr_io.copy(args.read2, outs.read2s)

    if args.bcs is not None:
        out_path, _ = cr_utils.splitexts(outs.corrected_bcs)
        _, in_ext = cr_utils.splitexts(args.bcs)
        outs.corrected_bcs = out_path + in_ext
        cr_io.copy(args.bcs, outs.corrected_bcs)

def join(args, outs, chunk_defs, chunk_outs):
    summary_files = [
        args.extract_reads_summary,
        args.correct_barcodes_summary,
        args.trim_reads_summary,
    ]

    summary_files = [sum_file for sum_file in summary_files if not sum_file is None]

    cr_report.merge_jsons(summary_files, outs.summary)

    cr_io.copy(args.raw_barcode_counts, outs.raw_barcode_counts)
    cr_io.copy(args.corrected_barcode_counts, outs.corrected_barcode_counts)
    cr_io.copy(args.barcode_summary, outs.barcode_summary)
    outs.gem_groups = args.gem_groups
    outs.read_groups = args.read_groups
    outs.align = args.align
    outs.bam_comments = args.bam_comments

    outs.read1s = [co.read1s for co in chunk_outs]
    outs.read2s = [co.read2s for co in chunk_outs]
    outs.corrected_bcs = [co.corrected_bcs for co in chunk_outs]
