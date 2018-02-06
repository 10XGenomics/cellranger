#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import itertools
import tenkit.fasta as tk_fasta
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.report as cr_report
import cellranger.stats as cr_stats
import cellranger.utils as cr_utils
import cellranger.vdj.report as vdj_report

__MRO__ = """
stage CORRECT_BARCODES(
    in  fastq[] read1s,
    in  fastq[] read2s,
    in  int[]   gem_groups,
    in  json    barcode_counts,
    in  float   barcode_confidence_threshold,
    in  string  barcode_whitelist,
    in  int     initial_reads,
    out tsv[]   corrected_bcs,
    out json    corrected_barcode_counts,
    out pickle  chunked_reporter,
    out json    summary,
    out h5      barcode_summary,
    src py      "stages/vdj/correct_barcodes",
) split using (
    in  fastq   read1_chunk,
    in  fastq   read2_chunk,
    in  int     gem_group,
)
"""

def split(args):
    chunks = []

    if cr_chem.is_paired_end(args.chemistry_def):
        assert len(args.read1s) == len(args.read2s)

    assert len(args.gem_groups) == len(args.read1s)

    for read1, read2, gem_group in itertools.izip_longest(args.read1s, args.read2s, args.gem_groups):
        chunks.append({
            'read1_chunk': read1,
            'read2_chunk': read2,
            'gem_group': gem_group,
        })
    return {'chunks': chunks}

def main(args, outs):
    # Load barcode whitelist
    if args.barcode_whitelist is not None:
        barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)

    reporter = vdj_report.VdjReporter()

    # Load barcode count distribution
    barcode_dist = cr_utils.load_barcode_dist(args.barcode_counts,
                                              barcode_whitelist,
                                              args.gem_group)

    if args.barcode_whitelist is not None:
        barcode_whitelist_set = set(barcode_whitelist)
    else:
        barcode_whitelist_set = None

    in_read1_fastq = cr_utils.open_maybe_gzip(args.read1_chunk)
    in_read2_fastq = cr_utils.open_maybe_gzip(args.read2_chunk) if args.read2_chunk else []

    outs.corrected_bcs += cr_constants.LZ4_SUFFIX
    out_file = cr_utils.open_maybe_gzip(outs.corrected_bcs, 'w')

    bc_counter = cr_fastq.BarcodeCounter(args.barcode_whitelist, outs.corrected_barcode_counts)

    # Correct barcodes, add processed bc tag to fastq
    read_pair_iter = itertools.izip_longest(tk_fasta.read_generator_fastq(in_read1_fastq), \
                                            tk_fasta.read_generator_fastq(in_read2_fastq))
    for read1, read2 in itertools.islice(read_pair_iter, args.initial_reads):
        read1_header = cr_fastq.AugmentedFastqHeader(read1[0])

        raw_bc = read1_header.get_tag(cr_constants.RAW_BARCODE_TAG)
        bc_qual = read1_header.get_tag(cr_constants.RAW_BARCODE_QUAL_TAG)

        processed_bc = None

        if raw_bc:
            if barcode_whitelist_set is not None and raw_bc not in barcode_whitelist_set:
                processed_bc = cr_stats.correct_bc_error(args.barcode_confidence_threshold,
                                                         raw_bc, bc_qual, barcode_dist)
            else:
                # Disallow Ns in no-whitelist case
                if 'N' in raw_bc:
                    processed_bc = None
                else:
                    processed_bc = raw_bc

            if processed_bc:
                bc_counter.count(None, processed_bc, None)

                # Add gem group to barcode sequence
                processed_bc = cr_utils.format_barcode_seq(processed_bc, gem_group=args.gem_group)

            reporter.vdj_barcode_cb(raw_bc, processed_bc)

        out_file.write('%s\n' % (processed_bc if processed_bc is not None else ''))

    in_read1_fastq.close()
    if in_read2_fastq:
        in_read2_fastq.close()
    out_file.close()

    bc_counter.close()

    reporter.save(outs.chunked_reporter)

def join(args, outs, chunk_defs, chunk_outs):
    outs.corrected_bcs = [co.corrected_bcs for co in chunk_outs]

    bc_counter = cr_fastq.BarcodeCounter(args.barcode_whitelist, outs.corrected_barcode_counts, gem_groups=args.gem_groups)
    for chunk_def, chunk_out in zip(chunk_defs, chunk_outs):
        bc_counter.merge(chunk_def.gem_group, chunk_out.corrected_barcode_counts)
    bc_counter.close()

    outs.chunked_reporter = None
    reporter = cr_report.merge_reporters([chunk_out.chunked_reporter for chunk_out in chunk_outs])

    reporter.report_summary_json(outs.summary)

    reporter.report_barcodes_h5(outs.barcode_summary)
