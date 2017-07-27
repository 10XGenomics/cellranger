#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

import itertools
import numpy as np
import tenkit.bam as tk_bam
import cellranger.chemistry as cr_chem
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.report as vdj_report
import cellranger.vdj.utils as vdj_utils

__MRO__ = '''
stage REPORT_CONTIG_ALIGNMENTS(
    in  bam        contig_bam,
    in  string[][] barcodes_in_chunks,
    in  map        chemistry_def,
    out pickle     chunked_reporter,
    out json       summary,
    src py         "stages/vdj/report_contig_alignments",
) split using (
    in  string[]   contigs,
)
'''

# Number of chunks to generate for bulk (non-barcoded) data
BULK_NCHUNKS = 100

def split(args):
    # Use a dummy chunk to appease Martian if there was no input
    if args.contig_bam is None or not vdj_utils.bam_has_seqs(args.contig_bam):
        return {
            'chunks': [
                {
                    'contigs': None,
                },
            ],
        }

    # Reuse the chunk structure from assemble_vdj_pd to chunk sets of contigs by barcode chunks
    chunks = []
    bam_iter = tk_bam.create_bam_infile(args.contig_bam)

    all_contigs = list(bam_iter.references)

    if args.barcodes_in_chunks is None:
        # Contigs are not grouped by barcode. Just split into a few chunks.
        contigs_per_chunk = int(np.ceil(len(all_contigs) / float(BULK_NCHUNKS)))
        for ch in range(BULK_NCHUNKS):
            start = ch * contigs_per_chunk
            stop = min(len(all_contigs), (ch + 1) * contigs_per_chunk)
            contig_names = list(itertools.islice(all_contigs, start, stop))
            chunks.append({'contigs': contig_names})
    else:
        for barcode_chunk in args.barcodes_in_chunks:
            barcode_chunk = set(barcode_chunk)
            contig_names = [contig for contig in all_contigs if vdj_utils.get_barcode_from_contig_name(contig) in barcode_chunk]
            chunks.append({'contigs': contig_names, '__mem_gb':6})

    bam_iter.close()
    return {'chunks': chunks}


def main(args, outs):
    # Handle the dummy chunk
    if args.contigs is None:
        outs.chunked_reporter = None
        return

    # Calculate metrics on assigned contigs
    reporter = vdj_report.VdjReporter()

    bam = tk_bam.create_bam_infile(args.contig_bam)
    contigs = [contig for contig in bam.header['SQ'] if contig['SN'] in set(args.contigs)]

    for contig in contigs:
        contig_name = contig['SN']

        # Fetch indexed portion of BAM.
        read_iter = bam.fetch(str(contig_name))

        reporter.contig_mapping_frac_statistics_cb(read_iter, contig['LN'],
                                                   cr_chem.get_strandedness(args.chemistry_def))

    reporter.save(outs.chunked_reporter)


def normalize_metric(reporter, metric_name, denominator):
    metric = reporter._get_metric_attr(metric_name)
    for cutoff in metric.d:
        metric.d[cutoff] = int(metric.d[cutoff] / float(denominator))

def join(args, outs, chunk_defs, chunk_outs):
    if all(chunk_out.chunked_reporter is None for chunk_out in chunk_outs):
        cr_utils.write_empty_json(outs.summary)
        return

    # Merge reporters and save
    reporters = [chunk_out.chunked_reporter for chunk_out in chunk_outs]
    bam = tk_bam.create_bam_infile(args.contig_bam)
    ncontigs = len(bam.references)
    final_report = cr_report.merge_reporters(reporters)
    # These metrics have accumulated the total read coverage on all contigs.
    # Normalize by the number of contigs.
    normalize_metric(final_report, 'vdj_contig_depth_at_contig_percentiles', ncontigs)
    normalize_metric(final_report, 'vdj_contig_q40_depth_at_contig_percentiles', ncontigs)

    # Add unmapped reads to reporter to get correct alignment rates
    unmapped_read_count = cr_utils.get_unmapped_read_count_from_indexed_bam(args.contig_bam)
    final_report._get_metric_attr('vdj_contig_mapping_frac').add(unmapped_read_count, filter=False)

    # Save final report
    final_report.report_summary_json(outs.summary)
