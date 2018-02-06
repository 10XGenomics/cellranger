#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import json
import os
import pandas as pd
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annotations
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.report as vdj_report
import cellranger.vdj.utils as vdj_utils
import cellranger.vdj.reference as vdj_ref

__MRO__ = """
stage REPORT_CONTIGS(
    in  path  vdj_reference_path,
    in  json  cell_barcodes,
    in  fasta contigs,
    in  json  annotations,
    in  csv   filter_summary,
    in  tsv   contig_summary,
    in  tsv   umi_summary,
    in  json  reads_summary,
    in  json  assemble_metrics_summary,
    out json  summary,
    src py    "stages/vdj/report_contigs",
) split using (
)
"""

MIN_CHAIN_TYPE_CONTIG_FRAC = 0.05

MEM_GB_PER_UMI_SUMMARY_GB = 4.5

def split(args):
    mem_gb_annot = vdj_utils.get_mem_gb_from_annotations_json(args.annotations)

    umi_summary_bytes = os.path.getsize(args.umi_summary) if args.umi_summary else 0
    mem_gb_umi = MEM_GB_PER_UMI_SUMMARY_GB * float(umi_summary_bytes)/1e9

    mem_gb = max(mem_gb_annot, mem_gb_umi)

    print 'requested %d' % mem_gb
    return {
        'chunks': [{
            '__mem_gb': max(cr_constants.MIN_MEM_GB, mem_gb),
        }]
    }

def main(args, outs):
    reporter = vdj_report.VdjReporter()

    barcode_contigs = defaultdict(list)
    contig_annotations = {}

    # Get annotations for each contig
    for annotation in iter(json.load(open(args.annotations))):
        contig_annotations[annotation['contig_name']] = annotation

    if args.contig_summary and os.path.isfile(args.contig_summary):
        contig_summary = pd.read_csv(args.contig_summary, header=0, index_col=None, sep='\t',
                                     dtype={'component': int, 'num_reads': int,
                                            'num_pairs': int, 'num_umis': int,
                                            'umi_list': str,
                                     })
        contig_summary = contig_summary.groupby('barcode')
    else:
        contig_summary = None

    if args.umi_summary and os.path.isfile(args.umi_summary):
        umi_summary = pd.read_csv(args.umi_summary, header=0, index_col=None, sep='\t')
        umi_summary = umi_summary.groupby('barcode')
    else:
        umi_summary = None

    if args.filter_summary:
        filter_summary = vdj_utils.load_contig_summary_table(args.filter_summary)
    else:
        filter_summary = None

    # Get contigs for each barcode
    for contig_hdr, contig_seq in cr_utils.get_fasta_iter(open(args.contigs)):
        contig_name = contig_hdr.split(' ')[0]
        if not filter_summary is None and not vdj_utils.is_contig_filtered(filter_summary, contig_name):
            continue

        barcode = vdj_utils.get_barcode_from_contig_name(contig_name)
        barcode_contigs[barcode].append((contig_name, contig_seq))

    # Compute metrics for each barcode
    if args.cell_barcodes:
        barcodes = vdj_utils.load_cell_barcodes_json(args.cell_barcodes)
    else:
        # Pass an empty barcode JSON for bulk
        barcodes = {''}


    reference = vdj_ref.VdjReference(args.vdj_reference_path)

    for barcode in barcodes:
        contigs = barcode_contigs[barcode]
        annotations = [contig_annotations[contig[0]] for contig in contigs]

        reporter.vdj_barcode_contig_cb(barcode, contigs, annotations, reference)

        if not contig_summary is None and barcode in contig_summary.groups:
            bc_contig_summary = contig_summary.get_group(barcode)
        else:
            bc_contig_summary = None

        if not umi_summary is None and barcode in umi_summary.groups:
            bc_umi_summary = umi_summary.get_group(barcode)
        else:
            bc_umi_summary = None

        reporter.vdj_assembly_cb(bc_contig_summary, bc_umi_summary, annotations, reference)

    ## Compute post-assembly per-cell metrics
    # Load the assembly metrics summary to get the total assemblable reads
    if args.assemble_metrics_summary and args.reads_summary:
        assemblable_read_pairs_by_bc = cr_utils.get_metric_from_json(args.assemble_metrics_summary, 'assemblable_read_pairs_by_bc')
        assemblable_read_pairs = sum(assemblable_read_pairs_by_bc.get(bc, 0) for bc in barcodes)

        total_read_pairs = cr_utils.get_metric_from_json(args.reads_summary, 'total_read_pairs')

        reporter._get_metric_attr('vdj_assemblable_read_pairs_per_filtered_bc').set_value(assemblable_read_pairs, len(barcodes))
        reporter._get_metric_attr('vdj_sequencing_efficiency').set_value(assemblable_read_pairs, total_read_pairs)

    ## Try to autodetect the chain type
    # Find all chains w/ a significant presence.
    # If there's exactly one, set the chain type filter to that.
    # Otherwise, show all chain types.

    chain_count = defaultdict(int)
    for anno_dict in contig_annotations.itervalues():
        contig = vdj_annotations.AnnotatedContig.from_dict(anno_dict, reference)
        if contig.is_cell and contig.high_confidence and contig.productive:
            for anno in contig.annotations:
                if anno.feature.chain_type in vdj_constants.VDJ_CHAIN_TYPES:
                    chain_count[anno.feature.chain_type] += 1

    outs.chain_type = vdj_constants.ALL_CHAIN_TYPES

    print chain_count

    if len(chain_count) > 0:
        n_contigs = sum(chain_count.itervalues())
        sig_chains = [ct for ct, count in chain_count.iteritems() if tk_stats.robust_divide(count, n_contigs) >= MIN_CHAIN_TYPE_CONTIG_FRAC]
        if len(sig_chains) == 1:
            outs.chain_type = sig_chains[0]

    reporter.report_summary_json(outs.summary)

def join(args, outs, chunk_defs, chunk_outs):
    outs.chain_type = chunk_outs[0].chain_type
    cr_utils.copy(chunk_outs[0].summary, outs.summary)
