#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import collections
import json
import math
import os
import pandas as pd
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
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
    out json  summary,
    src py    "stages/vdj/report_contigs",
) split using (
)
"""

MEM_GB_PER_ANNOTATIONS_JSON_GB = 20

def split(args):
    annotation_json_gb = float(os.path.getsize(args.annotations))/1e9
    mem_gb = int(math.ceil(float(MEM_GB_PER_ANNOTATIONS_JSON_GB) * annotation_json_gb))
    print 'requested %d' % mem_gb
    return {
        'chunks': [{
            '__mem_gb': max(cr_constants.MIN_MEM_GB, mem_gb),
        }]
    }

def main(args, outs):
    reporter = vdj_report.VdjReporter()

    barcode_contigs = collections.defaultdict(list)
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

    reporter.report_summary_json(outs.summary)

def join(args, outs, chunk_defs, chunk_outs):
    cr_utils.copy(chunk_outs[0].summary, outs.summary)
