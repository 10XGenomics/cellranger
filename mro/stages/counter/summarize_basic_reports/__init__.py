#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import cellranger.report as cr_report
import cellranger.utils as cr_utils

__MRO__ = """
stage SUMMARIZE_BASIC_REPORTS(
    in  json   extract_reads_summary,
    in  path   reference_path,
    in  map    align,
    in  json   attach_bcs_and_umis_summary,
    in  h5     attach_bcs_and_umis_barcode_summary,
    in  json   mark_duplicates_summary,
    in  json   count_genes_reporter_summary,
    in  h5     count_genes_barcode_summary,
    in  json   filter_barcodes_summary,
    in  json   subsample_molecules_summary,
    in  h5     raw_gene_bc_matrices_h5,
    in  h5     filtered_gene_bc_matrices_h5,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    out h5     barcode_summary,
    out json   summary,
    src py     "stages/counter/summarize_basic_reports",
) split using (
)
"""

def split(args):
    mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups)
    chunks = [{
        '__mem_gb': mem_gb,
    }]
    return {'chunks': chunks}

def main(args, outs):
    summary_files = [
        args.extract_reads_summary,
        args.attach_bcs_and_umis_summary,
        args.mark_duplicates_summary,
        args.count_genes_reporter_summary,
        args.filter_barcodes_summary,
        args.subsample_molecules_summary,
    ]

    barcode_summary_files = [
        args.attach_bcs_and_umis_barcode_summary,
        args.count_genes_barcode_summary,
    ]

    cr_report.merge_jsons(summary_files, outs.summary, [cr_utils.build_alignment_param_metrics(args.align)])
    cr_report.merge_h5(barcode_summary_files, outs.barcode_summary)

def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]

    cr_utils.copy(chunk_out.summary, outs.summary)
    cr_utils.copy(chunk_out.barcode_summary, outs.barcode_summary)
