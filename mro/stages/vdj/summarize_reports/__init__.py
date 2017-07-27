#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import martian
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.data as cr_webshim_data
from cellranger.webshim.constants.shared import PIPELINE_VDJ

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  string sample_id,
    in  string sample_desc,
    in  path   vdj_reference_path,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    in  json   reads_summary,
    in  json   filter_umis_summary,
    in  json   filter_barcodes_summary,
    in  json   trim_reads_summary,
    in  json   filter_reads_summary,
    in  json   filter_contigs_summary,
    in  json   report_contigs_summary,
    in  json   report_contig_alignments_summary,
    in  json   group_clonotypes_summary,
    in  json   raw_consensus_summary,
    in  h5     barcode_summary,
    in  json   cell_barcodes,
    in  csv    barcode_umi_summary,
    in  h5     umi_info,
    in  csv    clonotype_summary,
    in  csv    barcode_support,
    out json   metrics_summary_json,
    out csv    metrics_summary_csv,
    out html   web_summary,
    out json   alerts,
    out h5     barcode_summary,
    out json   cell_barcodes,
    out csv    barcode_support,
    out csv    barcode_umi_summary,
    out h5     umi_info,
    src py     "stages/vdj/summarize_reports",
) split using (
)
"""

def split(args):
    mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups)
    return {
        'chunks': [{}],
        'join': {
            '__mem_gb': mem_gb,
        },
    }

def main(args, outs):
    pass

def join(args, outs, chunk_defs, chunk_outs):
    summary_files = [
        args.reads_summary,
        args.filter_umis_summary,
        args.filter_barcodes_summary,
        args.trim_reads_summary,
        args.filter_reads_summary,
        args.filter_contigs_summary,
        args.report_contigs_summary,
        args.report_contig_alignments_summary,
        args.raw_consensus_summary,
        args.group_clonotypes_summary,
    ]

    summary_files = [sum_file for sum_file in summary_files if not sum_file is None]

    cr_report.merge_jsons(summary_files, outs.metrics_summary_json)

    # Copy barcode summary h5
    if args.barcode_summary:
        cr_utils.copy(args.barcode_summary, outs.barcode_summary)

    # Copy cell barcodes
    if args.cell_barcodes:
        cr_utils.copy(args.cell_barcodes, outs.cell_barcodes)

    # Copy barcode support
    if args.barcode_support:
        cr_utils.copy(args.barcode_support, outs.barcode_support)

    # Copy barcode umi summary
    if args.barcode_umi_summary:
        cr_utils.copy(args.barcode_umi_summary, outs.barcode_umi_summary)

    # Copy umi info
    if args.umi_info:
        cr_utils.copy(args.umi_info, outs.umi_info)

    sample_data_paths = cr_webshim_data.SampleDataPaths(
        summary_path=outs.metrics_summary_json,
        barcode_summary_path=args.barcode_summary,
        vdj_clonotype_summary_path=args.clonotype_summary,
        vdj_barcode_support_path=args.barcode_support,
    )

    sample_properties = cr_webshim.get_sample_properties(args.sample_id, args.sample_desc, [], version=martian.get_pipelines_version())

    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)

    if args.barcode_whitelist is not None:
        cr_webshim.build_web_summary_html(outs.web_summary, sample_properties, sample_data, PIPELINE_VDJ,
                                          alerts_output_filename=outs.alerts)
        cr_webshim.build_metrics_summary_csv(outs.metrics_summary_csv, sample_properties, sample_data, PIPELINE_VDJ)
