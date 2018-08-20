#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import martian
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.io as cr_io
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.data as cr_webshim_data
from cellranger.webshim.constants.gex import CountSampleProperties
from cellranger.webshim.constants.shared import PIPELINE_COUNT

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  json[] summaries,
    in  string sample_id,
    in  string sample_desc,
    in  path   reference_path,
    in  path   analysis,
    in  h5     barcode_summary_h5,
    in  h5     filtered_gene_bc_matrices_h5,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    out json   metrics_summary_json,
    out csv    metrics_summary_csv,
    out html   web_summary,
    out json   alerts,
    src py     "stages/counter/summarize_reports",
) split using (
)
"""

def split(args):
    chunks = [{
        '__mem_gb': 2,
    }]
    return {'chunks': chunks, 'join': {'__mem_gb': 1}}

def main(args, outs):
    cr_report.merge_jsons(args.summaries, outs.metrics_summary_json)

    sample_data_paths = cr_webshim_data.SampleDataPaths(
        summary_path=outs.metrics_summary_json,
        barcode_summary_path=args.barcode_summary_h5,
        analysis_path=args.analysis,
    )

    genomes = cr_utils.get_reference_genomes(args.reference_path)
    sample_properties = CountSampleProperties(sample_id=args.sample_id,
                                              sample_desc=args.sample_desc,
                                              genomes=genomes,
                                              version=martian.get_pipelines_version())
    sample_properties = dict(sample_properties._asdict())

    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)

    cr_webshim.build_web_summary_html(outs.web_summary, sample_properties, sample_data, PIPELINE_COUNT,
                                      alerts_output_filename=outs.alerts)
    cr_webshim.build_metrics_summary_csv(outs.metrics_summary_csv, sample_properties, sample_data, PIPELINE_COUNT)

def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]

    cr_io.copy(chunk_out.web_summary, outs.web_summary)
    cr_io.copy(chunk_out.alerts, outs.alerts)
    cr_io.copy(chunk_out.metrics_summary_json, outs.metrics_summary_json)
    cr_io.copy(chunk_out.metrics_summary_csv, outs.metrics_summary_csv)
