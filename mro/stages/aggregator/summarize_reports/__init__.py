#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import json
import martian

import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.data as cr_webshim_data
from cellranger.webshim.constants.gex import AggrSampleProperties
from cellranger.webshim.constants.shared import PIPELINE_AGGR

__MRO__ = """
stage SUMMARIZE_AGGREGATED_REPORTS(
    in  string aggregation_id,
    in  string aggregation_desc,
    in  map    gem_group_index,
    in  h5     filtered_matrices_h5,
    in  h5     barcode_summary_h5,
    in  path   analysis,
    in  json   normalize_depth_summary,
    in  json   count_genes_summary,
    in  json   analyze_matrices_summary,
    out json   summary,
    out html   web_summary,
    src py     "stages/aggregator/summarize_reports",
) split using (
)
"""

def split(args):
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrices_h5)
    chunks = [{
        '__mem_gb': max(matrix_mem_gb, h5_constants.MIN_MEM_GB),
    }]
    return {'chunks': chunks}

def main(args, outs):
    summary = {}

    filtered_mat = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrices_h5)
    genomes = filtered_mat.get_genomes()

    # get metrics from other summaries
    if args.analyze_matrices_summary:
        with open(args.analyze_matrices_summary) as reader:
            analysis_summary = json.load(reader)
        summary.update(analysis_summary)

    with open(args.normalize_depth_summary, 'r') as reader:
        summary.update(json.load(reader))
        agg_batches = summary['batches']

    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)

    # build web summary
    sample_properties = AggrSampleProperties(sample_id=args.aggregation_id,
                                             sample_desc=args.aggregation_desc,
                                             genomes=genomes,
                                             version=martian.get_pipelines_version(),
                                             agg_batches=agg_batches)
    sample_properties = dict(sample_properties._asdict())

    sample_data_paths = cr_webshim_data.SampleDataPaths(
        summary_path=outs.summary,
        barcode_summary_path=args.barcode_summary_h5,
        analysis_path=args.analysis,
    )

    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)
    cr_webshim.build_web_summary_html(outs.web_summary, sample_properties, sample_data, PIPELINE_AGGR)

def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]
    cr_io.copy(chunk_out.summary, outs.summary)
    cr_io.copy(chunk_out.web_summary, outs.web_summary)
