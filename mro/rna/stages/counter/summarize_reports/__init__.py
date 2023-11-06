#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
import json
import os
import shutil

import cellranger.report as cr_report
import cellranger.webshim.common as cr_webshim
import cellranger.websummary.cmdline as cml
import cellranger.websummary.sample_properties as sp
from cellranger.reference_paths import get_reference_genomes
from cellranger.webshim.constants.shared import PIPELINE_COUNT
from cellranger.websummary.react_components import ReactComponentEncoder
from cellranger.websummary.web_summary_builder import build_web_summary_html_sc

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  json[]   summaries,
    in  string   sample_id,
    in  string   sample_desc,
    in  path     reference_path,
    in  path     analysis,
    in  h5       barcode_summary_h5,
    in  h5       filtered_gene_bc_matrices_h5,
    in  csv      filtered_barcodes,
    in  tps.json target_panel_summary,
    in  string   barcode_whitelist,
    in  json     antibody_histograms,
    in  json     antibody_treemap,
    in  json     antigen_histograms,
    in  json     antigen_treemap,
    in  csv      feature_reference,
    in  string   target_set_name,
    in  csv      per_feature_metrics_csv,
    in  bool     include_introns,
    in  string   throughput,
    out json     metrics_summary_json,
    out csv      metrics_summary_csv,
    out html     web_summary,
    out csv      feature_reference,
    out json     ws_data,
    src py       "stages/counter/summarize_reports",
) using (
    mem_gb   = 16,
    volatile = strict,
) retain (
    metrics_summary_json,
)
"""


def main(args, outs):
    id_dict = {"sample_id": args.sample_id, "sample_desc": args.sample_desc}
    if args.target_panel_summary is not None:
        args.summaries.append(args.target_panel_summary)
    cr_report.merge_jsons(args.summaries, outs.metrics_summary_json, [id_dict])

    sample_data_paths = sp.SampleDataPaths(
        summary_path=outs.metrics_summary_json,
        barcode_summary_path=args.barcode_summary_h5,
        analysis_path=args.analysis,
        filtered_barcodes_path=args.filtered_barcodes,
        feature_metrics_path=args.per_feature_metrics_csv,
        antibody_histograms_path=args.antibody_histograms,
        antibody_treemap_path=args.antibody_treemap,
        antigen_histograms_path=args.antigen_histograms,
        antigen_treemap_path=args.antigen_treemap,
    )

    cmdline = os.environ.get("CMDLINE")
    if cmdline:
        cmdline_parsed = cml.parse_cmdline_basename(cmdline)
    else:
        cmdline_parsed = "NA"
    genomes = get_reference_genomes(args.reference_path)
    sample_properties = sp.ExtendedCountSampleProperties(
        sample_id=args.sample_id,
        sample_desc=args.sample_desc,
        barcode_whitelist=args.barcode_whitelist,
        reference_path=args.reference_path,
        include_introns=args.include_introns,
        throughput=args.throughput,
        target_set=args.target_set_name,
        target_panel_summary=args.target_panel_summary,
        genomes=genomes,
        cmdline=cmdline_parsed,
    )

    # TODO: Move metrics CSV somewhere else
    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)
    cr_webshim.build_metrics_summary_csv(
        outs.metrics_summary_csv, sample_properties, sample_data, PIPELINE_COUNT
    )

    websummary_data = build_web_summary_html_sc(
        outs.web_summary, sample_properties, sample_data_paths, PIPELINE_COUNT
    )

    if args.feature_reference is not None:
        shutil.copy(args.feature_reference, outs.feature_reference)

    with open(outs.ws_data, "w") as out:
        json.dump(websummary_data, out, sort_keys=True, indent=4, cls=ReactComponentEncoder)
