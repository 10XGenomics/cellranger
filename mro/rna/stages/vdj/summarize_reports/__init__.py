#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import json

import cellranger.webshim.common as cr_webshim
import cellranger.websummary.vdj as vdj_web
from cellranger.analysis.singlegenome import TSNE_NAME
from cellranger.webshim.constants.shared import PIPELINE_VDJ
from cellranger.websummary.sample_properties import SampleDataPaths, VdjSampleProperties

__MRO__ = """
stage SUMMARIZE_VDJ_REPORTS(
    in  json         metrics_summary_json,
    in  ChemistryDef vdj_chemistry_def,
    in  json         cell_barcodes,
    in  csv          clonotype_summary,
    out csv          metrics_summary_csv,
    out html         web_summary,
    out json         web_summary_data,
    src py           "stages/vdj/summarize_reports",
) split (
) 
"""


def split(args):
    return {
        "chunks": [{}],
        "join": {
            "__mem_gb": 2,
        },
    }


def main(args, outs):
    pass


def join(args, outs, chunk_defs, chunk_outs):

    sample_data_paths = SampleDataPaths(
        summary_path=args.metrics_summary_json,
        vdj_clonotype_summary_path=args.clonotype_summary,
        vdj_cell_barcodes_path=args.cell_barcodes,
        vdj_all_contig_annotations_csv_path=args.all_contig_annotations_csv,
    )

    with open(args.metrics_summary_json) as f:
        metrics_summary = json.load(f)

    sample_properties = VdjSampleProperties(
        sample_id=metrics_summary["sample_id"],
        sample_desc=metrics_summary["sample_desc"],
        chemistry_def=args.vdj_chemistry_def,
        chain_type=metrics_summary["chain_type"],
    )
    sample_data = cr_webshim.load_sample_data(
        sample_properties, sample_data_paths, projections=TSNE_NAME
    )

    ws_data = vdj_web.build_vdj_web_summary_html(outs.web_summary, sample_properties, sample_data)
    with open(outs.web_summary_data, "w") as f:
        json.dump(ws_data, f, indent=4)
    cr_webshim.build_metrics_summary_csv(
        outs.metrics_summary_csv, sample_properties, sample_data, PIPELINE_VDJ
    )
