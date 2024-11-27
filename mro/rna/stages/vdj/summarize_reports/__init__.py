#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import json

import cellranger.report as cr_report
import cellranger.webshim.common as cr_webshim
import cellranger.websummary.vdj as vdj_web
from cellranger.analysis.singlegenome import TSNE_NAME
from cellranger.webshim.constants.shared import PIPELINE_VDJ
from cellranger.websummary.sample_properties import SampleDataPaths, VdjSampleProperties

__MRO__ = """
stage SUMMARIZE_VDJ_REPORTS(
    in  string       sample_id,
    in  string       sample_desc,
    in  ChemistryDef vdj_chemistry_def,
    in  json[]       summaries,
    in  int          total_read_pairs,
    in  json         cell_barcodes,
    in  csv          clonotype_summary,
    in  csv          barcode_support,
    in  string       receptor,
    in  int          n50_n50_rpu,
    out string       receptor,
    out json         metrics_summary_json,
    out csv          metrics_summary_csv,
    out html         web_summary,
    out json         web_summary_data,
    src py           "stages/vdj/summarize_reports",
) split (
) retain (
    metrics_summary_json,
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
    sample_info = {
        "sample_id": args.sample_id,
        "sample_desc": args.sample_desc,
        "chain_type": args.receptor,
        # Hack to pass this metric at the per-sample level, CELLRANGER-7783
        "VDJ_total_read_pairs": args.total_read_pairs,
        "n50_n50_rpu": args.n50_n50_rpu,
    }

    cr_report.merge_jsons(args.summaries, outs.metrics_summary_json, dicts=[sample_info])

    sample_data_paths = SampleDataPaths(
        summary_path=outs.metrics_summary_json,
        vdj_clonotype_summary_path=args.clonotype_summary,
        vdj_barcode_support_path=args.barcode_support,
        vdj_cell_barcodes_path=args.cell_barcodes,
    )

    outs.receptor = args.receptor

    sample_properties = VdjSampleProperties(
        sample_id=args.sample_id,
        sample_desc=args.sample_desc,
        chemistry_def=args.vdj_chemistry_def,
        chain_type=outs.receptor,
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
