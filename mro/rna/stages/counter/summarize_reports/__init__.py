#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
import json
import os
import shutil

import cellranger.report as cr_report
import cellranger.webshim.common as cr_webshim
import cellranger.websummary.sample_properties as sp
from cellranger.analysis.singlegenome import TSNE_NAME
from cellranger.matrix import CountMatrix
from cellranger.reference_paths import get_reference_genomes
from cellranger.webshim.constants.shared import PIPELINE_COUNT
from cellranger.websummary.react_components import ReactComponentEncoder
from cellranger.websummary.web_summary_builder import build_web_summary_html_sc

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  map<ChemistryDef> chemistry_defs,
    in  json[]            summaries,
    in  string            sample_id,
    in  string            sample_desc,
    in  path              reference_path,
    in  path              analysis,
    in  h5                barcode_summary_h5,
    in  h5                filtered_gene_bc_matrices_h5,
    in  csv               filtered_barcodes,
    in  tps.json          target_panel_summary,
    in  json              antibody_histograms,
    in  json              antibody_treemap,
    in  json              antigen_histograms,
    in  json              antigen_treemap,
    in  csv               feature_reference,
    in  string            target_set_name,
    in  csv               per_feature_metrics_csv,
    in  bool              include_introns,
    out json              metrics_summary_json,
    out csv               metrics_summary_csv,
    out html              web_summary,
    out csv               feature_reference,
    out json              ws_data,
    src py                "stages/counter/summarize_reports",
) split (
) using (
    volatile = strict,
) retain (
    metrics_summary_json,
)
"""


def split(args):
    _num_features, num_barcodes, nnz = CountMatrix.load_dims_from_h5(
        args.filtered_gene_bc_matrices_h5
    )
    mem_gib = 7 + CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=1.0)
    print(f"{num_barcodes=},{nnz=},{mem_gib=}")
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(args, outs, _chunk_defs, _chunk_outs):
    id_dict = {"sample_id": args.sample_id, "sample_desc": args.sample_desc}
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

    genomes = get_reference_genomes(args.reference_path)
    sample_properties = sp.ExtendedCountSampleProperties(
        sample_id=args.sample_id,
        sample_desc=args.sample_desc,
        genomes=genomes,
        reference_path=args.reference_path,
        chemistry_defs=args.chemistry_defs,
        include_introns=args.include_introns,
        target_set=args.target_set_name,
        target_panel_summary=args.target_panel_summary,
        cmdline=os.environ.get("CMDLINE", "NA"),
    )

    # TODO: Move metrics CSV somewhere else
    sample_data = cr_webshim.load_sample_data(
        sample_properties, sample_data_paths, projections=(TSNE_NAME,)
    )
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
