#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import json
import os
import shutil
from typing import TYPE_CHECKING

import martian

import cellranger.report as cr_report
import cellranger.webshim.common as cr_webshim
import cellranger.websummary.sample_properties as sp
from cellranger.analysis.singlegenome import UMAP_NAME
from cellranger.matrix import CountMatrix
from cellranger.webshim.constants.shared import PIPELINE_COUNT
from cellranger.websummary.react_components import ReactComponentEncoder
from cellranger.websummary.web_summary_builder import build_web_summary_html_sc

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  map<ChemistryDef> chemistry_defs,
    in  json[]            summaries,
    in  string            sample_id,
    in  string            sample_desc,
    in  ReferenceInfo     reference_info,
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
    in  json[]            cell_type_websummary_bundle_value,
    in  bool              cloud_cell_annotation_viable_but_not_requested,
    in  bool              skip_cell_annotation,
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
    ws_data,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_REPORTS stage inputs."""

    chemistry_defs: dict[str, mro_structs.ChemistryDef]
    summaries: list[mro_filetypes.json]
    sample_id: str
    sample_desc: str
    reference_info: mro_structs.ReferenceInfo
    analysis: None | mro_filetypes.path
    barcode_summary_h5: mro_filetypes.h5
    filtered_gene_bc_matrices_h5: mro_filetypes.h5
    filtered_barcodes: mro_filetypes.csv
    target_panel_summary: None | mro_filetypes.tps_json
    antibody_histograms: None | mro_filetypes.json
    antibody_treemap: None | mro_filetypes.json
    antigen_histograms: None | mro_filetypes.json
    antigen_treemap: None | mro_filetypes.json
    feature_reference: None | mro_filetypes.csv
    target_set_name: None | str
    per_feature_metrics_csv: None | mro_filetypes.csv
    include_introns: bool
    cell_type_websummary_bundle_value: None | list[None | mro_filetypes.json]
    cloud_cell_annotation_viable_but_not_requested: bool
    skip_cell_annotation: bool


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_REPORTS stage outputs."""

    metrics_summary_json: mro_filetypes.json
    metrics_summary_csv: mro_filetypes.csv
    web_summary: mro_filetypes.html
    feature_reference: None | mro_filetypes.csv
    ws_data: mro_filetypes.json


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_REPORTS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_REPORTS chunk outputs."""


def split(args: StageInputs):
    _num_features, num_barcodes, nnz = CountMatrix.load_dims_from_h5(
        args.filtered_gene_bc_matrices_h5
    )
    mem_gib = 7 + CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=1.0)
    print(f"{num_barcodes=},{nnz=},{mem_gib=}")
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(
    args: StageInputs,
    outs: StageOutputs,
    _chunk_defs: list[ChunkInputs],
    _chunk_outs: list[ChunkOutputs],
):
    id_dict = {"sample_id": args.sample_id, "sample_desc": args.sample_desc}
    cr_report.merge_jsons(
        args.summaries,
        outs.metrics_summary_json,
        [id_dict],
    )

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

    sample_properties = sp.ExtendedCountSampleProperties(
        sample_id=args.sample_id,
        sample_desc=args.sample_desc,
        genomes=CountMatrix.get_genomes_from_h5(args.filtered_gene_bc_matrices_h5),
        reference_path=(
            args.reference_info["transcriptome_info"]["reference_path"]
            if args.reference_info["transcriptome_info"]
            else None
        ),
        chemistry_defs=args.chemistry_defs,
        include_introns=args.include_introns,
        target_set=args.target_set_name,
        target_panel_summary=args.target_panel_summary,
        cmdline=os.environ.get("CMDLINE", "NA"),
        cloud_cell_annotation_viable_but_not_requested=args.cloud_cell_annotation_viable_but_not_requested,
        skip_cell_annotation=args.skip_cell_annotation,
    )

    # TODO: Move metrics CSV somewhere else
    sample_data = cr_webshim.load_sample_data(
        sample_properties, sample_data_paths, projections=(UMAP_NAME,)
    )
    cr_webshim.build_metrics_summary_csv(
        outs.metrics_summary_csv, sample_properties, sample_data, PIPELINE_COUNT
    )

    bundle_values = None
    if args.cell_type_websummary_bundle_value:
        bundle_values = []
        for path in args.cell_type_websummary_bundle_value:
            if path is not None:
                with open(path) as f:
                    bundle_values.append(json.load(f))

    websummary_data = build_web_summary_html_sc(
        outs.web_summary,
        sample_properties,
        sample_data_paths,
        PIPELINE_COUNT,
        cell_type_websummary_bundle_value=bundle_values,
    )

    if args.feature_reference is not None:
        shutil.copy(args.feature_reference, outs.feature_reference)

    with open(outs.ws_data, "w") as out:
        json.dump(websummary_data, out, sort_keys=True, indent=4, cls=ReactComponentEncoder)
