#!/usr/bin/env python
#
# Copyright (c) 2022 10x Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import martian

import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.websummary.sample_properties as wsp
from cellranger.analysis.singlegenome import TSNE_NAME, UMAP_NAME
from cellranger.websummary.aggr_websummary_builder import (
    build_web_summary_data_aggr,
    build_web_summary_html_aggr,
)
from cellranger.websummary.react_components import ReactComponentEncoder

if TYPE_CHECKING:

    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage SUMMARIZE_AGGREGATED_REPORTS(
    in  string sample_id,
    in  string sample_desc,
    in  map    gem_group_index,
    in  h5     filtered_matrices_h5,
    in  path   analysis,
    in  map[]  sample_defs,
    in  json   normalize_depth_summary,
    in  json   analyze_matrices_summary,
    in  json   antibody_histograms,
    in  json   antibody_treemap,
    in  json   crispr_analysis_metrics,
    in  string product_type,
    in  bool   skip_tsne,
    in  json[] cell_type_websummary_bundle_value,
    in  bool   cloud_cell_annotation_viable_but_not_requested,
    out json   summary,
    out html   web_summary,
    out json   web_summary_data,
    src py     "stages/aggregator/summarize_aggregated_reports",
) split (
) using (
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_AGGREGATED_REPORTS stage inputs."""

    sample_id: str
    sample_desc: str
    gem_group_index: dict
    filtered_matrices_h5: mro_filetypes.h5
    analysis: mro_filetypes.path
    sample_defs: list[dict]  # TODO: be more specific with the type here
    normalize_depth_summary: mro_filetypes.json
    analyze_matrices_summary: mro_filetypes.json | None
    antibody_histograms: mro_filetypes.json | None
    antibody_treemap: mro_filetypes.json | None
    crispr_analysis_metrics: mro_filetypes.json | None
    product_type: str
    skip_tsne: bool
    cell_type_websummary_bundle_value: None | list[None | mro_filetypes.json]
    cloud_cell_annotation_viable_but_not_requested: bool


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_AGGREGATED_REPORTS stage outputs."""

    summary: mro_filetypes.json
    web_summary: mro_filetypes.html
    web_summary_data: mro_filetypes.json


ChunkInputs = StageInputs
ChunkOutputs = StageOutputs


def split(args: StageInputs):
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrices_h5)
    join_mem_gb = max(matrix_mem_gb, h5_constants.MIN_MEM_GB)
    return {"chunks": [], "join": {"__mem_gb": join_mem_gb}}


def join(
    args: StageInputs,
    outs: StageOutputs,
    _chunk_defs: list[ChunkInputs],
    _chunk_outs: list[ChunkOutputs],
):
    summary = {}

    genomes = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrices_h5)
    target_set = cr_matrix.CountMatrix.load_feature_ref_from_h5_file(
        args.filtered_matrices_h5
    ).get_target_feature_ids()
    assert args.product_type in [
        cr_constants.SPATIAL_PRODUCT_TYPE,
        cr_constants.SINGLE_CELL_PRODUCT_TYPE,
    ]
    is_spatial = args.product_type == cr_constants.SPATIAL_PRODUCT_TYPE

    # get metrics from other summaries
    if args.analyze_matrices_summary:
        with open(args.analyze_matrices_summary) as reader:
            summary.update(json.load(reader))
    with open(args.normalize_depth_summary) as reader:
        summary.update(json.load(reader))
        agg_batches = summary["batches"]
    if args.crispr_analysis_metrics:
        with open(args.crispr_analysis_metrics) as reader:
            summary.update(json.load(reader))
    summary["samples_count"] = len(args.sample_defs)
    with open(outs.summary, "w") as f:
        json.dump(summary, f, indent=4, sort_keys=True)

    # build web summary
    # First, get the properties of the data we want to run through the builder
    # In this case, an aggr run
    sample_properties = wsp.AggrCountSampleProperties(
        sample_id=args.sample_id,
        sample_desc=args.sample_desc,
        genomes=genomes,
        agg_batches=agg_batches,
        is_spatial=is_spatial,
        target_set=target_set,
        cloud_cell_annotation_viable_but_not_requested=args.cloud_cell_annotation_viable_but_not_requested,
    )
    # Second, build a SampleDataPaths object listing the necessary file paths to produce the websummary
    sample_data_paths = wsp.SampleDataPaths(
        summary_path=outs.summary,
        barcode_summary_path=None,
        analysis_path=args.analysis,
        antibody_histograms_path=args.antibody_histograms,
        antibody_treemap_path=args.antibody_treemap,
    )
    # Call the websummary builder.
    gg_id_to_name = {int(id): name[0] for id, name in args.gem_group_index.items()}

    projection = UMAP_NAME if args.skip_tsne else TSNE_NAME

    bundle_values = None
    if args.cell_type_websummary_bundle_value:
        bundle_values = []
        for path in args.cell_type_websummary_bundle_value:
            if path is not None:
                with open(path) as f:
                    bundle_values.append(json.load(f))

    build_web_summary_html_aggr(
        filename=outs.web_summary,
        sample_properties=sample_properties,
        gg_id_to_name_map=gg_id_to_name,
        sample_data_paths=sample_data_paths,
        sample_defs=args.sample_defs,
        projection=projection,
        cell_type_websummary_bundle_value=bundle_values,
    )

    # Do it again because ReactComponentEncoder transforms and deletes data while encoding
    ws_data = build_web_summary_data_aggr(
        sample_properties=sample_properties,
        gg_id_to_name_map=gg_id_to_name,
        sample_data_paths=sample_data_paths,
        sample_defs=args.sample_defs,
        projection=projection,
        cell_type_websummary_bundle_value=bundle_values,
    )
    with open(outs.web_summary_data, "w") as f:
        json.dump(ws_data, f, indent=4, cls=ReactComponentEncoder)
