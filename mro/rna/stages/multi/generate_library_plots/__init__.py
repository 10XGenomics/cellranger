#!/usr/bin/env python3
#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Generate library plots for multi web summary.

For a single library, generate the barcode rank, jibes plots as plotly JSON, and barnyard count
biplots, so that they can be passed forward to the WRITE_MULTI_WEB_SUMMARY stage and inserted into
the sample web summary
"""

from __future__ import annotations

import json
import pickle
from typing import TYPE_CHECKING

import h5py
import martian

import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
from cellranger.analysis.multigenome import MultiGenomeAnalysis
from cellranger.analysis.singlegenome import UMAP_NAME
from cellranger.feature.feature_assignments import CellsPerFeature
from cellranger.feature_ref import HASHTAG_TAG
from cellranger.multi.barcode_rank_plots import write_library_to_barcode_rank_json
from cellranger.rna.library import BarcodeMultiplexingType, CellLevel
from cellranger.webshim.common import populate_biplot_chart
from cellranger.webshim.jibes_web import make_jibes_biplot_histogram
from cellranger.websummary.analysis_tab_aux import (
    cmo_tags_on_umap_from_path,
    initialize_barnyard_biplot_chart,
)
from cellranger.websummary.analysis_tab_core import umi_on_projection_from_path
from cellranger.websummary.react_components import ReactComponentEncoder

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage GENERATE_LIBRARY_PLOTS(
    in  bool               disable_count,
    # for barcode rank
    in  h5                 barcode_summary_h5,
    in  csv                filtered_barcodes,
    in  ReferenceInfo      reference_info,
    # for jibes biplot
    in  pickle             tag_assigner_pickle,
    # for barnyard biplot and cmo UMAP plot
    in  path               analysis,
    # cmo UMAP plot
    in  json               cells_per_tag,
    in  json               non_singlet_barcodes,
    in  BarcodeAssignments force_sample_barcodes,
    in  string             multiplexing_method,
    out json               library_to_barcode_rank,
    out json               jibes_biplot_histogram,
    out json               barnyard_biplot,
    out json               cmo_projection_plot,
    src py                 "stages/multi/generate_library_plots",
) split (
) using (
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """GENERATE_LIBRARY_PLOTS stage inputs."""

    disable_count: bool
    barcode_summary_h5: mro_filetypes.h5
    filtered_barcodes: mro_filetypes.csv
    reference_info: mro_structs.ReferenceInfo
    tag_assigner_pickle: None | mro_filetypes.pickle
    analysis: None | mro_filetypes.path
    cells_per_tag: None | mro_filetypes.json
    non_singlet_barcodes: None | mro_filetypes.json
    force_sample_barcodes: None | mro_structs.BarcodeAssignments
    multiplexing_method: None | str


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """GENERATE_LIBRARY_PLOTS stage outputs."""

    library_to_barcode_rank: None | mro_filetypes.json
    jibes_biplot_histogram: None | mro_filetypes.json
    barnyard_biplot: None | mro_filetypes.json
    cmo_projection_plot: None | mro_filetypes.json


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """GENERATE_LIBRARY_PLOTS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """GENERATE_LIBRARY_PLOTS chunk outputs."""


def split(args: StageInputs):
    if args.multiplexing_method is not None:
        multiplexing_method = BarcodeMultiplexingType(args.multiplexing_method)
        hashtag_gib = 4 if multiplexing_method.is_cell_multiplexed() else 0
    else:
        multiplexing_method = None
        hashtag_gib = 0

    with open(args.filtered_barcodes) as f:
        filtered_barcodes = sum(1 for _ in f)
    with h5py.File(args.barcode_summary_h5, "r") as f:
        raw_barcodes = len(f["bc_sequence"])
    barcodes_gib = round(108 * raw_barcodes / (1024**3), 1)

    mem_gib = 2 + hashtag_gib + barcodes_gib
    print(
        f"{multiplexing_method=},{filtered_barcodes=},{raw_barcodes=},"
        f"{barcodes_gib=},{hashtag_gib=},{mem_gib=}"
    )
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(
    args: StageInputs,
    outs: StageOutputs,
    _chunk_defs: list[ChunkInputs],
    _chunk_outs: list[ChunkOutputs],
):
    if args.disable_count:
        outs.library_to_barcode_rank = None
        outs.jibes_biplot_histogram = None
        outs.barnyard_biplot = None
        outs.cmo_projection_plot = None
        return

    # BARCODE RANK PLOTS FOR MULTI WEBSUMMARY

    # read in data
    barcode_summary = h5py.File(args.barcode_summary_h5, "r")
    genomes: list[str] = args.reference_info["genomes"]
    cell_barcodes = (
        cr_utils.get_cell_associated_barcode_set(args.filtered_barcodes)
        if args.filtered_barcodes
        else None
    )

    write_library_to_barcode_rank_json(
        cell_barcodes,
        barcode_summary,
        genomes,
        outs.library_to_barcode_rank,
        restrict_barcodes=set(),
    )

    # JIBES PLOTS FOR MULTI WEBSUMMARY
    # make the jibes biplot
    force_sample_barcodes = args.force_sample_barcodes["sample_barcodes"]
    cells_per_tag = args.force_sample_barcodes["cells_per_tag"] or args.cells_per_tag

    if (
        args.multiplexing_method is not None
        and BarcodeMultiplexingType(args.multiplexing_method).is_cell_multiplexed()
        and args.analysis is not None
        and cells_per_tag is not None
        and args.non_singlet_barcodes is not None
    ):
        multiplexing_method = BarcodeMultiplexingType(args.multiplexing_method)

        if multiplexing_method.type == CellLevel.Hashtag:
            cmo_umi_projection_plot = umi_on_projection_from_path(
                args.analysis,
                projection=UMAP_NAME,
                library_type=multiplexing_method.multiplexing_library_type(),
                tag_type=HASHTAG_TAG,
            )
        else:
            cmo_umi_projection_plot = umi_on_projection_from_path(
                args.analysis,
                projection=UMAP_NAME,
                library_type=multiplexing_method.multiplexing_library_type(),
            )
        if cmo_umi_projection_plot is not None:
            cmo_umi_projection_plot["layout"]["title"] = None

        cells_per_tag = CellsPerFeature.load_from_file(cells_per_tag)
        non_singlet_barcodes = CellsPerFeature.load_from_file(args.non_singlet_barcodes)

        cmo_tags_projection_plot = cmo_tags_on_umap_from_path(
            args.analysis,
            cells_per_tag,
            non_singlet_barcodes,
            multiplexing_method,
        )
        if cmo_tags_projection_plot is not None:
            cmo_tags_projection_plot["layout"]["title"] = None

        cmo_projection_plot = {
            "cmo_umi_projection_plot": cmo_umi_projection_plot,
            "cmo_tags_projection_plot": cmo_tags_projection_plot,
        }

        with open(outs.cmo_projection_plot, "w") as out:
            json.dump(
                tk_safe_json.json_sanitize(cmo_projection_plot),
                out,
                sort_keys=True,
                indent=4,
                cls=ReactComponentEncoder,
            )
    else:
        outs.cmo_projection_plot = None

    if args.tag_assigner_pickle is not None and force_sample_barcodes is None:
        with open(args.tag_assigner_pickle, "rb") as p:
            jibes_biplot_histogram = make_jibes_biplot_histogram(pickle.load(p))

        with open(outs.jibes_biplot_histogram, "w") as out:
            json.dump(
                tk_safe_json.json_sanitize(jibes_biplot_histogram),
                out,
                sort_keys=True,
                indent=4,
                cls=ReactComponentEncoder,
            )

    else:
        outs.jibes_biplot_histogram = None

    # The barnyard biplot is derived from the multi-genome secondary analysis, which is
    # absent when secondary analysis is disabled (`--no-secondary`). In that case fall
    # through to the `else` branch and emit no biplot, mirroring the single-genome path.
    if len(genomes) >= 2 and args.analysis is not None:
        barnyard_biplot = initialize_barnyard_biplot_chart()
        multi_genome_analysis = MultiGenomeAnalysis.load_default_format(args.analysis)
        populate_biplot_chart(barnyard_biplot, multi_genome_analysis)
        with open(outs.barnyard_biplot, "w") as out:
            json.dump(
                tk_safe_json.json_sanitize(barnyard_biplot),
                out,
                sort_keys=True,
                cls=ReactComponentEncoder,
            )

    else:
        outs.barnyard_biplot = None
