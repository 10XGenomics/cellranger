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

import h5py

import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
from cellranger.analysis.multigenome import MultiGenomeAnalysis
from cellranger.analysis.singlegenome import UMAP_NAME
from cellranger.feature.feature_assignments import CellsPerFeature
from cellranger.feature_ref import HASHTAG_TAG
from cellranger.hdf5 import estimate_mem_gib_for_dataset
from cellranger.multi.barcode_rank_plots import write_library_to_barcode_rank_json
from cellranger.webshim.common import populate_biplot_chart
from cellranger.webshim.jibes_web import make_jibes_biplot_histogram
from cellranger.websummary.analysis_tab_aux import (
    cmo_tags_on_umap_from_path,
    initialize_barnyard_biplot_chart,
)
from cellranger.websummary.analysis_tab_core import umi_on_projection_from_path
from cellranger.websummary.react_components import ReactComponentEncoder

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


def split(args):
    with h5py.File(args.barcode_summary_h5, "r") as f:
        barcodes_gib = estimate_mem_gib_for_dataset(f["bc_sequence"])
    # memory usage of the stage scales approximately with the number of barcodes
    mem_gib = 6 + round(2.5 * barcodes_gib, 1)
    print(f"{barcodes_gib=},{mem_gib=}")
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(args, outs, _chunk_defs, _chunk_outs):
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
        and rna_library.BarcodeMultiplexingType(args.multiplexing_method).is_cell_multiplexed()
        and args.analysis is not None
        and cells_per_tag is not None
        and args.non_singlet_barcodes is not None
    ):
        multiplexing_method = rna_library.BarcodeMultiplexingType(args.multiplexing_method)

        if multiplexing_method.type == rna_library.CellLevel.Hashtag:
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

    if len(genomes) >= 2:
        assert args.analysis is not None

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
