#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""For a single sample, generate the gene expression clustering/umap plots, antibody clustering/umap plots,.

and sample-level rank plots as plotly JSON, so that they can be passed forward to the
WRITE_MULTI_WEB_SUMMARY_JSON stage and inserted into the sample web summary
"""

import json
from math import ceil

import h5py

import cellranger.rna.library as rna_library
from cellranger.analysis.singlegenome import UMAP_NAME
from cellranger.matrix import CountMatrix
from cellranger.multi.barcode_rank_plots import write_library_to_barcode_rank_json
from cellranger.websummary.analysis_tab_core import (
    projection_diffexp_plots_from_path,
    umi_on_projection_from_path,
)
from cellranger.websummary.react_components import ReactComponentEncoder
from cellranger.websummary.treemaps import MIN_ANTIBODY_UMI, make_antibody_treemap_plot

__MRO__ = """
stage GENERATE_SAMPLE_PLOTS(
    in  h5   matrices_h5,
    in  h5   raw_matrices_h5,
    in  path analysis,
    in  h5   barcode_summary,
    out json sample_projection_plots,
    out json sample_library_to_barcode_rank,
    out json sample_treemap_plots,
    src py   "stages/multi/generate_sample_plots",
) split (
) using (
    volatile = strict,
)
"""

# Differential expression clustering plot names and library types.
DIFFEXP_CLUSTERING_PLOT_NAMES = [
    ("gex_diffexp_clustering_plots", rna_library.GENE_EXPRESSION_LIBRARY_TYPE),
    ("antibody_diffexp_clustering_plots", rna_library.ANTIBODY_LIBRARY_TYPE),
]

# pairs colored-by-UMI UMAP plot names with library types
# names should be kept in sync with SampleUmapPlots field names in write_websummary_json.rs
UMI_ON_UMAP_PLOT_NAMES = [
    ("antibody_umi_on_umap", rna_library.ANTIBODY_LIBRARY_TYPE),
    ("crispr_umi_on_umap", rna_library.CRISPR_LIBRARY_TYPE),
    ("custom_umi_on_umap", rna_library.CUSTOM_LIBRARY_TYPE),
]


def split(args):
    _features, _cells, nnz = CountMatrix.load_dims_from_h5(args.matrices_h5)
    mem_gib = ceil(
        CountMatrix.get_mem_gb_from_matrix_h5(args.matrices_h5) + (3 + 31 * nnz / 1024**3)
    )

    return {
        "chunks": [],
        "join": {"__mem_gb": mem_gib},
    }


def write_umap_plots(args, outs):
    """Write UMAP plots."""
    if args.analysis is None:
        outs.sample_projection_plots = None
        return

    # build umap/clustering, differential expression plots for gene expression
    library_types = CountMatrix.load_library_types_from_h5_file(args.matrices_h5)
    umap_plots = {
        plot_name: (
            projection_diffexp_plots_from_path(
                args.analysis, projection=UMAP_NAME, library_type=library_type
            )
            if library_type in library_types
            else None
        )
        for plot_name, library_type in DIFFEXP_CLUSTERING_PLOT_NAMES
    }

    umap_plots.update(
        (
            plot_name,
            umi_on_projection_from_path(
                args.analysis, projection=UMAP_NAME, library_type=library_type
            ),
        )
        for plot_name, library_type in UMI_ON_UMAP_PLOT_NAMES
    )

    with open(outs.sample_projection_plots, "w") as out:
        json.dump(umap_plots, out, sort_keys=True, cls=ReactComponentEncoder)


def write_barcode_rank_plots(args, outs):
    """Write per-sample barcode-rank plots."""
    genomes = CountMatrix.get_genomes_from_h5(args.matrices_h5)
    cell_barcodes = set(CountMatrix.load_bcs_from_h5(args.matrices_h5))
    barcode_summary = h5py.File(args.barcode_summary, "r")

    # If raw matrix for this specific sample is available, restrict to those barcodes in rank plot
    restrict_barcodes = (
        None if args.raw_matrices_h5 is None else CountMatrix.load_bcs_from_h5(args.raw_matrices_h5)
    )

    write_library_to_barcode_rank_json(
        cell_barcodes,
        barcode_summary,
        genomes,
        outs.sample_library_to_barcode_rank,
        restrict_barcodes=restrict_barcodes,
    )


def call_make_treemap(filtered_matrix, lib_type):
    """Generates inputs for make_antibody_treemap_plot function."""
    filtered_counts_matrix = CountMatrix.load_h5_file(filtered_matrix).select_features_by_type(
        lib_type
    )
    return make_antibody_treemap_plot(filtered_counts_matrix, lib_type, MIN_ANTIBODY_UMI, False)


def write_treemap_plots(args, outs):
    """Write per-sample treemap feature UMI distribution plots."""
    # build UMAP/clustering, differential expression plots for antigen and antibody
    sample_library_types = CountMatrix.load_library_types_from_h5_file(args.matrices_h5)
    treemap_plots = {
        library_type: call_make_treemap(args.matrices_h5, library_type)
        for library_type in sample_library_types
        if library_type in [rna_library.ANTIBODY_LIBRARY_TYPE, rna_library.ANTIGEN_LIBRARY_TYPE]
    }
    with open(outs.sample_treemap_plots, "w") as out:
        json.dump(treemap_plots, out, sort_keys=True, cls=ReactComponentEncoder)


def join(args, outs, _chunk_defs, _chunk_outs):
    write_umap_plots(args, outs)
    write_barcode_rank_plots(args, outs)
    write_treemap_plots(args, outs)
