#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""For a single sample, generate the gene expression clustering/tsne plots, antibody clustering/tsne plots,.

and sample-level rank plots as plotly JSON, so that they can be passed forward to the
WRITE_MULTI_WEB_SUMMARY_JSON stage and inserted into the sample web summary
"""

import json

import h5py

import cellranger.rna.library as rna_library
from cellranger import matrix as cr_matrix
from cellranger.multi.barcode_rank_plots import write_library_to_barcode_rank_json
from cellranger.reference_paths import get_reference_genomes
from cellranger.websummary.analysis_tab_core import (
    tsne_diffexp_plots_from_path,
    umi_on_tsne_from_path,
)
from cellranger.websummary.react_components import ReactComponentEncoder
from cellranger.websummary.treemaps import MIN_ANTIBODY_UMI, make_antibody_treemap_plot

__MRO__ = """
stage GENERATE_SAMPLE_PLOTS(
    in  h5   matrices_h5,
    in  h5   raw_matrices_h5,
    in  path analysis,
    in  h5   barcode_summary,
    in  path reference_path,
    out json sample_tsne_plots,
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

# pairs colored-by-UMI TSNE plot names with library types
# names should be kept in sync with SampleTsnePlots field names in write_websummary_json.rs
UMI_ON_TSNE_PLOT_NAMES = [
    ("antibody_umi_on_tsne", rna_library.ANTIBODY_LIBRARY_TYPE),
    ("crispr_umi_on_tsne", rna_library.CRISPR_LIBRARY_TYPE),
    ("custom_umi_on_tsne", rna_library.CUSTOM_LIBRARY_TYPE),
]


def split(args):
    _features, _cells, nnz = cr_matrix.CountMatrix.load_dims_from_h5(args.matrices_h5)
    return {
        "chunks": [],
        "join": {"__mem_gb": 3 + 31 * nnz / 1024**3},
    }


def write_tsne_plots(args, outs):
    """Write t-SNE plots."""
    if args.analysis is None:
        outs.sample_tsne_plots = None
        return

    # build tsne/clustering, differential expression plots for gene expression
    library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(args.matrices_h5)
    tsne_plots = {
        plot_name: tsne_diffexp_plots_from_path(args.analysis, library_type=library_type)
        if library_type in library_types
        else None
        for plot_name, library_type in DIFFEXP_CLUSTERING_PLOT_NAMES
    }

    tsne_plots.update(
        (plot_name, umi_on_tsne_from_path(args.analysis, library_type=library_type))
        for plot_name, library_type in UMI_ON_TSNE_PLOT_NAMES
    )

    with open(outs.sample_tsne_plots, "w") as out:
        json.dump(tsne_plots, out, sort_keys=True, cls=ReactComponentEncoder)


def write_barcode_rank_plots(args, outs):
    """Write per-sample barcode-rank plots."""
    # FIXME: would be better to load this directly into a set instead of
    # building a list first.
    cell_barcodes = set(cr_matrix.CountMatrix.load_bcs_from_h5(args.matrices_h5))
    barcode_summary = h5py.File(args.barcode_summary, "r")
    genomes = get_reference_genomes(args.reference_path)

    # If raw matrix for this specific sample is available, restrict to those barcodes in rank plot
    restrict_barcodes = (
        None
        if args.raw_matrices_h5 is None
        else cr_matrix.CountMatrix.load_bcs_from_h5(args.raw_matrices_h5)
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
    filtered_counts_matrix = cr_matrix.CountMatrix.load_h5_file(
        filtered_matrix
    ).select_features_by_type(lib_type)
    return make_antibody_treemap_plot(filtered_counts_matrix, lib_type, MIN_ANTIBODY_UMI, False)


def write_treemap_plots(args, outs):
    """Write per-sample treemap feature UMI distribution plots."""
    # build tsne/clustering, differential expression plots for antigen and antibody
    sample_library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(args.matrices_h5)

    treemap_plots = {
        library_type: call_make_treemap(args.matrices_h5, library_type)
        for library_type in sample_library_types
        if library_type in [rna_library.ANTIBODY_LIBRARY_TYPE, rna_library.ANTIGEN_LIBRARY_TYPE]
    }

    with open(outs.sample_treemap_plots, "w") as out:
        json.dump(treemap_plots, out, sort_keys=True, cls=ReactComponentEncoder)


def join(args, outs, _chunk_defs, _chunk_outs):
    write_tsne_plots(args, outs)
    write_barcode_rank_plots(args, outs)
    write_treemap_plots(args, outs)
