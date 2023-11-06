#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Functions to barcode rank plots for multi websummary."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import cellranger.rna.library as rna_library
import tenkit.safe_json as tk_safe_json
from cellranger.webshim.common import plot_basic_barcode_rank
from cellranger.websummary.react_components import ReactComponentEncoder
from cellranger.websummary.summary_tab import get_empty_rank_plot

if TYPE_CHECKING:
    import h5py

LIBRARY_TYPES = [
    rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.ANTIGEN_LIBRARY_TYPE,
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.CUSTOM_LIBRARY_TYPE,
    rna_library.MULTIPLEXING_LIBRARY_TYPE,
]


def get_barcode_rank_json(
    lib_type: str,
    cell_barcodes: set[bytes],
    barcode_summary: h5py.File,
    genomes: list[str],
    restrict_barcodes: list[bytes] | None = None,
):
    """Creates a barcode rank plot for a library type, dumps it to a json file."""
    prefix = rna_library.get_library_type_metric_prefix(lib_type)
    chart = get_empty_rank_plot()
    del chart["layout"]["title"]  # remove duplicate title
    plot = plot_basic_barcode_rank(
        chart, cell_barcodes, barcode_summary, genomes, prefix, restrict_barcodes=restrict_barcodes
    )
    return tk_safe_json.json_sanitize(plot) if plot else None


def make_all_lib_type_rank_plots(
    cell_barcodes: set[bytes],
    barcode_summary: h5py.File,
    genomes: list[str],
    restrict_barcodes: list[bytes] | None = None,
):
    """Returns a dictionary mapping each library type to a JSON file containing its rank plot.

    only non-empty entries are retained.
    """
    library_to_barcode_rank = {
        lib_type: get_barcode_rank_json(
            lib_type, cell_barcodes, barcode_summary, genomes, restrict_barcodes=restrict_barcodes
        )
        for lib_type in LIBRARY_TYPES
    }

    # get rid of empties and return
    return {k: v for k, v in library_to_barcode_rank.items() if v}


def write_library_to_barcode_rank_json(
    cell_barcodes: set[bytes],
    barcode_summary: h5py.File,
    genomes: list[str],
    output_file: str,
    restrict_barcodes: list[bytes] | None = None,
):
    """Convenience function to call make_all_lib_type_rank_plots and write output to JSON."""
    library_to_barcode_rank = make_all_lib_type_rank_plots(
        cell_barcodes, barcode_summary, genomes, restrict_barcodes=restrict_barcodes
    )

    with open(output_file, "w") as output:
        json.dump(
            library_to_barcode_rank,
            output,
            indent=4,
            sort_keys=True,
            cls=ReactComponentEncoder,
        )
