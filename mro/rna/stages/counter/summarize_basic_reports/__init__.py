#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Compute summary metrics from the filtered barcode matrix h5 and combine all basic summaries."""

from __future__ import annotations

import h5py

import cellranger.report as cr_report
import cellranger.rna.library as rna_library
import cellranger.rna.refactor_report_matrix as refactor_rna_report_mat
import cellranger.rna.report_matrix as rna_report_mat
import cellranger.utils as cr_utils
from cellranger.cell_calling_helpers import CellCallingParam, get_recovered_cells
from cellranger.fast_utils import FilteredBarcodes
from cellranger.matrix import CountMatrix

__MRO__ = """
stage SUMMARIZE_BASIC_REPORTS(
    in  string           sample,
    in  h5               matrices_h5,
    in  csv              filtered_barcodes,
    in  csv              per_barcode_metrics,
    in  json             matrix_computer_summary,
    in  h5               barcode_summary,
    in  CellCallingParam recovered_cells,
    in  json[]           summary_jsons,
    in  tps.json         target_panel_summary,
    in  bool             sample_bcs_only,
    out json             summary,
    src py               "stages/counter/summarize_basic_reports",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    _num_features, num_barcodes, nnz = CountMatrix.load_dims_from_h5(args.matrices_h5)
    matrix_gib = CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=0.7)
    with h5py.File(args.barcode_summary) as f:
        len_barcode_summary = len(f["bc_sequence"])
    barcode_summary_gib = round(272 * len_barcode_summary / 1024**3, 1)
    mem_gib = 3 + matrix_gib + barcode_summary_gib
    print(
        f"{len_barcode_summary=},{num_barcodes=},{nnz=},{matrix_gib=},{barcode_summary_gib=},{mem_gib=}"
    )
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(args, outs, chunk_defs, chunk_outs):
    if args.target_panel_summary is not None:
        args.summary_jsons.append(args.target_panel_summary)
    raw_matrix = CountMatrix.load_h5_file(args.matrices_h5)
    genomes = raw_matrix.get_genomes()
    filtered_barcodes = FilteredBarcodes(args.filtered_barcodes)
    genome_filtered_bcs = {
        genome: set(barcodes or ())
        for genome, barcodes in (
            dict.fromkeys(genomes) | filtered_barcodes.per_genome_barcodes()
        ).items()
    }

    # Re-compute various metrics on the filtered matrix
    reads_summary = cr_utils.merge_jsons_as_dict([args.matrix_computer_summary])

    matrix_summary = rna_report_mat.report_genomes(
        raw_matrix,
        reads_summary=reads_summary,
        barcode_summary_h5_path=args.barcode_summary,
        per_barcode_metrics_path=args.per_barcode_metrics,
        recovered_cells=get_recovered_cells(CellCallingParam(args.recovered_cells), args.sample),
        cell_bc_seqs=genome_filtered_bcs,
        sample_bc_seqs=raw_matrix.bcs if args.sample_bcs_only else None,
    )

    is_targeted = raw_matrix.feature_ref.has_target_features()

    del genome_filtered_bcs
    del raw_matrix

    library_types = CountMatrix.load_library_types_from_h5_file(args.matrices_h5)
    antibody_present = rna_library.ANTIBODY_LIBRARY_TYPE in library_types
    crispr_present = rna_library.CRISPR_LIBRARY_TYPE in library_types
    custom_present = rna_library.CUSTOM_LIBRARY_TYPE in library_types
    per_cell_metrics = refactor_rna_report_mat.compute_per_cell_metrics(
        filtered_barcodes=filtered_barcodes,
        per_barcode_metrics_path=args.per_barcode_metrics,
        genomes=genomes,
        antibody_present=antibody_present,
        crispr_present=crispr_present,
        custom_present=custom_present,
        is_targeted=is_targeted,
    )
    cr_report.merge_jsons(args.summary_jsons, outs.summary, [matrix_summary, per_cell_metrics])
