#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Functions for summarizing a CountMatrix."""

from __future__ import annotations

from collections import OrderedDict
from functools import reduce
from typing import AnyStr

import h5py as h5
import numpy as np
from six import ensure_binary

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import tenkit.stats as tk_stats
from cellranger import csv_utils


def _get_total_reads(summary, library_type):
    prefix = rna_library.get_library_type_metric_prefix(library_type)
    return int(summary.get(prefix + "total_read_pairs", 0))


def _get_conf_mapped_reads(summary, genomes, library_type):
    prefix = rna_library.get_library_type_metric_prefix(library_type)
    conf_mapped_metrics = [
        prefix
        + "_".join(
            [
                ref,
                cr_constants.TRANSCRIPTOME_REGION,
                cr_constants.CONF_MAPPED_READ_TYPE,
                "reads_frac",
            ]
        )
        for ref in genomes
    ]
    total_reads = _get_total_reads(summary, library_type)
    return sum(float(summary.get(metric, 0)) * float(total_reads) for metric in conf_mapped_metrics)


def _get_barcode_summary_h5_indices(summary_h5: h5.File, barcodes: list[bytes]):
    """Get the indices into a barcode summary HDF5 file for the given barcodes.

    Does NOT preserve order of the barcodes argument.

    Args:
      summary_h5 (h5py.File): HDF5 file with datasets
      barcodes (iterable of bytes): List of barcode strings (e.g., "ACGT-1") to get.

    Returns:
      np.array(int): List of indices into the HDF5 arrays
    """
    if len(barcodes) == 0:
        return np.zeros(0, dtype=int)
    max_len = max(len(bc) for bc in barcodes)
    bc_arr = np.fromiter(barcodes, count=len(barcodes), dtype="S%d" % max_len)
    return np.flatnonzero(np.isin(summary_h5["bc_sequence"][:], bc_arr))


MIN_COUNTS_PER_BARCODE = 2
MATRIX_REPORT_READ_TYPES = [
    cr_constants.CONF_MAPPED_BC_READ_TYPE,
    cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
]
MATRIX_USE_MATRIX_FOR_READ_TYPE = [cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE]


def _report_genome_agnostic_metrics(
    ftype,
    matrix,
    barcode_summary_h5,
    recovered_cells,
    sample_bc_indices: list[int] | None,
    cell_bcs_union: set[bytes],
    matrix_cell_bc_indices: list[int],
    bc_summary_cell_bc_indices: list[int],
    total_reads,
    gex_reads_per_barcode,
    total_conf_mapped_reads,
    library_prefix,
):
    """Report metrics that are computed across all barcodes and all genomes."""
    d = {}

    matrix = matrix.view()

    genomes = matrix.get_genomes()
    if len(genomes) == 0:
        # For genomeless features, use "multi" in place of the genome prefix
        genomes = [rna_library.MULTI_REFS_PREFIX]

    for g in genomes:
        assert isinstance(g, str)
    # Get number of cell bcs across all genomes
    n_cell_bcs_union = len(cell_bcs_union)

    for bc in cell_bcs_union:
        assert isinstance(bc, bytes)

    d["filtered_bcs_transcriptome_union"] = n_cell_bcs_union
    # TODO: This doesn't seem to be used anywhere anymore
    d[rna_library.add_multi_prefix("filtered_bcs")] = n_cell_bcs_union

    # Report reads/cell across all genomes
    mean_reads_per_cell = tk_stats.robust_divide(total_reads, n_cell_bcs_union)
    d[
        f"{rna_library.MULTI_REFS_PREFIX}_{cr_constants.TRANSCRIPTOME_REGION}_total_raw_reads_per_filtered_bc"
    ] = mean_reads_per_cell
    # Create a feature-barcode dual whose name makes sense
    d["reads_per_cell"] = mean_reads_per_cell

    if ftype == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
        assert isinstance(next(iter(gex_reads_per_barcode.keys())), bytes)
        # Report mean reads in cell per cell across all genomes.
        n_filtered_reads = sum(
            raw_reads
            for barcode, raw_reads in gex_reads_per_barcode.items()
            if barcode in cell_bcs_union
        )
        d["filtered_reads_per_filtered_bc"] = tk_stats.robust_divide(
            n_filtered_reads, n_cell_bcs_union
        )

    d[
        f"{rna_library.MULTI_REFS_PREFIX}_{cr_constants.TRANSCRIPTOME_REGION}_total_conf_mapped_reads_per_filtered_bc"
    ] = tk_stats.robust_divide(total_conf_mapped_reads, n_cell_bcs_union)

    # Split the matrix by genome
    if genomes[0] != rna_library.MULTI_REFS_PREFIX:
        genome_matrices = OrderedDict((g, matrix.select_features_by_genome(g)) for g in genomes)
    else:
        # Genomeless feature types
        genome_matrices = OrderedDict(((rna_library.MULTI_REFS_PREFIX, matrix),))

    # Total UMI counts across all matrices and all filtered barcodes
    total_umi_counts = 0
    for mat in genome_matrices.values():
        total_umi_counts += mat.select_barcodes(matrix_cell_bc_indices).sum()

    # Deviation from cell load
    if recovered_cells is None:
        d[f"{rna_library.MULTI_REFS_PREFIX}_filtered_bcs_difference_from_recovered_cells"] = 0
        d[
            f"{rna_library.MULTI_REFS_PREFIX}_filtered_bcs_relative_difference_from_recovered_cells"
        ] = 0
    else:
        d[f"{rna_library.MULTI_REFS_PREFIX}_filtered_bcs_difference_from_recovered_cells"] = int(
            n_cell_bcs_union
        ) - int(recovered_cells)
        d[
            f"{rna_library.MULTI_REFS_PREFIX}_filtered_bcs_relative_difference_from_recovered_cells"
        ] = tk_stats.robust_divide(n_cell_bcs_union - recovered_cells, recovered_cells)

    # Duplicate these metrics across genomes for backwards-compat
    for genome in genomes:
        d[f"{genome}_total_raw_reads_per_filtered_bc"] = tk_stats.robust_divide(
            total_reads, n_cell_bcs_union
        )
        d[f"{genome}_total_conf_mapped_reads_per_filtered_bc"] = tk_stats.robust_divide(
            total_conf_mapped_reads, n_cell_bcs_union
        )

        for read_type in MATRIX_REPORT_READ_TYPES:
            metric = f"{genome}_total_{read_type}_reads_per_filtered_bc"
            if read_type in MATRIX_USE_MATRIX_FOR_READ_TYPE:
                n_reads = total_umi_counts
            else:
                # Note the extra underscore after the library_type prefix.
                # This is induced by the Reporter framework.
                h5_keys = [
                    cr_utils.format_barcode_summary_h5_key(
                        library_prefix, g, cr_constants.TRANSCRIPTOME_REGION, read_type
                    )
                    for g in genomes
                ]
                h5_keys = [x for x in h5_keys if x in barcode_summary_h5]
                if sample_bc_indices is not None:
                    n_reads = sum(
                        np.array(barcode_summary_h5[h5_key][:][sample_bc_indices]).sum()
                        for h5_key in h5_keys
                    )
                else:
                    n_reads = sum(np.array(barcode_summary_h5[h5_key]).sum() for h5_key in h5_keys)
            d[metric] = tk_stats.robust_divide(n_reads, n_cell_bcs_union)

    # Report frac reads in cells across all genomes
    total_conf_mapped_reads_in_cells = 0
    total_conf_mapped_barcoded_reads = 0

    for genome in genome_matrices:
        h5_key = f"{library_prefix}_{genome}_{cr_constants.TRANSCRIPTOME_REGION}_{cr_constants.CONF_MAPPED_BC_READ_TYPE}_reads"
        cmb_reads = barcode_summary_h5[h5_key][:]
        total_conf_mapped_reads_in_cells += cmb_reads[bc_summary_cell_bc_indices].sum()
        if sample_bc_indices is not None:
            total_conf_mapped_barcoded_reads += cmb_reads[sample_bc_indices].sum()
        else:
            total_conf_mapped_barcoded_reads += cmb_reads.sum()
    frac_reads_in_cells = tk_stats.robust_divide(
        total_conf_mapped_reads_in_cells, total_conf_mapped_barcoded_reads
    )
    d["multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac"] = frac_reads_in_cells
    # Create a feature-barcode dual whose name makes sense
    d["feature_reads_in_cells"] = frac_reads_in_cells

    # Compute fraction of reads usable (conf mapped, barcoded, filtered barcode)
    unique_barcodes = set(cell_bcs_union)
    in_unique_barcodes_vectorized = np.vectorize(lambda x: x in unique_barcodes)
    filtered_bc_h5_row = in_unique_barcodes_vectorized(np.array(barcode_summary_h5["bc_sequence"]))

    usable_reads = 0

    for genome in genomes:
        h5_key = cr_utils.format_barcode_summary_h5_key(
            library_prefix,
            genome,
            cr_constants.TRANSCRIPTOME_REGION,
            cr_constants.CONF_MAPPED_BC_READ_TYPE,
        )

        if h5_key not in barcode_summary_h5:
            continue

        usable_reads += (filtered_bc_h5_row * np.array(barcode_summary_h5[h5_key])).sum()

    # Fraction reads usable
    d[f"{rna_library.MULTI_REFS_PREFIX}_transcriptome_usable_reads_frac"] = tk_stats.robust_divide(
        usable_reads, total_reads
    )
    # Create a feature barcoding dual whose name makes sense
    d["frac_feature_reads_usable"] = tk_stats.robust_divide(usable_reads, total_reads)

    # Usable reads
    d[f"{rna_library.MULTI_REFS_PREFIX}_usable_reads"] = usable_reads

    # Usable reads per cell
    reads_usable_per_cell = tk_stats.robust_divide(usable_reads, n_cell_bcs_union)
    d[f"{rna_library.MULTI_REFS_PREFIX}_usable_reads_per_filtered_bc"] = reads_usable_per_cell
    # Create a feature barcoding dual whose name makes sense
    d["feature_reads_usable_per_cell"] = reads_usable_per_cell

    # Compute matrix density
    filtered_mat = matrix.select_barcodes(matrix_cell_bc_indices)
    total_nonzero_entries = filtered_mat.get_num_nonzero()
    filtered_shape = filtered_mat.get_shape()
    total_entries = filtered_shape[0] * filtered_shape[1]
    print(total_entries, total_nonzero_entries, filtered_shape)
    d[f"{rna_library.MULTI_REFS_PREFIX}_filtered_gene_bc_matrix_density"] = tk_stats.robust_divide(
        total_nonzero_entries, total_entries
    )

    return d


def _report(
    matrix: cr_matrix.CountMatrixView,
    genome: str,
    barcode_summary_h5: h5.File,
    matrix_cell_bc_indices: list[int],
    bc_summary_cell_bc_indices: list[int],
    sample_bc_indices: list[int] | None,
    library_prefix: str,
):
    d = {}
    matrix = matrix.view()

    filtered_mat = matrix.select_barcodes(matrix_cell_bc_indices)
    filtered_mat_shape = filtered_mat.get_shape()

    n_cell_bcs = len(matrix_cell_bc_indices)

    # Don't compute metrics if no cells detected
    if n_cell_bcs == 0:
        return d

    # Compute matrix density
    d["filtered_gene_bc_matrix_density"] = tk_stats.robust_divide(
        filtered_mat.get_num_nonzero(), filtered_mat_shape[0] * filtered_mat_shape[1]
    )

    counts_per_gene: np.ndarray = filtered_mat.sum(axis=1)
    genes_top_n = min(cr_constants.TOP_N, len(counts_per_gene))
    top_genes_with_counts = {
        filtered_mat.int_to_feature_id(i): int(count)
        for i, count in cr_matrix.top_n(counts_per_gene, genes_top_n)
    }
    d["filtered_bcs_top_genes_with_reads"] = top_genes_with_counts

    unique_bcs_per_gene = filtered_mat.count_ge(axis=1, threshold=MIN_COUNTS_PER_BARCODE)
    top_genes_with_unique_bcs = {
        filtered_mat.int_to_feature_id(i): int(count)
        for i, count in cr_matrix.top_n(unique_bcs_per_gene, genes_top_n)
    }
    d["filtered_bcs_top_genes_with_unique_bcs"] = top_genes_with_unique_bcs

    # Total genes and counts
    total_genes_detected = np.count_nonzero(counts_per_gene)
    total_counts = int(counts_per_gene.sum())
    d["filtered_bcs_total_unique_genes_detected"] = total_genes_detected
    d["filtered_bcs_total_counts"] = total_counts

    def _summarize_per_barcode(a):
        mean = np.mean(a)
        stddev = np.std(a)
        return {
            "mean": mean,
            "median": np.median(a),
            "cv": tk_stats.robust_divide(float(stddev), float(mean)),
            "iqr": np.percentile(a, 75) - np.percentile(a, 25),
        }

    # Unique genes per bc
    unique_genes_per_bc = filtered_mat.count_ge(axis=0, threshold=cr_constants.MIN_COUNTS_PER_GENE)
    unique_genes_stats = _summarize_per_barcode(unique_genes_per_bc)
    for stat, value in unique_genes_stats.items():
        d[f"filtered_bcs_{stat}_unique_genes_detected"] = value

    # Counts per bc
    counts_per_bc = filtered_mat.sum(axis=0)
    counts_per_bc_stats = _summarize_per_barcode(counts_per_bc)
    for stat, value in counts_per_bc_stats.items():
        d[f"filtered_bcs_{stat}_counts"] = value

    # Cumulative fraction of counts going to top bcs
    filt_total_umis = filtered_mat.sum()
    raw_total_umis = matrix.sum()
    d["filtered_bcs_cum_frac"] = tk_stats.robust_divide(filt_total_umis, raw_total_umis)

    # cDNA PCR duplication in top bcs
    dupe_candidate_h5_key = cr_utils.format_barcode_summary_h5_key(
        library_prefix,
        genome,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_BC_READ_TYPE,
    )
    if dupe_candidate_h5_key in barcode_summary_h5:
        n_reads = barcode_summary_h5[dupe_candidate_h5_key][:][bc_summary_cell_bc_indices].sum()
        n_deduped_reads = filt_total_umis
    else:
        n_reads = 0
        n_deduped_reads = 0
    d[f"filtered_bcs_{cr_constants.CDNA_PCR_DUPE_TYPE}_dupe_reads_frac"] = (
        1 - tk_stats.robust_divide(n_deduped_reads, n_reads)
    )

    # Reads per top bc for the various read types (computed over top bcs)
    for read_type in MATRIX_REPORT_READ_TYPES:
        # Compute (n_reads)/(n_bcs) over all bcs and over top bcs
        per_bc_metric = f"filtered_bcs_{read_type}_reads_per_filtered_bc"

        # Cumulative fraction of reads going to top bcs
        frac_metric = f"filtered_bcs_{read_type}_reads_cum_frac"

        if read_type in MATRIX_USE_MATRIX_FOR_READ_TYPE:
            n_reads = filt_total_umis
            n_all_reads = raw_total_umis
        else:
            h5_key = cr_utils.format_barcode_summary_h5_key(
                library_prefix, genome, cr_constants.TRANSCRIPTOME_REGION, read_type
            )
            if h5_key in barcode_summary_h5:
                counts = barcode_summary_h5[h5_key][:]
                n_reads = counts[bc_summary_cell_bc_indices].sum()
                if sample_bc_indices is not None:
                    n_all_reads = counts[sample_bc_indices].sum()
                else:
                    n_all_reads = counts.sum()
            else:
                n_reads = 0
                n_all_reads = 0
        d[per_bc_metric] = tk_stats.robust_divide(n_reads, n_cell_bcs)
        d[frac_metric] = tk_stats.robust_divide(n_reads, n_all_reads)
    return d


def load_per_barcode_metrics(per_barcode_metrics_path: AnyStr) -> dict[bytes, int] | None:
    """Return a dictionary of barcode to raw reads."""
    if not per_barcode_metrics_path:
        return None

    reader = csv_utils.load_csv_filter_comments(
        per_barcode_metrics_path, "per barcode metrics", ["barcode", "raw_reads"]
    )
    return {ensure_binary(row["barcode"]): int(row["raw_reads"]) for row in reader}


def report_genomes(
    matrix: cr_matrix.CountMatrix,
    reads_summary,
    barcode_summary_h5_path: AnyStr,
    per_barcode_metrics_path: AnyStr,
    recovered_cells,
    sample_bc_seqs,
    cell_bc_seqs,
):
    """Report on all genomes in this matrix."""
    barcode_summary_h5 = h5.File(ensure_binary(barcode_summary_h5_path), "r")
    gex_reads_per_barcode = load_per_barcode_metrics(per_barcode_metrics_path)

    sample_bc_indices = (
        _get_barcode_summary_h5_indices(barcode_summary_h5, sample_bc_seqs)
        if sample_bc_seqs is not None
        else None
    )
    cell_bcs_union = reduce(lambda a, x: a | set(x), cell_bc_seqs.values(), set())
    cell_bcs_union_indices = matrix.bcs_to_ints(list(cell_bcs_union))
    bc_summary_cell_bc_union_indices = (
        cell_bcs_union_indices
        if sample_bc_seqs is None
        else _get_barcode_summary_h5_indices(barcode_summary_h5, cell_bcs_union)
    )

    metrics = {}

    genomes = matrix.get_genomes()
    # assert len(cell_bc_seqs) == len(genomes)

    # Compute genome-agnostic metrics
    feature_types = sorted(list({f.feature_type for f in matrix.feature_ref.feature_defs}))
    for ftype in feature_types:
        total_reads = _get_total_reads(reads_summary, ftype)
        if total_reads == 0:
            continue

        conf_mapped_reads = _get_conf_mapped_reads(reads_summary, genomes, ftype)

        submatrix = matrix.view().select_features_by_type(ftype)
        prefix = rna_library.get_library_type_metric_prefix(ftype)

        m = _report_genome_agnostic_metrics(
            ftype,
            submatrix,
            barcode_summary_h5,
            recovered_cells,
            sample_bc_indices,
            cell_bcs_union,
            cell_bcs_union_indices,
            bc_summary_cell_bc_union_indices,
            total_reads,
            gex_reads_per_barcode,
            conf_mapped_reads,
            prefix,
        )

        if rna_library.has_genomes(ftype):
            for genome in genomes:
                # Compute genome-specific metrics
                genome_matrix = matrix.view().select_features_by_genome(genome)
                cell_bc_indices = matrix.bcs_to_ints(cell_bc_seqs[genome])
                bc_summary_cell_bc_indices = (
                    cell_bc_indices
                    if sample_bc_seqs is None
                    else _get_barcode_summary_h5_indices(barcode_summary_h5, cell_bc_seqs[genome])
                )
                genome_summary = _report(
                    genome_matrix,
                    genome,
                    barcode_summary_h5,
                    cell_bc_indices,
                    bc_summary_cell_bc_indices,
                    sample_bc_indices,
                    prefix,
                )

                for key, value in genome_summary.items():
                    key = "_".join([genome, key])
                    m[key] = value
        else:
            # This feature has no genomes
            genome_summary = _report(
                submatrix,
                rna_library.MULTI_REFS_PREFIX,
                barcode_summary_h5,
                cell_bcs_union_indices,
                bc_summary_cell_bc_union_indices,
                sample_bc_indices,
                prefix,
            )
            for key, value in genome_summary.items():
                key = "_".join([rna_library.MULTI_REFS_PREFIX, key])
                m[key] = value

        # Prepend feature type to metric keys
        m_prefixed = {(prefix + k): v for k, v in m.items()}
        metrics.update(m_prefixed)

    return metrics
