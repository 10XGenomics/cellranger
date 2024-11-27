#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Functions for refatoring report_matrix.py.

right now, computes only median confidently mapped reads per singlet per genome
"""
from __future__ import annotations

from typing import TYPE_CHECKING

# pylint: disable=too-many-format-args
import pandas as pd

from cellranger.constants import NO_BARCODE

if TYPE_CHECKING:
    from cellranger.fast_utils import (  # pylint: disable=no-name-in-module,unused-import
        FilteredBarcodes,
    )


def set_default_no_cells(genomes, antibody_present, crispr_present, custom_present):
    """Set default values for metrics."""
    metrics = {}
    metrics["median_total_reads_per_singlet"] = 0

    # rna/report_matrix is a convoluted mess, so handling edge case here
    # for some metrics computed there

    for genome in genomes:
        metrics[f"{genome}_filtered_bcs_median_unique_genes_detected"] = 0
        metrics[f"{genome}_filtered_bcs_mean_unique_genes_detected"] = 0
        metrics[f"{genome}_filtered_bcs_total_unique_genes_detected"] = 0
        metrics[f"{genome}_filtered_bcs_median_counts"] = 0
        metrics[f"{genome}_filtered_bcs_mean_counts"] = 0

    if antibody_present:
        metrics["ANTIBODY_multi_filtered_bcs"] = 0
        metrics["ANTIBODY_filtered_bcs_transcriptome_union"] = 0
        metrics["ANTIBODY_multi_filtered_bcs_median_counts"] = 0
        metrics["ANTIBODY_multi_filtered_bcs_mean_counts"] = 0
        metrics["ANTIBODY_multi_usable_reads_per_filtered_bc"] = 0

    if crispr_present:
        metrics["CRISPR_multi_filtered_bcs"] = 0
        metrics["CRISPR_multi_filtered_bcs_median_counts"] = 0
        metrics["CRISPR_multi_filtered_bcs_mean_counts"] = 0
        metrics["CRISPR_multi_usable_reads_per_filtered_bc"] = 0

    if custom_present:
        metrics["Custom_multi_filtered_bcs"] = 0
        metrics["Custom_multi_filtered_bcs_median_counts"] = 0
        metrics["Custom_multi_filtered_bcs_mean_counts"] = 0
        metrics["Custom_multi_usable_reads_per_filtered_bc"] = 0
    return metrics


def is_valid_barcode(bc: str | bytes):
    if isinstance(bc, str):
        bc = bc.encode()
    return bc != NO_BARCODE


def compute_per_cell_metrics(
    filtered_barcodes: FilteredBarcodes,
    per_barcode_metrics_path,
    genomes,
    antibody_present=False,
    crispr_present=False,
    custom_present=False,
    is_targeted=False,
):
    """Right now, computes only median confidently mapped reads per singlet per genome."""
    # input validation
    if not per_barcode_metrics_path:
        return {}
    per_barcode_metrics = pd.read_csv(per_barcode_metrics_path)

    num_cells = filtered_barcodes.num_cells()

    if num_cells == 0:
        return set_default_no_cells(
            genomes=genomes,
            antibody_present=antibody_present,
            crispr_present=crispr_present,
            custom_present=custom_present,
        )

    metrics = {}

    filtered_barcode_indices = per_barcode_metrics["barcode"].apply(
        lambda bc: is_valid_barcode(bc) and filtered_barcodes.contains(bc)
    )
    metrics["median_total_reads_per_singlet"] = per_barcode_metrics.loc[
        filtered_barcode_indices, "raw_reads"
    ].median()

    metrics.update(
        _compute_per_genome_metrics(filtered_barcodes, genomes, is_targeted, per_barcode_metrics)
    )

    return metrics


def _compute_per_genome_metrics(
    filtered_barcodes: FilteredBarcodes, genomes, is_targeted, per_barcode_metrics
):
    metrics = {}

    # Dictionary from column names in the per_barcode_metrics to metric json keys
    genome_metrics = {
        "mapped_reads": "median_reads_per_singlet",
        "conf_reads": "median_conf_reads_per_singlet",
    }
    if is_targeted:
        genome_metrics.update(
            {
                "ontarget_reads": "median_reads_per_singlet_ontarget",
                "offtarget_reads": "median_reads_per_singlet_offtarget",
                "conf_ontarget_reads": "median_conf_reads_per_singlet_ontarget",
                "conf_offtarget_reads": "median_conf_reads_per_singlet_offtarget",
            }
        )

    for genome in genomes:
        # must get the per-barcode metrics for only barcodes from this genome
        barcode_metric_keys = [
            f"{barcode_metric_key}_{genome}" for barcode_metric_key in genome_metrics
        ]

        filtered_rows = per_barcode_metrics["barcode"].apply(
            lambda bc, genome=genome: is_valid_barcode(bc)
            and filtered_barcodes.contains(bc, genome)
        )
        metric_values = per_barcode_metrics.loc[
            filtered_rows,
            barcode_metric_keys,
        ].median()

        for barcode_metric_key, metric_key in genome_metrics.items():
            metric_key = f"{genome}_{metric_key}"
            metrics[metric_key] = metric_values.get(f"{barcode_metric_key}_{genome}", 0)

    return metrics
