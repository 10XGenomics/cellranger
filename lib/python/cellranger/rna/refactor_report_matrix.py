#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Functions for refatoring report_matrix.py.

right now, computes only median confidently mapped reads per singlet per genome
"""
from __future__ import annotations

# pylint: disable=too-many-format-args
import pandas as pd

import cellranger.feature.utils as feature_utils
import cellranger.utils as cr_utils


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


def set_default_no_barcode_metrics():
    """Set default values for metrics."""
    metrics = {}
    # metrics["median_total_reads_per_singlet"] = 0
    return metrics


def compute_per_cell_metrics(
    filtered_barcodes_path,
    per_barcode_metrics_path,
    genomes,
    antibody_present=False,
    crispr_present=False,
    custom_present=False,
    is_targeted=False,
):
    """Right now, computes only median confidently mapped reads per singlet per genome."""
    # input validation
    input_files = [per_barcode_metrics_path, filtered_barcodes_path]
    input_files_present = feature_utils.all_files_present(input_files)
    if not input_files_present:
        raise ValueError("Per barcode metrics or filtered_barcodes CSV is not present")

    filtered_barcodes = [
        x.decode("utf8") for x in feature_utils.get_gex_cell_list(filtered_barcodes_path)
    ]
    num_cells = len(filtered_barcodes)

    try:
        per_barcode_metrics = pd.read_csv(per_barcode_metrics_path)
    except pd.errors.EmptyDataError:
        return set_default_no_barcode_metrics()

    if num_cells == 0:
        return set_default_no_cells(
            genomes=genomes,
            antibody_present=antibody_present,
            crispr_present=crispr_present,
            custom_present=custom_present,
        )

    metrics = {}

    filtered_barcode_indices = per_barcode_metrics["barcode"].isin(filtered_barcodes)
    metrics["median_total_reads_per_singlet"] = per_barcode_metrics.loc[
        filtered_barcode_indices, "raw_reads"
    ].median()
    del filtered_barcodes

    metrics.update(
        _compute_per_genome_metrics(
            filtered_barcodes_path, genomes, is_targeted, per_barcode_metrics
        )
    )

    return metrics


def _compute_per_genome_metrics(filtered_barcodes_path, genomes, is_targeted, per_barcode_metrics):
    metrics = {}
    genome_to_barcodes = cr_utils.load_barcode_csv(filtered_barcodes_path)

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
        genome_barcodes = [x.decode("utf8") for x in genome_to_barcodes[genome.encode()]]
        barcode_metric_keys = [
            f"{barcode_metric_key}_{genome}" for barcode_metric_key in genome_metrics
        ]
        metric_values = per_barcode_metrics.loc[
            per_barcode_metrics["barcode"].isin(genome_barcodes), barcode_metric_keys
        ].median()

        for barcode_metric_key, metric_key in genome_metrics.items():
            metric_key = f"{genome}_{metric_key}"
            metrics[metric_key] = metric_values[f"{barcode_metric_key}_{genome}"]

    return metrics
