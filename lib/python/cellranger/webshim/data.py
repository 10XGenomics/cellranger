#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Manage data used by the webshim


from __future__ import annotations

import collections
import json
import os
from collections.abc import Sequence
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.analysis.multigenome as cr_mg_analysis
import cellranger.analysis.singlegenome as cr_sg_analysis
import cellranger.utils as cr_utils
import cellranger.vdj.utils as vdj_utils
import cellranger.webshim.constants.gex as gex_constants

if TYPE_CHECKING:
    from cellranger.analysis.singlegenome import Projection

from cellranger.websummary.sample_properties import (
    AggrCountSampleProperties,
    CountSampleProperties,
    ExtendedCountSampleProperties,
    SampleProperties,
)

FILTERED_BCS_TRANSCRIPTOME_UNION = "filtered_bcs_transcriptome_union"

BarcodeRankPlotSegment = collections.namedtuple(
    "BarcodeRankPlotSegment", ["start", "end", "cell_density", "legend"]
)


def get_plot_segment(start_index, end_index, sorted_bc, cell_barcodes: set[bytes], legend=False):
    """Helper function to build a plot segment."""
    assert end_index > start_index
    num_cells = sum(1 for i in range(start_index, end_index) if sorted_bc[i] in cell_barcodes)
    density = float(num_cells) / float(end_index - start_index)
    return BarcodeRankPlotSegment(
        start=start_index, end=end_index, cell_density=density, legend=legend
    )


def generate_counter_barcode_rank_plot_data(
    cell_barcodes: set[bytes], barcode_summary, key, restrict_barcodes=None
):
    """A helper function to generate the data required to generate the barcode rank.

    plot for RNA counter. The barcode rank plot consists of
    multiple plot segments

    If the barcodes are ordered by UMI counts
        - All the barcodes until the first non cell barcodes go into
        the `Cell` segment.
        - All the barcodes beyond the last cell barcode goes into the
        `Background` segment.
        - The remaining barcodes are further divided into multiple plot
        segments
    Input:
        - cell_barcodes: A set of bytes with the cell barcodes
        - barcode_summary: The barcode summary from the h5.
        - key: Specifies the entry to look up in barcode summary
            for getting the UMI counts
        - restrict_barcodes: Optional list of cell barcodes to restrict to
    Output:
        - sorted_counts: UMI counts sorted in descending order
        - plot_segments: List of BarcodeRankPlotSegment
    """
    if restrict_barcodes:
        restrict_indices = np.nonzero(np.isin(barcode_summary["bc_sequence"][:], restrict_barcodes))
        counts_per_bc = barcode_summary[key][:][restrict_indices]
        barcode_sequences = barcode_summary["bc_sequence"][:][restrict_indices]
    else:
        counts_per_bc = barcode_summary[key][:]
        barcode_sequences = barcode_summary["bc_sequence"][:]

    srt_order = compute_sort_order(counts_per_bc, barcode_sequences, cell_barcodes)
    sorted_bc = barcode_sequences[srt_order]
    sorted_counts = counts_per_bc[srt_order]
    del srt_order

    plot_segments = compute_plot_segments(sorted_bc, sorted_counts, cell_barcodes)

    return sorted_counts, plot_segments


class SampleData:
    def __init__(
        self,
        sample_properties,
        sample_data_paths,
        plot_preprocess_func,
        projections: Sequence[Projection],
    ):
        if sample_properties:  # Can be None
            assert isinstance(sample_properties, SampleProperties)
        self.summary = load_metrics_summary(sample_data_paths.summary_path)
        # Not guaranteed to exist for all pipelines
        self.num_cells = self.summary.get("filtered_bcs_transcriptome_union")

        self.analyses, self.original_cluster_sizes = load_analyses(
            sample_data_paths.analysis_path, plot_preprocess_func, sample_properties, projections
        )

        self.barcode_summary = load_barcode_summary(sample_data_paths.barcode_summary_path)

        self.has_target_set = (
            hasattr(sample_properties, "target_set") and sample_properties.target_set is not None
        )

        self.targeting_method = self.summary.get("targeting_method")

        self.is_visium_hd = (
            sample_properties.is_visium_hd
            if isinstance(sample_properties, ExtendedCountSampleProperties)
            else False
        )

        self.filter_probes = (
            sample_properties.filter_probes
            if isinstance(sample_properties, ExtendedCountSampleProperties)
            else None
        )

        self.vdj_clonotype_summary = (
            pd.read_csv(ensure_str(sample_data_paths.vdj_clonotype_summary_path))
            if sample_data_paths.vdj_clonotype_summary_path
            else None
        )

        self.vdj_barcode_support = (
            pd.read_csv(
                ensure_str(sample_data_paths.vdj_barcode_support_path),
                converters={"barcode": ensure_binary},
            )
            if sample_data_paths.vdj_barcode_support_path
            and os.path.getsize(sample_data_paths.vdj_barcode_support_path) != 0
            else None
        )

        self.cell_barcodes = (
            cr_utils.get_cell_associated_barcode_set(sample_data_paths.filtered_barcodes_path)
            if sample_data_paths.filtered_barcodes_path
            else None
        )

        if sample_data_paths.vdj_cell_barcodes_path:
            assert self.cell_barcodes is None
            self.cell_barcodes = {
                ensure_binary(bc)
                for bc in vdj_utils.load_cell_barcodes_json(
                    sample_data_paths.vdj_cell_barcodes_path
                )
            }

        self.feature_metrics = load_feature_counts(sample_data_paths.feature_metrics_path)
        self.antibody_histograms = load_antibody_data(sample_data_paths.antibody_histograms_path)
        self.antibody_treemap = load_treemap_data(sample_data_paths.antibody_treemap_path)
        self.antigen_histograms = load_antibody_data(sample_data_paths.antigen_histograms_path)
        self.raw_normalized_heatmap = load_antibody_data(
            sample_data_paths.raw_normalized_heatmap_path
        )
        self.isotype_scatter = load_antibody_data(sample_data_paths.isotype_scatter_path)
        self.gex_fbc_correlation_heatmap = load_antibody_data(
            sample_data_paths.gex_fbc_correlation_heatmap_path
        )

    def get_analysis(self, analysis_type):
        if self.analyses is None:
            return None
        r = [a for a in self.analyses if isinstance(a, analysis_type)]
        assert len(r) == 0 or len(r) == 1
        if len(r) == 0:
            return None
        else:
            return r[0]

    def counter_barcode_rank_plot_data(self, key):
        """Generate the data required to generate the barcode rank.

        plot for RNA counter. The barcode rank plot consists of
        multiple plot segments

        If the barcodes are ordered by UMI counts
            - All the barcodes until the first non cell barcodes go into
              the `Cell` segment.
            - All the barcodes beyond the last cell barcode goes into the
              `Background` segment.
            - The remaining barcodes are further divided into multiple plot
              segments
        Input:
            - key: Specifies the entry to look up in barcode summary
                   for getting the UMI counts
        Output:
            - sorted_counts: UMI counts sorted in descending order
            - plot_segments: List of BarcodeRankPlotSegment
        """
        assert self.cell_barcodes is not None
        return generate_counter_barcode_rank_plot_data(
            self.cell_barcodes, self.barcode_summary, key
        )

    def vdj_barcode_rank_plot_data(self):
        """Generate the data required to generate the barcode rank plot for VDJ."""
        assert self.cell_barcodes is not None
        counts_per_bc = self.vdj_barcode_support["count"].to_numpy()
        srt_order = compute_sort_order(
            counts_per_bc, self.vdj_barcode_support["barcode"], self.cell_barcodes
        )
        sorted_bc = self.vdj_barcode_support["barcode"].to_numpy()[srt_order]
        sorted_counts = counts_per_bc[srt_order]
        del srt_order

        plot_segments = compute_plot_segments(sorted_bc, sorted_counts, self.cell_barcodes)
        return sorted_counts, plot_segments

    def is_targeted(self):
        return self.has_target_set

    def is_antibody_only(self):
        """Returns true if data is antibody only, currently determined by the presence of a metric.

        Returns:
            bool: Whether or not the data is antibody-only.
        """
        return self.summary is not None and FILTERED_BCS_TRANSCRIPTOME_UNION not in self.summary


def compute_sort_order(counts_per_bc, bc_sequences, cell_barcodes: set[bytes]):
    assert len(counts_per_bc) == len(bc_sequences)
    is_cell = np.full(len(counts_per_bc), False)
    for i, bc in enumerate(bc_sequences):
        is_cell[i] = bc in cell_barcodes
    srt_order = sorted(
        range(len(counts_per_bc)), key=lambda x: (counts_per_bc[x], is_cell[x]), reverse=True
    )
    return srt_order


def compute_plot_segments(sorted_bc, sorted_counts, cell_barcodes: set[bytes]):
    # find the first barcode which is not a cell
    first_non_cell = len(sorted_bc)
    for i, bc in enumerate(sorted_bc):
        if bc not in cell_barcodes:
            first_non_cell = i
            break

    # find the last barcode which is a cell
    last_cell = 0
    for i in reversed(range(len(sorted_bc))):
        if sorted_bc[i] in cell_barcodes:
            last_cell = i
            break

    ranges = [0, first_non_cell, last_cell + 1, len(sorted_bc)]

    plot_segments = []
    plot_segments.append(
        BarcodeRankPlotSegment(start=0, end=ranges[1], cell_density=1.0, legend=True)
    )
    plot_segments.append(
        BarcodeRankPlotSegment(start=ranges[2], end=ranges[3], cell_density=0.0, legend=True)
    )

    # Subdivide the mixed section
    mixed_segments = segment_log_plot_by_length(sorted_counts, ranges[1], ranges[2])
    for i in range(len(mixed_segments) - 1):
        plot_segments.append(
            get_plot_segment(
                mixed_segments[i], mixed_segments[i + 1], sorted_bc, cell_barcodes, legend=False
            )
        )

    return plot_segments


def segment_log_plot_by_length(y_data, x_start, x_end):
    """Given the extends of the mixed region [x_start, x_end), compute.

    the x-indices that would divide the plot into segments of a
    prescribed length (in pixel coordinates) with a minimum number
    of barcodes in each segment
    """
    if x_end <= x_start:
        return []

    SEGMENT_NORMALIZED_MAX_LEN = 0.02
    MIN_X_SPAN = 20

    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))

    this_segment_len = 0.0
    segment_idx = [x_start]

    for i in range(x_start, x_end):
        if i == 0:
            continue
        last_i = max(x_start, i - 1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        this_segment_len += np.linalg.norm([dx, dy])
        if this_segment_len >= SEGMENT_NORMALIZED_MAX_LEN and i > (segment_idx[-1] + MIN_X_SPAN):
            segment_idx.append(i + 1)
            this_segment_len = 0.0

    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)

    return segment_idx


def load_metrics_summary(summary_path):
    if summary_path is None:
        return None
    with open(summary_path) as f:
        return json.load(f)


def load_feature_counts(feature_metrics_path):
    if feature_metrics_path is None:
        return None
    return SummaryFeatureCounts.from_file(feature_metrics_path)


def load_antibody_data(antibody_histograms_path):
    if antibody_histograms_path is None:
        return None
    with open(antibody_histograms_path) as f:
        return json.load(f)


def load_treemap_data(antibody_treemap_path):
    if antibody_treemap_path is None:
        return None
    with open(antibody_treemap_path) as f:
        return json.load(f)


def load_analyses(
    base_dir,
    plot_preprocess_func,
    sample_properties,
    projections: Sequence[Projection],
):
    """Returns (analysis_object, original_cluster_sizes)."""
    if base_dir is None:
        return None, None

    if isinstance(sample_properties, CountSampleProperties):
        if len(sample_properties.genomes) == 1:
            analyses = [
                cr_sg_analysis.SingleGenomeAnalysis.load_default_format(
                    base_dir, method="pca", projections=projections
                )
            ]
        else:
            analyses = [
                cr_sg_analysis.SingleGenomeAnalysis.load_default_format(
                    base_dir, method="pca", projections=projections
                ),
                cr_mg_analysis.MultiGenomeAnalysis.load_default_format(base_dir),
            ]

    if analyses[0] is None:
        return None, None

    # Subsample barcodes for plotting purposes
    original_cluster_sizes = None
    analysis = analyses[0]
    if analysis is not None and isinstance(analysis, cr_sg_analysis.SingleGenomeAnalysis):
        original_cluster_sizes = analysis.get_cluster_sizes()
        is_single_cell = (
            isinstance(sample_properties, CountSampleProperties)
            and not sample_properties.is_spatial
        )
        if analysis.matrix.bcs_dim > gex_constants.MAX_WEBSHIM_BCS_DIM and (
            is_single_cell or isinstance(sample_properties, AggrCountSampleProperties)
        ):
            # Downsample the barcodes to 10k if the sample is count single-cell or aggr even if spatial.
            # Don't want to do this in spatial because if affects the spatial plots in XL slides.
            analysis.subsample_bcs(gex_constants.MAX_WEBSHIM_BCS_DIM)

    return plot_preprocess_func(analyses), original_cluster_sizes


def load_barcode_summary(h5_path):
    if h5_path is None:
        return None
    return h5py.File(h5_path, "r")


class SummaryFeatureCounts:
    """Stores and queries a dataframe that summarizes relevant features per gene."""

    def __init__(self, df: pd.DataFrame, filepath=None):
        """Stores and queries a dataframe that summarizes relevant features per gene.

        There should be one row per gene and any number of columns.

        Args:
            df (pandas dataframe): Dataframe with n_genes rows and columns with features per genes.
        """
        self.df = df
        # for easy querying
        self.feature_dict = df.set_index("feature_id").to_dict(orient="index")
        self.filepath = filepath

    @classmethod
    def from_file(cls, path) -> SummaryFeatureCounts:
        """Loads the features from a csv file."""
        path = ensure_str(path)
        assert os.path.isfile(path)
        return cls(pd.read_csv(ensure_str(path)), path)

    def to_csv(self, path):
        path = ensure_str(path)
        self.df.to_csv(path, index=False)

    def get_df(self) -> pd.DataFrame:
        return self.df

    def get_value_for_feature(self, feature_id, key):
        """Get the given key for the given feature."""
        if feature_id not in self.feature_dict:
            return None
        if key not in self.feature_dict[feature_id]:
            return None
        return self.feature_dict[feature_id][key]
