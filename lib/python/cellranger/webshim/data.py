#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Manage data used by the webshim

import h5py
import json
import pandas as pd

import cellranger.analysis.singlegenome as cr_sg_analysis
import cellranger.analysis.multigenome as cr_mg_analysis
import cellranger.webshim.constants.gex as gex_constants
import cellranger.utils as cr_utils

import collections
import numpy as np

BarcodeRankPlotSegment = collections.namedtuple('BarcodeRankPlotSegment', ['start', 'end', 'cell_density', 'legend'])

def get_plot_segment(start_index, end_index, sorted_bc, cell_barcodes, legend=False):
    """
    Helper function to build a plot segment.
    """
    assert end_index > start_index
    num_cells = sum([1 for i in xrange(start_index, end_index) if sorted_bc[i] in cell_barcodes])
    density = float(num_cells)/float(end_index-start_index)
    return BarcodeRankPlotSegment(start=start_index, end=end_index, cell_density=density, legend=legend)

class SampleDataPaths:
    def __init__(self, summary_path=None, barcode_summary_path=None, analysis_path=None,
                 vdj_clonotype_summary_path=None,
                 vdj_barcode_support_path=None,
                 filtered_barcodes_path=None):
        self.summary_path = summary_path
        self.barcode_summary_path = barcode_summary_path
        self.analysis_path = analysis_path
        self.vdj_clonotype_summary_path = vdj_clonotype_summary_path
        self.vdj_barcode_support_path = vdj_barcode_support_path
        self.filtered_barcodes_path = filtered_barcodes_path

class SampleData:
    def __init__(self, sample_properties, sample_data_paths, plot_preprocess_func):
        self.summary = load_metrics_summary(sample_data_paths.summary_path)

        # Not guaranteed to exist for all pipelines
        self.num_cells = self.summary.get('filtered_bcs_transcriptome_union')

        self.analyses, self.original_cluster_sizes = load_analyses(sample_data_paths.analysis_path,
                                                                   plot_preprocess_func,
                                                                   sample_properties)

        self.barcode_summary = load_barcode_summary(sample_data_paths.barcode_summary_path)

        self.vdj_clonotype_summary = pd.read_csv(sample_data_paths.vdj_clonotype_summary_path) if sample_data_paths.vdj_clonotype_summary_path else None

        self.vdj_barcode_support = pd.read_csv(sample_data_paths.vdj_barcode_support_path, usecols=[1]) if sample_data_paths.vdj_barcode_support_path else None

        self.cell_barcodes = cr_utils.get_cell_associated_barcode_set(sample_data_paths.filtered_barcodes_path) if sample_data_paths.filtered_barcodes_path else None

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
        """
        Generate the data required to generate the barcode rank
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
        counts_per_bc = self.barcode_summary[key][:]
        srt_order = sorted(range(len(counts_per_bc)), key=lambda x: counts_per_bc[x], reverse=True)
        sorted_bc = self.barcode_summary['bc_sequence'][:][srt_order]
        sorted_counts = counts_per_bc[srt_order]
        del srt_order
 
        # find the first barcode which is not a cell
        first_non_cell = len(sorted_bc)
        for i, bc in enumerate(sorted_bc):
            if bc not in self.cell_barcodes:
                first_non_cell = i
                break

        # find the last barcode which is a cell
        last_cell = 0 
        for i in reversed(xrange(len(sorted_bc))):
            if sorted_bc[i] in self.cell_barcodes:
                last_cell = i
                break

        ranges = [0, first_non_cell, last_cell+1, len(sorted_bc)]

        plot_segments = []
        plot_segments.append(BarcodeRankPlotSegment(start=0, end=ranges[1], cell_density=1.0, legend=True))
        plot_segments.append(BarcodeRankPlotSegment(start=ranges[2], end=ranges[3], cell_density=0.0, legend=True))

        # Subdivide the mixed section
        mixed_segments = segment_log_plot_by_length(sorted_counts, ranges[1], ranges[2])
        for i in xrange(len(mixed_segments)-1):
            plot_segments.append(get_plot_segment(mixed_segments[i], mixed_segments[i+1], sorted_bc, self.cell_barcodes, legend=False))

        return sorted_counts, plot_segments

def segment_log_plot_by_length(y_data, x_start, x_end):
    """
    Given the extends of the mixed region [x_start, x_end), compute
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

    for i in xrange(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        this_segment_len += np.linalg.norm([dx, dy])
        if this_segment_len >= SEGMENT_NORMALIZED_MAX_LEN and i > (segment_idx[-1] + MIN_X_SPAN):
            segment_idx.append(i+1)
            this_segment_len = 0.0

    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)

    return segment_idx

def load_metrics_summary(summary_path):
    if summary_path is None:
        return None
    with open(summary_path) as f:
        return json.load(f)

def load_analyses(base_dir, plot_preprocess_func, sample_properties):
    """ Returns (analysis_object, original_cluster_sizes) """
    if base_dir is None:
        return None, None

    if len(sample_properties['genomes']) == 1:
        analyses = [cr_sg_analysis.SingleGenomeAnalysis.load_default_format(base_dir, 'pca')]
    else:
        analyses = [cr_sg_analysis.SingleGenomeAnalysis.load_default_format(base_dir, 'pca'),
                    cr_mg_analysis.MultiGenomeAnalysis.load_default_format(base_dir)]

    if analyses[0] is None:
        return None, None

    # Subsample barcodes for plotting purposes
    original_cluster_sizes = None
    analysis = analyses[0]
    if analysis is not None and isinstance(analysis, cr_sg_analysis.SingleGenomeAnalysis):
        original_cluster_sizes = analysis.get_cluster_sizes()
        if analysis.matrix.bcs_dim > gex_constants.MAX_WEBSHIM_BCS_DIM:
            analysis.subsample_bcs(gex_constants.MAX_WEBSHIM_BCS_DIM)

    return plot_preprocess_func(analyses, sample_properties), original_cluster_sizes

def load_barcode_summary(h5_path):
    if h5_path is None:
        return None
    return h5py.File(h5_path, 'r')
