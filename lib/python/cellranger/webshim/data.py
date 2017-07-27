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

class SampleDataPaths:
    def __init__(self, summary_path=None, barcode_summary_path=None, analysis_path=None,
                 vdj_clonotype_summary_path=None,
                 vdj_barcode_support_path=None):
        self.summary_path = summary_path
        self.barcode_summary_path = barcode_summary_path
        self.analysis_path = analysis_path
        self.vdj_clonotype_summary_path = vdj_clonotype_summary_path
        self.vdj_barcode_support_path = vdj_barcode_support_path

class SampleData:
    def __init__(self, sample_properties, sample_data_paths, plot_preprocess_func):
        self.summary = load_metrics_summary(sample_data_paths.summary_path)

        # Not guaranteed to exist for all pipelines
        self.num_cells = self.summary.get('filtered_bcs_transcriptome_union')

        self.analysis, self.original_cluster_sizes = load_analysis(sample_data_paths.analysis_path,
                                                                   plot_preprocess_func,
                                                                   sample_properties)

        self.barcode_summary = load_barcode_summary(sample_data_paths.barcode_summary_path)

        self.vdj_clonotype_summary = pd.read_csv(sample_data_paths.vdj_clonotype_summary_path) if sample_data_paths.vdj_clonotype_summary_path else None

        self.vdj_barcode_support = pd.read_csv(sample_data_paths.vdj_barcode_support_path, usecols=[1]) if sample_data_paths.vdj_barcode_support_path else None


def load_metrics_summary(summary_path):
    if summary_path is None:
        return None
    with open(summary_path) as f:
        return json.load(f)

def load_analysis(base_dir, plot_preprocess_func, sample_properties):
    """ Returns (analysis_object, original_cluster_sizes) """
    if base_dir is None:
        return None, None

    if len(sample_properties['genomes']) == 1:
        analysis = cr_sg_analysis.SingleGenomeAnalysis.load_default_format(base_dir)
    else:
        analysis = cr_mg_analysis.MultiGenomeAnalysis.load_default_format(base_dir)

    if analysis is None:
        return None, None

    # Subsample barcodes for plotting purposes
    original_cluster_sizes = None
    if analysis is not None and isinstance(analysis, cr_sg_analysis.SingleGenomeAnalysis):
        original_cluster_sizes = analysis.get_cluster_sizes()
        if analysis.matrix.bcs_dim > gex_constants.MAX_WEBSHIM_BCS_DIM:
            analysis.subsample_bcs(gex_constants.MAX_WEBSHIM_BCS_DIM)

    return plot_preprocess_func(analysis, sample_properties), original_cluster_sizes

def load_barcode_summary(h5_path):
    if h5_path is None:
        return None
    return h5py.File(h5_path, 'r')
