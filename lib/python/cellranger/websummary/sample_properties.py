# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""These define information about the sample required to generate a web summary."""
# pylint: disable=too-few-public-methods,missing-docstring,too-many-arguments


from __future__ import annotations

import martian

from cellranger.version import get_version


class SampleProperties:
    def __init__(self, sample_id, sample_desc, version_from_git=False, throughput=None):
        self.sample_id = sample_id
        self.sample_desc = sample_desc
        self.throughput = throughput
        if not version_from_git:
            self.version = martian.get_pipelines_version()
        else:
            self.version = get_version()


class CountSampleProperties(SampleProperties):
    """Various versions of this class are passed around for Count, Aggr, Reanalyze, Spatial.

    web summaries, etc.
    """

    def __init__(
        self,
        sample_id,
        sample_desc,
        genomes,
        version_from_git=False,
        is_spatial=False,
        target_set=None,
        target_panel_summary=None,
        feature_ref_path=None,
        include_introns=False,
        reorientation_mode=None,
        filter_probes=None,
        throughput=None,
        aligner=None,
        redundant_loupe_alignment=False,
        loupe_alignment_file=None,
        v1_filtered_fbm=None,
        cmdline=None,
    ):
        super().__init__(
            sample_id, sample_desc, version_from_git=version_from_git, throughput=throughput
        )
        self.genomes = genomes
        self.is_spatial = is_spatial
        self.target_set = target_set
        self.target_panel_summary = target_panel_summary
        self.feature_ref_path = feature_ref_path
        self.include_introns = include_introns
        self.throughput = throughput
        self.reorientation_mode = reorientation_mode
        self.filter_probes = filter_probes
        self.aligner = aligner
        self.redundant_loupe_alignment = redundant_loupe_alignment
        self.loupe_alignment_file = loupe_alignment_file
        self.v1_filtered_fbm = v1_filtered_fbm
        self.cmdline = cmdline

    @property
    def is_targeted(self):
        return self.target_set is not None

    @property
    def is_lt(self):
        return False


class ExtendedCountSampleProperties(CountSampleProperties):
    """Properties for a count run."""

    def __init__(
        self,
        sample_id,
        sample_desc,
        genomes,
        barcode_whitelist,
        reference_path,
        target_set=None,
        target_panel_summary=None,
        feature_ref_path=None,
        version_from_git=False,
        is_spatial=False,
        include_introns=False,
        throughput=None,
        reorientation_mode=None,
        filter_probes=None,
        aligner=None,
        redundant_loupe_alignment=False,
        loupe_alignment_file=None,
        v1_filtered_fbm=None,
        cmdline=None,
    ):
        super().__init__(
            sample_id,
            sample_desc,
            genomes,
            version_from_git=version_from_git,
            is_spatial=is_spatial,
            target_set=target_set,
            target_panel_summary=target_panel_summary,
            feature_ref_path=feature_ref_path,
            include_introns=include_introns,
            throughput=throughput,
            reorientation_mode=reorientation_mode,
            filter_probes=filter_probes,
            aligner=aligner,
            redundant_loupe_alignment=redundant_loupe_alignment,
            loupe_alignment_file=loupe_alignment_file,
            v1_filtered_fbm=v1_filtered_fbm,
            cmdline=cmdline,
        )
        self.barcode_whitelist = barcode_whitelist
        self.reference_path = reference_path

    @property
    def is_lt(self):
        return self.barcode_whitelist == "9K-LT-march-2021"


class AggrCountSampleProperties(CountSampleProperties):
    """Properties from an Aggr Run."""

    def __init__(
        self,
        sample_id,
        sample_desc,
        genomes,
        agg_batches,
        is_spatial,
        target_set=None,
        target_panel_summary=None,
        feature_ref_path=None,
        version_from_git=False,
    ):
        super().__init__(
            sample_id,
            sample_desc,
            genomes,
            version_from_git=version_from_git,
            is_spatial=is_spatial,
            target_set=target_set,
            target_panel_summary=target_panel_summary,
            feature_ref_path=feature_ref_path,
            throughput=None,
        )
        self.agg_batches = agg_batches


class VdjSampleProperties(SampleProperties):
    def __init__(self, sample_id, sample_desc, chain_type, version_from_git=False):
        super().__init__(
            sample_id,
            sample_desc,
            version_from_git=version_from_git,
            throughput=None,  # TODO need to set throughput in summarize_vdj_reports if we support VDJ HT
        )
        self.chain_type = chain_type


class SampleDataPaths:  # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        summary_path=None,
        barcode_summary_path=None,
        analysis_path=None,
        filtered_barcodes_path=None,
        feature_metrics_path=None,
        antibody_histograms_path=None,
        antibody_treemap_path=None,
        raw_normalized_heatmap_path=None,
        isotype_scatter_path=None,
        gex_fbc_correlation_heatmap_path=None,
        antigen_histograms_path=None,
        antigen_treemap_path=None,
        vdj_clonotype_summary_path=None,
        vdj_barcode_support_path=None,
        vdj_cell_barcodes_path=None,
    ):
        assert filtered_barcodes_path is None or vdj_cell_barcodes_path is None
        self.summary_path = summary_path
        self.barcode_summary_path = barcode_summary_path
        self.analysis_path = analysis_path
        self.filtered_barcodes_path = filtered_barcodes_path
        self.feature_metrics_path = feature_metrics_path
        self.antibody_histograms_path = antibody_histograms_path
        self.antibody_treemap_path = antibody_treemap_path
        self.raw_normalized_heatmap_path = raw_normalized_heatmap_path
        self.isotype_scatter_path = isotype_scatter_path
        self.gex_fbc_correlation_heatmap_path = gex_fbc_correlation_heatmap_path
        self.antigen_histograms_path = antigen_histograms_path
        self.antigen_treemap_path = antigen_treemap_path
        self.vdj_clonotype_summary_path = vdj_clonotype_summary_path
        self.vdj_barcode_support_path = vdj_barcode_support_path
        self.vdj_cell_barcodes_path = vdj_cell_barcodes_path


class CellTypeDataPaths:  # pylint: disable=too-many-instance-attributes
    def __init__(self, summary_path=None, cell_type_bar_chart=None, cell_type_umap_plot=None):
        self.summary_path = summary_path
        self.cell_type_bar_chart = cell_type_bar_chart
        self.cell_type_umap_plot = cell_type_umap_plot
