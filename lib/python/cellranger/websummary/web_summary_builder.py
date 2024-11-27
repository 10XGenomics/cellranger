#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
# pylint: disable=too-few-public-methods


from __future__ import annotations

import copy
from typing import TYPE_CHECKING

import cellranger.rna.library as rna_library
import cellranger.websummary.sample_properties as wsp
import cellranger.websummary.violin_plots as violin_plots
from cellranger.analysis.singlegenome import TSNE_NAME, UMAP_NAME, Projection
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC, TARGETING_METHOD_TL
from cellranger.webshim.common import load_sample_data
from cellranger.webshim.constants.shared import CELLRANGER_COMMAND_NAME, PIPELINE_AGGR
from cellranger.websummary.analysis_tab_aux import (
    antibody_histogram_plot,
    antibody_qc_plot,
    antibody_treemap_plot,
    antigen_histogram_plot,
    barnyard_table,
    gdna_table,
    median_gene_plot,
    seq_saturation_plot,
    targeted_table,
)
from cellranger.websummary.analysis_tab_core import analysis_by_clustering, umi_on_projection_plot
from cellranger.websummary.metrics import (
    INFO_THRESHOLD,
    LTMetricAnnotations,
    MetricAnnotations,
    TargetedMetricAnnotations,
)
from cellranger.websummary.react_components import WebSummaryData
from cellranger.websummary.react_summarize import write_html_file
from cellranger.websummary.sample_properties import CountSampleProperties, SampleDataPaths
from cellranger.websummary.summary_tab import (
    ANTIBODY_CELL_CALLING_METRIC_KEYS,
    CELL_CALLING_ALARM_KEYS,
    CELL_CALLING_METRIC_KEYS,
    FEATURE_BARCODINGS,
    FILTERED_BCS_TRANSCRIPTOME_UNION,
    TARGETED_CELL_CALLING_ALARM_KEYS_HC,
    TARGETED_CELL_CALLING_METRIC_KEYS_HC,
    TARGETED_CELL_CALLING_METRIC_KEYS_TL,
    add_data,
    aggregation_by_sample_table,
    batch_correction_table,
    cell_or_spot_calling_table,
    feature_barcode_aggregation_table,
    feature_barcode_application_table,
    feature_barcode_sequencing_table,
    hero_metrics,
    mapping_table,
    normalization_summary_table,
    pipeline_info_table,
    sequencing_table,
    summary_image_table,
)

if TYPE_CHECKING:
    from cellranger.webshim.data import SampleData
    from cellranger.websummary.react_components import BarnyardPanel


def _build_summary_tab_common(
    metadata,
    sample_data: SampleData,
    species_list,
    sample_properties,
    pipeline,
    ws_data,
    zoom_images,
    regist_images,
    sample_defs=None,
):
    """Code to create summary tab shared by spatial and single cell, updates ws_data."""
    alarm_list = ws_data.alarms
    summary_tab = ws_data.summary_tab

    add_data(summary_tab, alarm_list, hero_metrics(metadata, sample_data, species_list))
    add_data(
        summary_tab,
        alarm_list,
        pipeline_info_table(
            sample_data,
            sample_properties,
            pipeline,
            metadata,
            species_list,
            sample_defs=sample_defs,
        ),
    )
    if pipeline != PIPELINE_AGGR:
        add_data(
            summary_tab,
            alarm_list,
            sequencing_table(
                metadata,
                sample_data,
                species_list,
                is_targeted=sample_data.is_targeted(),
                is_hd=sample_data.is_visium_hd,
            ),
        )
        add_data(
            summary_tab,
            alarm_list,
            mapping_table(metadata, sample_data, species_list),
        )
    # TODO: Detecting antibody-only/GEX-Less was done by two different conditions in two different
    # stages before.  In `count` the check was if 'cells' not in web_sum_data.summary_tab
    # and in the check was `aggr` was if len(species_list) == 0:
    # currently trying the condition below which should count for both.

    # Default for Single Cell or Spatial with or without antibody.
    metric_keys = (
        ANTIBODY_CELL_CALLING_METRIC_KEYS
        if sample_data.is_antibody_only()
        else CELL_CALLING_METRIC_KEYS
    )
    alarm_keys = CELL_CALLING_ALARM_KEYS

    # Targeted hybridcap or RTL with or without antibody and not aggr.
    # The aggr pipeline doesn't generate targeted metrics.
    if sample_data.is_targeted() and not isinstance(
        sample_properties, wsp.AggrCountSampleProperties
    ):
        if sample_data.targeting_method == TARGETING_METHOD_HC:
            metric_keys = TARGETED_CELL_CALLING_METRIC_KEYS_HC + ANTIBODY_CELL_CALLING_METRIC_KEYS
            alarm_keys = TARGETED_CELL_CALLING_ALARM_KEYS_HC
        elif sample_data.targeting_method == TARGETING_METHOD_TL:
            metric_keys = TARGETED_CELL_CALLING_METRIC_KEYS_TL + ANTIBODY_CELL_CALLING_METRIC_KEYS
            alarm_keys = CELL_CALLING_ALARM_KEYS

    add_data(
        summary_tab,
        alarm_list,
        cell_or_spot_calling_table(
            metadata,
            sample_data,
            sample_properties,
            species_list,
            metric_keys,
            alarm_keys,
        ),
    )

    if sample_properties.is_spatial and not isinstance(
        sample_properties, wsp.AggrCountSampleProperties
    ):
        add_data(
            summary_tab,
            alarm_list,
            summary_image_table(
                metadata,
                sample_data,
                species_list,
                metric_keys=[],
                alarm_keys=[],
                zoom_images=zoom_images,
                regist_images=regist_images,
            ),
        )
    # Add the command line args to the summary tab
    if sample_properties.cmdline:
        summary_tab["cmdline"] = {
            "help": {
                "data": [
                    [
                        "",
                        [
                            f"<span style='font-size: 18px;'><code><pre style='white-space: pre-wrap;'>{sample_properties.cmdline!s}</pre></code></span>"
                        ],
                    ]
                ],
                "title": "Command Line Arguments",
            },
        }


def _add_gdna_to_websummary(ws_data: WebSummaryData, g_table: BarnyardPanel) -> None:
    """Add the gDNA to the websummary, storing the arrays as shared resources to be unpacked.

    Args:
        ws_data: web summary data
        g_table: a gDNA table

    Returns:
        None
    """
    data = g_table.plot["data"]
    # Convert the x,y arrays for the plots into shared resources to save disk space.
    for i in range(2):
        cdata = data[i]
        for axis in ["x", "y"]:
            cdata[axis] = ws_data.shared_resources.add_shared_resource(cdata[axis])
    ws_data.analysis_tab["gdna"] = g_table


def _build_analysis_tab_common(
    metadata,
    sample_data,
    pipeline,
    species_list,
    sample_properties,
    ws_data,
    *,
    projection: Projection,
):
    """Code to create analysis tab shared by spatial and single cell."""
    alarm_list = ws_data.alarms
    analysis_tab = ws_data.analysis_tab
    add_data(
        analysis_tab,
        alarm_list,
        seq_saturation_plot(sample_data, sample_properties),
    )
    add_data(
        analysis_tab,
        alarm_list,
        median_gene_plot(sample_data, sample_properties, species_list),
    )
    ws_data.clustering_selector = analysis_by_clustering(
        sample_data, spatial=sample_properties.is_spatial, projection=projection
    )
    # gDNA
    g_table = gdna_table(metadata, sample_data, species_list)
    if g_table is not None:
        _add_gdna_to_websummary(ws_data, g_table)

    violin_plot = None
    if sample_properties.is_spatial:
        violin_plot = violin_plots.make_gene_umi_violin_plots(
            sample_data,
            library_type=GENE_EXPRESSION_LIBRARY_TYPE,
            is_spatial=sample_properties.is_spatial,
        )
    if violin_plot:
        add_data(
            analysis_tab,
            alarm_list,
            violin_plot,
        )


def _build_antibody_tab_common(sample_data, sample_properties, ws_data, *, projection: Projection):
    """Code to create antibody tab shared by spatial and single cell."""
    alarm_list = ws_data.alarms
    antibody_tab = ws_data.antibody_tab
    add_data(
        antibody_tab,
        alarm_list,
        antibody_histogram_plot(sample_data),
    )
    add_data(
        antibody_tab,
        alarm_list,
        antibody_treemap_plot(sample_data),
    )
    add_data(
        antibody_tab,
        alarm_list,
        antibody_qc_plot(sample_data, "gex_fbc_correlation_heatmap"),
    )
    add_data(
        antibody_tab,
        alarm_list,
        antibody_qc_plot(sample_data, "isotype_scatter"),
    )
    # Raw and isotype normalized count ab correlation matrix plot
    add_data(
        antibody_tab,
        alarm_list,
        antibody_qc_plot(sample_data, "raw_normalized_heatmap"),
    )

    ws_data.antibody_clustering_selector = analysis_by_clustering(
        sample_data,
        spatial=sample_properties.is_spatial,
        library_type=rna_library.ANTIBODY_LIBRARY_TYPE,
        projection=projection,
    )


def _build_antigen_tab_common(sample_data, sample_properties, ws_data):
    """Code to create antibody tab shared by spatial and single cell."""
    alarm_list = ws_data.alarms
    antigen_tab = ws_data.antigen_tab
    add_data(
        antigen_tab,
        alarm_list,
        antigen_histogram_plot(sample_data),
    )


def _build_diagnostic_values(sample_data, ws_data):
    """Some metrics are useful for debugging but not yet ready to be displayed to the customer,.

    we hide these in the web summary json inside a "diagnostics" field if we encounter them.
    """
    # assert isinstance(sample_data, SampleData)
    assert isinstance(ws_data, WebSummaryData)
    diagnostic_metrics = [
        "tso_frac",
        "i1_bases_with_q30_frac",
        "i2_bases_with_q30_frac",
        "low_support_umi_reads_frac",
        "corrected_bc_frac",
    ]
    mets = {}
    for met in diagnostic_metrics:
        if met in sample_data.summary:
            mets[met] = sample_data.summary[met]
    ws_data.diagnostics = mets


def build_web_summary_data_common(
    sample_properties,
    sample_data: SampleData,
    pipeline,
    metadata,
    command,
    zoom_images=None,
    regist_images=None,
    sample_defs=None,
    *,
    projection: Projection,
):
    """Produce common data shared by both spatial and cell ranger, VDJ is currently independent."""
    is_spatial = command != CELLRANGER_COMMAND_NAME
    wsd = WebSummaryData(sample_properties, command, pipeline)

    species_list = sample_properties.genomes
    _build_summary_tab_common(
        metadata,
        sample_data,
        species_list,
        sample_properties,
        pipeline,
        wsd,
        zoom_images,
        regist_images,
        sample_defs,
    )
    if FILTERED_BCS_TRANSCRIPTOME_UNION in sample_data.summary:
        # if the sample is not Ab only build the GEX tab
        _build_analysis_tab_common(
            metadata,
            sample_data,
            pipeline,
            species_list,
            sample_properties,
            wsd,
            projection=projection,
        )
    # use presence/absence of the histograms as a proxy for whether to build the antibody/antigen tabs
    if sample_data.antibody_histograms is not None:
        _build_antibody_tab_common(sample_data, sample_properties, wsd, projection=projection)
    if sample_data.antigen_histograms is not None:
        _build_antigen_tab_common(sample_data, sample_properties, wsd)
    _build_diagnostic_values(sample_data, wsd)

    for fb in FEATURE_BARCODINGS:
        # Not all three of theses will be present

        # feature barcoding sequencing info (Count)
        add_data(
            wsd.summary_tab,
            wsd.alarms,
            feature_barcode_sequencing_table(metadata, sample_data, species_list, fb),
        )
        # feature barcoding application metric (Count)
        if sample_properties.is_spatial or pipeline != PIPELINE_AGGR:
            add_data(
                wsd.summary_tab,
                wsd.alarms,
                feature_barcode_application_table(metadata, sample_data, species_list, fb),
            )
        # aggregation metrics (AGGR)
        add_data(
            wsd.summary_tab,
            wsd.alarms,
            feature_barcode_aggregation_table(metadata, sample_data, sample_properties, fb),
        )

    # Targeting
    add_data(
        wsd.analysis_tab,
        wsd.alarms,
        targeted_table(metadata, sample_data, species_list, is_spatial=is_spatial),
    )
    return wsd


def build_web_summary_html_sc_and_aggr(
    sample_properties,
    sample_data_paths,
    pipeline,
    metadata,
    command_name,
    *,
    projection: Projection,
    sample_defs=None,
):
    """Produce the main websummary."""
    assert isinstance(sample_properties, CountSampleProperties)
    assert isinstance(sample_data_paths, SampleDataPaths)
    sample_data = load_sample_data(sample_properties, sample_data_paths, (projection,))

    species_list = sample_properties.genomes
    web_sum_data = build_web_summary_data_common(
        sample_properties,
        sample_data,
        pipeline,
        metadata,
        command_name,
        sample_defs=sample_defs,
        projection=projection,
    )

    # Single cell specific stuff
    add_data(
        web_sum_data.summary_tab,
        web_sum_data.alarms,
        normalization_summary_table(metadata, sample_data, sample_properties),
    )
    add_data(
        web_sum_data.summary_tab,
        web_sum_data.alarms,
        aggregation_by_sample_table(metadata, sample_data, sample_properties),
    )
    add_data(
        web_sum_data.summary_tab,
        web_sum_data.alarms,
        batch_correction_table(metadata, sample_data, species_list),
    )

    # t-SNE/UMAP plot appears as a constant left plot in the Cell Ranger
    # clustering selector
    if web_sum_data.clustering_selector:
        web_sum_data.clustering_selector.left_plots = umi_on_projection_plot(
            sample_data,
            spatial=sample_properties.is_spatial,
            projection=projection,
            library_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
            tag_type=None,
        )
        if pipeline == PIPELINE_AGGR and projection == UMAP_NAME:
            web_sum_data.alarms.extend(
                [
                    {
                        "formatted_value": None,
                        "title": "UMAP Projection",
                        "message": """UMAP projection is now the default for Cellranger Aggr. Run cellranger aggr with --enable-tsne=true to enable t-SNE projection.""",
                        "level": INFO_THRESHOLD,
                    }
                ]
            )
    if web_sum_data.antibody_clustering_selector:
        web_sum_data.antibody_clustering_selector.left_plots = umi_on_projection_plot(
            sample_data,
            projection=projection,
            spatial=sample_properties.is_spatial,
            library_type=rna_library.ANTIBODY_LIBRARY_TYPE,
            tag_type=None,
        )

    # Barnyard
    b_table = barnyard_table(metadata, sample_data, sample_properties, species_list)
    if b_table:
        web_sum_data.analysis_tab["barnyard"] = b_table
    return web_sum_data


def build_web_summary_html_sc(
    filename: str,
    sample_properties: CountSampleProperties,
    sample_data_paths: SampleDataPaths,
    pipeline: str,
    sample_defs=None,
):
    if sample_properties.is_lt:
        metadata = LTMetricAnnotations(intron_mode_alerts=sample_properties.include_introns)
    elif sample_properties.is_targeted:
        metadata = TargetedMetricAnnotations()
    else:
        metadata = MetricAnnotations(intron_mode_alerts=sample_properties.include_introns)
    command = CELLRANGER_COMMAND_NAME
    web_sum_data = build_web_summary_html_sc_and_aggr(
        sample_properties,
        sample_data_paths,
        pipeline,
        metadata,
        command,
        projection=TSNE_NAME,
        sample_defs=sample_defs,
    )
    # to circumvent self._check_valid
    data = copy.deepcopy(web_sum_data)
    data = data.to_dict_for_json()

    write_html_file(filename, web_sum_data)
    return data
