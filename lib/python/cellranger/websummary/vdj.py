# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

from six import ensure_binary

import cellranger.vdj.chain_types as chain_types
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.constants.shared as shared
from cellranger.vdj.clonotype import extract_clonotype_id_from_name
from cellranger.webshim.constants.shared import CELLRANGER_COMMAND_NAME, PIPELINE_VDJ
from cellranger.websummary.metrics import VDJMetricAnnotations
from cellranger.websummary.react_components import WebSummaryData
from cellranger.websummary.react_summarize import write_html_file
from cellranger.websummary.sample_properties import VdjSampleProperties
from cellranger.websummary.summary_tab import add_data, pipeline_info_table

if TYPE_CHECKING:
    from cellranger.webshim.data import SampleData

REFERENCE_GENOME_KEY = "vdj_reference_genomes"
VDJ_RANK_PLOT_HELP = [
    [
        "Barcode Rank Plot",
        [
            "The plot shows the count of filtered UMIs mapped to each barcode. "
            "A barcode must have a contig that aligns to a V segment to be identified as a targeted cell. "
            "(In the denovo case, the only requirement is a contig's presence.) "
            "There must also be at least three filtered UMIs with at least two read pairs each. "
            "It is possible that a barcode with at least as many filtered UMIs as another cell-associated barcode is not identified as a targeted cell. "
            "The color of the graph is based on the local density of cell-associated barcodes. "
            "Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts."
        ],
    ]
]


def vdj_hero_metrics(metadata, sample_data: SampleData):
    """Make metrics.

    Args:
        metadata (VdjMetricAnnotation): object
        sample_data: object

    Returns:
        (dict): Dictionary that will populate the hero metrics in the web summary
    """
    if sample_data is None or sample_data.summary is None:
        return None

    # TODO: This is a HACK because hero metrics are fixed keys
    fake_name = {
        "vdj_filtered_bcs": "filtered_bcs_transcriptome_union",
        "vdj_total_raw_read_pairs_per_filtered_bc": "multi_transcriptome_total_raw_reads_per_filtered_bc",
        "multi_vdj_assembly_contig_pair_productive_full_len_bc_count": "filtered_bcs_median_unique_genes_detected",
    }

    HERO_METRIC_KEYS = ["vdj_filtered_bcs", "vdj_total_raw_read_pairs_per_filtered_bc"]
    if sample_data.summary.get(REFERENCE_GENOME_KEY, None):
        HERO_METRIC_KEYS.append("multi_vdj_assembly_contig_pair_productive_full_len_bc_count")

    data = {}
    for key in HERO_METRIC_KEYS:
        metrics = metadata.gen_metric_list(sample_data.summary, [key])
        for metric in metrics:
            data[fake_name[metric.key]] = metric.gen_metric_dict()

    # No need to raise alert on vdj_filtered_bcs since that will be done in the
    # cell calling table

    return data


def create_table_with_alarms(
    table_key, title, metric_keys, alarm_keys, metadata, sample_data, chain_type=None
):
    """Construct an instance.

    Args:
        table_key: Return te result in a dictionary with this key
        title: The title for the table in the web summary
        metric_keys: Use these keys to generate metric rows
        alarm_keys: Use these keys for raising alerts
        metadata: VdjMetricAnnotation object
        sample_data: webshim.data.SampleData object
        chain_type: (Optional) One of chain_types.CHAIN_TYPE_SPECS defining if it's a T or B or both
    """
    if sample_data is None or sample_data.summary is None:
        return None

    data_dict = {}

    metrics = metadata.gen_metric_list(sample_data.summary, metric_keys, chain_type=chain_type)
    if metrics:
        data_dict["help"] = {
            "title": title,
            "data": metadata.gen_metric_helptext(metric_keys, chain_type=chain_type),
        }
        data_dict["table"] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}

    if not data_dict:
        return None

    result = {table_key: data_dict}

    # Alerts.
    if alarm_keys:
        alarms = metadata.gen_metric_list(sample_data.summary, alarm_keys, chain_type=chain_type)
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]

        if new_alarms:
            result["alarms"] = new_alarms

    return result


def vdj_enrichment_table(metadata, sample_data, chain_type):
    """Make table.

    Args:
        metadata (VdjMetricAnnotation): object
        sample_data (webshim.data.SampleData): object
        chain_type: One of chain_types.CHAIN_TYPE_SPECS defining if it's a T or B or both

    Returns:
        (dict): Dictionary that will populate the enrichment table in the web
            summary which contains the reads mapped to any vdj gene and reads
            mapped to individual chains.
    """
    VDJ_ENRICHMENT_KEYS = [
        "multi_vdj_recombinome_mapped_reads_frac",
        "{chain}_vdj_recombinome_mapped_reads_frac",
    ]
    return create_table_with_alarms(
        "vdj_enrichment",
        "Enrichment",
        VDJ_ENRICHMENT_KEYS,
        VDJ_ENRICHMENT_KEYS,
        metadata,
        sample_data,
        chain_type=chain_type,
    )


def vdj_expression_table(metadata, sample_data, chain_type):
    """Make table.

    Args:
        metadata (VdjMetricAnnotation): object
        sample_data (webshim.data.SampleData): object
        chain_type: One of chain_types.CHAIN_TYPE_SPECS defining if it's a T or B or both

    Returns:
        (dict): Dictionary that will populate the expression table in the web
            summary which contains the median UMIs for each chain
    """
    VDJ_EXPRESSION_KEYS = [
        "{chain}_vdj_assembly_umis_per_cell_median",
    ]
    return create_table_with_alarms(
        "vdj_expression",
        "V(D)J Expression",
        VDJ_EXPRESSION_KEYS,
        VDJ_EXPRESSION_KEYS,
        metadata,
        sample_data,
        chain_type=chain_type,
    )


def vdj_sequencing_table(metadata, sample_data):
    """Make table.

    Args:
        metadata (VdjMetricAnnotation): object
        sample_data (webshim.data.SampleData): object

    Returns:
        (dict): Dictionary that will populate the sequencing table in the web
            summary which contains the number of read pairs, valid barcodes and
            various Q30 metrics
    """
    VDJ_SEQUENCING_KEYS = [
        "VDJ_total_read_pairs",
        "VDJ_unprocessed_read_pairs",
        "vdj_good_bc_frac",
        "VDJ_bc_bases_with_q30_frac",
        "VDJ_read_bases_with_q30_frac",
        "VDJ_read2_bases_with_q30_frac",
        "VDJ_umi_bases_with_q30_frac",
    ]
    return create_table_with_alarms(
        "vdj_sequencing",
        "Sequencing",
        VDJ_SEQUENCING_KEYS,
        VDJ_SEQUENCING_KEYS,
        metadata,
        sample_data,
    )


def load_paired_metrics_from_chain_type(chain_type):
    chain_type = ensure_binary(chain_type)
    assert chain_type in chain_types.CHAIN_TYPE_SPECS
    if chain_type == chain_types.TR_CHAIN_TYPE:
        return ["TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac"]
    elif chain_type == chain_types.IG_CHAIN_TYPE:
        return [
            "IGK_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
            "IGL_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
        ]
    else:
        return [
            "TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac",
            "IGK_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
            "IGL_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
        ]


def vdj_annotation_table(metadata, sample_data, chain_type):
    """Make table.

    Args:
        metadata (VdjMetricAnnotation): object
        sample_data (webshim.data.SampleData): object
        chain_type: One of chain_types.CHAIN_TYPE_SPECS defining if it's a T or B or both

    Returns:
        (dict): Dictionary that will populate the annotation table in the web summary
    """
    VDJ_ANNOTATION_KEYS = [
        "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
    ]
    VDJ_ANNOTATION_KEYS.extend(load_paired_metrics_from_chain_type(chain_type))
    VDJ_ANNOTATION_KEYS.extend(
        [
            "multi_raw_vdj_paired_clonotype_diversity",
            "{chain}_vdj_assembly_contig_bc_frac",
            "{chain}_vdj_assembly_cdr_detected_bc_frac",
            "{chain}_vdj_assembly_contig_full_len_bc_frac",
            "{chain}_vdj_assembly_prod_cdr_bc_frac",
        ]
    )
    return create_table_with_alarms(
        "vdj_annotation",
        "V(D)J Annotation",
        VDJ_ANNOTATION_KEYS,
        VDJ_ANNOTATION_KEYS,
        metadata,
        sample_data,
        chain_type=chain_type,
    )


def vdj_cell_call_table(metadata, sample_data):
    """Make table.

    Args:
        metadata (VdjMetricAnnotation): object
        sample_properties (webshim.constants.vdj.VdjSampleProperties): objet
        sample_data (webshim.data.SampleData): object

    Returns:
        (dict): Dictionary that will populate the cell calling table in the web
            summary and plots the barcode rank plot
    """
    VDJ_CELL_CALL_KEYS = [
        "vdj_filtered_bcs",
        "vdj_total_raw_read_pairs_per_filtered_bc",
        "vdj_assemblable_read_pairs_per_filtered_bc",
        "vdj_filtered_bcs_cum_frac",
    ]

    chart = {
        "layout": {
            "title": "Barcode Rank",
            "width": 470,
            "height": 313,
            "margin": {"l": 60, "r": 0, "t": 30, "b": 40},
            "hovermode": "closest",
            "xaxis": {
                "title": "Barcodes",
                "type": "log",
                "fixedrange": False,
                "showline": True,
            },
            "yaxis": {
                "title": "UMI counts",
                "type": "log",
                "fixedrange": False,
                "showline": True,
            },
        },
        "data": [
            {
                "x": [],
                "y": [],
                "name": "Cells",
                "hoverinfo": "name",
                "type": "scattergl",
                "mode": "lines",
                "line": {
                    "color": shared.BC_PLOT_CMAP(1.0),
                    "width": shared.BC_RANK_PLOT_LINE_WIDTH,
                },
            },
            {
                "x": [],
                "y": [],
                "name": "Background",
                "hoverinfo": "name",
                "type": "scattergl",
                "mode": "lines",
                "line": {
                    "color": shared.BC_PLOT_CMAP(0.0),
                    "width": shared.BC_RANK_PLOT_LINE_WIDTH,
                },
            },
        ],
        "config": {
            "staticPlot": False,
            "displayModeBar": True,
            "modeBarButtons": [["toImage", "resetScale2d"]],
            "dragmode": "zoom",
            "showAxisDragHandles": True,
            "scrollZoom": False,
        },
    }

    data_dict = create_table_with_alarms(
        "cells", "Cells", VDJ_CELL_CALL_KEYS, VDJ_CELL_CALL_KEYS, metadata, sample_data
    )

    barcode_rank_plot = cr_webshim.plot_vdj_barcode_rank(chart, None, sample_data)
    if barcode_rank_plot:
        data_dict["cells"]["barcode_knee_plot"] = barcode_rank_plot
        data_dict["cells"]["help"]["data"] = VDJ_RANK_PLOT_HELP + data_dict["cells"]["help"]["data"]

    return data_dict


def vdj_clonotype_table(sample_data):
    """Make table.

    Args:
        sample_data (webshim.data.SampleData): the clonotype summary

    Returns:
        (dict): Dictionary that will populate the table containing top 10
            clonotypes in the analysis tab
    """
    assert sample_data.vdj_clonotype_summary is not None
    table_heading = ["Clonotype ID", "CDR3s", "Number of cells", "Fraction of cells"]
    table_help = {
        "title": "Top 10 Clonotype CDR3 Sequences",
        "data": [
            [
                "",
                [
                    'This table lists the CDR3 sequence of the first exact subclonotype of the 10 \
                    most abundant clonotypes in this sample. For each of the top 10 clonotypes, \
                    the constant region, number of cells, and what percentage of the \
                    dataset those cells occupy (Fraction of cells) are also displayed. For the full table \
                    and more details, please refer to the "clonotypes.csv" and \
                    "consensus_annotations.csv" files produced by the pipeline.',
                ],
            ],
        ],
    }
    table_rows = []
    for _, row in sample_data.vdj_clonotype_summary.iloc[0:10].iterrows():
        table_rows.append(
            [
                str(extract_clonotype_id_from_name(row["clonotype_id"])),
                row["cdr3s_aa"].replace(";", "<br>"),
                "{}".format(row["frequency"]),
                "{:.2%}".format(row["proportion"]),
            ]
        )
    result = {"help": table_help, "table": {"header": table_heading, "rows": table_rows}}
    return result


def vdj_clonotype_chart(sample_properties, sample_data):
    """Make chart.

    Args:
        sample_properties (webshim.constants.vdj.VdjSampleProperties): object
        sample_data (webshim.data.SampleData): the clonotype summary

    Returns:
        (dict): Dictionary that will be used to plot the clonotype histogram in
            the analysis tab
    """
    assert sample_data.vdj_clonotype_summary is not None
    chart = {
        "layout": {
            "showlegend": False,
            "xaxis": {
                "type": "category",
                "title": "Clonotype ID",
            },
            "yaxis": {
                "title": "Fraction of Cells",
            },
            "margin": {"l": 60, "t": 0, "r": 40},
            "hovermode": "closest",
        },
        "data": [{"type": "bar"}],
        "config": {
            "staticPlot": False,
            "displayModeBar": True,
            "modeBarButtons": [["toImage"]],
            "dragmode": "zoom",
        },
    }

    x = []
    y = []
    for _, row in sample_data.vdj_clonotype_summary.iloc[0:10].iterrows():
        x.append(extract_clonotype_id_from_name(row["clonotype_id"]))
        y.append(row["proportion"])

    chart["data"][0]["x"] = x
    chart["data"][0]["y"] = y

    result = {
        "help": {
            "helpText": "This histogram displays the fraction of cells (percentage of cells) \
            occupied by the 10 most abundant clonotypes in this sample. The clonotype IDs on the X \
            axis correspond to the clonotype IDs listed in the Top 10 Clonotype CDR3 Sequences \
            below.",
            "title": "Top 10 Clonotype Frequencies",
        },
        "plot": chart,
    }
    return result


def vdj_summary_tab(web_sum_data, metadata, sample_properties, sample_data: SampleData):
    """Make tab.

    Args:
        metadata: VdjMetricAnnotation object
        sample_properties: VdjSampleProperties
        sample_data: webshim.data.SampleData object

    Returns:
        (dict): Dictionary that defines the contents in the "Summary" tab
        Alarms: Alarms that need to be raised
    """
    assert isinstance(sample_properties, VdjSampleProperties)
    assert isinstance(web_sum_data, WebSummaryData)
    summary_tab_json = web_sum_data.summary_tab
    alarm_list = web_sum_data.alarms
    add_data(
        summary_tab_json,
        alarm_list,
        pipeline_info_table(sample_data, sample_properties, PIPELINE_VDJ),
    )
    if sample_data.summary.get(REFERENCE_GENOME_KEY, None):
        add_data(
            summary_tab_json,
            alarm_list,
            vdj_enrichment_table(metadata, sample_data, chain_type=sample_properties.chain_type),
        )
        add_data(
            summary_tab_json,
            alarm_list,
            vdj_expression_table(metadata, sample_data, chain_type=sample_properties.chain_type),
        )
        add_data(
            summary_tab_json,
            alarm_list,
            vdj_annotation_table(metadata, sample_data, chain_type=sample_properties.chain_type),
        )
    add_data(summary_tab_json, alarm_list, vdj_cell_call_table(metadata, sample_data))
    add_data(summary_tab_json, alarm_list, vdj_sequencing_table(metadata, sample_data))
    add_data(summary_tab_json, alarm_list, vdj_hero_metrics(metadata, sample_data))


def vdj_analysis_tab(web_sum_data, sample_properties, sample_data):
    """Make tab.

    Args:
        sample_properties: VdjSampleProperties
        sample_data: webshim.data.SampleData object

    Returns:
        Dictionary that defines the contents in the "Analysis" tab
    """
    if sample_data.vdj_clonotype_summary is None:
        return

    web_sum_data.vdj_tab = {
        "vdj_clonotype": vdj_clonotype_table(sample_data),
        "vdj_clonotype_hist": vdj_clonotype_chart(sample_properties, sample_data),
    }


def vdj_web_summary_json(metadata, sample_properties, sample_data):
    """Make json.

    Args:
        metadata: VdjMetricAnnotation object
        sample_properties: VdjSampleProperties
        sample_data: webshim.data.SampleData object

    Returns:
        (dict): Dictionary that is used to generate the web summary
    """
    web_sum_data = WebSummaryData(sample_properties, CELLRANGER_COMMAND_NAME, PIPELINE_VDJ)
    vdj_summary_tab(web_sum_data, metadata, sample_properties, sample_data)
    vdj_analysis_tab(web_sum_data, sample_properties, sample_data)
    return web_sum_data


def build_vdj_web_summary_html(filename, sample_properties, sample_data):
    """Make html.

    Args:
        filename: File to write the web summary html to
        sample_properties: VdjSampleProperties
        sample_data: webshim.data.SampleData object

    Returns:
        (dict): Dictionary that is used to generate the web summary
    """
    metadata = VDJMetricAnnotations()
    web_sum_data = vdj_web_summary_json(metadata, sample_properties, sample_data)
    # to circumvent self._check_valid
    data = copy.deepcopy(web_sum_data)
    data = data.to_dict_for_json()
    write_html_file(filename, web_sum_data)
    return data
