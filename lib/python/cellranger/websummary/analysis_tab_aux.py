#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
"""This file shows plots in the analysis tab that aren't related to the outputs of the SC_RNA_ANALYZER pipeline.

e.g. targeted UMI plots, sequencing saturation plots, etc
splitting things up in this way allows turing repo to utilize UMAP/clustering plots without pulling in a lot of transitive dependencies
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

import cellranger.analysis.jibes_constants as jibes
import cellranger.rna.library as rna_library
import cellranger.targeted.utils as cr_tgt_utils
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.constants.shared as shared_constants
import cellranger.websummary.plotly_tools as pltly
from cellranger.analysis.multigenome import MultiGenomeAnalysis
from cellranger.analysis.singlegenome import UMAP_NAME, SingleGenomeAnalysis
from cellranger.targeted.targeted_constants import (
    GDNA_CONTENT_METRIC,
    GDNA_PLOT_NAME,
    GDNA_UNSPLICED_THRESHOLD,
    TARGETING_METHOD_HC,
    TARGETING_METHOD_TL,
)
from cellranger.targeted.utils import OFFTARGET_WS_LABEL, TARGETED_WS_LABEL
from cellranger.webshim.jibes_plotting import make_color_map
from cellranger.websummary.analysis_tab_core import projection_layout_config
from cellranger.websummary.helpers import get_projection_key
from cellranger.websummary.metrics import SpatialAggrMetricAnnotations
from cellranger.websummary.numeric_converters import round_floats_in_list
from cellranger.websummary.react_components import BarnyardPanel

if TYPE_CHECKING:
    from cellranger.feature.feature_assignments import CellsPerFeature

READS_PER_UMI_LAYOUT_CONFIG = {
    "xaxis": {
        "title": "UMIs per gene, log10",
        "showline": False,
        "zeroline": False,
        "fixedrange": False,
    },
    "yaxis": {
        "title": "Reads per gene, log10",
        "showline": False,
        "zeroline": False,
        "fixedrange": False,
    },
    "legend": {
        "orientation": "h",
        "y": -0.2,
    },
    "hovermode": "closest",
}

HELP_TEXT_KEY = "helpText"
TITLE_KEY = "title"
CANONICAL_VISIUM_HD_BIN_NAME = "square_008um"

SEQ_SATURATION_PLOT_HELP = {
    "helpText": "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per cell), up to the observed sequencing depth. "
    "Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. ",
    "title": "Sequencing Saturation",
}

TARGETED_RPU_PLOT_HELP = [
    [
        "Reads vs UMIs per Gene by Targeting Status",
        [
            "This plot shows the number of reads per gene (log10) in cell-associated barcodes as "
            "a function of the number of UMIs per gene (log10) in cell-associated barcodes, with "
            "genes colored by whether or not they were targeted. "
            "The yellow line represents the optimal threshold "
            "(specifically, its y-intercept = log10(Reads per UMI threshold)) "
            "in a mixture model that classifies genes into two classes based on Reads per UMI values. "
            "If sequencing saturation is low, this line will not be shown. "
            "Ideally, targeted genes will lie above the dotted line while non-targeted genes will be below it.",
        ],
    ]
]

SPATIAL_TARGETED_RPU_PLOT_HELP = [
    [
        "Reads vs UMIs per Gene by Targeting Status",
        [
            "This plot shows the number of reads per gene (log10) in tissue-covered spots as a "
            "function of the number of UMIs per gene (log10) in tissue-covered spots, with genes "
            "colored by whether or not they were targeted. "
            "The yellow line represents the optimal threshold "
            "(specifically, its y-intercept = log10(Reads per UMI threshold)) "
            "in a mixture model that classifies genes into two classes based on Reads per UMI values. "
            "If sequencing saturation is low, this line will not be shown. "
            "Ideally, targeted genes will lie above the dotted line while non-targeted genes will be below it.",
        ],
    ]
]

SPATIAL_SEQ_SATURATION_PLOT_HELP = {
    "helpText": "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per spot), up to the observed sequencing depth. "
    "Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. ",
    "title": "Sequencing Saturation",
}

RTL_SEQ_SATURATION_PLOT_HELP = {
    HELP_TEXT_KEY: "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per spot), up to the observed sequencing depth. "
    "Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted probe ligation products have been sequenced. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. ",
    TITLE_KEY: "Sequencing Saturation",
}

VISIUM_HD_RTL_SEQ_SATURATION_PLOT_HELP = {
    HELP_TEXT_KEY: "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in total number of read pairs), up to the observed sequencing depth. "
    "Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted probe ligation products have been sequenced. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. ",
    TITLE_KEY: "Sequencing Saturation",
}

MEDIAN_GENE_PLOT_HELP = {
    "helpText": "This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",
    "title": "Median Genes per Cell",
}

SPATIAL_MEDIAN_GENE_PLOT_HELP = {
    HELP_TEXT_KEY: "This plot shows the Median Genes per Spot as a function of downsampled sequencing depth in mean reads per spot, up to the observed sequencing depth. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",
    TITLE_KEY: "Median Genes per Spot",
}

SPATIAL_HD_MEAN_GENE_PLOT_HELP_JSON_STRING = (
    f'"{HELP_TEXT_KEY}": "This plot shows the Mean Genes '
    "per {bin_size} µm bin as a function of downsampled sequencing depth in mean reads per {bin_size} µm bin, up to the observed sequencing depth. "
    'The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",'
    f'"{TITLE_KEY}": '
    '"Median Genes per {bin_size} µm bin"'
)

BARNYARD_PLOT_HELP = [
    [
        "Cell UMI Counts Plot",
        [
            "Each point represents a cell-barcode. "
            "The axes measure the total UMI counts in each cell-barcode that mapped to each transcriptome reference. "
            "The points are colored by the number of inferred cells in the GEM associated with each barcode. "
            "A multiplet represents either a GEM inferred to have encapsulated >1 cell or a barcode sequence that was shared by multiple single-cell GEMs."
        ],
    ]
]


GDNA_PLOT_HELP = [
    [
        "Segmented Linear Model Plot",
        [
            "Each point represents a gene that has probes targeting both exon-junction-spanning and non-exon-junction-spanning regions, 'spliced' and 'unspliced', respectively. "
            "Unspliced probes can stem from open gDNA and from RNA. "
            "Spliced probes are expected to stem only from RNA. "
            "A segmented linear model is used to estimate where the unspliced and spliced counts begin to deviate. "
            "The mean of unspliced counts in purple estimates the UMI background level per unspliced probe. "
            "Counts less than this have a high probability of stemming from gDNA."
        ],
    ]
]


def gdna_table(metadata, sample_data, species_list) -> BarnyardPanel | None:
    """Make the equivalent of a barnyard table but with plots and metrics.

    related to gDNA
    """
    if (
        sample_data is None
        or sample_data.summary is None
        or GDNA_PLOT_NAME not in sample_data.summary
        or GDNA_CONTENT_METRIC not in sample_data.summary
        or GDNA_UNSPLICED_THRESHOLD not in sample_data.summary
    ):
        return None
    metric_keys = [GDNA_CONTENT_METRIC, GDNA_UNSPLICED_THRESHOLD]
    metrics = metadata.gen_metric_list(sample_data.summary, metric_keys, species_list)
    chart = sample_data.summary[GDNA_PLOT_NAME]
    data = {"table": {"rows": [[metric.name, metric.value_string] for metric in metrics]}}
    help_text = GDNA_PLOT_HELP + metadata.gen_metric_helptext(metric_keys)
    data["help"] = {"title": "UMIs from Genomic DNA", "data": help_text}
    return BarnyardPanel(chart, data)


def barnyard_table(metadata, sample_data, sample_properties, species_list):
    """Barnyard table and barnyard plot."""
    if (
        len(species_list) <= 1
        or sample_data is None
        or sample_data.summary is None
        or sample_data.get_analysis(MultiGenomeAnalysis) is None
    ):
        return None

    metric_keys = [
        "filtered_bcs_observed_all",
        "filtered_bcs_inferred_multiplets",
        "filtered_bcs_inferred_multiplet_rate",
    ]

    gems = {}
    metrics = metadata.gen_metric_list(sample_data.summary, metric_keys, species_list)
    gems["table"] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}

    helptext = metadata.gen_metric_helptext(metric_keys) + BARNYARD_PLOT_HELP
    gems["help"] = {"title": "GEM Partitions", "data": helptext}
    chart = {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "updatemenus": [
                {
                    "buttons": [
                        {
                            "label": "Linear Scale",
                            "method": "update",
                            "args": [
                                {"visible": [True, True]},
                                {"yaxis": {"type": ""}, "xaxis": {"type": ""}},
                            ],
                        },
                        {
                            "label": "Log Scale",
                            "method": "update",
                            "args": [
                                {"visible": [True, True]},
                                {"yaxis": {"type": "log"}, "xaxis": {"type": "log"}},
                            ],
                        },
                    ],
                    "xanchor": "center",
                    "yanchor": "bottom",
                    "pad": {"b": 20},
                }
            ],
            "title": "Cell UMI Counts",
            "showlegend": True,
            "hovermode": "closest",
        },
        "data": [
            {
                "x": [],
                "y": [],
                "mode": "markers",
                "type": "scattergl",
            },
        ],
    }
    plot = cr_webshim.plot_barnyard_barcode_counts(chart, sample_properties, sample_data)
    return BarnyardPanel(plot, gems)


def targeted_table(metadata, sample_data, species_list, is_spatial=False):
    # pylint: disable=missing-function-docstring,too-many-locals
    if (
        not sample_data.is_targeted()
        or sample_data.targeting_method != TARGETING_METHOD_HC
        or isinstance(metadata, SpatialAggrMetricAnnotations)
        or sample_data.feature_metrics is None
    ):
        return None
    if is_spatial:
        metric_keys = [
            "multi_frac_conf_transcriptomic_reads",
            "num_genes",
            "spatial_num_genes_quantifiable",
            "spatial_num_rpu_enriched_genes",
            "spatial_mean_reads_per_umi_per_gene_cells",
        ]
    else:
        metric_keys = [
            "multi_frac_conf_transcriptomic_reads",
            "num_genes",
            "num_genes_quantifiable",
            "num_rpu_enriched_genes",
            "mean_reads_per_umi_per_gene_cells",
        ]

    # umi filtering keys -- set to NA if it was disabled and/or inactive
    umi_filtering_active = sample_data.summary.get("filtered_target_umi_count_threshold", 0) > 1
    # if enrichment calculation is disabled, log_rpu_threshold set to str(np.nan)
    enrichment_calculations_disabled = not isinstance(
        sample_data.summary.get("log_rpu_threshold"), float
    )

    table_rows = []
    helptext = []
    for metric in metric_keys:
        on_target_metric = metadata.gen_metric_list(
            sample_data.summary, [metric + TARGETED_WS_LABEL.metric_suffix], species_list
        )
        off_target_metric = metadata.gen_metric_list(
            sample_data.summary, [metric + OFFTARGET_WS_LABEL.metric_suffix], species_list
        )
        value_on_target = on_target_metric[0].value_string
        value_off_target = off_target_metric[0].value_string

        if is_spatial:
            if enrichment_calculations_disabled and metric == "spatial_num_rpu_enriched_genes":
                value_on_target = value_off_target = "N/A"
        elif enrichment_calculations_disabled and metric == "num_rpu_enriched_genes":
            value_on_target = value_off_target = "N/A"

        metric_name = on_target_metric[0].name
        metric_helptext = metadata.gen_metric_helptext([metric + TARGETED_WS_LABEL.metric_suffix])
        table_rows.append([metric_name, value_on_target, value_off_target])
        helptext.extend(metric_helptext)

    if not is_spatial:
        for metric in ["filtered_target_umi_count_threshold", "filtered_target_umi_reads_frac"]:
            on_target_metric = metadata.gen_metric_list(sample_data.summary, [metric], species_list)
            value_on_target = on_target_metric[0].value_string if umi_filtering_active else "N/A"
            # get metric name from helptext in case metric doesn't exist
            metric_helptext = metadata.gen_metric_helptext([metric])
            helptext.extend(metric_helptext)
            metric_name = metric_helptext[0][0]
            table_rows.append([metric_name, value_on_target, ""])

    # make Off-Target header stay on one line instead of wrapping
    # &#8209; is a non-breaking hyphen character
    table_heading = ["", "Targeted", "Non&#8209;Targeted"]

    helptext.extend(SPATIAL_TARGETED_RPU_PLOT_HELP if is_spatial else TARGETED_RPU_PLOT_HELP)
    table_help = {
        "title": "Targeted Enrichment",
        "data": helptext,
    }
    data_dict = {"help": table_help, "table": {"header": table_heading, "rows": table_rows}}

    targeted = {
        "plot": reads_per_umi_plot(sample_data),
        "gene_targeting_table": data_dict,
    }

    # add alarm keys if enrichment calculation not disabled
    if not enrichment_calculations_disabled:
        alarm_keys = ["frac_on_target_genes_enriched", "frac_off_target_genes_enriched"]
        alarms = metadata.gen_metric_list(sample_data.summary, alarm_keys, species_list)
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
    else:
        new_alarms = []

    return {"targeted": targeted, shared_constants.ALARMS: new_alarms}


def seq_saturation_plot(sample_data, sample_properties):
    # pylint: disable=missing-function-docstring
    is_targeted = sample_data.is_targeted() and sample_data.targeting_method == TARGETING_METHOD_HC
    chart = {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "showlegend": bool(is_targeted),
            "hovermode": "closest",
            "xaxis": {
                "title": "Mean Reads per {}".format(
                    "Cell" if sample_properties.is_spatial is False else "Spot"
                ),
                "fixedrange": False,
            },
            "yaxis": {
                "title": "Sequencing Saturation",
                "range": [0, 1],
                "fixedrange": False,
            },
        },
        "data": [],  # data entries are built in the function
    }

    plot = cr_webshim.plot_subsampled_scatterplot_metric(
        chart,
        sample_properties,
        sample_data,
        metric_suffix="subsampled_duplication_frac",
        show_multi_genome_only=True,
        is_targeted=is_targeted,
        show_targeted_only=False,
    )

    if plot:
        # keep targeting color scheme consistent
        if is_targeted:
            plot = _add_targeting_line_colors(plot)
        return {
            "seq_saturation_plot": {
                "plot": plot,
                "help": (
                    SEQ_SATURATION_PLOT_HELP
                    if sample_properties.is_spatial is False
                    else (
                        RTL_SEQ_SATURATION_PLOT_HELP
                        if sample_properties.is_spatial
                        and sample_data.targeting_method == TARGETING_METHOD_TL
                        else SPATIAL_SEQ_SATURATION_PLOT_HELP
                    )
                ),
            }
        }
    else:
        return None


def reads_per_umi_plot(sample_data):
    # pylint: disable=invalid-name
    """Reads vs umi plot with separation boundary."""
    per_feature_metrics = sample_data.feature_metrics.get_df()

    data = []

    # plot off-target genes first so they are in the background
    groups = [
        (OFFTARGET_WS_LABEL.label, 0, OFFTARGET_WS_LABEL.color),
        (TARGETED_WS_LABEL.label, 1, TARGETED_WS_LABEL.color),
    ]
    for group in groups:
        group_name = group[0]
        group_value = group[1]
        group_color = group[2]
        per_feature_metrics_subset = per_feature_metrics[
            per_feature_metrics[cr_tgt_utils.TARGETING_COLNAME] == group_value
        ]

        # avoid plotting totally overlapping points, esp at (0,0) -- this is useless (can
        # only hover/see the last one plotted), slows down the webpage with tons of points,
        # and triggers a plotly bug where traces are not toggled on and off correctly
        per_feature_metrics_subset = per_feature_metrics_subset[
            per_feature_metrics_subset[cr_tgt_utils.READS_IN_CELLS_COLNAME] > 0
        ]
        per_feature_metrics_subset.drop_duplicates(
            subset=[
                cr_tgt_utils.UMIS_IN_CELLS_COLNAME,
                cr_tgt_utils.READS_IN_CELLS_COLNAME,
            ],
            inplace=True,
        )

        gene_names = per_feature_metrics_subset[cr_tgt_utils.FEATURE_NAME_COLNAME].tolist()
        x = np.log10(per_feature_metrics_subset[cr_tgt_utils.UMIS_IN_CELLS_COLNAME]).tolist()
        y = np.log10(per_feature_metrics_subset[cr_tgt_utils.READS_IN_CELLS_COLNAME]).tolist()
        data.append(
            {
                "name": group_name,
                "x": x,
                "y": y,
                "type": "scattergl",
                "mode": "markers",
                "marker": {
                    "opacity": 0.5,
                    "size": 5,
                    "color": group_color,
                },
                "text": [f"Gene: {gene_name}" for gene_name in gene_names],
            }
        )

    # plot separation boundary on top, if it was determined
    a = 1.0
    b = sample_data.summary["log_rpu_threshold"]
    if isinstance(b, float):
        x0 = 0
        x1 = np.log10(np.nanmax(per_feature_metrics[cr_tgt_utils.UMIS_IN_CELLS_COLNAME]) + 1)
        data.append(
            {
                "name": f"Reads per UMI threshold ({math.pow(10, b):.2f})",
                "x": [x0, x1],
                "y": [x0 * a + b, x1 * a + b],
                "mode": "lines",
                "line": {
                    "color": "#E6B72B",
                    "width": 3,
                },
            }
        )

    layout = READS_PER_UMI_LAYOUT_CONFIG.copy()
    layout["title"] = "Reads vs UMIs per Gene by Targeting Status"
    reads_vs_umis_plot = {"config": pltly.PLOT_CONFIG, "layout": layout, "data": data}
    return reads_vs_umis_plot


def median_gene_plot(sample_data, sample_properties, species_list):
    # pylint: disable=missing-function-docstring
    unit = "Cell" if sample_properties.is_spatial is False else "Spot"
    is_targeted = sample_data.is_targeted() and sample_data.targeting_method == TARGETING_METHOD_HC
    chart = {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "showlegend": True,
            "hovermode": "closest",
            "xaxis": {
                "title": f"Mean Reads per {unit}",
                "fixedrange": False,
            },
            "yaxis": {
                "title": "Median {}Genes per {}".format("Targeted " if is_targeted else "", unit),
                "fixedrange": False,
            },
        },
        "data": [],  # data entries are built in the function
    }

    plot = cr_webshim.plot_subsampled_scatterplot_metric(
        chart,
        sample_properties,
        sample_data,
        metric_suffix="subsampled_filtered_bcs_median_unique_genes_detected",
        references=species_list,
        is_targeted=is_targeted,
        show_targeted_only=True,
    )

    if plot:
        # keep targeting color scheme consistent
        if is_targeted:
            plot = _add_targeting_line_colors(plot)
        return {
            "median_gene_plot": {
                "plot": plot,
                "help": (
                    MEDIAN_GENE_PLOT_HELP
                    if sample_properties.is_spatial is False
                    else SPATIAL_MEDIAN_GENE_PLOT_HELP
                ),
            }
        }
    else:
        return None


def _add_targeting_line_colors(plot):
    """Change colors from default colors for traces in targeting plots in order to keep colors consistent."""
    on_target_trace_indices = [
        i
        for i in range(len(plot["data"]))
        if TARGETED_WS_LABEL.label.lower() in plot["data"][i]["name"].lower()
    ]
    off_target_trace_indices = [
        i
        for i in range(len(plot["data"]))
        if OFFTARGET_WS_LABEL.label.lower() in plot["data"][i]["name"].lower()
    ]

    for on_target_trace_idx in on_target_trace_indices:
        plot["data"][on_target_trace_idx]["line"]["color"] = TARGETED_WS_LABEL.color
    for off_target_trace_idx in off_target_trace_indices:
        plot["data"][off_target_trace_idx]["line"]["color"] = OFFTARGET_WS_LABEL.color

    return plot


def cmo_tags_on_umap_from_path(
    analysis_dir,
    cells_per_tag: CellsPerFeature,
    non_singlet_barcodes: CellsPerFeature,
    multiplexing_method: rna_library.BarcodeMultiplexingType,
):
    """Given the base analysis directory for a single genome analysis h5,.

    Get the CMO UMAP colored by CMO tag labels
    """
    if analysis_dir is None or cells_per_tag is None or non_singlet_barcodes is None:
        return None

    analysis = SingleGenomeAnalysis.load_default_format(
        analysis_dir, method="pca", projections=(UMAP_NAME,)
    )

    if analysis is None:
        return None

    return cmo_tags_on_umap_helper(
        analysis, cells_per_tag, non_singlet_barcodes, multiplexing_method
    )


def _filter_out_missing_barcodes(
    analysis: SingleGenomeAnalysis, barcode_groups: list[set[bytes]]
) -> set[bytes]:
    """With manually specified CMO it is possible that a Barcode will be in the filtered matrix, but not in the.

    analysis matrix produced by PREPROCESS_MATRIX.  This can happen if the barcode has zero umi counts across all
    features for example.  To avoid issues this causes, before plotting we remove any barcodes that aren't in the
    matrix.

    Args:
        analysis: the SingleGenomeAnalysis
        barcode_groups: a list of barcode sets that we'll filter over

    Returns:
       The set of barcodes in the matrix
    """
    matrix_bcs = set(analysis.matrix.bcs)
    for group in barcode_groups:
        to_remove = []
        for barcode in group:
            if barcode not in matrix_bcs:
                print(f"Warning: Barcode {barcode} is not in the analysis matrix.")
                to_remove.append(barcode)
        if to_remove:
            for bc in to_remove:
                group.remove(bc)
    return matrix_bcs


# whole-library UMAP for CMO tag data, colored by CMO tag.
def cmo_tags_on_umap_helper(
    analysis: SingleGenomeAnalysis,
    cells_per_tag: CellsPerFeature,
    non_singlet_barcodes: CellsPerFeature,
    multiplexing_method: rna_library.BarcodeMultiplexingType,
):
    # pylint: disable=too-many-locals
    """CMO tag labels on CMO UMAP helper function."""
    assert multiplexing_method.is_cell_multiplexed(), "Unsupported multiplexing method!"
    library_type = multiplexing_method.multiplexing_library_type()
    key = get_projection_key(library_type, 2)
    if key not in analysis.umap:
        return None
    umap_coordinates = analysis.get_umap(key=key).transformed_umap_matrix

    title = "UMAP Projection of Cells by CMO"
    data = []

    total_bcs = float(len(analysis.matrix.bcs))

    # get list of tags
    # sanity check that the info in cells_per_tag and non_singlet_barcodes match up
    assigned_bcs = set()
    tags = []
    for tag in cells_per_tag:
        if len(cells_per_tag[tag]) == 0:
            continue
        tags.append(tag)
        assigned_bcs = assigned_bcs.union(set(cells_per_tag[tag]))

    tags.sort()
    color_map = make_color_map(tags, jibes_plot=True)
    unassigned_and_blanks = set(list(analysis.matrix.bcs)) - assigned_bcs

    unassigned = set(non_singlet_barcodes[jibes.UNASSIGNED_FACTOR_NAME.encode()])
    blanks = set(non_singlet_barcodes[jibes.BLANK_FACTOR_NAME.encode()])
    multiplets = set(non_singlet_barcodes[jibes.MULTIPLETS_FACTOR_NAME.encode()])

    assert unassigned_and_blanks == unassigned.union(blanks)

    def add_data(tag, tag_prop, tag_indices):
        tag = tag.decode("utf8")
        data.append(
            {
                "name": f"{tag} ({tag_prop:.1%})",
                "x": round_floats_in_list(umap_coordinates[tag_indices, 0]),
                "y": round_floats_in_list(umap_coordinates[tag_indices, 1]),
                "type": "scattergl",
                "mode": "markers",
                "marker": {"opacity": 0.9, "size": 4, "color": color_map[tag]},
                "text": f"{tag}: {tag_prop:.1%}",
            }
        )

    # Calculate proportions before any filtering so numbers match what is expected
    blanks_prop = float(len(blanks)) / total_bcs
    multiplet_prop = float(len(multiplets)) / total_bcs
    unassigned_prop = float(len(unassigned)) / total_bcs

    # Filter out missing Barcodes and get a set to filter out the Tag barcodes furth down
    matrix_bcs = _filter_out_missing_barcodes(analysis, [blanks, multiplets, unassigned])
    # the order here (blanks, tags, multiplets, unassigned) is meant to match that of the jibes biplot
    # add trace for blanks
    blanks_indices = analysis.matrix.bcs_to_ints(blanks)
    add_data(b"Blank", blanks_prop, blanks_indices)

    # add traces for non-multiplet barcodes for each CMO tag
    for tag in tags:
        tag_barcodes = {bc for bc in cells_per_tag[tag] if bc not in multiplets}
        tag_prop = float(len(tag_barcodes)) / total_bcs
        tag_barcodes = tag_barcodes.intersection(matrix_bcs)
        tag_indices = analysis.matrix.bcs_to_ints(tag_barcodes)
        add_data(tag, tag_prop, tag_indices)

    # add trace for multiplets
    multiplet_indices = analysis.matrix.bcs_to_ints(multiplets)
    add_data(b"Multiplet", multiplet_prop, multiplet_indices)

    # add trace for unassigned
    unassigned_indices = analysis.matrix.bcs_to_ints(unassigned)
    add_data(b"Unassigned", unassigned_prop, unassigned_indices)

    # Note: the help text has been included in umap_cluster plot
    layout = projection_layout_config(projection=UMAP_NAME).copy()
    layout["title"] = title
    cmo_tag_umap_plot = {
        "config": pltly.PLOT_CONFIG,
        "layout": layout,
        "data": data,
    }
    return cmo_tag_umap_plot


def antibody_histogram_plot(sample_data):
    """Generate the json to render antibody hsitograms."""
    ab_histograms = sample_data.antibody_histograms
    if ab_histograms:
        antibody_histograms = {"antibody_histogram_plot": ab_histograms}
        return antibody_histograms
    else:
        return None


def antibody_treemap_plot(sample_data):
    """Generate the json to render antibody treemp."""
    ab_treemap = sample_data.antibody_treemap
    if ab_treemap:
        antibody_treemap = {"antibody_treemap_plot": ab_treemap}
        return antibody_treemap
    else:
        return None


def antibody_qc_plot(sample_data, plot_name):
    """Generate the json to render antibody QC plot."""
    if plot_name == "raw_normalized_heatmap":
        plot_key = "raw_normalized_heatmap"
    elif plot_name == "isotype_scatter":
        plot_key = "isotype_scatter"
    elif plot_name == "gex_fbc_correlation_heatmap":
        plot_key = "gex_fbc_correlation_heatmap"
    else:
        return None

    plot = getattr(sample_data, plot_key, None)
    if plot:
        return plot
    else:
        return None


def antigen_histogram_plot(sample_data):
    """Generate the json to render antigen hsitograms."""
    ag_histograms = sample_data.antigen_histograms
    if ag_histograms:
        antigen_histograms = {"antigen_histogram_plot": ag_histograms}
        return antigen_histograms
    else:
        return None
