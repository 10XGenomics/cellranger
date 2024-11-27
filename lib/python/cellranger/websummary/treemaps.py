#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Code to produce plotly treemap plots."""

from __future__ import annotations

import json

import numpy as np
import plotly.graph_objects as go

from cellranger.rna.library import ANTIBODY_LIBRARY_TYPE
from cellranger.websummary.isotypes import LINK_HELP_TEXT
from cellranger.websummary.plotly_tools import MODE_BAR_BUTTONS, PLOT_CONFIG

MIN_ANTIBODY_UMI = 1


def make_treemap_help(is_antibody: bool, is_spatial: bool) -> dict:
    """Create a help text for the antibody treemap plots.

    Args:
        is_antibody (bool): A boolean indicating whether the features are antibodies.
        is_spatial (bool): A boolean indicating whether the sample is spatial or not.

    Returns:
        dict: A dictionary containing the help text and title for the treemap plot.

    """
    singular = "Antibody" if is_antibody else "Antigen"
    barcode_text = "spot" if is_spatial else "cell"
    cell_calls = f"{barcode_text} barcodes" if is_antibody else f"GEX {barcode_text} barcodes"
    link_text = LINK_HELP_TEXT if is_spatial else ""

    return {
        "helpText": f"Relative composition of {link_text} {singular.lower()} counts for features with at least {MIN_ANTIBODY_UMI} UMI.\n"
        f"Box size represents fraction of total normalized antibody UMIs from {cell_calls} that are derived from this {singular.lower()}. "
        f"Hover over a box to view more information on a particular {singular.lower()}, including number of associated barcodes.",
        "title": f"Distribution of {singular} Counts",
    }


# pylint: disable=too-many-locals
def make_antibody_treemap_plot(ab_matrix, lib_type, subset_features, is_spatial=True):
    """Create a "TreeMap" plot from an antibody or antigen filtered matrix.

    Args:
        ab_matrix (str): a preloaded filtered antibody count matrix
        lib_type (str): antibody or antigen library type
        subset_features (int): number of total umis per feature to filter on
        is_spatial(bool): if the sample is spatial or not
    Returns:
        A dict with the JSON data needed to make the bar plot representing this histogram.
    """
    if ab_matrix.get_shape()[0] == 0:
        return None
    is_antibody = ANTIBODY_LIBRARY_TYPE == lib_type

    # Generate the plot text
    ## The Feature Name label in the box will use the feature ref "id"
    antibody_label = [feature.id.decode("utf-8") for feature in ab_matrix.feature_ref.feature_defs]
    tags = [feature.tags for feature in ab_matrix.feature_ref.feature_defs]
    ## Calculate the fraction total umi for each feature and convert them to strings with parentheses for plotting
    ab_frac_values = ab_matrix.get_frac_counts_per_feature()
    ab_percent_strings = [f"{x:.2%}" for x in ab_frac_values]
    ab_percent_strings = list(map("({})".format, ab_percent_strings))
    ## generate secondary name labels for spatial samples
    ### Check if there is a gene_symbol or secondary_name in the feature def tag dict
    has_secondary_name = any("gene_symbol" in d for d in tags) or any(
        "secondary_name" in d for d in tags
    )
    if has_secondary_name:
        # These account for different versions of the feature reference
        keys = ["secondary_name", "gene_symbol"]
        secondary_names = [tag[key] for tag in tags for key in keys if key in tag]
    else:
        secondary_names = None
    ## Make the label for the box
    antibody_label = list(map(" ".join, zip(antibody_label, ab_percent_strings)))
    ## Calculate the number of tissue/cell associated barcodes per feature
    barcodes_per_feature = list(map(str, ab_matrix.get_numbcs_per_feature()))

    # Subset Features
    ## build an index based on the minimum number UMIs
    ab_index = np.where(np.array(ab_matrix.get_counts_per_feature()) >= subset_features)[0]
    ## Based on a minimum number of UMIs subset all lists to features that comply with that minimum
    ab_frac_values = np.array(ab_frac_values)[ab_index].tolist()
    antibody_label = np.array(antibody_label)[ab_index].tolist()
    ab_percent_strings = np.array(ab_percent_strings)[ab_index].tolist()
    barcodes_per_feature = np.array(barcodes_per_feature)[ab_index].tolist()
    if has_secondary_name:
        secondary_names = np.array(secondary_names)[ab_index].tolist()

    # Hover Lables
    ## Create the lables and hover lables for the plots depending on the type of assay it is
    if is_spatial:
        cell_or_tissue = "Tissue"
        barcodes_or_cells = " barcodes"
    else:
        cell_or_tissue = "Cell"
        barcodes_or_cells = " cells"

    hover_antibody_labs = [
        "<b>Feature ID (% Total Normalized Antibody UMI):</b><br>" + ab for ab in antibody_label
    ]
    hover_barcodes_per_feature_labs = [
        f"<b>{cell_or_tissue} Associated Barcodes per Feature:</b><br>" + name
        for name in barcodes_per_feature
    ]
    barcodes_per_feature_labs = [name + barcodes_or_cells for name in barcodes_per_feature]
    if has_secondary_name:
        secondary_labs = ["<b>Secondary Name:</b><br>" + name for name in secondary_names]
        hover_text = list(
            map(
                "<br>".join,
                zip(hover_antibody_labs, hover_barcodes_per_feature_labs, secondary_labs),
            )
        )
    else:
        hover_text = list(
            map("<br>".join, zip(hover_antibody_labs, hover_barcodes_per_feature_labs))
        )

    ## This adds html to the end of the string to remove the ugly "trace 0" box
    hover_text = [text + "<extra></extra>" for text in hover_text]

    # Make the Plot
    treemap_plot = go.Figure(
        go.Treemap(
            labels=["<b>" + label + "</b>" for label in antibody_label],
            values=ab_frac_values,
            parents=[""] * len(antibody_label),
            text=(
                list(
                    map(
                        "<br>".join,
                        zip(
                            [
                                "<b>" + name + barcodes_or_cells + "</b>"
                                for name in barcodes_per_feature
                            ],
                            ["<b>" + name + "</b>" for name in secondary_names],
                        ),
                    )
                )
                if has_secondary_name
                else barcodes_per_feature_labs
            ),
            marker_colorscale="Blackbody",
            hovertemplate=hover_text,
        )
    )

    # Convert the figure to a json string and then to a dict
    treemap_plot = json.loads(treemap_plot.to_json())
    # set the config and remove the reset axis
    PLOT_CONFIG[MODE_BAR_BUTTONS] = [["toImage"]]
    treemap_plot["config"] = PLOT_CONFIG
    treemap_plot_data = {"plot": treemap_plot, "help": make_treemap_help(is_antibody, is_spatial)}
    return treemap_plot_data
