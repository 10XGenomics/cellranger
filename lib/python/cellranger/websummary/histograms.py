#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Code to produce plotly histogram plots in the form of bar-plots."""


from __future__ import annotations

from math import ceil, floor, log10
from typing import NamedTuple

import numpy as np
from six import ensure_str

from cellranger.feature.feature_assigner import MIN_COUNTS_PER_ANTIBODY
from cellranger.rna.library import ANTIBODY_LIBRARY_TYPE
from cellranger.webshim.jibes_plotting import make_color_map
from cellranger.websummary.isotypes import LINK_HELP_TEXT
from cellranger.websummary.plotly_tools import MODE_BAR_BUTTONS, PLOT_CONFIG

# pylint: disable=invalid-name

DEFAULT_MAX_HISTOGRAMS = 120
MHC_ALLELE = "mhc_allele"
TARGETING_ANTIGEN = "targeting_antigen"
HASHTAG_LABEL = " (Hashtag)"


class HistogramData(NamedTuple):
    """Data class used to make histograms."""

    data: np.array
    name: str
    group: str
    name_suffix: str


class BinningSummary(NamedTuple):
    """Data for summary stats used for a histogram trace."""

    min: float
    max: float
    bin_width: float
    total: float


def _get_bin_range(data: HistogramData, min_bin_width=0.05):
    assert min_bin_width > 0.0, "Must have positive bin width"
    min_val = np.amin(data.data)
    max_val = np.amax(data.data)
    std_dev = np.std(data.data)
    n = data.data.shape[0]
    total = data.data.sum()
    ## Attempt to do a bin size calculation similar to plotly, whose
    # logic is here:
    # https://github.com/plotly/plotly.js/blob/7588c9e057f5317332621e56933736eaece217b7/src/plots/cartesian/axes.js#L390-L408
    # I am not sure where the choice of 2 * sigma / n^(2/5) came from, but is similar to n^(1/3) rules
    # that are commonly used: https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width
    plotly_size = 2 * std_dev / pow(n, 0.4)
    bin_width = max(min_bin_width, plotly_size)
    return BinningSummary(min_val, max_val, bin_width, total)


def _round_to_1(x):
    """Returns the number rounded to one significant digit."""
    return round(x, -int(floor(log10(abs(x)))))


def make_histogram_help(is_antibody, is_spatial):
    """Returns help for antibody/antigen umi count histogram."""
    spot_or_cell = "spots" if is_spatial else "cells"
    link_text = LINK_HELP_TEXT if is_spatial else ""

    if is_antibody:
        return {
            "helpText": f"Only {link_text} antibodies with total UMI counts over {MIN_COUNTS_PER_ANTIBODY:,} are plotted, and only up to the top {DEFAULT_MAX_HISTOGRAMS} Antibodies by UMI count are shown. "
            f"The X-axis is the UMI counts in the log scale, while the Y-axis is the number of {spot_or_cell} observed with that UMI count.",
            "title": "Histogram of Antibody UMI Counts",
        }
    else:
        return {
            "helpText": f"Only Antigens with non-zero total UMI counts in the library are shown. The X-axis is the UMI counts in log scale while the Y-axis is the number of {spot_or_cell} observed with that UMI count. The histogram excludes the 0 data point i.e. the number of barcodes with 0 UMI counts for a given Antigen.\nDouble-click on the name of any feature in the plot legend to view the UMI counts only for that feature.",
            "title": "Histogram of Antigen UMI Counts",
        }


def make_histogram_as_bar_plot(
    counts: list[HistogramData],
    color_map,
    x_axis_title: str,
    y_axis_title: str,
    max_histograms=DEFAULT_MAX_HISTOGRAMS,
):
    """Create a "Histogram" in the form of a bar plot.

     Rather than passing an array of
    data to plotly with plot type histogram, to save space we pre-compute the bins
    and heights in order to make a more memory efficient `bar` plot to show the histogram data.

    Args:
        counts: A list of HistogramData objects
        x_axis_title: title for the x-axis
        y_axis_title: title for the y-axis
        color_map: A dictionary of names to colors
        max_histograms: If more that max_histograms are specified, only the top 24 based on sum counts will be used

    Returns:
        A dict with the JSON data needed to make the bar plot representing this histogram.
    """
    # pylint: disable=too-many-locals
    data = []
    if len(counts) > 0:
        bin_info = [_get_bin_range(x) for x in counts]

        ## Trim down if we have too many histograms
        cutoff = len(bin_info) > max_histograms
        if cutoff:
            totals = np.fromiter((x.total for x in bin_info), count=len(bin_info), dtype=float)
            # Note with argsort tied values can be arbitrarily excluded or included to cap at the max
            indices_over_cutoff = set(np.argsort(totals)[::-1][:DEFAULT_MAX_HISTOGRAMS])
            old_bin = bin_info
            old_counts = counts
            counts = []
            bin_info = []
            for i, (cnts, bininfo) in enumerate(zip(old_counts, old_bin)):
                if i in indices_over_cutoff:
                    counts.append(cnts)
                    bin_info.append(bininfo)

        ## Of all datasets to be plotted, find the minimum bin size, and min/max of range
        max_range = max(x.max for x in bin_info)
        min_range = min(x.min for x in bin_info)
        bin_width = min(x.bin_width for x in bin_info)
        bin_width = _round_to_1(bin_width)  # Use one significant digit

        low_end = floor(min_range / bin_width) * bin_width
        high_end = ceil(max_range / bin_width) * bin_width
        bin_cutoffs = np.arange(low_end, high_end + 0.01, bin_width)
        mid_points = (bin_width / 2.0) + bin_cutoffs[:-1]
        mid_points = list(mid_points)

        for cnts in counts:
            hist, _ = np.histogram(cnts.data, bin_cutoffs)
            trace = {
                "type": "bar",
                "x": mid_points,
                "y": list(hist),
                "width": bin_width,
                "opacity": 0.6,
                "name": cnts.name,
                "marker": {"color": color_map[cnts.name]},
                "showlegend": True,
            }
            if cnts.group is not None and cnts.group != "":
                trace["legendgroup"] = cnts.group
                trace["legendgrouptitle"] = {"text": ensure_str(cnts.group)}
            if cnts.name_suffix is not None and cnts.name_suffix != "":
                trace["name"] = " ".join([trace["name"], cnts.name_suffix])
            data.append(trace)
    plot = {
        "data": data,
        "layout": {
            "barmode": "overlay",
            "xaxis": {"title": x_axis_title},
            "yaxis": {"title": y_axis_title},
            "legend": {"groupclick": "toggleitem"},
        },
    }
    return plot


def make_antibody_histograms(ab_matrix, is_spatial, hashtags=None):
    """Take an antibody or antigen barcode matrix and generate json for web summary histograms.

    This is also used to define whether the antibody tab should be inserted into the websummary.
    If present the tab gets inserted. If not it doesn't.
    """
    if ab_matrix.get_shape()[0] == 0:
        return None
    is_antibody = ANTIBODY_LIBRARY_TYPE in ab_matrix.get_library_types()
    if is_antibody:
        over_threshold_features = np.flatnonzero(
            ab_matrix.get_counts_per_feature() >= MIN_COUNTS_PER_ANTIBODY
        )
    else:
        over_threshold_features = np.flatnonzero(ab_matrix.get_counts_per_feature() > 0)

    if ab_matrix.features_dim > len(over_threshold_features):
        ab_matrix = ab_matrix.select_features(over_threshold_features)

    ab_counts = np.log10(ab_matrix.m.toarray() + 1.0)

    # get allele information (if available) for legend grouping
    if not is_antibody and MHC_ALLELE in ab_matrix.feature_ref.all_tag_keys:
        legend_group = [ensure_str(f.tags[MHC_ALLELE]) for f in ab_matrix.feature_ref.feature_defs]
    else:
        legend_group = [None] * ab_matrix.features_dim
    # get control information (if available)
    if not is_antibody and TARGETING_ANTIGEN in ab_matrix.feature_ref.all_tag_keys:
        control_suffix = [
            "(Negative Control)" if f.tags[TARGETING_ANTIGEN] == "False" else ""
            for f in ab_matrix.feature_ref.feature_defs
        ]
    else:
        control_suffix = [""] * ab_matrix.features_dim
    vector_names = [ensure_str(x) for x in ab_matrix.feature_ids_map.keys()]

    # Specify which features are hashtags.
    hashtags = hashtags or []
    vector_names = [f"{x}{HASHTAG_LABEL}" if x in hashtags else x for x in vector_names]
    color_map = make_color_map(vector_names, jibes_plot=False)
    data = []
    for (
        feature_id,
        index,
    ) in ab_matrix.feature_ids_map.items():
        cnts = ab_counts[index, :]
        if not is_antibody:
            cnts = cnts[cnts != 0]
        hashtag_label_value = HASHTAG_LABEL if feature_id.decode("utf-8") in hashtags else ""
        data.append(
            HistogramData(
                cnts,
                (f"{ensure_str(feature_id)}{hashtag_label_value}"),
                legend_group[index],
                control_suffix[index],
            )
        )

    spot_or_cell = "Spots" if is_spatial else "Cells"

    histogram_plot = make_histogram_as_bar_plot(
        data, color_map, "Log10(1+Count)", f"Number of {spot_or_cell}"
    )
    PLOT_CONFIG[MODE_BAR_BUTTONS] = [["toImage"]]
    histogram_plot["config"] = PLOT_CONFIG
    histogram_data = {
        "plot": histogram_plot,
        "help": make_histogram_help(is_antibody, is_spatial),
    }
    return histogram_data
