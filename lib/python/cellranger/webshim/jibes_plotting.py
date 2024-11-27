#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""File meant to provide some functions used by others, without a dependency on websummary."""
from __future__ import annotations

from six import ensure_str

import cellranger.websummary.numeric_converters
from cellranger.analysis import jibes_constants as jibes

# https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/economist_pal/

# Add three colors for multiplets, blanks and unassigned

_JIBES_COLOR_PALETTE_TAGS = [
    "#3E647D",
    "#7B92AB",
    "#82C0E9",
    "#2D6D66",
    "#BFA19C",
    "#008BBC",
    "#97B6B0",
    "#D7D29E",
    "#1A476F",
    "#90353B",
    "#9C8847",
    "#938DD2",
]
_JIBES_MULTIPLET_COLOR = "#C10534"
_JIBES_BLANK_COLOR = "#6E8E84"
_JIBES_UNASSIGNED_COLOR = "#CAC27E"
JIBES_COLOR_PALETTE = _JIBES_COLOR_PALETTE_TAGS + [
    _JIBES_MULTIPLET_COLOR,
    _JIBES_BLANK_COLOR,
    _JIBES_UNASSIGNED_COLOR,
]


def make_color_map(features, jibes_plot=True):
    """Make an economist style color map for the websummary.

    :param features: list of tags
    :return: dict of {features:color}
    """
    color_map = {}
    for i, feature in enumerate(sorted(features)):
        color_map[ensure_str(feature)] = JIBES_COLOR_PALETTE[i % len(_JIBES_COLOR_PALETTE_TAGS)]
    if jibes_plot:
        color_map[ensure_str(jibes.MULTIPLETS_FACTOR_NAME)] = _JIBES_MULTIPLET_COLOR
        color_map[ensure_str(jibes.BLANK_FACTOR_NAME)] = _JIBES_BLANK_COLOR
        color_map[ensure_str(jibes.UNASSIGNED_FACTOR_NAME)] = _JIBES_UNASSIGNED_COLOR
    return color_map


def make_histogram_plot(
    feature_names, feature_counts, color_map, x_axis_lab: str = "Log 10 (1 + UMI Counts)"
):
    """Make a histogram plot.

    Args:
        feature_names: labels for features
        feature_counts: the data to plot
        color_map: the colors to plot

    Returns:
        histogram data for websummary JSON
    """
    traces = []
    for i, tag in enumerate(feature_names):
        trace = {
            "x": feature_counts[i],
            "type": "histogram",
            "opacity": 0.6,
            "name": ensure_str(tag),
            "marker": {
                "color": color_map[ensure_str(tag)],
            },
        }
        traces.append(trace)
    plot = {
        "data": traces,
        "layout": {"barmode": "overlay", "xaxis": {"title": x_axis_lab}},
    }
    return plot


JIBES_SHARED_VECTORS_RESOURCE_KEY_PREFIX = "shared_jibes_data_vector"
RESOURCE_PREFIX = "_resources"


def _make_rc_vectors(original_vectors):
    shared_resources = {}
    # this is sort of duplicating functionality that should be encapsulated by ResourceCollection in react_components.py
    # but we are sort of doing this in an ad-hoc way and passing it off to the Rust websummary code to
    # put it in the right place
    vectors = []
    for i, v in enumerate(original_vectors):
        store_key = f"{JIBES_SHARED_VECTORS_RESOURCE_KEY_PREFIX}_{i}"
        use_key = f"{RESOURCE_PREFIX}_{store_key}"
        shared_resources[store_key] = (
            cellranger.websummary.numeric_converters.array_to_float32_base64(v)
        )
        vectors.append(use_key)
    return shared_resources, vectors
