#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Websummary code for JIBES models."""
from __future__ import annotations

# sample 100k points for jibes biplot
import numpy as np
from six import ensure_str

from cellranger.analysis import jibes_constants as jibes
from cellranger.analysis.jibes import format_count_column_name
from cellranger.websummary import react_components as rc

JIBES_BIPLOT_N_POINTS = 100000


MIN_QV_VAL = 1.0
MAX_QV_VAL = 100.0

# https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/economist_pal/

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
# Add three colors for multiplets, blanks and unassigned

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
    feature_names, feature_counts, color_map, x_axis_lab: str = "Log 10 (1 + Count)"
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
        shared_resources[store_key] = rc.array_to_float32_base64(v)
        vectors.append(use_key)
    return shared_resources, vectors


def get_jibes_threshold_values(assignments):
    """Gets clipped assignment values for plotting in the websummary.

    Args:
        assignments: A pandas data frame used by jibes assigner

    Returns:
        A list of threshold values
    """
    # Clamp confidence values to avoid numeric issues
    assn_probs = assignments[jibes.ASSIGNMENT_PROB_COL_NAME].values
    ### Numeric issues values close to 1 by a small delta, sanity check we're close to
    # the expected range, then clamp the values to the QV min/max.
    over_1 = np.sum(assn_probs > 1.0001).sum() > 0
    under_0 = np.sum(assn_probs < -0.0001).sum() > 0
    if over_1 or under_0:
        # Should never ever happen
        raise Exception(
            "Normalized probabilities were outside of tolerance ranges.  This likely indicates a bug, please contact support@10xgenomics.com."
        )
    max_confidence = 1 - np.power(10, -(MAX_QV_VAL / 10.0))
    min_confidence = 1 - np.power(10, -(MIN_QV_VAL / 10.0))
    assn_probs = assn_probs.clip(min_confidence, max_confidence)
    phred_values = -10.0 * np.log10(1.0 - assn_probs)
    phred_values = phred_values.tolist()
    return phred_values


def make_jibes_biplot_histogram(assigner):
    """Make a biplot.

    Args:
        fitter (JibesTagAssigner): instance

    Returns:
        dict: jibes biplot and histogram data for websummary JSON
    """
    # pylint: disable=too-many-locals
    # isinstance(assigner, JibesTagAssigner)
    fitter = assigner.fitter
    if fitter is None:
        return None
    assignments = assigner.jibes_assignments

    # we want roughly JIBES_BIPLOT_N_POINTS so find out how many rows to subsample
    n_rows_to_downsample = int(JIBES_BIPLOT_N_POINTS)
    # create a random mask for the rows with that many points
    downsample_indices = np.full(assignments.shape[0], False)
    downsample_indices[:n_rows_to_downsample] = True
    np.random.shuffle(downsample_indices)

    # downsample the assignments
    assignments = assignments[downsample_indices]

    # downsample the data vectors

    original_vectors = [
        np.asarray(assignments[format_count_column_name(name)]).ravel().tolist()
        for name in fitter.data.column_names
    ]

    for vector in original_vectors:
        assert len(vector) == len(assignments)

    shared_resources, vectors = _make_rc_vectors(original_vectors)
    vector_names = fitter.data.column_names
    color_map = make_color_map(vector_names, jibes_plot=True)
    histogram_data = make_histogram_plot(vector_names, vectors, color_map)
    phred_values = get_jibes_threshold_values(assignments)
    slider_values = [-10.0 * np.log10(1.0 - jibes.JIBES_MIN_CONFIDENCE)]

    biplot_data = rc.BiplotsWithThresholding(
        vectors,
        vector_names,
        True,
        False,
        phred_values,
        slider_values,
        "Assignment Phred Score Threshold",
        "Singlets",
        "Filtered Singlets",
        "% Filtered",
        [
            ensure_str(jibes.BLANK_FACTOR_NAME),
            ensure_str(jibes.MULTIPLETS_FACTOR_NAME),
            ensure_str(jibes.UNASSIGNED_FACTOR_NAME),
        ],
        [ensure_str(x) for x in assignments[jibes.ASSIGNMENT_COL_NAME].values.tolist()],
        color_map,
    )

    jibes_biplot_histogram = {
        "biplot": biplot_data,
        "histogram": histogram_data,
        "_resources": shared_resources,
    }

    return jibes_biplot_histogram
