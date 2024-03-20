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
from cellranger.webshim.jibes_plotting import _make_rc_vectors, make_color_map, make_histogram_plot
from cellranger.websummary import react_components as rc

JIBES_BIPLOT_N_POINTS = 100000


MIN_QV_VAL = 1.0
MAX_QV_VAL = 100.0


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
