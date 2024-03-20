#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Gets the sufficient statistics from rust stage and generate plotly plot."""


import json

import numpy as np

import tenkit.safe_json as tk_safe_json
from cellranger.targeted.targeted_constants import GDNA_PLOT_NAME
from cellranger.websummary.numeric_converters import array_to_float32_base64
from cellranger.websummary.plotly_tools import PLOT_CONFIG

__MRO__ = """
stage GET_GDNA_PLOT(
    in  json gdna_plot_sufficient_stats,
    out json summary,
    src py   "stages/targeted/get_gdna_plot",
)
"""


def main(args, outs):
    with open(args.gdna_plot_sufficient_stats) as f:
        data_in = json.load(f)
    data_in["unspliced_counts"] = np.array(data_in["unspliced_counts"])
    data_in["spliced_counts"] = np.array(data_in["spliced_counts"])

    # Plotting code is from
    # lib/python/cellranger/analysis/segment_model_fitter.py
    alpha = 0.02
    x_crit = data_in["model_crit_point"]
    a_1 = float(data_in["model_constant"])
    cis = data_in["spliced_counts"] < x_crit
    constant_trace = {
        "x": array_to_float32_base64(data_in["spliced_counts"][cis].tolist()),
        "y": array_to_float32_base64(data_in["unspliced_counts"][cis].tolist()),
        "mode": "markers",
        "type": "scatter",
        "marker": {"color": "#71439A", "opacity": alpha},
        "name": "Constant",
    }
    cis = np.invert(cis)
    linear_trace = {
        "x": array_to_float32_base64(data_in["spliced_counts"][cis].tolist()),
        "y": array_to_float32_base64(data_in["unspliced_counts"][cis].tolist()),
        "type": "scatter",
        "mode": "markers",
        "marker": {"color": "#00A9A1", "opacity": alpha},
        "name": "Linear",
    }
    constant_line = {
        "x": [0, float(x_crit)],
        "y": [a_1, a_1],
        "type": "line",
        "marker": {"color": "#D413EB"},
        "name": "Mean",
        "mode": "lines",
    }
    chart = {
        "config": PLOT_CONFIG,
        "layout": {
            "xaxis": {"title": "Spliced Probe Counts (Ln(1+Count))"},
            "yaxis": {"title": "Unspliced Probe Counts (Ln(1+Count))"},
        },
        "data": [constant_trace, linear_trace, constant_line],
    }
    gdna_summary = {GDNA_PLOT_NAME: chart}
    with open(outs.summary, "w") as outf:
        tk_safe_json.dump_numpy(gdna_summary, outf)
