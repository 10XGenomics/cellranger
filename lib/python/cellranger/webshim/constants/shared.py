#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import math

PIPELINE_AGGR = "aggr"
PIPELINE_COUNT = "count"
PIPELINE_REANALYZE = "reanalyze"
PIPELINE_VDJ = "vdj"
PIPELINE_META_COUNT = "meta-count"
CELLRANGER_COMMAND_NAME = "Cell Ranger"

# TODO delete this if it isn't being used anywhere
WEBSHIM_PIPELINES = [
    PIPELINE_AGGR,
    PIPELINE_COUNT,
    PIPELINE_REANALYZE,
    PIPELINE_VDJ,
    PIPELINE_META_COUNT,
]

ALARM_WARN = "warn"
ALARM_ERROR = "error"
DEFAULT_CHART_WIDTH = 600
DEFAULT_CHART_HEIGHT = 400
ALARMS = "alarms"

# Formatting of xy coordinates in javascript lists going to plotly
DATA_VALUE_FORMAT = "%.4f"

# these are the options for making the views fixed;
# add these to the modeBarButtonsToRemove list to
# disable movement & zoom
#
# src: src/components/modebar/buttons.js in Plotly
CHARTS_PLOTLY_MODEBAR_TRANSFORM_BUTTONS = [
    "zoom2d",
    "pan2d",
    "zoomIn2d",
    "zoomOut2d",
    "autoScale2d",
    # 'resetScale2d'  can't totally disable interaction, it seems-- keep reset option
]

CHARTS_PLOTLY_EXPORT_BUTTONS = [
    "toImage",
    "sendDataToCloud",
]

CHARTS_PLOTLY_FIXED_CONFIG = {
    "modeBarButtonsToRemove": CHARTS_PLOTLY_MODEBAR_TRANSFORM_BUTTONS
    + CHARTS_PLOTLY_EXPORT_BUTTONS,
    "displaylogo": False,
    "showLink": False,
}

CHARTS_PLOTLY_MOVABLE_CONFIG = {
    "modeBarButtonsToRemove": CHARTS_PLOTLY_EXPORT_BUTTONS,
    "displaylogo": False,
    "showLink": False,
}

BC_RANK_PLOT_LINE_WIDTH = 3
# Gradient scheme used in the barcode rank plot
BC_PLOT_COLORS = [
    "#dddddd",
    "#d1d8dc",
    "#c6d3dc",
    "#bacfdb",
    "#aecada",
    "#a3c5d9",
    "#97c0d9",
    "#8cbbd8",
    "#80b7d7",
    "#74b2d7",
    "#6aadd6",
    "#66abd4",
    "#62a8d2",
    "#5ea5d1",
    "#59a2cf",
    "#559fce",
    "#519ccc",
    "#4d99ca",
    "#4997c9",
    "#4594c7",
    "#4191c5",
    "#3d8dc4",
    "#3a8ac2",
    "#3787c0",
    "#3383be",
    "#3080bd",
    "#2c7cbb",
    "#2979b9",
    "#2676b7",
    "#2272b6",
    "#1f6eb3",
    "#1d6ab0",
    "#1a65ac",
    "#1861a9",
    "#155ca6",
    "#1358a2",
    "#10539f",
    "#0e4f9b",
    "#0b4a98",
    "#094695",
    "#09438f",
    "#0a4189",
    "#0c3f83",
    "#0d3d7c",
    "#0e3b76",
    "#103970",
    "#11366a",
    "#123463",
    "#14325d",
    "#153057",
]


def BC_PLOT_CMAP(density):
    """Colormap utility fn to map a number to one of the colors in the gradient.

    color scheme defined above
    Input
    - density : A real number in the range [0,1]
    """
    assert density >= 0.0
    assert density <= 1.0
    levels = len(BC_PLOT_COLORS)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return BC_PLOT_COLORS[ind]


HISTOGRAM_METRIC_ORDER_INTEGER_BIN = "order_integer"

# Freq expects ints, prop expects floats
HISTOGRAM_METRIC_ORDER_DECREASING_FREQUENCY = "order_frequency"
HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION = "order_proportion"

HISTOGRAM_METRIC_DEFAULT_ORDERING = HISTOGRAM_METRIC_ORDER_INTEGER_BIN

HISTOGRAM_METRIC_ORDER = [
    HISTOGRAM_METRIC_ORDER_INTEGER_BIN,
    HISTOGRAM_METRIC_ORDER_DECREASING_FREQUENCY,
    HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION,
]
