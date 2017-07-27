#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
PIPELINE_AGGR = 'aggr'
PIPELINE_COUNT = 'count'
PIPELINE_REANALYZE = 'reanalyze'
PIPELINE_VDJ = 'vdj'
WEBSHIM_PIPELINES = [PIPELINE_AGGR, PIPELINE_COUNT, PIPELINE_REANALYZE, PIPELINE_VDJ]

ALARM_WARN = 'warn'
ALARM_ERROR = 'error'
DEFAULT_CHART_WIDTH = 600
DEFAULT_CHART_HEIGHT = 400

# Formatting of xy coordinates in javascript lists going to plotly
DATA_VALUE_FORMAT = '%.4f'

# these are the options for making the views fixed;
# add these to the modeBarButtonsToRemove list to
# disable movement & zoom
#
# src: src/components/modebar/buttons.js in Plotly
CHARTS_PLOTLY_MODEBAR_TRANSFORM_BUTTONS = [
    'zoom2d',
    'pan2d',
    'zoomIn2d',
    'zoomOut2d',
    'autoScale2d',
    #'resetScale2d'  can't totally disable interaction, it seems-- keep reset option
]

CHARTS_PLOTLY_EXPORT_BUTTONS = [
    'toImage',
    'sendDataToCloud',
]

CHARTS_PLOTLY_FIXED_CONFIG = {
    'modeBarButtonsToRemove': CHARTS_PLOTLY_MODEBAR_TRANSFORM_BUTTONS+CHARTS_PLOTLY_EXPORT_BUTTONS,
    'displaylogo': False,
    'showLink': False
}

CHARTS_PLOTLY_MOVABLE_CONFIG = {
    'modeBarButtonsToRemove': CHARTS_PLOTLY_EXPORT_BUTTONS,
    'displaylogo': False,
    'showLink': False
}

HISTOGRAM_METRIC_ORDER_INTEGER_BIN = 'order_integer'

# Freq expects ints, prop expects floats
HISTOGRAM_METRIC_ORDER_DECREASING_FREQUENCY = 'order_frequency'
HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION = 'order_proportion'

HISTOGRAM_METRIC_DEFAULT_ORDERING = HISTOGRAM_METRIC_ORDER_INTEGER_BIN

HISTOGRAM_METRIC_ORDER = [
    HISTOGRAM_METRIC_ORDER_INTEGER_BIN,
    HISTOGRAM_METRIC_ORDER_DECREASING_FREQUENCY,
    HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION,
]
