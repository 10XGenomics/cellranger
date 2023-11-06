# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Shared ploting code for the web summary."""


from __future__ import annotations

import copy

from websummary.summarize import DEFAULT_FONT as DEFAULT_WEB_FONT

BUTTON_RESET_SCALED_2D = "resetScale2d"
BUTTON_TO_IMAGE = "toImage"
TO_IMAGE_BUTTON_OPTIONS = "toImageButtonOptions"
MODE_BAR_BUTTONS = "modeBarButtons"
CONFIG = "config"
SHOW_AXIS_DRAG_HANDLES = "showAxisDragHandles"

PLOT_CONFIG = {
    "staticPlot": False,
    "displayModeBar": True,
    MODE_BAR_BUTTONS: [[BUTTON_TO_IMAGE, BUTTON_RESET_SCALED_2D]],
    SHOW_AXIS_DRAG_HANDLES: True,
    TO_IMAGE_BUTTON_OPTIONS: {"width": 730, "height": 650},  # ensure that this is square
    "scrollZoom": False,
}

# Make a spatial config that is similar but with no zoom or reset axis button
SPATIAL_PLOT_CONFIG = copy.deepcopy(PLOT_CONFIG)
SPATIAL_PLOT_CONFIG[MODE_BAR_BUTTONS] = [[]]
del SPATIAL_PLOT_CONFIG[TO_IMAGE_BUTTON_OPTIONS]
SPATIAL_PLOT_CONFIG[SHOW_AXIS_DRAG_HANDLES] = False


def get_layout_from_chart(chart) -> dict:
    """Gets the layout from a plotly chart that's either a naked dict or an object.

    Args:
        chart:

    Returns:
        The layout dictionary.
    """
    if isinstance(chart, dict):
        layout = chart["layout"]
    elif hasattr(chart, "layout"):
        layout = chart.layout
    else:
        raise TypeError("Must pass an object or dictionary to update layout font.")
    return layout


def add_font_to_layout(chart, font=DEFAULT_WEB_FONT):
    """Updates a charts `layout` dictionary to use a specific font."""
    fnt = {"family": font}
    layout = get_layout_from_chart(chart)
    font_value = layout.get("font", {})
    if isinstance(font_value, dict) and len(font_value) == 0:
        font_value.update(fnt)
        layout["font"] = font_value
