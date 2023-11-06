#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Utility functions for dealing with altair plots."""

import altair as alt
import pandas as pd

LAYER = "layer"
HCONCAT = "hconcat"
VCONCAT = "vconcat"
DATA = "data"
DATA_NAME = "name"
DATA_FORMAT = "format"
DATA_FORMAT_CSV_DICT = {"type": "csv"}


def sanitise_data_names(chart_dict):
    """Sanitise data names in chart dict."""
    if LAYER in chart_dict:
        for layer_dict in chart_dict[LAYER]:
            sanitise_data_names(layer_dict)
    if HCONCAT in chart_dict:
        for hconcat_dict in chart_dict[HCONCAT]:
            sanitise_data_names(hconcat_dict)
    if VCONCAT in chart_dict:
        for vconcat_dict in chart_dict[VCONCAT]:
            sanitise_data_names(vconcat_dict)
    if DATA in chart_dict:
        if DATA_NAME not in chart_dict[DATA]:
            raise ValueError("Found nameless data in Altair JSON. Not expected VEGA out.")
        else:
            chart_dict[DATA][DATA_NAME] += ".csv"
            chart_dict[DATA][DATA_FORMAT] = DATA_FORMAT_CSV_DICT


def sanitise_chart_dict(chart_dict):
    """Clean up the chart dict."""
    # Convert all datasets into CSVs from JSONs
    raw_dataset_names = list(chart_dict["datasets"].keys())
    for raw_dataset in raw_dataset_names:
        csv_string = pd.DataFrame.from_dict(chart_dict["datasets"][raw_dataset]).to_csv(index=False)
        new_dataset_name = f"{raw_dataset}.csv"
        chart_dict["datasets"][new_dataset_name] = csv_string
        chart_dict["datasets"].pop(raw_dataset)
    sanitise_data_names(chart_dict)


def chart_to_json(
    chart: alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart,
) -> dict:
    """Convert an altair chart to a Dict.

    We store data in CSV format rather than as a JSON. This avoids having to
    write out field names for all rows and usually saves us a lot in
    storage.

    Args:
        chart (alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart):
            Altair Chart of various flavours.

    Returns:
        dict: Dict corresponding to the chart.
    """
    chart_dict = chart.to_dict()
    sanitise_chart_dict(chart_dict)
    return chart_dict
