#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Utility functions for dealing with altair plots."""

import altair as alt
import pandas as pd
from pandas.api.types import is_bool_dtype, is_datetime64_any_dtype, is_numeric_dtype

LAYER = "layer"
HCONCAT = "hconcat"
VCONCAT = "vconcat"
DATA = "data"
DATA_NAME = "name"
DATA_FORMAT = "format"
DATASETS_KEY = "datasets"
CSV_PARSE_KEY = "parse"
DATA_FORMAT_CSV_DICT = {"type": "csv"}


def get_vega_dtype(dtype: pd.DataFrame.dtypes):
    """Convert pandas dtype to vega-lite dtype.

    Vega-lite CSV parse
    dtypes are - https://vega.github.io/vega-lite-api/api/csv.html#parse .
    """
    if is_bool_dtype(dtype):
        return "boolean"
    elif is_datetime64_any_dtype(dtype):
        return "date"
    elif is_numeric_dtype(dtype):  # A little tricky as bools are also considered numeric.
        return "number"
    else:
        raise ValueError("Expected code to be unreachable")


def sanitise_data_names(chart_dict, name_to_format_dict):
    """Sanitise data names in chart dict.

    Uses name_to_format_dict to get column types to parse.
    """
    if LAYER in chart_dict:
        for layer_dict in chart_dict[LAYER]:
            sanitise_data_names(layer_dict, name_to_format_dict)
    if HCONCAT in chart_dict:
        for hconcat_dict in chart_dict[HCONCAT]:
            sanitise_data_names(hconcat_dict, name_to_format_dict)
    if VCONCAT in chart_dict:
        for vconcat_dict in chart_dict[VCONCAT]:
            sanitise_data_names(vconcat_dict, name_to_format_dict)
    if DATA in chart_dict:
        if DATA_NAME not in chart_dict[DATA]:
            raise ValueError("Found nameless data in Altair JSON. Not expected VEGA out.")
        else:
            old_data_name = chart_dict[DATA][DATA_NAME]
            chart_dict[DATA][DATA_FORMAT] = name_to_format_dict[old_data_name]
            chart_dict[DATA][DATA_NAME] += ".csv"


def sanitise_chart_dict(chart_dict):
    """Clean up the chart dict."""
    # Convert all datasets into CSVs from JSONs
    raw_dataset_names = list(chart_dict[DATASETS_KEY].keys())
    name_to_format_dict = {}
    for raw_dataset in raw_dataset_names:
        name_to_format_dict[raw_dataset] = DATA_FORMAT_CSV_DICT.copy()
        df = pd.DataFrame.from_dict(chart_dict[DATASETS_KEY][raw_dataset])
        csv_string = df.to_csv(index=False)
        # Getting columns with special parse
        name_to_format_dict[raw_dataset][CSV_PARSE_KEY] = {
            x: get_vega_dtype(y)
            for x, y in df.dtypes.items()
            if (is_numeric_dtype(y) or is_datetime64_any_dtype(y) or is_bool_dtype(y))
        }
        new_dataset_name = f"{raw_dataset}.csv"
        chart_dict[DATASETS_KEY][new_dataset_name] = csv_string
        chart_dict[DATASETS_KEY].pop(raw_dataset)
    sanitise_data_names(chart_dict, name_to_format_dict)


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
