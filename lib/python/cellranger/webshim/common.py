#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import collections
import copy
import csv
import math
import random
import re
import sys
from collections.abc import Sequence
from typing import Any  # pylint: disable=unused-import

import numpy as np

import cellranger.analysis.clustering as cr_clustering
import cellranger.constants as cr_constants
import cellranger.report as cr_report  # pylint: disable=no-name-in-module
import cellranger.rna.library as rna_library
import cellranger.targeted.utils as cr_tgt_utils
import cellranger.utils as cr_utils
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.report as vdj_report
import cellranger.webshim.constants.gex as ws_gex_constants
import cellranger.webshim.constants.shared as shared_constants
import cellranger.webshim.constants.vdj as ws_vdj_constants
import tenkit.safe_json as tk_safe_json
from cellranger.analysis.multigenome import MultiGenomeAnalysis
from cellranger.analysis.singlegenome import Projection, SingleGenomeAnalysis
from cellranger.reference_paths import get_ref_name_from_genomes
from cellranger.webshim.data import SampleData, generate_counter_barcode_rank_plot_data
from cellranger.webshim.jibes_plotting import _make_rc_vectors, make_color_map, make_histogram_plot
from cellranger.websummary.helpers import get_projection_key
from cellranger.websummary.sample_properties import CountSampleProperties, SampleProperties


def add_prefix(prefix, name):
    if "%s" in name:
        return name % prefix
    else:
        return "_".join([prefix, name]) if prefix else name


def format_name(display_name, prefix, prefixes, prefix_format_func=None):
    if prefix is not None and prefix_format_func is not None:
        prefix = prefix_format_func(prefix)

    # Default multi -> '' if no format func given
    if prefix_format_func is None and prefix == rna_library.MULTI_REFS_PREFIX:
        prefix = ""

    if len(prefixes) > 1 or "%s" in display_name:
        display_name = add_prefix(prefix, display_name)

    # Replace underscores w/ spaces
    display_name = display_name.replace("_", " ")

    # Collapse whitespace
    display_name = re.sub(r"\s+", " ", display_name)

    return display_name


def format_description(description, prefix, prefixes, prefix_format_func=None):
    if prefix is not None and prefix_format_func is not None:
        prefix = prefix_format_func(prefix)

    # Default multi -> '' if no format func given
    if prefix_format_func is None and prefix == rna_library.MULTI_REFS_PREFIX:
        prefix = ""

    if "%s" in description:
        # Escape stray percents
        description = re.sub("%([^s])", "%%\\1", description)

        # Only add the prefix if there are multiple possibilities
        s = str(prefix) if len(prefixes) > 1 and prefix else ""
        description = description % s

    # Collapse whitespace
    description = re.sub(r"\s+", " ", description)

    return description


def format_value(value, format_type):
    if format_type == "string":
        return str(value)

    value = value if not math.isnan(float(value)) else 0

    if format_type == "percent":
        return f"{float(value):.1%}"
    elif format_type == "integer":
        return f"{value:,.0f}"
    elif format_type[0] == "%":
        return format_type % float(value)
    raise Exception(f"Invalid format type: {format_type}")


def lookup_name(data, name):
    name_parts = name.split("/", 1)
    data = data.get(name_parts[0])
    if len(name_parts) == 1 or data is None:
        return data
    else:
        return lookup_name(data, name_parts[1])


def add_alarm(value, formatted_value, level, alarm_dict, alarms):
    # The numeric value may have become a string via JSON ('NaN', 'Infinity', '-Infinity')
    # Note that by default python's json encoder outputs a bare NaN symbol which is even worse
    test = 'float("{}") {}'.format(str(value), alarm_dict["test"])
    raised = eval(test, {}, {})
    if raised:
        title = alarm_dict["title"]
        message = alarm_dict["message"]
        alarms[title] = {
            "title": title,
            "message": message,
            "value": formatted_value,
            "level": level,
        }
    return raised


def add_alarms(data, alarms_dict, alarms):
    # Prefixes on alarm metrics not implemented
    full_name = alarms_dict["name"]
    value = lookup_name(data, full_name)
    if value is None:
        return

    formatted_value = format_value(value, alarms_dict["format"])
    error_dict = alarms_dict.get(shared_constants.ALARM_ERROR)
    if error_dict is not None and add_alarm(
        value, formatted_value, shared_constants.ALARM_ERROR, error_dict, alarms
    ):
        return
    warn_dict = alarms_dict.get(shared_constants.ALARM_WARN)
    if warn_dict is not None:
        add_alarm(value, formatted_value, shared_constants.ALARM_WARN, warn_dict, alarms)


def add_table_rows(
    data: dict, name: str, metric_dict: dict, rows, style_func, target_func, prefixes=None
):
    """Add rows to a table.

    Args:
        data: Dict containing summary metrics
        name: Metric name (key in data dict)
        metric_dict: Metric display definition dict
    """
    format_type = metric_dict["format"]
    display_name = metric_dict.get("display_name", name)
    description = metric_dict.get("description", display_name)

    if not prefixes:
        prefixes = [None]

    values = []
    for prefix in prefixes:
        # Find the metric and construct a table row
        full_name = add_prefix(prefix, name)
        value = lookup_name(data, full_name)

        if value is not None:
            prefix_format_func = metric_dict.get("prefix_format")
            formatted_name = format_name(
                display_name, prefix, prefixes, prefix_format_func=prefix_format_func
            )
            formatted_description = format_description(
                description, prefix, prefixes, prefix_format_func=prefix_format_func
            )
            formatted_value = format_value(value, format_type)

            style = style_func(metric_dict, value)

            rows.append(
                [
                    {
                        "v": formatted_name,
                        "f": formatted_description,
                        "s": style,
                    },
                    {
                        "v": target_func(metric_dict),
                        "s": style,
                    },
                    {
                        "v": formatted_value,
                        "s": style,
                        "r": value,
                    },
                ]
            )
            values.append(value)


def _is_metric_hidden(sample_properties, metric_dict):
    if "hidden" in metric_dict:
        return eval(metric_dict["hidden"], sample_properties, sample_properties)
    return False


def build_tables(
    sample_properties,
    table_dicts,
    alarm_table_dicts,
    sample_data,
    style_func=lambda *args: "",
    target_func=lambda *args: "",
    metric_headers=[],
    all_prefixes={},
):
    tables = []
    for table_dict in table_dicts:
        table_name = table_dict["name"]
        table_metrics = table_dict["metrics"]

        rows = []
        for metric_dict in table_metrics:
            name = metric_dict["name"]
            prefix = metric_dict.get("prefix")

            if _is_metric_hidden(sample_properties, metric_dict):
                continue

            if prefix is not None:
                # Search sample properties for prefix values
                prefixes = sample_properties.get(prefix)

                # Otherwise search the given prefix values
                if not prefixes:
                    prefixes = all_prefixes.get(prefix)

                # Apply the prefix filter for this metric
                if "prefix_filter" in metric_dict:
                    if prefixes is not None:
                        prefixes = [x for x in prefixes if x in metric_dict["prefix_filter"]]
                    else:
                        sys.stderr.write(
                            f"Warning: no metric prefix values to filter for metric {name}. The prefix name {prefix} is probably wrong.\n"
                        )

                add_table_rows(
                    sample_data.summary,
                    name,
                    metric_dict,
                    rows,
                    style_func,
                    target_func,
                    prefixes=prefixes,
                )
            else:
                add_table_rows(
                    sample_data.summary, name, metric_dict, rows, style_func, target_func
                )

        if len(rows) > 0:
            tables.append(
                {
                    "name": table_name,
                    "rows": rows,
                    "headers": metric_headers,
                }
            )

        # TODO: Temporary hack to view antibody rank plot, until the new websummary arrives
        if (len(rows) == 0 and table_name == "Mapping") or (
            len(rows) == 0 and table_name == "Cells"
        ):
            tables.append(
                {
                    "name": table_name,
                    "rows": rows,
                    "headers": metric_headers,
                }
            )

    # Build alarm table
    alarms = collections.OrderedDict()
    for alarm_dict in alarm_table_dicts:
        add_alarms(sample_data.summary, alarm_dict, alarms)

    return tables, list(alarms.values())


def convert_numpy_array_to_line_chart(array, ntype):
    """Given an array of values, sort them in descending order.

    Args:
        array: a numpy array
        ntype: a datatype to return the array values as

    Returns:
        the x,y coordinates for a line that could represent the sorted order by
        collapsing identical Y values into a single range of ranks.

    Example:
        .. code-block:: python

            array = np.array([1,1,1,2,3,3,3,3,3,4,5,9,20])
            convert_numpy_array_to_line_chart(array, float)
            # Result
            [[0, 20.0],
            [1, 9.0],
            [2, 5.0],
            [3, 4.0],
            [4, 3.0],
            [8, 3.0],
            [9, 2.0],
            [10, 1.0],
            [12, 1.0]]
    """
    array = np.sort(array)[::-1]
    rows = []
    previous_count = None
    for (index,), count in np.ndenumerate(array):
        if index in (0, len(array) - 1):
            rows.append([index, ntype(count)])
        elif previous_count != count:
            assert previous_count is not None
            previous_index = rows[-1][0]
            if previous_index != index - 1:
                rows.append([index - 1, ntype(previous_count)])
            rows.append([index, ntype(count)])
        previous_count = count
    return rows


def _plot_barcode_rank(chart, counts, num_cells):
    """Generate a generic barcode rank plot."""
    rows = convert_numpy_array_to_line_chart(counts, int)

    for row in rows:
        index, count = row[0], row[1]
        if index < num_cells:
            series_list = [chart["data"][0]]
        elif index == num_cells:
            # Connect the two lines
            series_list = [chart["data"][0], chart["data"][1]]
        else:
            series_list = [chart["data"][1]]

        for series in series_list:
            series["x"].append(index)
            series["y"].append(count)

    # Handle case where there is no background
    bg_series = chart["data"][1]
    if len(bg_series["x"]) == 0:
        bg_series["x"].append(0)
        bg_series["y"].append(0)

    return chart


def build_plot_data_dict(plot_segment, counts, show_name_and_hover=True, color=None):
    """Construct the data for a plot segment by appropriately slicing the counts.

    Args:
        plot_segment: BarcodeRankPlotSegment containing [start, end)
          of the segment, the cell density and legend visibility option
        counts: Reverse sorted UMI counts for all barcodes.
        show_name_and_hover: boolean whether to add hover and name text specific
            to the classic CR rank plot (the user might want to set these
            themselves for alternate rank plots)
    """
    # -1 for continuity between two charts
    start = max(0, plot_segment.start - 1)
    end = plot_segment.end
    plot_rows = convert_numpy_array_to_line_chart(counts[start:end], int)

    if color is None:
        color = shared_constants.BC_PLOT_CMAP(plot_segment.cell_density)

    data_dict = {
        "x": [],
        "y": [],
        "type": "scattergl",
        "mode": "lines",
        "line": {
            "color": color,
            "width": shared_constants.BC_RANK_PLOT_LINE_WIDTH,
        },
        "showlegend": plot_segment.legend,
    }
    offset = 1 + start  # it's a log-log plot, hence the 1
    for index, count in plot_rows:
        data_dict["x"].append(index + offset)
        data_dict["y"].append(count)

    # Handle case where the data is empty
    if len(data_dict["x"]) == 0:
        data_dict["x"].append(0)
        data_dict["y"].append(0)

    if show_name_and_hover:
        name = "Cells" if plot_segment.cell_density > 0 else "Background"
        data_dict["name"] = name
        if plot_segment.cell_density > 0.0:
            n_barcodes = plot_segment.end - plot_segment.start
            n_cells = int(round(plot_segment.cell_density * n_barcodes))
            hoverbase = f"{100 * plot_segment.cell_density:.0f}% Cells<br>({n_cells}/{n_barcodes})"
        else:
            hoverbase = "Background"
        hover = []
        for index, count in plot_rows:
            rank = index + offset
            hover.append(f"{hoverbase}<br>Rank {rank:,}<br>UMIs {count:,}")
        data_dict["hoverinfo"] = "text"
        data_dict["text"] = hover

    return data_dict


def _plot_segmented_barcode_rank(chart, counts, plot_segments):
    """Generate the RNA counter barcode rank plot.

    Inputs:
        - chart: chart element to populate data
        - counts: UMI counts reverse sorted
        - plot_segments: A list of BarcodeRankPlotSegments
    """
    for segment in plot_segments:
        chart["data"].append(build_plot_data_dict(segment, counts))

    return chart


def plot_basic_barcode_rank(
    chart, cell_barcodes: set[bytes], barcode_summary, genomes, lib_prefix, restrict_barcodes=None
):
    """Generate a basic RNA counter barcode rank plot without depending on SampleData/SampleProperties.

    Args:
        chart: chart element to populate data
        cell_barcodes: set of cell barcodes as bytes
        barcode_summary: barcode summary from the barcode_summary.h5
        lib_prefix: The library prefix to create the plot for
        restrict_barcodes: Optional list of cell barcodes to restrict to
    """
    genome = rna_library.MULTI_REFS_PREFIX

    if (
        lib_prefix
        == rna_library.get_library_type_metric_prefix(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
        and len(genomes) == 1
    ):
        genome = genomes[0]

    key = cr_utils.format_barcode_summary_h5_key(
        lib_prefix,
        genome,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
    )

    if key not in barcode_summary:
        return None

    counts_per_bc, plot_segments = generate_counter_barcode_rank_plot_data(
        cell_barcodes, barcode_summary, key, restrict_barcodes=restrict_barcodes
    )
    return _plot_segmented_barcode_rank(chart, counts_per_bc, plot_segments)


def plot_barcode_rank(chart, sample_properties, sample_data):
    """Generate the RNA counter barcode rank plot."""
    # TODO: The old PD plots use this function and pass things around as dictionaries,
    # so we convert here if this function was called from the PD code
    if isinstance(sample_properties, dict):
        sample_properties = CountSampleProperties(
            sample_id=sample_properties.get("sample_id", None),
            sample_desc=sample_properties.get("sample_desc", None),
            genomes=sample_properties.get("genomes", None),
        )

    assert isinstance(sample_data, SampleData)
    assert isinstance(sample_properties, SampleProperties)
    if (
        not isinstance(sample_properties, CountSampleProperties)
        or sample_data.barcode_summary is None
        or sample_data.cell_barcodes is None
    ):
        return None

    if len(sample_properties.genomes) == 0:
        return None

    # UMI counts per BC across all genomes present
    if len(sample_properties.genomes) > 1:
        genome = rna_library.MULTI_REFS_PREFIX
    else:
        genome = sample_properties.genomes[0]

    gex_prefix = rna_library.get_library_type_metric_prefix(
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE
    )
    ab_prefix = rna_library.get_library_type_metric_prefix(rna_library.ANTIBODY_LIBRARY_TYPE)
    gex_key = cr_utils.format_barcode_summary_h5_key(
        gex_prefix,
        genome,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
    )
    ab_key = cr_utils.format_barcode_summary_h5_key(
        ab_prefix,
        rna_library.MULTI_REFS_PREFIX,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
    )

    if gex_key in sample_data.barcode_summary:
        counts_per_bc, plot_segments = sample_data.counter_barcode_rank_plot_data(gex_key)
        return _plot_segmented_barcode_rank(chart, counts_per_bc, plot_segments)
    elif ab_key in sample_data.barcode_summary:  # in case there was only Antibody library
        counts_per_bc, plot_segments = sample_data.counter_barcode_rank_plot_data(ab_key)
        return _plot_segmented_barcode_rank(chart, counts_per_bc, plot_segments)
    else:
        # Not guaranteed to exist, depending on pipeline
        return None


# TODO: This unneccesary argument is needed because of how the PD code is structured.
def plot_vdj_barcode_rank(chart, sample_properties, sample_data):
    """Generate the VDJ barcode rank plot."""
    if not sample_data.cell_barcodes or sample_data.vdj_barcode_support is None:
        return None

    counts_per_bc, plot_segments = sample_data.vdj_barcode_rank_plot_data()
    return _plot_segmented_barcode_rank(chart, counts_per_bc, plot_segments)


def plot_clonotype_table(chart, sample_properties, sample_data):
    if sample_data.vdj_clonotype_summary is None:
        return None

    clonotypes = sample_data.vdj_clonotype_summary.iloc[0:10]

    # This column used to be called 'cdr3s'; allow the webshim to work on older data
    cdr3_aa_col = "cdr3s_aa"
    if cdr3_aa_col not in clonotypes:
        cdr3_aa_col = "cdr3s"

    col_defs = collections.OrderedDict(
        [
            (
                "clonotype_id",
                {
                    "label": "Clonotype ID",
                    "format": "string",
                    "title": "Clonotype ID",
                    "style": "text-align: left",
                },
            ),
            (
                cdr3_aa_col,
                {
                    "label": "CDR3s",
                    "format": "string",
                    "title": "CDR3s in clonotype",
                    "style": "text-align: left",
                },
            ),
            (
                "frequency",
                {
                    "label": "Frequency",
                    "format": "integer",
                    "title": "Number of cells with clonotype",
                    "style": "text-align: right",
                },
            ),
            (
                "proportion",
                {
                    "label": "Proportion",
                    "format": "%0.4f",
                    "title": "Fraction of cell with clonotype",
                    "style": "text-align: right",
                },
            ),
        ]
    )

    cols = []
    for name, col_def in col_defs.items():
        if name not in clonotypes:
            raise ValueError(f"Column not found in clonotype summary: {name}")
        cols.append(
            {
                "label": col_def["label"],
                "title": col_def["title"],
            }
        )

    rows = []
    for _, cl_row in clonotypes.iterrows():
        row = []
        for col_name, col_def in col_defs.items():
            value = cl_row[col_name]
            formatted_value = format_value(value, col_def["format"])

            # Make the CDR3 list bit more readable
            formatted_value = formatted_value.replace(";", "; ")

            row.append(
                {
                    "v": tk_safe_json.json_sanitize(value),
                    "f": formatted_value,
                    "s": col_def["style"],
                }
            )
        rows.append(row)

    chart["table"].update({"rows": rows, "cols": cols})

    return chart


def plot_histogram_metric(
    chart,
    sample_properties,
    sample_data,
    *,
    metric_name,
    order_by=shared_constants.HISTOGRAM_METRIC_DEFAULT_ORDERING,
    **kwargs,
):
    """Plot a HistogramMetric from the summary json."""
    summary_data = sample_data.summary
    items = list(summary_data.get(metric_name, {}).items())
    if len(items) < 1:
        return None

    ordering = order_by

    if ordering == shared_constants.HISTOGRAM_METRIC_ORDER_INTEGER_BIN:
        items.sort(key=lambda x: nullable_int_sort_key(x[0]))

    elif ordering == shared_constants.HISTOGRAM_METRIC_ORDER_DECREASING_FREQUENCY:
        items.sort(key=lambda x: nullable_int_sort_key(x[1], negate=True))

    elif ordering == shared_constants.HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION:
        items.sort(key=lambda x: -convert_to_float_gracefully(x[1]))

    x, y = zip(*items)
    chart["data"][0].update({"x": x, "y": y})

    return chart


def plot_barnyard_barcode_counts(chart, sample_properties, sample_data):
    analysis = sample_data.get_analysis(MultiGenomeAnalysis)
    if analysis is None:
        return None

    chart["data"] = []
    for label_info in ws_gex_constants.GEM_CALL_LABELS:
        name = label_info["label"]
        name = name.replace("genome0", analysis.result["genome0"])
        name = name.replace("genome1", analysis.result["genome1"])
        chart["data"].append(
            {
                "x": [],
                "y": [],
                "name": name,
                "mode": "markers",
                "marker": {"color": label_info["color"], "opacity": 0.15},
            }
        )

    call_to_series = {v["key"]: i for i, v in enumerate(ws_gex_constants.GEM_CALL_LABELS)}
    for count0, count1, call in zip(
        analysis.result["count0"], analysis.result["count1"], analysis.result["call"]
    ):
        series = chart["data"][call_to_series[call]]
        series["x"].append(int(count0))
        series["y"].append(int(count1))

    chart["layout"]["xaxis"] = {
        "title": "{} UMI counts".format(analysis.result["genome0"]),
        "rangemode": "tozero",
        "autorange": True,
        "fixedrange": False,
    }
    chart["layout"]["yaxis"] = {
        "title": "{} UMI counts".format(analysis.result["genome1"]),
        "rangemode": "tozero",
        "autorange": True,
        "fixedrange": False,
    }

    return chart


def plot_preprocess(analyses):
    if analyses is None or len(analyses) == 0:
        return None

    sg_analyses = [an for an in analyses if isinstance(an, SingleGenomeAnalysis)]
    sg_analysis = sg_analyses[0] if len(sg_analyses) > 0 else None
    if sg_analysis is None or sg_analysis.is_zero_matrix():
        return None

    # Limit the number of K-means clusterings displayed to limit the HTML filesize
    new_clusterings = {}
    for key, clu in sg_analysis.clusterings.items():
        if not (
            clu.clustering_type == cr_clustering.CLUSTER_TYPE_KMEANS
            and clu.num_clusters > ws_gex_constants.MAX_WEBSHIM_KMEANS_K
        ):
            new_clusterings[key] = clu
    sg_analysis.clusterings = new_clusterings

    return analyses


def load_sample_data(sample_properties, sample_data_paths, projections: Sequence[Projection]):
    return SampleData(sample_properties, sample_data_paths, plot_preprocess, projections)


def plot_tsne(chart, sample_properties, sample_data):
    """Plot cells in t-SNE space, colored by clustering label."""
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None

    matrix = sample_data.get_analysis(SingleGenomeAnalysis).matrix
    library_types = matrix.get_library_types()
    if rna_library.GENE_EXPRESSION_LIBRARY_TYPE in library_types:
        key = get_projection_key(rna_library.GENE_EXPRESSION_LIBRARY_TYPE, 2)
    elif rna_library.ANTIBODY_LIBRARY_TYPE in library_types:
        key = get_projection_key(rna_library.ANTIBODY_LIBRARY_TYPE, 2)

    args = [
        analysis.get_tsne(key=key).transformed_tsne_matrix,
        ws_gex_constants.TSNE_CLUSTER_DESCRIPTION,
        1,
        2,
    ]
    return clustering_plot_func(chart, sample_properties, sample_data, plot_dimensions, args)


def plot_dimensions(
    chart,
    transformed_matrix,
    description,
    pc1,
    pc2,
    clip=None,
    clustering=None,
    diff_expr=None,
    values=None,
    original_cluster_sizes=None,
):
    """Plot cells in a 2-d space, colored by clustering label."""
    assert (values is not None and clustering is None) or (
        clustering is not None and values is None and original_cluster_sizes is not None
    )

    if clustering is not None:
        values = clustering.clusters

    if original_cluster_sizes is None:
        value_freqs = np.bincount(values)[1:].tolist()
    else:
        # Cluster size array starts w/ cluster1 at index0; make it equivalent to np.bincount on arbitrary values
        value_freqs = [0] + original_cluster_sizes.tolist()

    n, m = transformed_matrix.shape
    if m < max(pc1, pc2):
        return None

    max_value = values.max()

    chart["data"] = [None] * max_value
    for value in np.unique(values):
        chart["data"][value - 1] = {
            "x": [],
            "y": [],
            "name": "{} {} - {} cells".format(
                description,
                format_value(value, "integer"),
                format_value(value_freqs[value], "integer"),
            ),
            "hoverinfo": "name",
            "mode": "markers",
            "type": "scattergl",
            "marker": {
                "size": 4,
            },
        }

    for i in range(n):
        r1 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc1 - 1]
        r2 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc2 - 1]
        value = values[i]

        series = chart["data"][value - 1]
        series["x"].append(r1)
        series["y"].append(r2)

    if clip is not None:
        xmin, xmax = np.percentile(transformed_matrix[:, pc1 - 1], clip)
        ymin, ymax = np.percentile(transformed_matrix[:, pc2 - 1], clip)
        chart["layout"]["xaxis"] = {
            "range": [xmin, xmax],
        }
        chart["layout"]["yaxis"] = {
            "range": [ymin, ymax],
        }

    chart["config"] = shared_constants.CHARTS_PLOTLY_MOVABLE_CONFIG
    return chart


def plot_dimensions_color(chart, transformed_matrix, values, description, vmin, vmax, pc1, pc2):
    _, m = transformed_matrix.shape
    if m < max(pc1, pc2):
        return None

    series = chart["data"][0]

    index_order = list(range(transformed_matrix.shape[0]))
    random.shuffle(index_order)
    for i in index_order:
        r1 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc1 - 1]
        r2 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc2 - 1]
        value = int(values[i])
        text = "{}: {}".format(description, format_value(value, "integer"))

        value = min(max(value, vmin), vmax)

        series["x"].append(r1)
        series["y"].append(r2)
        series["marker"]["color"].append(value)
        series["text"].append(text)

    return chart


def plot_tsne_totalcounts(chart, sample_properties, sample_data):
    """Plot cells colored by total counts."""
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if not analysis:
        return None

    matrix = sample_data.get_analysis(SingleGenomeAnalysis).matrix
    library_types = matrix.get_library_types()
    if rna_library.GENE_EXPRESSION_LIBRARY_TYPE in library_types:
        library_type = rna_library.GENE_EXPRESSION_LIBRARY_TYPE
        key = get_projection_key(library_type, 2)
    elif rna_library.ANTIBODY_LIBRARY_TYPE in library_types:
        library_type = rna_library.ANTIBODY_LIBRARY_TYPE
        key = get_projection_key(library_type, 2)

    matrix = analysis.matrix.select_features_by_type(library_type)
    reads_per_bc = matrix.get_counts_per_bc()
    vmin, vmax = np.percentile(reads_per_bc, ws_gex_constants.TSNE_TOTALCOUNTS_PRCT_CLIP)

    return plot_dimensions_color(
        chart,
        analysis.get_tsne(key=key).transformed_tsne_matrix,
        reads_per_bc,
        ws_gex_constants.TSNE_TOTALCOUNTS_DESCRIPTION,
        vmin,
        vmax,
        1,
        2,
    )


def _plot_differential_expression(
    chart, analysis, clustering=None, diff_expr=None, original_cluster_sizes=None
):
    n_clusters = clustering.clusters.max()

    # Get the union of top DE genes
    top_genes = set()

    # Limit the number of entries in the DE table
    n_genes = int(np.floor(float(ws_gex_constants.MAX_DE_TABLE_ENTRIES) / (n_clusters**2)))
    if n_genes < 1:
        n_genes = 1
    elif n_genes > ws_gex_constants.MAX_TOP_N_GENES:
        n_genes = ws_gex_constants.MAX_TOP_N_GENES

    cols = [
        {"type": "string", "label": "Gene ID"},
        {"type": "string", "label": "Gene name"},
    ]

    for i in range(n_clusters):
        # Filter genes by mean count and sort by log2 fold-change, descending
        means = diff_expr.data[:, 0 + 3 * i]
        log2fcs = diff_expr.data[:, 1 + 3 * i]

        keep_indices = np.flatnonzero(means >= ws_gex_constants.TOP_DE_GENES_MIN_MEAN)
        top_gene_indices = keep_indices[log2fcs[keep_indices].argsort()[::-1]][:n_genes]

        for j in top_gene_indices:
            top_genes.add(analysis.matrix.int_to_feature_id(j))

        cols.append(
            {
                "type": "number",
                "label": "L2FC",
                "title": "Log2 fold-change in cluster %d vs other cells" % (i + 1),
            }
        )
        cols.append(
            {
                "type": "number",
                "label": "p-value",
                "title": "Adjusted p-value of differential expression in cluster %d" % (i + 1),
            }
        )

    rows = []
    for gene_id in top_genes:
        i = analysis.matrix.feature_id_to_int(gene_id)
        gene_name = analysis.matrix.feature_id_to_name(gene_id)

        row = [gene_id, gene_name]
        for j in range(n_clusters):
            log2fc = diff_expr.data[i, 1 + (3 * j)]
            adj_p_value = diff_expr.data[i, 2 + (3 * j)]

            if log2fc <= 0 or adj_p_value >= ws_gex_constants.PVALUE_DEEMPHASIS_CUTOFF:
                style = "#DDD"
            else:
                style = "#000"

            row.append(
                {
                    "v": tk_safe_json.json_sanitize(log2fc),
                    "f": format_value(log2fc, "%.2f"),
                    "s": style,
                }
            )
            row.append(
                {
                    "v": tk_safe_json.json_sanitize(adj_p_value),
                    "f": format_value(adj_p_value, "%.0e"),
                    "s": style,
                }
            )

        rows.append(row)

    # Sort by log2fc, descending, in first cluster
    if n_clusters > 0:
        rows = sorted(rows, key=lambda row: row[2]["v"], reverse=True)

    chart["table"].update({"rows": rows, "cols": cols})
    return chart


def plot_differential_expression(chart, sample_properties, sample_data):
    sg_analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    return clustering_plot_func(
        chart, sample_properties, sample_data, _plot_differential_expression, [sg_analysis]
    )


def clustering_plot_func(chart, sample_properties, sample_data, plot_func, args=[], kwargs={}):
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None

    new_charts = []
    for clustering_key, clustering in analysis.clusterings.items():
        kwargs["clustering"] = clustering
        kwargs["original_cluster_sizes"] = sample_data.original_cluster_sizes[clustering_key]
        kwargs["diff_expr"] = analysis.differential_expression[clustering_key]
        new_chart = plot_func(copy.deepcopy(chart), *args, **kwargs)
        if new_chart is not None:
            new_chart["filters"] = {ws_gex_constants.CLUSTERS_FILTER_TITLE: clustering.description}
            new_charts.append(new_chart)
    return new_charts


def make_chart_filters(analyses):
    if analyses is None:
        return {}
    sg_analyses = [an for an in analyses if isinstance(an, SingleGenomeAnalysis)]
    if len(sg_analyses) == 0 or analyses[0] is None:
        return {}
    assert len(sg_analyses) == 1
    analysis = sg_analyses[0]

    filter_values = [
        x.description for x in cr_clustering.sort_clusterings(analysis.clusterings.values())
    ]

    return {
        ws_gex_constants.CLUSTERS_FILTER_TITLE: {
            "values": filter_values,
            "selected": filter_values[0],
        },
    }


def plot_subsampled_scatterplot_metric(
    chart,
    sample_properties,
    sample_data,
    *,
    subsample_type="raw_rpc",
    ref_prefix=None,
    references=None,
    is_targeted=False,
    show_targeted_only=False,
    metric_suffix=None,
    multi_genome_only=False,
    show_multi_genome_only=False,
):
    """Modifies chart data entry to add traces for subsampled metrics with specified suffixes.

    Metrics must take the form `<reference>_<subsample_type>_<subsample_depth>_<metric_suffix>`,
    where metric suffix is specified as a kwarg

    Args:
        metric_suffix (str): suffix for the subsampled metric given a metric of
            the form `<reference>_<subsample_type>_<subsample_depth>_<metric_suffix>`
    """
    summary_data = sample_data.summary or {}

    # Just use ref_prefix if specified; otherwise try to construct it out of the set of references
    if ref_prefix is None:
        if references is not None:
            ref_prefix = "(" + "|".join(references) + ")_"
        else:
            ref_prefix = "(.+)_"

    # if targeted, determine which suffix we're expecting
    if is_targeted:
        targeting_groups_to_plot = [cr_constants.ON_TARGET_SUBSAMPLE]
        if not show_targeted_only:
            targeting_groups_to_plot.append(cr_constants.OFF_TARGET_SUBSAMPLE)
        targeting_group_suffix = "_({})".format("|".join(targeting_groups_to_plot))
    else:
        targeting_group_suffix = ""

    # Regular expression to match <reference>_<subsample_type>_<subsample_depth>_<metric_suffix> pattern
    metric_pattern = (
        rf"^{ref_prefix}({subsample_type})_([0-9]+)_{metric_suffix}{targeting_group_suffix}$"
    )

    metric_search_results = [re.search(metric_pattern, key) for key in summary_data.keys()]

    # Only return chart when multiple references are present for multi-genome only metrics
    if multi_genome_only:
        references = {result.group(1) for result in metric_search_results if result is not None}
        if len(references) <= 2:
            return None

    # Find relevant metrics and extract data
    points = []

    for val, search_result in zip(summary_data.values(), metric_search_results):
        if search_result and not isinstance(val, float):
            val = float(val)
        # Only count corresponding keys and non-zero values
        if search_result and val > 0:
            genome = search_result.group(1)
            subsample_type = search_result.group(2).replace("_", " ")
            if show_multi_genome_only and (genome != "multi"):
                continue
            subsample_depth = int(search_result.group(3))

            if is_targeted:
                targeting_group = search_result.group(4)
                if targeting_group not in targeting_groups_to_plot:
                    continue
                subsample_type = f"{subsample_type}_{targeting_group}"

            if genome == "multi":
                if is_targeted:
                    trace_label = cr_tgt_utils.reformat_targeted_label(targeting_group)
                else:
                    trace_label = ""
            elif is_targeted:
                trace_label = f"{genome}_{cr_tgt_utils.reformat_targeted_label(targeting_group)}"
            else:
                trace_label = genome

            points.append(
                (
                    subsample_depth,
                    val,
                    trace_label,
                    genome,
                    subsample_type,
                )
            )

    # Don't produce chart if there's no data
    if len(points) == 0:
        return None

    # Sort extracted points by subsample depth (makes line plot connect properly)
    sorted_points = sorted(points, key=lambda x: (x[2], x[0]))

    # Build up a set of traces for <reference>_<subsample_type> pairs
    traces = {}

    for point in sorted_points:
        depth, value, trace_label, genome, subsample_type = point
        trace_key = (genome, subsample_type)

        if trace_key not in traces:
            # add (0,0)
            traces[trace_key] = {
                "x": [0, depth],
                "y": [0, value],
                "name": trace_label,
                "mode": "lines",
                "line": {
                    "width": 3,
                },
            }
        else:
            traces[trace_key]["x"].append(depth)
            traces[trace_key]["y"].append(value)

    # Set extents of saturation line
    if (
        chart.get("layout")
        and chart["layout"].get("shapes")
        and chart["layout"]["shapes"][0]["type"] == "line"
    ):
        chart["layout"]["shapes"][0]["x1"] = max(max(trace["x"]) for trace in traces.values())

    # Add new data into chart
    chart["data"] = [v for _, v in sorted(traces.items())]
    return chart


def build_charts(sample_properties, chart_dicts, sample_data, all_prefixes={}, module=None):
    modules = [module, globals()] if module else [globals()]

    filters = make_chart_filters(sample_data.analyses)

    charts = []
    for chart_dict in chart_dicts:
        chart_dict = copy.deepcopy(chart_dict)
        function = chart_dict.pop("function")
        f = None
        for mod in modules:
            f = mod.get(function)
            if f is not None:
                break
        assert f

        if function is not None and f is None:
            raise ValueError(f'Could not find webshim chart function "{function}"')

        kwargs: dict[str, Any] = chart_dict.pop("kwargs", {})

        kwargs_prefixes: list[str] = chart_dict.pop("kwargs_prefixes", [])
        for prefix in kwargs_prefixes:
            kwargs[prefix] = all_prefixes[prefix]

        new_chart_obj = f(chart_dict, sample_properties, sample_data, **kwargs)
        if new_chart_obj is None:
            continue

        new_charts = new_chart_obj if isinstance(new_chart_obj, list) else [new_chart_obj]
        charts.extend(new_charts)

    return charts, filters


def filter_vdj_prefixes(all_prefixes, sample_properties):
    """Only get subset of metric prefix values."""
    chain_filter = sample_properties.chain_type
    if chain_filter is None:
        return all_prefixes

    # NOTE: Assumes chains (TRA, TRB, etc) are prefixed with the chain_type (TR or IG)
    result = {}
    for key, values in all_prefixes.items():
        # Only filter any prefix that is a candidate for filtering
        # (i.e., contains some values prefixed by the selected chain type)
        if values is not None and any(v.startswith(chain_filter) for v in values):
            result[key] = [
                v
                for v in values
                if v.startswith(chain_filter) or v == rna_library.MULTI_REFS_PREFIX
            ]
        else:
            result[key] = values
    return result


def filter_vdj_alarms(all_alarms, sample_properties):
    """Only get subset of metric alarms."""
    chain_filter = sample_properties.chain_type
    if chain_filter is None:
        return all_alarms

    result = []
    for alarm in all_alarms:
        # No filters specified; don't filter
        if "filters" not in alarm:
            result.append(alarm)
            continue

        for f in alarm["filters"]:
            if f.get("chain_type") == chain_filter:
                result.append(alarm)

    return result


def get_constants_for_pipeline(pipeline, sample_properties):
    """Get the appropriate metrics/alarms/charts for a pipeline."""
    if pipeline == shared_constants.PIPELINE_VDJ:
        metrics, alarms, charts = (
            ws_vdj_constants.METRICS,
            ws_vdj_constants.METRIC_ALARMS,
            ws_vdj_constants.CHARTS,
        )

        metric_prefixes = filter_vdj_prefixes(
            vdj_report.VdjReporter().get_all_prefixes(), sample_properties
        )

        alarms = filter_vdj_alarms(alarms, sample_properties)

    else:
        metrics, alarms, charts = (
            ws_gex_constants.METRICS,
            ws_gex_constants.METRIC_ALARMS,
            ws_gex_constants.CHARTS,
        )

        metric_prefixes = cr_report.Reporter().get_all_prefixes()

    return metrics, alarms, charts, metric_prefixes


def get_custom_features(sample_data):
    """Infer the set of distinct custom feature types present in a dataset."""
    # Add the "custom feature type" prefix using the sample data
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)

    # Adding a dummy value here suppresses extraneous web summary dashboards
    custom_features = ["dummy"]

    if (
        analysis
        and analysis.matrix.feature_ref.get_count_of_feature_type(rna_library.CUSTOM_LIBRARY_TYPE)
        > 0
    ):
        custom_features.append(rna_library.CUSTOM_LIBRARY_TYPE)
    return custom_features


def get_genomes(sample_data):
    """Infer the set of genomes present in a dataset."""
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis:
        return [x for x in analysis.matrix.get_genomes() if x != ""]
    else:
        return ["dummy"]


def build_web_summary_json(sample_properties, sample_data, pipeline):
    """Build a websummary json.

    Args:
        sample_properties (SampleProperties): object
        sample_data (SampleData): class
        pipeline (str): ?
    """
    view = copy.deepcopy(sample_properties.__dict__)
    sample_properties_as_dict = copy.deepcopy(view)

    metrics, alarms, charts, all_prefixes = get_constants_for_pipeline(pipeline, sample_properties)

    all_prefixes["custom_features"] = get_custom_features(sample_data)

    tables, alarms = build_tables(
        sample_properties_as_dict, metrics, alarms, sample_data, all_prefixes=all_prefixes
    )
    if tables:
        view["tables"] = tables
    if alarms:
        view["alarms"] = alarms

    all_prefixes["references"] = list(set(all_prefixes["references"] + get_genomes(sample_data)))

    charts, filters = build_charts(
        sample_properties_as_dict, charts, sample_data, all_prefixes=all_prefixes
    )
    if charts:
        view["charts"] = charts
    if filters:
        view["filters"] = filters

    # Selected metrics that the web summary template needs
    info = build_info_dict(sample_properties, sample_data, pipeline)
    if info:
        view["info"] = info

    return view


def build_info_dict(sample_properties, sample_data, pipeline):
    """Add miscellaneous metrics required by the web summary template."""
    if sample_data.summary is None:
        return None
    info = {}

    info["chemistry_description"] = sample_data.summary.get("chemistry_description")

    info["references"] = []

    if pipeline in [shared_constants.PIPELINE_AGGR, shared_constants.PIPELINE_REANALYZE]:
        genomes = sample_properties.genomes
        if genomes is not None:
            info["references"].append(
                {
                    "type": cr_constants.REFERENCE_TYPE,
                    "name": get_ref_name_from_genomes(genomes),
                }
            )

    else:
        # Find all references in the summary
        reference_metric_prefixes = [
            cr_constants.REFERENCE_METRIC_PREFIX,
            vdj_constants.REFERENCE_METRIC_PREFIX,
        ]

        for prefix in reference_metric_prefixes:
            type_metric = f"{prefix}{cr_constants.REFERENCE_TYPE_KEY}"

            if type_metric not in sample_data.summary:
                continue

            ref_type = sample_data.summary[type_metric]

            name_metric = f"{prefix}{cr_constants.REFERENCE_GENOMES_KEY}"
            if name_metric not in sample_data.summary:
                raise ValueError(
                    f"Reference metadata metric {name_metric} not found in metrics summary."
                )
            ref_name = sample_data.summary.get(name_metric)

            if isinstance(ref_name, list):
                ref_name = get_ref_name_from_genomes(ref_name)

            info["references"].append({"type": ref_type, "name": ref_name})

    return info


def build_metrics_summary_csv(filename, sample_properties, sample_data, pipeline):
    metrics, alarms, _, all_prefixes = get_constants_for_pipeline(pipeline, sample_properties)

    sample_properties_as_dict = copy.deepcopy(sample_properties.__dict__)

    all_prefixes["custom_features"] = get_custom_features(sample_data)

    tables, _ = build_tables(
        sample_properties_as_dict, metrics, alarms, sample_data, all_prefixes=all_prefixes
    )
    if not tables:
        sys.stderr.write("No metrics tables were generated, skipping CSV generation.\n")
        return

    csv_metrics = collections.OrderedDict()
    for table in tables:
        if not table:
            continue
        for metric, _, value in table["rows"]:
            if isinstance(metric, dict):
                metric = metric["v"]
            if isinstance(value, dict):
                value = value["v"]
            if metric not in csv_metrics:
                csv_metrics[metric] = value

    with open(filename, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(csv_metrics.keys())
        writer.writerow(csv_metrics.values())


def nullable_int_sort_key(s, negate=False):
    """Return a sort key for s that will sort by int(s).

    Values of s that cannot be converted to an integer will be sorted
    to the begining.
    """
    try:
        if not negate:
            return (True, int(s))
        else:
            return (True, -int(s))
    except ValueError:
        return (False, 0)


def convert_to_float_gracefully(s):
    try:
        return float(s)
    except ValueError:
        return sys.float_info.max


# the shared_resources implementation in python doesn't interop with Rust websummary at all
# need to just define a key to put the data under and use that


def make_antibody_histograms(ab_matrix):
    ab_counts = np.log10(ab_matrix.m.toarray() + 1.0).tolist()

    shared_resources, vectors = _make_rc_vectors(ab_counts)
    vector_names = ab_matrix.feature_ids_map.keys()
    color_map = make_color_map(vector_names, jibes_plot=False)

    histogram_data = make_histogram_plot(vector_names, vectors, color_map)

    antibody_histograms = {
        "histogram": histogram_data,
        "_resources": shared_resources,
    }

    return antibody_histograms


def plot_umi_depth(
    chart: dict,
    sample_properties: dict | CountSampleProperties,
    sample_data: SampleData,
):
    """Generate a UMI depth plot.

    Each point is a barcode. X = Reads/UMI for barcode, Y = UMIs per barcode (log). This is useful
    in differentiating between index hopping vs other assay issues when the barcode rank plot
    shows poor separation between cells and background.

    See CELLRANGER-4776 for more detail.
    """
    if sample_data.barcode_summary is None:
        return None

    # TODO: The old PD plots use this function and pass things around as dictionaries,
    # so we convert here if this function was called from the PD code
    if isinstance(sample_properties, dict):
        sample_properties = CountSampleProperties(
            sample_id=sample_properties.get("sample_id", None),
            sample_desc=sample_properties.get("sample_desc", None),
            genomes=sample_properties.get("genomes", None),
        )

    if len(sample_properties.genomes) == 0:
        return None

    # UMI counts per BC across all genomes present
    if len(sample_properties.genomes) > 1:
        genome = rna_library.MULTI_REFS_PREFIX
    else:
        genome = sample_properties.genomes[0]

    gex_prefix = rna_library.get_library_type_metric_prefix(
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE
    )
    gex_umi_key = cr_utils.format_barcode_summary_h5_key(
        gex_prefix,
        genome,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
    )
    gex_reads_key = cr_utils.format_barcode_summary_h5_key(
        gex_prefix,
        genome,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_BC_READ_TYPE,
    )

    # Not guaranteed to exist, depending on pipeline and libraries
    if (
        gex_umi_key not in sample_data.barcode_summary
        or gex_reads_key not in sample_data.barcode_summary
    ):
        return None

    # minimum number of umis per barcode to control the number of points being plotted
    MIN_COUNT = 10

    # load from bc summary
    umis = sample_data.barcode_summary[gex_umi_key][:]
    select = umis >= MIN_COUNT
    if select.sum() == 0:
        return None

    y: np.ndarray = umis[select].astype(np.float64)
    reads: np.ndarray = sample_data.barcode_summary[gex_reads_key][:]
    x = reads[select].astype(np.float64) / y
    data = chart["data"][0]
    data["x"] = x.tolist()
    data["y"] = y.tolist()

    return chart
