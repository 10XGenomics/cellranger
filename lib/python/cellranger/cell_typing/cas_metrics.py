#!/usr/bin/env python3
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#


"""Calculate Cell Typing Metrics."""

import csv
import json
import os
import textwrap
from collections import defaultdict

import altair as alt
import h5py
import numpy as np
import pandas as pd

import cellranger.altair_utils as alt_utils
import cellranger.analysis.constants as cr_analysis_constants
import cellranger.webshim.constants.gex as gex_constants
from cellranger.altair_plot_utils import make_alt_spatial_plot
from cellranger.cell_typing.cas_postprocessing import BARCODE_KEY, COARSE_CELL_TYPES_KEY

alt.data_transformers.disable_max_rows()


def get_cas_cluster_purity(cell_types: str, analysis: str) -> dict:
    """Calculate the cluster purity based on cell annotation cell types and analysis data.

    Args:
        cell_types (str): Path to the cell annotation cell types CSV file.
        analysis (str): Path to the analysis directory.

    Returns:
        dict: A dictionary representing the cluster purity. will sum to 1 for each cluster:
    """
    with open(cell_types) as cell_types_file:
        cell_types = dict(
            (row[BARCODE_KEY], row[COARSE_CELL_TYPES_KEY])
            for row in csv.DictReader(cell_types_file)
        )

    # Get clusters
    analysis_h5 = h5py.File(os.path.join(analysis, "analysis.h5"), "r")
    barcodes = (x.decode("utf-8") for x in analysis_h5["matrix/barcodes"])
    clusters = analysis_h5["/clustering/_gene_expression_graphclust/clusters"]
    clusters_data = dict(zip(barcodes, clusters))
    # Join clusters and cell types
    joined_data = (
        {"barcode": barcode, "clusters": clusters_data[barcode], COARSE_CELL_TYPES_KEY: cell_type}
        for barcode, cell_type in cell_types.items()
        if barcode in clusters_data
    )

    # Calculate the cluster "purity"
    cluster_break_up = defaultdict(lambda: defaultdict(int))
    for data in joined_data:
        cluster = data["clusters"]
        cell_type = data[COARSE_CELL_TYPES_KEY]
        cluster_break_up[cluster][cell_type] += 1

    # Calculate the fractions
    cluster_purity = {}
    for cluster, cluster_cell_types in cluster_break_up.items():
        total_count = sum(cluster_cell_types.values())
        for cell_type, count in cluster_cell_types.items():
            cluster_purity.setdefault(cluster, {})[cell_type] = count / total_count
    cluster_purity = {"cas_cluster_purity": dict(sorted(cluster_purity.items()))}

    return cluster_purity


def get_df_from_analysis(
    analysis: str,
    get_umap_coords: bool = False,
) -> pd.DataFrame:
    """Get a dataframe of stuff from the analysis.h5.

    Returns barcode and clusters from graph-clust if get_umap_coords is False.
    Also includes UMAP coords if get_umap_coords is True.
    """
    analysis_h5_path = os.path.join(analysis, "analysis.h5")
    analysis_h5 = h5py.File(analysis_h5_path, "r")
    cluster_labels = analysis_h5["/clustering/_gene_expression_graphclust/clusters"]
    barcodes = (x.decode("utf-8") for x in analysis_h5["matrix/barcodes"])
    dict_of_data = {"barcode": barcodes, "clusters": cluster_labels}

    if get_umap_coords:
        umap_results = np.array(
            analysis_h5[cr_analysis_constants.ANALYSIS_H5_UMAP_GROUP]["_gene_expression_2"][
                "transformed_umap_matrix"
            ]
        )
        umap_x, umap_y = zip(*umap_results)
        dict_of_data["umap_x"] = umap_x
        dict_of_data["umap_y"] = umap_y

    cluster_df = pd.DataFrame(dict_of_data)
    return cluster_df


def get_df_from_cloupe_projection(cloupe_projection: str) -> pd.DataFrame:
    """Get a dataframe of the cloupe projection and rename columns for the umap plot."""
    cloupe_projection_df = pd.read_csv(cloupe_projection)
    cloupe_projection_df.rename(
        columns={"Barcode": "barcode", "X Coordinate": "umap_x", "Y Coordinate": "umap_y"},
        inplace=True,
    )
    return cloupe_projection_df


def get_cell_types_df(cell_types: str, analysis=None):
    """Process cell type information and optionally merge it with analysis data.

    Args:
        cell_types (str): Path to the CSV file containing cell type information.
        analysis (optional): Path to the analysis directory.
            If provided, the cell type information is merged based on the 'barcode' column.

    Returns:
        cell_types (pd.DataFrame): Processed DataFrame containing cell type information.
            If 'analysis' is provided, cell type information is merged and wrapped cell type
            names are included in the 'cell_type_wrapped' column.
    """
    cell_types = pd.read_csv(cell_types)
    if analysis:
        cluster_df = get_df_from_analysis(analysis)
        cell_types = pd.merge(cluster_df, cell_types, on="barcode", how="left")
    cell_types = cell_types[[col for col in cell_types.columns if not col.endswith("distance")]]
    cell_types = cell_types.dropna()
    cell_types["cell_type_wrapped"] = cell_types["coarse_cell_type"].apply(
        lambda x: "|".join(textwrap.wrap(x, width=15))
    )
    return cell_types


def cell_type_bar_chart(
    cell_types: str,
    analysis: str,
    return_chart: bool = False,
    stacked: bool = True,
) -> dict:
    """Generate a stacked bar chart representing the cluster purity based on cell annotation cell types and clusters.

    Args:
        cell_types (str): Path to the cell annotation cell types CSV file.
        analysis (str): Path to the analysis directory.
        return_chart(bool): Return the chart for plotting in a notebook or a dict for plotting in the web summary
        stacked(bool): Make a stacked bar chart per cluster or a bar chart of cells types as percent of total cells

    Returns:
        dict: A JSON representation of the stacked bar chart created using Altair.
    """
    if stacked:
        cell_types = get_cas_cluster_purity(cell_types=cell_types, analysis=analysis)
        cell_types = pd.DataFrame(cell_types["cas_cluster_purity"]).fillna(0)
        cell_types = cell_types.reset_index(names=COARSE_CELL_TYPES_KEY)
        cell_types = pd.melt(cell_types, id_vars=COARSE_CELL_TYPES_KEY, value_name="Fraction")
        cell_types["cluster"] = cell_types["variable"]
        cell_types["cell_type_wrapped"] = cell_types[COARSE_CELL_TYPES_KEY].apply(
            lambda x: "|".join(textwrap.wrap(x, width=15))
        )
        chart = (
            alt.Chart(cell_types)
            .mark_bar()
            .encode(
                x=alt.X("cluster:N", axis=alt.Axis(title="Cluster", labelAngle=0)),
                y=alt.Y("sum(Fraction):Q", axis=alt.Axis(title="Percent Cell Type", format="%")),
                color=alt.Color(
                    "cell_type_wrapped:N",
                    legend=alt.Legend(
                        title="Cell Type",
                        labelExpr="split(datum.label, '|')",
                        orient="bottom",
                        columns=4,
                    ),
                ),
                tooltip=[
                    alt.Tooltip(COARSE_CELL_TYPES_KEY, title="Cell Type"),
                    alt.Tooltip("Fraction", title="Percent Cell Type", format=".0%"),
                ],
            )
            .properties(width=320, height=300)
        )
    else:
        cluster_df = get_df_from_analysis(analysis)
        cell_types = pd.read_csv(cell_types)
        cell_types = pd.merge(cluster_df, cell_types, on="barcode", how="left")
        cell_types = cell_types[[col for col in cell_types.columns if not col.endswith("distance")]]
        cell_types = cell_types.dropna()
        cell_types["cell_type_wrapped"] = cell_types[COARSE_CELL_TYPES_KEY].apply(
            lambda x: "|".join(textwrap.wrap(x, width=15))
        )

        total_count = len(cell_types)
        chart = (
            alt.Chart(cell_types)
            .transform_aggregate(
                count="count()", groupby=[COARSE_CELL_TYPES_KEY, "cell_type_wrapped"]
            )
            .transform_calculate(fraction=f"datum.count / {total_count}")
            .mark_bar()
            .encode(
                x=alt.X(
                    "fraction:Q",
                    axis=alt.Axis(title="Percentage", format="%"),
                    scale=alt.Scale(domain=(0, 1)),
                ),
                y=alt.Y("cell_type:N", sort="-x", axis=alt.Axis(title="Cell Type")),
                color=alt.Color(
                    "cell_type_wrapped:N",
                    legend=alt.Legend(
                        title="Cell Type",
                        labelExpr="split(datum.label, '|')",
                        orient="bottom",
                        columns=8,
                    ),
                ),
                tooltip=[
                    alt.Tooltip("cell_type:N", title="Cell Type"),
                    alt.Tooltip("fraction:Q", title="Percentage", format=".2%"),
                ],
            )
        )

    if return_chart:
        return chart
    chart = alt_utils.chart_to_json(chart)
    return chart


def cell_type_umap(
    cell_types: str, analysis: str, cloupe_projection: str, return_chart: bool = False
) -> dict:
    """Generate a UMAP scatter plot representing the UMAP projection colored by cell annotation cell types and labeled with cluster.

    Args:
        cell_types (str): Path to the cell annotation cell types CSV file.
        analysis (str): Path to the analysis directory.
        cloupe_projection (str): Path to the cloupe projection CSV file.
        return_chart(bool): Return the chart for plotting in a notebook or a dict for plotting in the web summary

    Returns:
        dict: A JSON representation of the UMAP scatter plot colored with cell type created using Altair.
    """
    if analysis:
        cluster_df = get_df_from_analysis(analysis, get_umap_coords=True)
    if cloupe_projection:
        cluster_df = get_df_from_cloupe_projection(cloupe_projection)
    cell_types = get_cell_types_df(cell_types)
    cell_types = pd.merge(cluster_df, cell_types, on="barcode", how="left")
    cell_types = cell_types.dropna()
    if len(cell_types) > gex_constants.MAX_WEBSHIM_BCS_DIM:
        cell_types = cell_types.sample(n=gex_constants.MAX_WEBSHIM_BCS_DIM, random_state=1)
    # interpolate the umap spot size
    num_spots = [1000, 5000]
    spot_sizes = [50, 10]
    num_barcodes = len(cell_types)
    spot_size = np.interp(num_barcodes, num_spots, spot_sizes)

    highlight = alt.selection_multi(fields=["cell_type_wrapped"], bind="legend")

    common_columns = [
        "umap_x",
        "umap_y",
        "barcode",
        "coarse_cell_type",
        "cell_type_wrapped",
        "fine_cell_type",
    ]
    tool_tip = [
        alt.Tooltip("coarse_cell_type", title="Cell Type"),
        alt.Tooltip("fine_cell_type", title="Sub-type"),
    ]
    if cloupe_projection:
        cell_types = cell_types[common_columns]
    else:
        cell_types = cell_types[common_columns + ["clusters"]]
        tool_tip.append(alt.Tooltip("clusters", title="Cluster"))
    # make the chart
    chart = (
        alt.Chart(cell_types)
        .mark_circle(size=spot_size)
        .encode(
            x=alt.X("umap_x", axis=alt.Axis(title="UMAP X")),
            y=alt.Y("umap_y", axis=alt.Axis(title="UMAP Y")),
            color=alt.condition(
                highlight,
                "cell_type_wrapped:N",
                alt.value("lightgray"),
                legend=alt.Legend(
                    title="Cell Type",
                    labelExpr="split(datum.label, '|')",
                    orient="bottom",
                    columns=4,
                ),
            ),
            opacity=alt.condition(highlight, alt.value(1), alt.value(0.2)),
            tooltip=tool_tip,
        )
        .add_selection(highlight)
        .configure_legend(titleFontSize=15, labelFontSize=12)
        .interactive()
        .properties(width=450)
    )
    if return_chart:
        return chart
    chart = alt_utils.chart_to_json(chart)
    return chart


def spatial_cell_types_plot(
    tissue_positions_path: str,
    tissue_lowres_image: str,
    cell_types: str,
    analysis: str,
    scale_factors: str,
    return_chart: bool = False,
) -> dict | alt.Chart:
    """Generate a spatial plot of cell annotation cell types over tissue.

    Args:
        tissue_positions_path (str): Path to the tissue positions file.
        tissue_lowres_image (str): Path to the low-resolution tissue image.
        cell_types (str): Cell types data for the tissue.
        analysis (str): Path to analysis information for the cell types.
        scale_factors(str): Path to the scalefactors_json.json.
        return_chart (bool, optional): If True, return the Altair chart object. If False (default),
            return the plot as a JSON representation.

    Returns:
        dict or alt.Chart: The generated spatial plot as a dictionary or an Altair chart object.
            If return_chart is True, returns the chart object; otherwise, returns the JSON representation of the chart.
    """
    # Make the tissue positions with cell type calls
    cell_types = get_cell_types_df(cell_types, analysis)

    with open(scale_factors) as f:
        scale_factor_dict = json.load(f)

    def transfer_to_lowres(x: float) -> float:
        return float(x) * scale_factor_dict["tissue_lowres_scalef"]

    tissue_positions_df = pd.read_csv(
        tissue_positions_path,
        usecols=["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"],
        converters={
            "pxl_row_in_fullres": transfer_to_lowres,
            "pxl_col_in_fullres": transfer_to_lowres,
        },
    ).rename(
        columns={
            "pxl_row_in_fullres": "pixl_row_in_lowres",
            "pxl_col_in_fullres": "pixl_col_in_lowres",
        }
    )
    tissue_positions_df = pd.merge(tissue_positions_df, cell_types, on="barcode", how="inner")
    tissue_positions_df = tissue_positions_df.dropna()

    final_plot = make_alt_spatial_plot(
        tissue_lowres_image=tissue_lowres_image,
        plotting_df=tissue_positions_df,
        col_to_plot="cell_type_wrapped",
        col_to_tooltip=COARSE_CELL_TYPES_KEY,
        tooltip_title="Cell Type",
        add_color_selector=True,
        data_is_continuous=False,
        return_chart=return_chart,
    )
    return final_plot
