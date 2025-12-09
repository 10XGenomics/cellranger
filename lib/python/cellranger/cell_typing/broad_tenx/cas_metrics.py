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
from cellranger.cell_typing.common_cell_typing import (
    BARCODE_KEY,
    COARSE_CELL_TYPES_KEY,
    FINE_CELL_TYPES_KEY,
)

COMMON_COLUMNS = [
    "umap_x",
    "umap_y",
    "barcode",
    "coarse_cell_type",
    "cell_type_wrapped",
    "fine_cell_type",
]
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
    barcodes = (x.decode("utf-8") for x in analysis_h5["matrix/barcodes"][:])
    clusters = analysis_h5["/clustering/_gene_expression_graphclust/clusters"][:]
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
    cluster_labels = analysis_h5["/clustering/_gene_expression_graphclust/clusters"][:]
    barcodes = (x.decode("utf-8") for x in analysis_h5["matrix/barcodes"][:])
    dict_of_data = {"barcode": barcodes, "clusters": cluster_labels}

    if get_umap_coords:
        umap_results = np.array(
            analysis_h5[cr_analysis_constants.ANALYSIS_H5_UMAP_GROUP]["_gene_expression_2"][
                "transformed_umap_matrix"
            ][:]
        )
        umap_x, umap_y = zip(*umap_results)
        dict_of_data["umap_x"] = umap_x
        dict_of_data["umap_y"] = umap_y

    cluster_df = pd.DataFrame(dict_of_data)
    return cluster_df


def sample_df(segmentation_umap_df, group_by_key: str, max_samples=20000):
    """Samples rows from a DataFrame to create a subset with a maximum number of rows.

    Args:
        segmentation_umap_df (pd.DataFrame): The input DataFrame containing data to be sampled.
        group_by_key (str): The column name used to group the data for mandatory sampling.
        max_samples (int, optional): The maximum number of rows to include in the sampled subset.
                                      Defaults to 20,000.

    Returns:
        pd.DataFrame: A DataFrame containing the sampled subset of rows. If the input DataFrame
                      has fewer rows than `max_samples`, the original DataFrame is returned.
    """
    if len(segmentation_umap_df) > max_samples:
        # Select one random row from each cluster to ensure each cluster is represented
        mandatory = segmentation_umap_df.groupby(f"{group_by_key}", group_keys=False).apply(
            lambda x: x.sample(1)
        )

        # Calculate how many additional rows are needed
        remaining_count = max_samples - len(mandatory)

        # Drop the already selected rows and sample the remaining needed rows
        remaining_df = segmentation_umap_df.drop(mandatory.index)
        additional = remaining_df.sample(remaining_count)

        # Concatenate the mandatory and additional rows and shuffle the result
        return pd.concat([mandatory, additional]).sample(frac=1).reset_index(drop=True)
    return segmentation_umap_df


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


def prepare_cell_types_data(cell_types: str, analysis: str, cloupe_projection: str) -> pd.DataFrame:
    """Prepare and process the cell types DataFrame for plotting."""
    if analysis:
        cluster_df = get_df_from_analysis(analysis, get_umap_coords=True)
    if cloupe_projection:
        cluster_df = get_df_from_cloupe_projection(cloupe_projection)
    cell_types = get_cell_types_df(cell_types)
    cell_types = pd.merge(cluster_df, cell_types, on="barcode", how="left")
    cell_types = cell_types.dropna()
    if len(cell_types) > gex_constants.MAX_WEBSHIM_BCS_DIM:
        cell_types = cell_types.sample(n=gex_constants.MAX_WEBSHIM_BCS_DIM, random_state=1)
    return cell_types


def interpolate_spot_size(num_barcodes: int) -> float:
    """Interpolate the UMAP spot size based on the number of barcodes."""
    num_spots = [1000, 5000]
    spot_sizes = [50, 10]
    return np.interp(num_barcodes, num_spots, spot_sizes)


def base_umap_plot(  # pylint: disable=too-many-arguments, too-many-positional-arguments
    cell_types: pd.DataFrame,
    spot_size: float,
    tool_tip: list,
    color_encoding: alt.Color,
    return_chart: bool,
    highlight=None,
    additional_params=None,
    configure_legend=False,
    chart_width=450,
) -> dict:
    """Create the UMAP chart using Altair."""
    chart = (
        alt.Chart(cell_types)
        .mark_circle(size=spot_size)
        .encode(
            x=alt.X("umap_x:Q", axis=alt.Axis(title="UMAP X")),
            y=alt.Y("umap_y:Q", axis=alt.Axis(title="UMAP Y")),
            color=color_encoding,
            tooltip=tool_tip,
        )
        .add_params(highlight, *(additional_params or []))
        .configure_legend(
            **(
                {"titleFontSize": 15, "labelFontSize": 12, "orient": "bottom", "columns": 4}
                if configure_legend
                else {}
            )
        )
        .interactive()
        .properties(width=chart_width)
    )
    if return_chart:
        return chart  # type: ignore
    return alt_utils.chart_to_json(chart)


def cell_type_umap(
    cell_types: str,
    analysis: str,
    cloupe_projection: str,
    configure_legend: bool,
    return_chart: bool = False,
) -> dict:
    """Generate a UMAP scatter plot colored by cell annotation cell types."""
    cell_types = prepare_cell_types_data(cell_types, analysis, cloupe_projection)
    col_names = cell_types.columns
    has_summary_score = "summary_score" in col_names
    spot_size = interpolate_spot_size(len(cell_types))
    highlight = alt.selection_point(fields=["cell_type_wrapped"], bind="legend")
    tool_tip = [
        alt.Tooltip("coarse_cell_type", title="Cell Type"),
        alt.Tooltip("fine_cell_type", title="Sub-type"),
    ]
    conditional_columns = []  # Additional columns may be added

    if not cloupe_projection:
        conditional_columns.append("clusters")
        tool_tip.append(alt.Tooltip("clusters", title="Cluster"))

    # Account for presence of summary score
    if has_summary_score:
        conditional_columns.append("summary_score")
        tool_tip.append(
            alt.Tooltip("summary_score", title="Summary Score"),
        )
        lower_slider = alt.binding_range(
            min=0, max=1.01, step=0.01, name="Lower Quality Score Cutoff "
        )
        upper_slider = alt.binding_range(
            min=0, max=1, step=0.01, name="Upper Quality Score Cutoff "
        )
        lower_cutoff = alt.param(name="LowerCutoff", value=0, bind=lower_slider)
        upper_cutoff = alt.param(name="UpperCutoff", value=1, bind=upper_slider)
        pd_predicate = (alt.datum.summary_score >= lower_cutoff) & (
            alt.datum.summary_score <= upper_cutoff
        )
        combined_predicate = highlight & pd_predicate
    color_encoding = alt.condition(
        combined_predicate if has_summary_score else highlight,
        alt.Color(
            "cell_type_wrapped:N",
            scale=alt.Scale(scheme="tableau20"),
            title="Cell Type",
            legend=alt.Legend(
                title="Cell Type",
                labelExpr="split(datum.label, '|')",
            ),
        ),
        alt.value("lightgray"),
    )
    cell_types = cell_types[COMMON_COLUMNS + conditional_columns]
    return base_umap_plot(
        cell_types,
        spot_size,
        tool_tip,
        color_encoding,
        return_chart,
        highlight,
        additional_params=[lower_cutoff, upper_cutoff] if has_summary_score else None,
        configure_legend=configure_legend,
    )


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


def get_cell_type_summary_mets(cell_types: str):
    """Get a summary of coarse and fine cell types. Number of fine and coarse categories.

    Args:
        cell_types (str): Path to the cell annotation cell types CSV file.

    Returns:
        dict: Summary of cell types. num_coarse_cell_categories, num_fine_cell_categories
    """
    coarse_cell_types = set()
    fine_cell_types = set()

    with open(cell_types) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            coarse_cell_types.add(row[COARSE_CELL_TYPES_KEY])
            fine_cell_types.add(row[FINE_CELL_TYPES_KEY])

    num_coarse_cell_categories = len(coarse_cell_types)
    num_fine_cell_categories = len(fine_cell_types)
    return {
        "num_coarse_cell_categories": num_coarse_cell_categories,
        "num_fine_cell_categories": num_fine_cell_categories,
    }
