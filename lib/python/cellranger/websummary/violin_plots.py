# !/usr/bin/env python

# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

"""Code to produce altair umi and Genes distribution plots for the websummary."""
import altair as alt
import pandas as pd
import polars as pl

import cellranger.altair_utils as alt_utils
import cellranger.matrix as cr_matrix
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
from cellranger.webshim.data import SampleData
from cellranger.websummary.numeric_converters import round_floats_in_list

alt.data_transformers.disable_max_rows()


def _make_boxplot(counts: str):
    """Helper function to make a base boxplot.

    Args:
        counts(string): the name of the colum with the counts for plotting
    """
    boxplot = alt.Chart().mark_boxplot(color="grey").encode(y=alt.Y(f"{counts}:Q", title=""))
    return boxplot


def _make_violin_plot(
    data, counts: str, group_by_name: str, y_label: str, set_scale: bool, logy: bool = False
):
    """Helper function to make a base violin plot.

    Args:
        data(DataFrame): data to use for the plot. can be pandas or polars
        counts(string): the name of the colum with the counts for plotting
        group_by_name (string): name of the column that data should be grouped by
        y_label(string): name of the y-axis label
        set_scale(bool): set a custom color scale. Useful for spatial ws plots
        logy(bool): log the y-axis
    """
    # Conditional color scaling
    if set_scale:
        color_scale = alt.Scale(domain=data.variable.unique(), range=["black", "red"])
        x_scale = alt.Scale(nice=False, zero=True, padding=100)
    else:
        color_scale = alt.Scale()
        x_scale = alt.Scale(nice=False, zero=False)

    y_axis = (
        alt.Y(f"{counts}:Q", title=y_label).scale(type="log")
        if logy
        else alt.Y(f"{counts}:Q", title=y_label)
    )

    violin = (
        alt.Chart()
        .transform_density(
            f"{counts}",
            as_=[f"{counts}", "density"],
            groupby=[f"{group_by_name}"],
            extent=[1, max(data[f"{counts}"]) + 0.5] if logy else [0, max(data[f"{counts}"]) + 0.5],
        )
        .mark_area(orient="horizontal", opacity=0.4)
        .encode(
            y=y_axis,
            color=alt.Color(
                f"{group_by_name}:N",
                legend=None,
                scale=color_scale,
            ),
            x=alt.X(
                "density:Q",
                stack="center",
                impute=None,
                title=None,
                scale=x_scale,
                axis=alt.Axis(labels=False, values=[0], grid=False, ticks=False),
            ),
        )
    )
    return violin


def _combine_box_violin_plots(
    boxplot, violin, data, group_by_name, plot_width, rotate_labels=False
):
    """Combine violin and box plots.

    Args:
        boxplot = altar base boxplot
        violin = altair base violin plot
        data = plotting data
        group_by_name (string): name of the column that data should be grouped by
        rotate_labels = rotate the x axis labels
    """
    violin_box_plot = (
        alt.layer(violin, boxplot, data=data)
        .properties(width=plot_width)
        .facet(
            column=alt.Column(
                f"{group_by_name}:N",
                header=(
                    alt.Header(
                        title=None,
                        labelOrient="bottom",
                        labelFontSize=15,
                        labelAngle=-45,
                        labelAlign="right",
                        labelPadding=20,
                    )
                    if rotate_labels
                    else alt.Header(
                        title=None,
                        labelOrient="bottom",
                        labelFontSize=15,
                    )
                ),
            ),
        )
        .configure_facet(spacing=0)
        .resolve_scale(x=alt.ResolveMode("independent"))
        .configure_view(stroke=None)
        .configure_axis(titleFontSize=15)
    )
    return violin_box_plot


def make_violin_plot_help(is_spatial: bool):
    """Produce the help text seen in the web summary.

    Args:
        is_spatial (bool): is the sample spatial?

    Returns:
        dict: dictionary of help text for the plot
    """
    barcode_key = "Spots" if is_spatial else "Cells"
    violin_plot_help = {
        "helpText": f"Distributions of genes and UMIs + 1. {barcode_key} marked as outliers shown as grey circles. "
        "Y-axis is log scale. "
        "Hover over the boxplot to see quartile values.",
        "title": "Gene and UMI Distribution",
    }
    return violin_plot_help


def make_gene_umi_violin_plots(
    sample_data: SampleData,
    library_type: str,
    is_spatial: bool,
    group_by_name: str = "variable",
):
    """Makes a plot with distributions from UMIs/Genes per barcode.

    Args:
        sample_data (SampleData): SampleData class object
        library_type (str): library type from from cellranger.rna.library
        is_spatial (bool): is the sample spatial?
        group_by_name (string): name of the column that data should be grouped by

    Returns:
        Dict: violin plot data and help for plotting in the websummary
    """
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if not isinstance(analysis, SingleGenomeAnalysis):
        return {}

    # Define the barcode
    barcode_key = "Spot" if is_spatial else "Cell"
    # Get the sum UMI for each barcode
    matrix = analysis.matrix.select_features_by_type(library_type)
    umi_per_bc = matrix.get_counts_per_bc()
    genes_per_bc = matrix.get_numfeatures_per_bc()

    # Tidy the data into a pandas data.frame
    plot_data = pd.DataFrame.from_dict(
        {
            f"UMIs per {barcode_key}": umi_per_bc,
            f"Genes per {barcode_key}": genes_per_bc,
        }
    )
    plot_data = pd.melt(
        plot_data,
        value_vars=[f"UMIs per {barcode_key}", f"Genes per {barcode_key}"],
        value_name="count",
    )
    plot_data["count"] = round_floats_in_list(plot_data["count"] + 1)
    # Make the box plot
    boxplot = _make_boxplot(counts="count")
    # Make violin plot
    violin = _make_violin_plot(
        data=plot_data,
        counts="count",
        group_by_name=group_by_name,
        y_label="(1+Count)",
        set_scale=True,
        logy=True,
    )
    # Combine the plots
    violin_box_plot = _combine_box_violin_plots(
        boxplot,
        violin,
        data=plot_data,
        group_by_name=group_by_name,
        plot_width=440,
        rotate_labels=False,
    )

    violin_box_plot = alt_utils.chart_to_json(violin_box_plot)
    violin_plot_help = make_violin_plot_help(is_spatial)
    violin_plot_data = {"violin_plots": {"help": violin_plot_help, "spec": violin_box_plot}}
    return violin_plot_data


def make_cell_types_df(
    h5_path: str,
    cell_types: str,
    barcode_lower_bound: int,
    group_by_name: str,
):
    """Makes a plot with distributions from UMIs cell types.

    Args:
        h5_path (str): path to filtered_feature_bc_matrix.h5
        cell_types (str): path to cell_types.csv
        barcode_lower_bound(int): minimum number of barcodes per cell type for plot
        group_by_name (string): name of the column that data should be grouped by

    Returns:
        tupel: [cell_types (DataFrame), cell_levels(int)]
    """
    mat = cr_matrix.CountMatrix.load_h5_file(h5_path)
    mat = mat.select_features_by_type("Gene Expression")
    counts = mat.get_counts_per_bc()
    bcs = cr_matrix.CountMatrix.load_bcs_from_h5_file_handle(h5_path)
    bcs = [b.decode("utf-8") for b in bcs]
    counts_per_bcs = pl.DataFrame({"barcode": bcs, "umi_counts": counts})
    cell_types = pl.read_csv(cell_types)
    cell_types = cell_types.join(counts_per_bcs, on="barcode")
    cell_types = cell_types.select([group_by_name, "umi_counts"])
    cell_types = cell_types.join(
        cell_types.group_by(group_by_name)
        .agg(pl.count().alias("count"))
        .filter(pl.col("count") >= barcode_lower_bound),
        on=group_by_name,
        how="inner",
    )
    cell_types = cell_types.with_columns(pl.col("umi_counts") + 1)
    cell_levels = cell_types[group_by_name].n_unique()
    return cell_types, cell_levels


def make_cell_types_violin_plot(
    h5_path: str,
    cell_types: str,
    barcode_lower_bound: int,
    final_plot_width: int,
    group_by_name: str = "coarse_cell_type",
    return_plot: bool = False,
):
    """Makes a plot with distributions from UMIs cell types.

    Args:
        h5_path (str): path to filtered_feature_bc_matrix.h5
        cell_types (str): path to cell_types.csv
        barcode_lower_bound(int): minimum number of barcodes per cell type for plot
        final_plot_width (int): plot width that should be displayed in the WS
        group_by_name (string): name of the column that data should be grouped by
        return_plot (bool): should the plot be returned (for notebooks) instead of the json for the WS

    Returns:
        Dict: violin plot data and help for plotting in the websummary
    """
    cell_types, cell_levels = make_cell_types_df(
        h5_path, cell_types, barcode_lower_bound, group_by_name
    )
    cell_types = cell_types.with_columns(pl.col("umi_counts") + 1)
    cell_levels = cell_types[f"{group_by_name}"].n_unique()

    # Make the box plot
    boxplot = _make_boxplot(counts="umi_counts")
    # Make violin plot
    violin = _make_violin_plot(
        data=cell_types,
        counts="umi_counts",
        group_by_name=group_by_name,
        y_label="1+UMI",
        set_scale=False,
        logy=True,
    )
    # Combine the plots
    violin_box_plot = _combine_box_violin_plots(
        boxplot=boxplot,
        violin=violin,
        data=cell_types,
        group_by_name=group_by_name,
        plot_width=final_plot_width / cell_levels,
        rotate_labels=True,
    )

    if return_plot:
        return violin_box_plot
    else:
        violin_box_plot = alt_utils.chart_to_json(violin_box_plot)
        return violin_box_plot


def make_cell_types_boxplot(
    h5_path: str,
    cell_types: str,
    barcode_lower_bound: int,
    final_plot_width: int,
    group_by_name: str = "coarse_cell_type",
    return_plot: bool = False,
):
    """Makes a plot with distributions from UMIs cell types.

    Args:
        h5_path (str): path to filtered_feature_bc_matrix.h5
        cell_types (str): path to cell_types.csv
        barcode_lower_bound(int): minimum number of barcodes per cell type for plot
        final_plot_width (int): plot width that should be displayed in the WS
        group_by_name (string): name of the column that data should be grouped by
        return_plot (bool): should the plot be returned (for notebooks) instead of the json for the WS

    Returns:
        Dict: violin plot data and help for plotting in the websummary
    """
    cell_types, cell_levels = make_cell_types_df(
        h5_path, cell_types, barcode_lower_bound, group_by_name
    )

    box_width = final_plot_width / (cell_levels + 1)
    box_plot = alt.Chart(cell_types).mark_boxplot(size=box_width).encode(
        alt.X(
            f"{group_by_name}:N",
            axis=alt.Axis(labelAngle=-45, labelFontSize=15, title=None),
        ),
        alt.Y("umi_counts:Q", title="1+UMI").scale(type="log"),
        alt.Color(f"{group_by_name}:N").legend(None),
        tooltip=[
            alt.Tooltip("umi_counts:Q", title="UMI", format=",.2f"),
        ],
    ) + alt.Chart(cell_types).transform_aggregate(
        min="min(umi_counts)",
        max="max(umi_counts)",
        mean="mean(umi_counts)",
        median="median(umi_counts)",
        q1="q1(umi_counts)",
        q3="q3(umi_counts)",
        count="max(count)",
        groupby=[f"{group_by_name}"],
    ).mark_bar(
        opacity=0
    ).encode(
        x=f"{group_by_name}:N",
        y="q1:Q",
        y2="q3:Q",
        tooltip=[
            alt.Tooltip("min:Q", title="Minimum UMI", format=".2f"),
            alt.Tooltip("q1:Q", title="Lower Quartile", format=".2f"),
            alt.Tooltip("mean:Q", title="Mean UMI", format=".2f"),
            alt.Tooltip("median:Q", title="Median UMI", format=".2f"),
            alt.Tooltip("q3:Q", title="Upper Quartile", format=".2f"),
            alt.Tooltip("max:Q", title="Maximum UMI", format=".2f"),
            alt.Tooltip("count:Q", title="# Barcodes"),
        ],
    ).properties(
        width=final_plot_width
    )
    if return_plot:
        return box_plot
    else:
        box_plot = alt_utils.chart_to_json(box_plot)
        return box_plot
