#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Code to produce altair umi and Genes distribution plots for the websummary."""
import altair as alt
import numpy as np
import pandas as pd

import cellranger.altair_utils as alt_utils
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
from cellranger.webshim.data import SampleData
from cellranger.websummary.numeric_converters import round_floats_in_list

alt.data_transformers.disable_max_rows()


def make_violin_plot_help(is_spatial: bool):
    """Produce the help text seen in the web summary.

    Args:
        is_spatial (bool): is the sample spatial?

    Returns:
        dict: dictionary of help text for the plot
    """
    barcode_key = "Spots" if is_spatial else "Cells"
    violin_plot_help = {
        "helpText": f"Distributions of genes and UMIs. {barcode_key} marked as outliers shown as grey circles. "
        "Hover over the boxplot to see quartile values.",
        "title": "Gene and UMI Distribution",
    }
    return violin_plot_help


def make_violin_plots(sample_data: SampleData, library_type: str, is_spatial: bool):
    """Makes a plot with distributions from UMIs/Genes per barcode.

    Args:
        sample_data (SampleData): SampleData class object
        library_type (str): library type from from cellranger.rna.library
        is_spatial (bool): is the sample spatial?

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
        {f"UMIs per {barcode_key}": umi_per_bc, f"Genes per {barcode_key}": genes_per_bc}
    )
    plot_data = pd.melt(
        plot_data,
        value_vars=[f"UMIs per {barcode_key}", f"Genes per {barcode_key}"],
        value_name="count",
    )
    plot_data["count"] = round_floats_in_list(np.log10(plot_data["count"] + 1))
    # plot_data["count"] = round_floats_in_list(plot_data["count"])

    # Make the plot
    boxplot = alt.Chart().mark_boxplot(color="grey").encode(y=alt.Y("count:Q", title=""))

    violin = (
        alt.Chart()
        .transform_density(
            "count",
            as_=["count", "density"],
            groupby=["variable"],
            extent=[0, max(plot_data["count"]) + 0.5],
        )
        .mark_area(orient="horizontal", opacity=0.4)
        .encode(
            y=alt.Y("count:Q", title="Log10(1+Count)"),
            color=alt.Color(
                "variable:N",
                legend=None,
                scale=alt.Scale(domain=plot_data.variable.unique(), range=["black", "red"]),
            ),
            x=alt.X(
                "density:Q",
                stack="center",
                impute=None,
                title=None,
                scale=alt.Scale(nice=False, zero=True, padding=100),
                axis=alt.Axis(labels=False, values=[0], grid=False, ticks=False),
            ),
        )
    )

    violin_box_plot = (
        alt.layer(violin, boxplot, data=plot_data)
        .properties(width=450)
        .facet(
            column=alt.Column(
                "variable:N",
                header=alt.Header(
                    title=None,
                    labelOrient="bottom",
                    labelFontSize=15,
                ),
            ),
        )
        .configure_facet(spacing=0)
        .configure_view(stroke=None)
        .configure_axis(titleFontSize=15)
        .resolve_scale(x=alt.ResolveMode("independent"))
    )

    violin_box_plot = alt_utils.chart_to_json(violin_box_plot)
    violin_plot_help = make_violin_plot_help(is_spatial)
    violin_plot_data = {"violin_plots": {"help": violin_plot_help, "spec": violin_box_plot}}
    return violin_plot_data
