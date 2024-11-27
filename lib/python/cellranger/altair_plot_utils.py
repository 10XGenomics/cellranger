#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Convenience functions for oft used altair plots."""

import altair as alt
import pandas as pd

from cellranger.altair_utils import chart_to_json
from cellranger.spatial.image import WebImage


# pylint: disable=too-many-locals
def make_alt_spatial_plot(
    tissue_lowres_image: str,
    plotting_df: pd.DataFrame,
    col_to_plot: str,
    col_to_tooltip: str,
    tooltip_title: str,
    add_color_selector: bool = False,
    data_is_continuous: bool = True,
    return_chart: bool = False,
) -> alt.Chart | str:
    """Generate an Altair spatial plot displaying data from a given DataFrame overlaid on a web image.

    Args:
        tissue_lowres_image (str): Path to the low-resolution tissue image.
        plotting_df (pd.DataFrame): DataFrame containing data to be plotted.
            Typically a tissue_positions.csv joined with other per-barcode info
        col_to_plot (str): Column name in plotting_df containing data for coloring the plot.
        col_to_tooltip (str): Column name in plotting_df containing data for tooltip information.
        tooltip_title (str): Title for the tooltip.
        add_color_selector(bool, optional): Should the plot have a color scheme dropdown?
        data_is_continuous(bool, optional): Is the data being plotted continuous? Necessary for selecting the color scheme if using the dropdown
        return_chart (bool, optional): If True, return the Altair chart object. If False (default),
            return the plot as a JSON representation.

    Returns:
        final_plot (alt.Chart or str): Altair chart object if return_chart is True, or JSON
        representation of the plot if return_chart is False.
    """
    img = WebImage(tissue_lowres_image)
    encoded_image = img.base64_encoded_str
    source = pd.DataFrame.from_records([{"x": img.width, "y": img.height, "img": encoded_image}])
    imgs = alt.Chart(source).mark_image().encode(url="img")

    slider_selector = alt.param(
        value=1, bind=alt.binding_range(min=0, max=1, step=0.1, name="Spot Opacity: ")
    )
    # Hover tooltip
    nearest = alt.selection_point(
        fields=[f"{col_to_tooltip}"],
        nearest=True,
        on="mouseover",
        empty="none",
    )
    # add cluster tooltip if it exists
    cluster_tooltip = (
        alt.Tooltip("clusters", title="Cluster")
        if plotting_df.columns.isin(["clusters"]).any()
        else None
    )
    plot_tooltip = [
        alt.Tooltip("barcode", title="Barcode"),
        alt.Tooltip(f"{col_to_tooltip}", title=f"{tooltip_title}"),
    ]
    # append the cluster tooltip if it exists
    if cluster_tooltip is not None:
        plot_tooltip.append(cluster_tooltip)
    array_plot = (
        alt.Chart(plotting_df)
        .encode(
            x=alt.X("pixl_col_in_lowres:Q", title="").scale(domain=(0, source["x"][0])),
            y=alt.Y("pixl_row_in_lowres:Q", scale=alt.Scale(reverse=True), title="").scale(
                domain=(source["y"][0], 0)
            ),
            color=alt.Color(
                f"{col_to_plot}:N",
                legend=alt.Legend(
                    title=f"{tooltip_title}", labelExpr="split(datum.label, '|')", symbolOpacity=1
                ),
            ),
            tooltip=plot_tooltip,
        )
        .mark_circle(size=10 if len(plotting_df) < 5000 else 5, opacity=slider_selector)
        .properties(
            width=img.width / 1.25,
            height=img.height / 1.25,
        )
        .add_selection(nearest)
        .add_selection(slider_selector)
    )

    final_plot = (
        alt.layer(imgs, array_plot)
        .configure_axis(grid=False, domain=False, tickSize=0, labels=False)
        .configure_view(strokeWidth=0)
    )

    if add_color_selector:
        continuous_colors = [
            "turbo",
            "cividis",
            "bluegreen",
            "yelloworangered",
            "viridis",
            "magma",
            "sinebow",
        ]
        discrete_colors = [
            "tableau20",
            "accent",
            "dark2",
            "set1",
            "category20",
            "category20b",
            "category20c",
        ]
        input_dropdown = alt.binding_select(
            options=continuous_colors if data_is_continuous else discrete_colors,
            name="Color Scheme: ",
        )
        color_select = alt.param(
            name="schemeselect",
            bind=input_dropdown,
            value="turbo" if data_is_continuous else "tableau20",
        )
        final_plot = final_plot.encode(
            color=alt.Color(scale=alt.Scale(scheme={"expr": "schemeselect"}))
        ).add_selection(color_select)
    if return_chart:
        return final_plot
    # only write out the plot to json
    final_plot = chart_to_json(final_plot)
    return final_plot
