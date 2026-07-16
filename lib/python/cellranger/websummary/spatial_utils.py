#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
"""Utilities designed to help create the web summary used in Spatial."""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

import martian
import pandas as pd

import cellranger.rna.library as rna_library
import cellranger.websummary.analysis_tab_core as cr_atc
import cellranger.websummary.plotly_tools as pltly
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
from cellranger.spatial.data_utils import IMAGEX_LOWRES, IMAGEY_LOWRES

# pylint: enable=no-name-in-module
from cellranger.spatial.image import WebImage, base64_encode_image
from cellranger.spatial.webp_image_utils import get_crop_bbox
from cellranger.websummary.numeric_converters import round_floats_in_list
from cellranger.websummary.react_components import (
    SharedCoordinatePlotCollection,
)

if TYPE_CHECKING:
    import numpy as np

SPATIAL_PIPELINE_NAME = "count"
SPATIAL_COMMAND_NAME = "Space Ranger"


class SpatialReportArgs(NamedTuple):
    sample_id: str
    sample_desc: str
    scalefactors: str
    matrix: str
    tissue_lowres_image: str
    detected_tissue_image: str
    qc_resampled_cyta_img: str
    qc_regist_target_img: str
    tissue_positions: str
    analysis: str
    target_set_name: str
    target_panel_summary: str
    feature_ref_path: str
    targeting_method: str
    reorientation_mode: str
    loupe_alignment_file: str
    filter_probes: bool
    aligner: str


def _make_array_plot_layout(title, encoded_image):
    xlo, xhi, ylo, yhi = encoded_image.cropbox

    layout = {
        "title": {
            "text": title,
            "yanchor": "top",
        },
        "autosize": True,
        "hovermode": "closest",
        "legend": {"y": 0.9, "yanchor": "top"},
        "xaxis": {
            "type": "linear",
            "fixedrange": True,
            "showgrid": False,
            "zeroline": False,
            "showticklabels": False,
            "range": [xlo, xhi],
            "scaleanchor": "y",
            "scaleratio": 1.0,
        },
        "yaxis": {
            "type": "linear",
            "fixedrange": True,
            "showgrid": False,
            "zeroline": False,
            "showticklabels": False,
            "range": [yhi, ylo],
            "domain": [0.0, 0.90],
        },
        "margin": {
            "b": 10,
            "t": 10,
        },
        "images": [
            {
                "source": encoded_image.base64_encoded_str,
                "xref": "x",
                "yref": "y",
                "x": 0,
                "y": 0,
                "sizex": encoded_image.width,
                "sizey": encoded_image.height,
                "sizing": "stretch",
                "opacity": 1,
                "layer": "below",
            },
        ],
        "sliders": [
            {
                "pad": {"t": 10, "r": 0, "b": 10},
                "yanchor": "center",
                "xanchor": "left",
                "active": 4,
                "showactive": True,
                "tickwidth": 1,
                "currentvalue": {
                    "xanchor": "center",
                    "prefix": "Spot Opacity: ",
                    "font": {"color": "#888", "size": 12},
                    "visible": False,
                },
                "steps": [
                    {"label": "0.0", "method": "restyle", "args": ["opacity", 0.00]},
                    {"label": "", "method": "restyle", "args": ["opacity", 0.25]},
                    {"label": "", "method": "restyle", "args": ["opacity", 0.50]},
                    {"label": "", "method": "restyle", "args": ["opacity", 0.75]},
                    {"label": "1.0", "method": "restyle", "args": ["opacity", 1.00]},
                ],
            }
        ],
    }
    return layout


def make_array_plot(title, encoded_image, data):
    """Generate plot layout for displaying spatial images overlaid with spots."""
    assert isinstance(encoded_image, WebImage)
    return {
        "config": pltly.SPATIAL_PLOT_CONFIG,
        "layout": _make_array_plot_layout(title, encoded_image),
        "data": data,
    }


def umi_on_spatial_plot(sample_data, coords, encoded_image, library_type):
    """:param sample_data: Instance of webshim.data.SampleData.

    :param encoded_image: Instance of WebImage
    :return:
    """
    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None

    # make a df with the barcodes as the index
    bc_umi_counts = pd.DataFrame(index=analysis.matrix.bcs)
    # Add the UMI column by getting the sum UMI for each barcode
    matrix = analysis.matrix.select_features_by_type(library_type)
    umi_per_bc = matrix.get_counts_per_bc()
    # This append operation converts the ints to floats to account for the nan
    # after dropping missing values we convert back
    org_type = umi_per_bc.dtype
    umi_col_name = "UMI"
    bc_umi_counts[umi_col_name] = umi_per_bc
    cbd = coords.join(bc_umi_counts)
    cbd = cbd.dropna()
    cbd = cbd.astype({umi_col_name: org_type})
    color = cr_atc.get_umi_color(cbd[umi_col_name].values)
    title = "Tissue Plot with Spots Colored by UMI Count"
    data = [
        {
            "name": "Spots",
            "x": round_floats_in_list(cbd[IMAGEX_LOWRES].values),
            "y": round_floats_in_list(cbd[IMAGEY_LOWRES].values),
            "type": "scatter",
            "mode": "markers",
            "marker": {
                "opacity": 0.9,
                "size": encoded_image.markersize,
                "sizemode": "diameter",
                "color": color,
                "colorscale": "Jet",
                "colorbar": {
                    "title": "UMI counts",
                    "yanchor": "top",
                    "y": 0.9,
                    "len": 0.9,
                    "xpad": 20,
                },
            },
            "hovertemplate": "(%{x:.3f}, %{y:.3f}) <br> " + "UMI counts: %{text}<extra></extra>",
            "text": [f"{umis:,d}" for umis in cbd[umi_col_name].values],
        }
    ]
    return make_array_plot(title, encoded_image, data)


def tissue_plots_by_clustering_spatial(
    sample_data, coords, encoded_image, library_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE
):
    """Get the data structure that represents the plots for all the tissues.

    :param sample_data: Instance of webshim.data.SampleData
    :param coords:
    :param encoded_image: Instance of WebImage
    :return: A SharedCoordinatePlotCollection object
    """
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(SingleGenomeAnalysis)

    return tissue_plots_by_clustering_spatial_helper(analysis, coords, encoded_image, library_type)


def tissue_plots_by_clustering_spatial_helper(
    analysis,
    coords,
    encoded_image,
    library_type,
    title="Tissue Plot with Spots Colored by Clustering",
):
    """Get the data structure that represents the plots for all the tissues.

    Args:
        analysis: the loaded analysis h5 object
        coords:
        encoded_image: Instance of WebImage
        library_type: GEX or Antibody

    Returns:
        SharedCoordinatePlotCollection: object
    """
    if analysis is None:
        return None
    # order clustering by order: graph, kmeans=2, =3 etc
    clusterings = cr_atc.sort_analysis_clusterings(analysis.clusterings)
    # align clustering label with image by barcode
    if library_type == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
        clusters_dict = {
            clustering_key: clustering.clusters
            for (clustering_key, clustering) in clusterings
            if clustering_key[0:8] != "antibody"
        }
    elif library_type == rna_library.ANTIBODY_LIBRARY_TYPE:
        clusters_dict = {
            clustering_key: clustering.clusters
            for (clustering_key, clustering) in clusterings
            if clustering_key[0:8] == "antibody"
        }
    clusters_dict["barcode"] = analysis.matrix.bcs
    clusters_df = pd.DataFrame.from_dict(clusters_dict)
    clusters_df.set_index("barcode", inplace=True)
    cluster_tissue_df = clusters_df.join(coords[[IMAGEX_LOWRES, IMAGEY_LOWRES]])
    # pylint: disable=too-many-function-args
    tissue_plot_data = SharedCoordinatePlotCollection(
        cluster_tissue_df[IMAGEX_LOWRES].to_list(),
        cluster_tissue_df[IMAGEY_LOWRES].to_list(),
        pltly.SPATIAL_PLOT_CONFIG,
        _make_array_plot_layout(title, encoded_image),
        "scatter",
        {"size": encoded_image.markersize, "sizemode": "diameter"},
    )
    for clustering_key, clustering in clusterings:
        if clustering_key in clusters_dict:
            tissue_plot_data.add_clustering(clustering_key, clustering)
    return tissue_plot_data


def make_path_and_encode(
    chart, filename: str, lossless=True, crop_box: np.ndarray | None = None, padding: int = 0
):
    """Make a path for the chart and encode it as a base64 image.

    Args:
        chart (PIL.Image.Image): An image object to be encoded.
        filename (str): path to the image file.
        lossless (bool): Whether to save the image using lossless compression.
        crop_box (np.ndarray): 2D array of coordinates of crop box.
        padding (int): Padding to use. Defaults to 0.

    Returns:
        str: base64 encoded image.
    """
    if crop_box is not None:
        crop_bbox = get_crop_bbox(chart, crop_box=crop_box, padding=padding)
        chart = chart.crop(crop_bbox)
    chart_path = martian.make_path(f"{filename}.webp").decode()
    chart.save(chart_path, format="WEBP", lossless=lossless)
    return base64_encode_image(chart_path)


def scale_geometry(geometry, scale_factor):
    """Scales the coordinates of a given geometry by a specified scale factor.

    Args:
        geometry (dict): A dictionary representing the geometry to be scaled.
                        It must have a "type" key with the value "Polygon" and a "coordinates" key
                        containing a list of rings, where each ring is a list of [x, y] coordinate pairs.
        scale_factor (float): The factor by which to scale the coordinates.

    Returns:
        dict: A new geometry dictionary with the scaled coordinates.
    """
    if geometry["type"] == "Polygon":
        scaled_coordinates = [
            [[x * scale_factor, y * scale_factor] for x, y in ring]
            for ring in geometry["coordinates"]
        ]
        return {"type": geometry["type"], "coordinates": scaled_coordinates}
    else:
        raise ValueError("Unsupported geometry type")
