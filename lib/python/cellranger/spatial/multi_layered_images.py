# Copyright (c) 2025 10x Genomics, Inc. All rights reserved.
"""Classes for multi-layered images."""
from dataclasses import dataclass
from typing import Self

import martian
import numpy as np
from PIL import Image

from cellranger.spatial.image import base64_encode_image

TISSUE_IMAGE_DISPLAY_WIDTH = 400


def rgb_to_hex(rgb_tuple: tuple[int, int, int]) -> str:
    """RGB to HEX converter."""
    return "#" + "".join(f"{val:02X}" for val in rgb_tuple)


@dataclass
class InitialFocus:
    """Data for initial focus."""

    x: int
    y: int
    width: int
    height: int


@dataclass
class LabeledImage:
    """Data for a labelled image."""

    label: str
    color: str | None
    css_transform: list[float] | None
    image: str | None = None
    legend_image: str | None = None  # different colormaps

    @classmethod
    def get_cell_labelled_image(
        cls, cell_annotation_matrix, cluster_num, cluster_name, clrmap_in_rgb, css_transform=None
    ) -> Self:
        """Get cell labelled image from celltypes."""
        nrows, ncols = cell_annotation_matrix.shape
        spatial_plot_fname = martian.make_path("spatial_plot.png").decode()
        color_rgb = clrmap_in_rgb[cluster_num]
        color_hex = rgb_to_hex(color_rgb)
        cluster_mask = cell_annotation_matrix == cluster_num

        img = np.zeros((nrows, ncols, 4), dtype=np.uint8)
        img[cluster_mask, -1] = 255
        img[cluster_mask, :-1] = color_rgb

        Image.fromarray(img).save(spatial_plot_fname)

        return cls(
            label=cluster_name,
            color=color_hex,
            image=base64_encode_image(spatial_plot_fname),
            css_transform=css_transform,
        )

    @classmethod
    def get_tissue_labelled_image(cls, tissue_image, css_transform=None) -> Self:
        return cls(
            label="Tissue image",
            color=None,
            image=base64_encode_image(tissue_image),
            css_transform=css_transform,
        )


@dataclass
class Layer:
    """A layer in the multi-layer Image."""

    name: str
    images: list[LabeledImage]
    initial_opacity: float | None = 1.0
    composite_indexed_image: str | None = None
    render_scale: float | None = None
    z_index: int | None = None


@dataclass
class MultiLayerImages:  # pylint: disable=invalid-name
    """The final multi-layer image."""

    focus: InitialFocus
    layers: list[Layer]
    fullScreen: bool
    checkAll: bool = True  # start the plot with all boxes checked?
    showFullScreenButton: bool = False
    legendTitle: str | None = None
    legendRightOffset: int | None = 20
    legendWidthPx: int | None = 120


@dataclass
class NamedMultiLayerImages:
    """Named multi-layer image."""

    name: str
    multi_layer_image: MultiLayerImages
