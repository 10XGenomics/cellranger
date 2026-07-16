#!/usr/bin/env python
#
# Copyright (c) 2025 10x Genomics, Inc. All rights reserved.
#
"""Compute segmentation plots for Visium HD data."""
import csv
import json
import os
from collections import Counter, defaultdict

import cv2
import martian
import numpy as np
import skimage
from PIL import Image

import cellranger.spatial.hd_cs_websummary_color_utils as hd_cs_websummary_color_utils
import cellranger.spatial.multi_layered_images as mli
from cellranger.fast_utils import CellId
from cellranger.spatial.image import base64_encode_image
from cellranger.websummary.spatial_utils import make_path_and_encode, scale_geometry

__MRO__ = """
stage SEGMENTATION_SPATIAL_PLOTS(
    in  path    analysis_csv,
    in  csv     cell_types,
    in  geojson cell_segmentations,
    in  bool    is_standalone_pipeline,
    in  webp    websummary_tissue_image,
    in  float   websummary_tissue_image_scale,
    out json    spatial_segmentation_chart,
    src py      "stages/spatial/segmentation_spatial_plots",
) split (
) using (
    volatile = strict,
)
"""

SEGMENTATION_MAX_RENDER_WIDTH = 6000
FOCUS_PADDING = 30

CLUSTER_KEY = "Cluster"
BARCODE_KEY = "Barcode"

ANALYSIS_CSV_CLUSTERING_NAME = "clustering"
ANALYSIS_CSV_GRAPHCLUST_CLUSTERING_DIR_NAME = "gene_expression_graphclust"
ANALYSIS_CSV_CLUSTERING_CSV_NAME = "clusters.csv"

CAS_BARCODE_KEY = "barcode"
CAS_COARSE_CELL_TYPES_KEY = "coarse_cell_type"


def get_cell_type_coloring(
    celltype_csv_path,
) -> tuple[dict[str, list[int]], dict[int, list[int]], dict[str, int]]:
    """Get cluster coloring from cell type CSV path."""
    cell_cluster_dict = {}
    with open(celltype_csv_path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            bc_num = CellId(row[CAS_BARCODE_KEY]).id
            cell_cluster_dict[bc_num] = row[CAS_COARSE_CELL_TYPES_KEY]
    cluster_order = {x: -y for x, y in Counter(cell_cluster_dict.values()).items()}
    # Apply 10x cluster colors
    cluster_list = sorted(set(cell_cluster_dict.values()))  # Ensure clusters are sorted
    num_clusters = len(cluster_list)
    cluster_colors = {
        cluster: hd_cs_websummary_color_utils.cluster_color_rgb(ind + 1, num_clusters)
        for ind, cluster in enumerate(cluster_list)
    }
    # Map each cell_id to its RGB cluster color
    cluster_color_dict = {
        cell_id: cluster_colors[cluster_idx] for cell_id, cluster_idx in cell_cluster_dict.items()
    }
    return cluster_colors, cluster_color_dict, cluster_order


def get_cluster_coloring(
    cluster_csv_path,
) -> tuple[dict[str, list[int]], dict[int, list[int]], dict[str, int]]:
    """Get cluster coloring from cluster CSV path."""
    cell_cluster_dict = {}
    with open(cluster_csv_path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            bc_num = CellId(row[BARCODE_KEY]).id
            cell_cluster_dict[bc_num] = int(row[CLUSTER_KEY])
    # Apply 10x cluster colors
    num_clusters = len(set(cell_cluster_dict.values()))  # Ensure clusters are sorted
    cluster_colors = {
        i + 1: hd_cs_websummary_color_utils.cluster_color_rgb(i + 1, num_clusters)
        for i in range(num_clusters)
    }
    # Map each cell_id to its RGB cluster color
    cluster_color_dict = {
        cell_id: cluster_colors[cluster_idx] for cell_id, cluster_idx in cell_cluster_dict.items()
    }
    final_cluster_colors = dict(
        (f"Cluster{cluster_idx}", color) for cluster_idx, color in cluster_colors.items()
    )
    cluster_order = {f"Cluster{cluster_idx}": cluster_idx for cluster_idx in cluster_colors}
    return final_cluster_colors, cluster_color_dict, cluster_order


def build_composite_indexed_image(
    sorted_clusters: list[tuple],
    coords_by_color: dict,
    height: int,
    width: int,
) -> str:
    """Return base64 data URI of an 8-bit grayscale+alpha WebP lossless.

    Pixel grayscale value = 1-based cluster index; alpha = 255 for cluster pixels,
    0 for background. Supports up to 255 clusters.
    """
    index_array = np.zeros((height, width), dtype=np.uint8)
    for cluster_idx, (_, edge_color, _) in enumerate(sorted_clusters, start=1):
        for coords in coords_by_color[tuple(edge_color)]:
            cv2.fillPoly(index_array, [coords], color=cluster_idx)

    # Build LA (grayscale + alpha) image: L = index, A = 255 if cluster else 0
    la_array = np.zeros((height, width, 2), dtype=np.uint8)
    la_array[:, :, 0] = index_array
    la_array[index_array > 0, 1] = 255

    img = Image.fromarray(la_array, mode="LA")
    path = martian.make_path("composite_indexed.webp").decode()
    img.save(path, format="WEBP", lossless=True)
    return base64_encode_image(path)


def split(args):
    if not all(
        (
            args.cell_segmentations,
            args.websummary_tissue_image,
            args.websummary_tissue_image_scale is not None,
        )
    ):
        return {"chunks": [], "join": {"__mem_gb": 1}}

    if not args.is_standalone_pipeline:
        return {"chunks": [], "join": {"__mem_gb": 16}}
    else:
        return {"chunks": [], "join": {"__mem_gb": 32, "__vmem_gb": 64}}


def join(
    args, outs, _chunk_defs, _chunk_outs
):  # pylint: disable=too-many-locals, too-many-statements
    if not all(
        (
            args.cell_segmentations,
            args.websummary_tissue_image,
            args.websummary_tissue_image_scale is not None,
        )
    ):
        martian.clear(outs)
        return

    with open(args.cell_segmentations) as f:
        cell_segmentations = json.load(f)

    tissue_hires_image_scale = args.websummary_tissue_image_scale

    for feature in cell_segmentations["features"]:
        feature["geometry"] = scale_geometry(feature["geometry"], tissue_hires_image_scale)

    tissue_image_height, tissue_image_width = skimage.io.imread(args.websummary_tissue_image).shape[
        :2
    ]

    # Compute bounding box of all polygon coordinates (websummary image space) for initial focus.
    poly_coords = [
        coord
        for feature in cell_segmentations["features"]
        if feature["geometry"]["type"] == "Polygon"
        for coord in feature["geometry"]["coordinates"][0]
    ]
    if poly_coords:
        coords_arr = np.array(poly_coords)
        focus_x = max(0, int(coords_arr[:, 0].min()) - FOCUS_PADDING)
        focus_y = max(0, int(coords_arr[:, 1].min()) - FOCUS_PADDING)
        focus_width = min(tissue_image_width, int(coords_arr[:, 0].max()) + FOCUS_PADDING) - focus_x
        focus_height = (
            min(tissue_image_height, int(coords_arr[:, 1].max()) + FOCUS_PADDING) - focus_y
        )
    else:
        focus_x, focus_y = 0, 0
        focus_width, focus_height = tissue_image_width, tissue_image_height

    # Downsample overlay rasters to at most SEGMENTATION_MAX_RENDER_WIDTH pixels wide.
    # Polygon coordinates are scaled before rendering so no post-processing resampling
    # is needed — nearest-neighbour accuracy is preserved by construction.
    render_scale = min(1.0, SEGMENTATION_MAX_RENDER_WIDTH / tissue_image_width)
    render_width = round(tissue_image_width * render_scale)
    render_height = round(tissue_image_height * render_scale)

    has_annotations = False
    # Extract and clean the cell_id
    if args.analysis_csv and os.path.exists(
        cluster_csv_path := os.path.join(
            args.analysis_csv,
            ANALYSIS_CSV_CLUSTERING_NAME,
            ANALYSIS_CSV_GRAPHCLUST_CLUSTERING_DIR_NAME,
            ANALYSIS_CSV_CLUSTERING_CSV_NAME,
        )
    ):
        cluster_colors, cluster_color_dict, cluster_order = get_cluster_coloring(cluster_csv_path)
        default_color = None
        has_annotations = True
    elif args.cell_types:
        cluster_colors, cluster_color_dict, cluster_order = get_cell_type_coloring(args.cell_types)
        default_color = None
        has_annotations = True
    else:
        cluster_colors, cluster_color_dict, cluster_order, default_color = (
            {},
            {},
            {},
            hd_cs_websummary_color_utils.cluster_color_rgb(1, 1),
        )

    # Group all polygon coordinates by edge_color
    coords_by_color = defaultdict(list)
    for feature in cell_segmentations["features"]:
        if feature["geometry"]["type"] != "Polygon":
            continue
        cell_id = feature["properties"].get("cell_id")
        edge_color = cluster_color_dict.get(cell_id, default_color)
        if edge_color is None:
            continue
        coords = (
            (np.array(feature["geometry"]["coordinates"][0], dtype=np.float64) * render_scale)
            .round()
            .astype(np.int32)
        )
        coords_by_color[tuple(edge_color)].append(coords)

    # initialize the rasters at render resolution
    geometry_raster = np.zeros((render_height, render_width, 4), dtype=np.uint8)

    # Get the segmentation layer — all boundaries in black so shared edges
    # between adjacent cells have a consistent, unambiguous appearance.
    # When there are no annotations, use bright green so cells are visible.
    boundary_color = (0, 255, 56, 255) if not has_annotations else (0, 0, 0, 255)
    for poly_list in coords_by_color.values():
        for coords in poly_list:
            cv2.polylines(
                geometry_raster,
                [coords],
                isClosed=True,
                color=boundary_color,
                thickness=1,
            )

    geometry_img = Image.fromarray(geometry_raster, mode="RGBA")
    cell_segmentations_encoded = make_path_and_encode(
        geometry_img,
        "cell_segmentation",
        lossless=True,
    )
    # If rendered at reduced resolution, upscale via css_transform so the <img>
    # displays at the same physical size as the full-resolution tissue image.
    inv_rs = 1.0 / render_scale
    segmentation_image = [
        mli.LabeledImage(
            label="Segmentations",
            image=cell_segmentations_encoded,
            color=None,
            css_transform=(
                [inv_rs, 0, 0, 0, 0, inv_rs, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
                if render_scale < 1.0
                else None
            ),
        )
    ]

    cluster_layers = None
    if has_annotations:
        saved_clusters = []
        cluster_color_mapping = {tuple(v): k for k, v in cluster_colors.items()}
        for edge_color in coords_by_color:
            cluster_name = cluster_color_mapping.get(edge_color)
            saved_clusters.append((cluster_name, edge_color, None))
        saved_clusters.sort(key=lambda x: cluster_order[x[0]])

        composite_uri = build_composite_indexed_image(
            saved_clusters,
            coords_by_color,
            render_height,
            render_width,
        )

        cluster_layers = mli.Layer(
            name="Cell Types" if args.cell_types else "Clusters",
            images=[
                mli.LabeledImage(
                    label=cluster_name,
                    image=None,
                    color=hd_cs_websummary_color_utils.rgb_to_hex(edge_color),
                    css_transform=None,
                )
                for cluster_name, edge_color, _ in saved_clusters
            ],
            initial_opacity=0.60,
            composite_indexed_image=composite_uri,
            render_scale=render_scale if render_scale < 1.0 else None,
        )

    tissue_image = [
        mli.LabeledImage(
            label="Tissue image",
            image=base64_encode_image(args.websummary_tissue_image),
            color=None,
            css_transform=None,
        )
    ]

    # Legend order: Tissue → Segmentations → Clusters/Cell Types.
    # z_index overrides DOM order so Segmentations (z=2) renders above
    # Clusters (z=1) regardless of their position in the array.
    layers = [
        mli.Layer(name="Tissue image", images=tissue_image, z_index=0),
        mli.Layer(name="Segmentations", images=segmentation_image, initial_opacity=0.60, z_index=2),
    ]
    if cluster_layers:
        cluster_layers.z_index = 1
        layers.append(cluster_layers)

    final_layers = mli.MultiLayerImages(
        focus=mli.InitialFocus(
            x=focus_x,
            y=focus_y,
            width=focus_width,
            height=focus_height,
        ),
        layers=layers,
        fullScreen=False,
        checkAll=True,
        showFullScreenButton=False,
        legendTitle=None,
        legendWidthPx=125,
    )

    with open(outs.spatial_segmentation_chart, "w") as f:
        json.dump(final_layers, f, indent=2, default=lambda x: x.__dict__)
