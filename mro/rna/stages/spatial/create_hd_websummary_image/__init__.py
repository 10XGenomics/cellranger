# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
"""Generate the tissue image for the HD websummary."""

import json

import cv2
from PIL import Image

from cellranger.spatial.data_utils import DARK_IMAGES_CHANNELS, HD_WS_MAX_DIM
from cellranger.spatial.image_util import normalized_image_from_counts, shrink_to_max

__MRO__ = """
stage CREATE_HD_WEBSUMMARY_IMAGE(
    in  png   tissue_hires_image,
    in  json  scalefactors,
    in  int   dark_images,
    in  bool  is_standalone_pipeline,
    out webp  websummary_tissue_image,
    out float websummary_tissue_image_scale,
    src py    "stages/spatial/create_hd_websummary_image",
) using (
    mem_gb   = 2,
    vmem_gb  = 10,
    volatile = strict,
)
"""


def main(args, outs):
    if args.dark_images == DARK_IMAGES_CHANNELS:
        tissue_hires_image = cv2.imread(args.tissue_hires_image, cv2.IMREAD_GRAYSCALE)
        tissue_hires_image = normalized_image_from_counts(
            tissue_hires_image, log1p=True, invert=True
        )
        hires_to_ws_scale = shrink_to_max(tissue_hires_image.shape, HD_WS_MAX_DIM)
        tissue_ws_image = cv2.resize(
            tissue_hires_image, (0, 0), fx=hires_to_ws_scale, fy=hires_to_ws_scale
        )
        Image.fromarray(tissue_ws_image).save(
            outs.websummary_tissue_image, format="WEBP", quality=20, method=6
        )
    else:
        quality = 80 if args.is_standalone_pipeline else 20
        with Image.open(args.tissue_hires_image) as img:
            hires_to_ws_scale = shrink_to_max((img.height, img.width), HD_WS_MAX_DIM)
            if hires_to_ws_scale < 1.0:
                img = img.resize(
                    (int(img.width * hires_to_ws_scale), int(img.height * hires_to_ws_scale)),
                    Image.Resampling.LANCZOS,
                )
            img.save(outs.websummary_tissue_image, format="WEBP", quality=quality, method=6)

    with open(args.scalefactors) as f:
        scalefactors = json.load(f)
    outs.websummary_tissue_image_scale = scalefactors["tissue_hires_scalef"] * hires_to_ws_scale
