# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
"""Webp image utils."""

import numpy as np
from PIL import Image


def compute_crop_bbox(width, height, crop_box, padding):
    """Compute crop bbox."""
    if padding is None:
        padding = 0
    left = max(int(np.min(crop_box[:, 0]) - padding), 0)
    upper = max(int(np.min(crop_box[:, 1]) - padding), 0)
    right = min(int(np.max(crop_box[:, 0]) + padding), width)
    lower = min(int(np.max(crop_box[:, 1]) + padding), height)
    return (left, upper, right, lower)


def get_crop_bbox(img: Image.Image, crop_box: np.ndarray, padding: int | None = None):
    """Get crop bbox."""
    return compute_crop_bbox(img.width, img.height, crop_box, padding)


def convert_img_to_webp(
    input_image_path, output_webp_path, quality=50, crop_box=None, padding: int = 0
):
    """Convert an image to a WebP image with the specified quality and optional cropping.

    Args:
        input_image_path (str): Path to the input PNG image.
        output_webp_path (str): Path to save the output WebP image.
        quality (int): Quality of the WebP image (default is 50).
        crop_box (np.ndarray): 2D array of coordinates of crop box.
        padding (int): padding around the crop box. Defaults to 0.

    Returns:
        str: Path to the saved WebP image.
    """
    with Image.open(input_image_path) as img:
        if crop_box is not None:
            crop_bbox = get_crop_bbox(img, crop_box=crop_box, padding=padding)
            img = img.crop(crop_bbox)
        img.save(output_webp_path, format="WEBP", quality=quality, method=6)
    return output_webp_path
