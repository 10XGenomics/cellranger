#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET

import cv2
import martian
import numpy as np
import skimage.io
from PIL import Image, TiffImagePlugin

import cellranger.spatial.tiffer as tiffer
from cellranger.spatial.data_utils import DARK_IMAGES_CHANNELS, DARK_IMAGES_NONE

CYT_IMG_PIXEL_SIZE = 4.625  # in microns
SD_SPOT_TISSUE_AREA_UM2 = 8660.0  # in sq.microns (100 um) * (100 * sqrt(3)/2 um)
HD_SPOT_TISSUE_AREA_UM2 = 4.0  # in sq.microns
XML_HEADER_PREFIXES = ["<?xml", "<OME xml"]  # prefix of an XML string
REGIST_TARGET_IMAGE_MAX_DIM = 6000  # Maximum dimension of the registration target image


def normalized_image_from_counts(
    counts: np.ndarray, maxcount=None, log1p=False, invert=False
) -> np.ndarray:
    """Given a matrix of UMI/Read counts, generate a normalized image that is optionally log transformed and/or inverted.

    Args:
        counts (np.ndarray): 2D UMI/Read counts
        maxcount (int, optional): Normalize using this max count
        log1p (bool, optional): Apply log1p transform. Defaults to False.
        invert (bool, optional): Invert the image. Defaults to False.

    Returns:
        np.ndarray: Normalized 8-bit grayscale image
    """
    image = counts.astype("float")
    if maxcount is None:
        maxcount = np.amax(image)

    if log1p:
        image = np.log1p(image)
        maxcount = np.log1p(maxcount)

    if invert:
        image = maxcount - image
    return (image * (255 / maxcount)).astype("uint8")


def shrink_to_max(shape, target):
    """Given an image shape (x,y), if the largest dimension is larger than.

    target, return a scalefactor to shrink the image to approximately target in that dimension,
    otherwise return 1.0
    """
    scale = 1.0 * target / max(shape)
    return min(1.0, scale)


def cv_read_image_standard(filename, maxval=255):
    """Canonical function to read image files."""
    img = cv2.imread(filename, cv2.IMREAD_UNCHANGED)
    if img is not None:
        if len(img.shape) == 3:
            img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
        cv2.normalize(img, img, maxval, 0, cv2.NORM_MINMAX)
        if img.dtype == "uint16":
            img = img.astype("uint8")

    return img


def cv_write_downsampled_image(cvimg, filename, maxres):
    """Scale, resize, and return the downsampled image."""
    scalef = shrink_to_max(cvimg.shape, maxres)
    rcvimg = cv2.resize(cvimg, (0, 0), fx=scalef, fy=scalef)
    params = [cv2.IMWRITE_JPEG_QUALITY, 80]
    cv2.imwrite(filename, rcvimg, params)

    return scalef


def cv_read_rgb_image(filepath) -> np.ndarray:
    """Read an RGB image."""
    img = cv2.imread(filepath, cv2.IMREAD_COLOR)  # red in BGR fromat
    return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)


def cv_composite_labelmap(background, labelmap, colormap, alpha):
    """Take a background image (generally grayscale) and composite with.

    a color image based on a labelmap

    alpha value is [0,1.0)  (note open - can't be 1.0)

    Example:
    cv_composite_labelmap(gray, mask, {0:(255,0,0), 1:(0,0,255), 0.75})

    takes a binary labelmap (mask) and composites it with the gray image
    to return a color image with 0.25 (25%) background and 0.75 (75%) labelmap
    where the zero mask values are red (255,0,0) and the one mask values are
    blue (0,0,255)
    """
    if background.shape[:2] != labelmap.shape[:2]:
        raise RuntimeError("cv_composite_labelmap called with incompatible images")

    if alpha < 0.0 or alpha >= 1.0:
        raise ValueError("cv_composite_labelmap passed illegal alpha value - should be [0.,1.)")

    if background.ndim < 3:
        output = cv2.cvtColor(background, cv2.COLOR_GRAY2RGB)
    else:
        output = background.copy()

    for label, color in colormap.items():
        mask = labelmap == label
        output[mask, :] = (1.0 - alpha) * output[mask, :] + np.multiply(alpha, color)

    return output


def write_detected_keypoints_image(img, keypoints, diameter, filename):
    """Draw the detected keypoints over a brightfield image."""
    output = img.copy()
    for kp in keypoints:
        col, row = int(kp.pt[0]), int(kp.pt[1])
        cv2.circle(output, (col, row), int(diameter / 2.0), (0, 255, 0), 3)
    cv2.imwrite(filename, output)


def write_aligned_fiducials_image(img, aligned_fiducials, diameter, filename):
    """Draw the aligned fiducials over a brightfield image."""
    output = img.copy()
    for spot in aligned_fiducials:
        x, y = spot
        col, row = int(x), int(y)
        embiggen = 1.2  # add 20% to spot size to be able to visualize underlying fiducial
        cv2.circle(output, (col, row), int(diameter / 2.0 * embiggen + 0.5), (0, 0, 255), 2)
    cv2.imwrite(filename, output)


def generate_fiducial_mask(
    mask_shape: tuple[int, int], fiducials_xy: np.ndarray, diameter: float, embiggen: float = 2.0
) -> np.ndarray:
    """Generate a mask that covers the fiducial marker locations."""
    fiducial_mask = np.zeros(mask_shape, dtype=np.uint8)
    for x, y in fiducials_xy:
        cv2.circle(
            fiducial_mask,
            (int(x), int(y)),
            int(diameter / 2.0 * embiggen + 0.5),
            255,
            -1,
        )
    return fiducial_mask


def _check_sizes(sizes, tmp_json):
    """Checks json response from tiffer against width and height pair in sizes."""
    width = int(tmp_json["width"])
    height = int(tmp_json["height"])
    if sizes is None:
        return (width, height)
    elif sizes[0] != width or sizes[1] != height:
        martian.throw(
            f"input images/pages have inconsistent sizes ({sizes[0]},{sizes[1]}) vs ({width},{height})"
        )

    return sizes


def select_multi_page_downsample(
    image_file_paths: list[str],
    target_size: int,
    save_to_path: str,
    skip_pages: set,
    selected_pages: list[int] | None = None,
) -> tuple[list[str], float]:
    """Downsample and save selected page from multipage tiff.

    Args:
        image_file_paths (str): image path of the multipage tiff
        target_size (int): downsample the page to this size
        save_to_path (str): save the selected page to this path
        skip_pages (Set): skip these pages
        selected_pages (List[int]): selected pages

    Raises:
        ValueError: The selected page number is out of the range

    Returns:
        List[str]: list of path of all the downsampled channels
        float: scale factor of all the downsampled channels
    """
    npages = tiffer.call_tiffer_get_num_pages(image_file_paths[0])
    downsample_img_paths = []
    ds_scalef = None
    target_basename, target_ext = os.path.splitext(save_to_path)

    if len(image_file_paths) == 1 and npages > 1:
        if selected_pages:
            for page in selected_pages:
                if page >= npages or page < 0:
                    raise ValueError(
                        f"Page {page+1} (1-indexed) does not exist for image with {npages} pages."
                    )
            pages_idx = selected_pages
        else:
            pages_idx = range(npages)

        sizes = None
        for index in pages_idx:
            if index not in skip_pages:
                tmp_json = tiffer.call_tiffer_resample(
                    image_file_paths[0],
                    target_size,
                    target_basename + "_" + str(index) + target_ext,
                    page=index,
                )
                downsample_img_paths.append(tmp_json["outpath"])
                ds_scalef = tmp_json["scale"]
                _check_sizes(sizes, tmp_json)
    else:
        sizes = None
        num_channels = len(image_file_paths)
        if selected_pages:
            for page in selected_pages:
                if page >= num_channels:
                    raise ValueError(
                        f"Page {page} does not exist for image with {len(image_file_paths)} pages."
                    )
            pages_idx = selected_pages
        else:
            pages_idx = range(num_channels)

        for index in pages_idx:
            tmp_json = tiffer.call_tiffer_resample(
                image_file_paths[index],
                target_size,
                target_basename + "_" + str(index) + target_ext,
            )
            downsample_img_paths.append(tmp_json["outpath"])
            ds_scalef = tmp_json["scale"]
            _check_sizes(sizes, tmp_json)

    return downsample_img_paths, ds_scalef


def prepare_registration_target_image(
    tissue_image_paths: list[str],
    dark_images: int,
    regist_target_im_max_dim: int,
    save_to_path: str,
    skip_pages: set,
    selected_pages: list[int] | None = None,
    read_metadata_from_tissue_image: bool | None = False,
) -> float:
    """Preprocess the tissue image as the target for registration.

    Args:
        tissue_image_paths (List[str]): list of tissue image paths
        dark_images (int): whether the image is dark image
        regist_target_im_max_dim (int): maximum dimension of the registration target
        save_to_path (str): path to save the registration target image
        skip_pages (Set): skip these pages.
        selected_pages (List[int]): selected pages for the registration target

    Returns:
        float: scale factor of the regitration target
    """
    if dark_images == DARK_IMAGES_CHANNELS:
        downsample_img_paths, scalef = select_multi_page_downsample(
            tissue_image_paths,
            regist_target_im_max_dim,
            "raw_channel.png",
            skip_pages,
            selected_pages,
        )
        process_combine_multi_channel(downsample_img_paths, save_to_path)
    else:
        regist_im_json = tiffer.call_tiffer_resample(
            tissue_image_paths[0],
            regist_target_im_max_dim,
            "raw_rgb.png",
        )
        img = skimage.io.imread(regist_im_json["outpath"])
        if dark_images == DARK_IMAGES_NONE:
            # BF H&E, save the green channel as registration target for RGB image.
            target_img = img[:, :, 1] if img.ndim == 3 else img
        else:
            # Colorized IF image
            processed_img_list = []
            for i in range(3):
                processed_img_list.append(preprocess_fl_channel(img[:, :, i], method="log"))
            target_img = np.max(np.stack(processed_img_list, axis=2), axis=2)
            target_img = 255 - target_img
        skimage.io.imsave(save_to_path, target_img)
        scalef = float(regist_im_json["scale"])

    if read_metadata_from_tissue_image:
        # read the tiff file back using PIL.Image and get the tag dictionary
        # so we can write the updated tiff info
        reg_target_img = Image.open(save_to_path)
        tiff_info = reg_target_img.tag_v2

        # get image description metadata for the registration target image
        #
        # the expectation is that all tiff images in tissue_image_paths
        # are of the same size, so we make a reasonable assumption that
        # the metadata from the first image in that list will apply to all
        # images
        image_desc = get_image_description_for_reg_target_image(
            tissue_image_paths[0], reg_target_img.width, reg_target_img.height, scalef
        )

        if image_desc:
            # add new image description
            tiff_info[TiffImagePlugin.IMAGEDESCRIPTION] = image_desc
        elif TiffImagePlugin.IMAGEDESCRIPTION in tiff_info:
            # tiffer writes out misformed metadata. Removing it for consistency
            tiff_info.pop(TiffImagePlugin.IMAGEDESCRIPTION)

        # now save the registration target image with updated image description
        reg_target_img.save(save_to_path, tiffinfo=tiff_info)

    return scalef


def process_combine_multi_channel(img_path_list: list[str], save_to_path: str) -> None:
    """Process the images and combine them into one.

    Args:
        img_path_list (List[str]): list of path of all the images
        save_to_path (str): path to save the final
    """
    processed_img_list = []
    for img_path in img_path_list:
        raw_img = skimage.io.imread(img_path)
        processed_img = preprocess_fl_channel(raw_img, method="log")
        processed_img_list.append(processed_img)
    if len(processed_img_list) > 1:
        final_img = np.max(np.stack(processed_img_list, axis=2), axis=2)
    else:
        final_img = processed_img_list[0]
    # save the reverted image as the registration target since it helps with
    # the initialization of the registration and won't affect the final registration
    reverted_img = 255 - final_img
    skimage.io.imsave(save_to_path, reverted_img)


def preprocess_fl_channel(
    image: np.ndarray[tuple[int, int], np.dtype[np.int16]], method: str = "log"
) -> np.ndarray[tuple[int, int], np.dtype[np.uint8]]:
    """Preprocess a single channel FL image.

    Preprocessing of the FL image can have multiple purposes. The "log" method is to
    bring up the weak autofluorescence of the tissue.

    Args:
        image (npt.NDArray[np.int16]): 2D array, a single channel of FL image.
        method (str, optional): preprocessing method. Defaults to "log".

    Raises:
        ValueError: preprocessing method is not support
        ValueError: the image is not a 2D image

    Returns:
        npt.NDArray[np.uint8]: the image after the preprocessing
    """
    if method not in ("log",):
        raise ValueError(f"method {method} is not suppported.")
    if len(image.shape) != 2:
        raise ValueError(f"Only support 2D image but image shape is {image.shape}")
    image[image == 0] += 1
    mask = image < 1
    log_img = np.log10(image)
    masked_log = np.ma.masked_array(log_img, mask=mask)
    masked_log = np.ma.filled(masked_log, np.nan)
    min_val = np.nanpercentile(masked_log, 1)
    max_val = np.nanpercentile(masked_log, 99)
    log_img = (log_img - min_val) / (max_val - min_val)
    log_img = np.clip(log_img, 0, 1)
    log_img = (log_img * 255).astype(np.uint8)
    return log_img


def get_xml_namespace(doc: ET.Element) -> str:
    """Get namespace from xml Element object."""
    # the namespace prefix is wrapped in curly braces
    m = re.match(r"\{.*\}", doc.tag)
    return m.group(0) if m else ""


def remove_xml_namespace(doc: ET.Element, namespace: str) -> None:
    """Remove namespace from xml Element object in place."""
    nsl = len(namespace)
    for elem in doc.iter():
        if elem.tag.startswith(namespace):
            elem.tag = elem.tag[nsl:]


def extract_metadata_xml_tree(microscope_img_path: str) -> ET.Element | None:
    """Extract the metadata XML from the microscope image path."""
    try:
        image_info = tiffer.try_call_tiffer_info(microscope_img_path)
    except RuntimeError as err:
        print(err)
        return None

    image_desc_xml_string = image_info.get(tiffer.TIFFER_INFO_IMAGE_DESCRIPTION_KEY)
    if not image_desc_xml_string:
        return None
    if not any(image_desc_xml_string.strip().startswith(x) for x in XML_HEADER_PREFIXES):
        print("Non-XML IMAGE_DESCRIPTION tag observed. " f"Observed tag :\n{image_desc_xml_string}")
        return None

    try:
        # parse the image description to get the root tag
        xml_tree = ET.fromstring(image_desc_xml_string)
    except ET.ParseError as exc:
        raise ValueError(
            "Microscope Image has OME compliant metadata with misformed metadata XML. "
            f"Observed XML :\n{image_desc_xml_string}"
        ) from exc

    # remove the namespace prefix from all the xml tags to
    # make searching easier and cleaner
    remove_xml_namespace(xml_tree, get_xml_namespace(xml_tree))

    return xml_tree


def get_image_description_for_reg_target_image(
    microscope_img_path: str, reg_target_width: int, reg_target_height: int, reg_target_scale: float
) -> str | None:
    """Generate a metadata xml string for the registration target image.

    The microscope image metadata is updated to get a new scale value based on
    the registration target image scale.
    Returns an XML string iff the metadata XML is well formed and contains an
    attribute named "Pixels" with the first such attribute having sub-attributes
    `"PhysicalSizeX"` and `"PhysicalSizeY"`.
    TODO: Rewrite the API with the OME's official library here:
    https://omero.readthedocs.io/en/v5.6.9/developers/Python.html
    Probably will need some dependency wrangling on our end.
    Note: the function assumes that the metadata is in the IMAGEDESCRIPTION tag
    and is in OME-compliant xml format.

    Args:
         microscope_img_path (str): path to the microscope image
         reg_target_width (int): width of the registration target image
         reg_target_height (int): height of the registration target image
         reg_target_scale (float): downsampling factor for registration
             target image

    Returns:
         str | None: updated image description. None is returned in case of
            errors.
    """
    # get xml tag
    xml_tree = extract_metadata_xml_tree(microscope_img_path)

    if xml_tree is None:
        return None

    # get the pixels tag
    pixels = xml_tree.findall(".//Pixels")

    if pixels:
        if "PhysicalSizeX" in pixels[0].attrib:
            new_physical_sizex = float(pixels[0].attrib["PhysicalSizeX"]) / reg_target_scale
            print(
                f"PhysicalSizeX read: {pixels[0].attrib['PhysicalSizeX']}.\n",
                f"registration target scale: {reg_target_scale}\n"
                f"written out: {new_physical_sizex}\n",
            )
            pixels[0].attrib["PhysicalSizeX"] = str(new_physical_sizex)
        else:
            return None

        if "PhysicalSizeY" in pixels[0].attrib:
            new_physical_sizey = float(pixels[0].attrib["PhysicalSizeY"]) / reg_target_scale
            print(
                f"PhysicalSizeY read: {pixels[0].attrib['PhysicalSizeY']}.\n",
                f"registration target scale: {reg_target_scale}\n"
                f"written out: {new_physical_sizey}\n",
            )
            pixels[0].attrib["PhysicalSizeY"] = str(new_physical_sizey)
        else:
            return None

        # add a new attribute to record the registration target scale factor
        pixels[0].attrib["PostCaptureDownsamplingScale"] = str(reg_target_scale)

        # registration target image is a single channel image
        pixels[0].attrib["SizeC"] = "1"

        if "SizeX" in pixels[0].attrib:
            pixels[0].attrib["SizeX"] = str(reg_target_width)

        if "SizeY" in pixels[0].attrib:
            pixels[0].attrib["SizeY"] = str(reg_target_height)

        # return the string version of the modified xml tree
        return ET.tostring(xml_tree, encoding="unicode")
    else:
        return None


def get_tiff_pixel_size(img_path: str) -> float | None:
    """Get pixel physical size in microns from tiff metadata.

    Note: the function assumes that the metadata is in the IMAGEDESCRIPTION tag
    and is in OME-compliant xml format.

    Args:
        img_path (str): path to tiff file

    Returns:
        float | None: physical pixel size in microns. None is returned in case
            of errors, or if pixel size cannot be reliably determined
    """
    # get xml tag
    xml_tree = extract_metadata_xml_tree(img_path)

    if xml_tree is None:
        return None

    # get the pixels tag
    pixels = xml_tree.findall(".//Pixels")

    if pixels:
        # physical size in X & Y should be very similar, so we average
        # the two to get a single value
        if "PhysicalSizeX" in pixels[0].attrib and "PhysicalSizeY" in pixels[0].attrib:
            return (
                float(pixels[0].attrib["PhysicalSizeX"]) + float(pixels[0].attrib["PhysicalSizeY"])
            ) / 2

    return None
