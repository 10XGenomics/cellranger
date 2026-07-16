#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import json
import os
import subprocess
from typing import TYPE_CHECKING

import martian
import numpy as np

import cellranger.cr_io as cr_io
import tenkit.log_subprocess as tk_subproc

if TYPE_CHECKING:
    from cellranger.spatial.bounding_box import BoundingBox

TIFFER_INFO_PAGE_KEY = "pages"
TIFFER_INFO_HEIGHT_KEY = "height"
TIFFER_INFO_WIDTH_KEY = "width"
TIFFER_INFO_COLORMODE_KEY = "colorMode"
TIFFER_INFO_CYTA_SW_KEY = "cytassistInstrumentSoftwareVersion"
TIFFER_INFO_FORMAT_KEY = "format"
TIFFER_INFO_IMAGE_DESCRIPTION_KEY = "imageDescription"
TIFFER_INFO_PIXEL_DEPTH_KEY = "depth"

TIFFER_RESAMPLE_SCALE_KEY = "scale"

TIFFER_JPEG_VALUE = "jpg"
TIFFER_GRAYSCALE_VALUE = "gray"

# Find lib/bin/tiffer by relative path.
_TIFFER_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
    "bin",
    "tiffer",
)


def try_call_tiffer_info(image_path):
    """Run tiffer info to get all information about a file and throw an exception if a failure occurs."""
    return call_tiffer_info(image_path, use_raise=True)


def call_tiffer_info(image_path, use_raise=False):
    """Run tiffer info to get all information about a file.

    if "use_raise" is true, throw an exception instead of calling martian.exit
    """
    call = [_TIFFER_PATH, "info", image_path]
    with tk_subproc.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
        (proc_stdout, proc_stderr) = proc.communicate(None)
        if proc.returncode != 0:
            msgs = [
                f"Error reading image to get image metadata from {image_path}:",
                proc_stdout.decode(),
                proc_stderr.decode(),
            ]
            if use_raise:
                raise RuntimeError("\n".join([x for x in msgs if x]))
            else:
                martian.exit("\n".join([x for x in msgs if x]))
    imageinfo = json.loads(proc_stdout)
    return imageinfo


def call_tiffer_checksum(image_path):
    """Runs "tiffer checksum" to get a checksum for an image file."""
    call = [_TIFFER_PATH, "checksum", image_path]
    proc = tk_subproc.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (proc_stdout, proc_stderr) = proc.communicate(None)
    if proc.returncode != 0:
        msgs = [
            f"Error getting checksum of file {image_path}:",
            proc_stdout.decode(),
            proc_stderr.decode(),
        ]
        martian.exit("\n".join([x for x in msgs if x]))
    return proc_stdout.strip().decode()


def try_call_tiffer_compatibilty_fixes(image_path, out_dir) -> str:
    """Perform channel correction on a cytassist image."""
    os.makedirs(out_dir, exist_ok=True)
    call = [_TIFFER_PATH, "compatibilityfixes", image_path, out_dir]
    with tk_subproc.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
        (proc_stdout, proc_stderr) = proc.communicate(None)
        decoded_proc_stdout = proc_stdout.decode().strip()
        if proc.returncode != 0:
            msgs = [
                f"Image correction was unsuccessful for file {image_path}. If the problem persists, contact support@10xgenomics.com.",
                decoded_proc_stdout,
                proc_stderr.decode(),
            ]
            raise RuntimeError("\n".join(msgs))
        elif not os.path.exists(decoded_proc_stdout):
            msgs = [
                f"Failed to produce corrected image for file {image_path}. If the problem persists, contact support@10xgenomics.com.",
                decoded_proc_stdout,
                proc_stderr.decode(),
            ]
            raise RuntimeError("\n".join(msgs))
    return decoded_proc_stdout


def call_tiffer_get_num_pages(image_path):
    """Run tiffer pages to get the number of pages in the image."""
    call = [_TIFFER_PATH, "pages", image_path]
    proc = tk_subproc.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (proc_stdout, proc_stderr) = proc.communicate(None)
    if proc.returncode != 0:
        msgs = [
            f"Error reading image to get number of channels/pages for {image_path}:",
            proc_stdout.decode(),
            proc_stderr.decode(),
        ]
        martian.exit("\n".join([x for x in msgs if x]))
    npages = int(proc_stdout.strip())
    return npages


def call_tiffer_mem_estimate(image_path, downsample_size):
    """Run tiffer memreq with given parameters."""
    call = [_TIFFER_PATH, "memreq", image_path, "--size", str(downsample_size)]
    proc = tk_subproc.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (proc_stdout, proc_stderr) = proc.communicate(None)
    if proc.returncode != 0:
        msgs = [
            f"Error reading image file {image_path} while estimating memory requirements:",
            proc_stdout.decode(),
            proc_stderr.decode(),
        ]
        martian.exit("\n".join([x for x in msgs if x]))
    output_bytes = int(proc_stdout.strip())
    return output_bytes


def call_tiffer_mem_estimate_gb(image_path, downsample_size):
    """Call_tiffer_mem_estimate but in GB."""
    output_bytes = call_tiffer_mem_estimate(image_path, downsample_size)
    bytes_gb = np.ceil(output_bytes / (1024.0 * 1024.0 * 1024.0))
    return bytes_gb


def call_tiffer_resample(image_path, downsample_size, full_output_path, page=None):
    """Run tiffer with given parameters."""
    full_output_path = os.path.abspath(full_output_path)
    output_path = os.path.dirname(full_output_path)
    output_file_name = os.path.basename(full_output_path)
    (output_file_name_noext, ext) = os.path.splitext(output_file_name)
    ext = ext[1:]

    call = [
        _TIFFER_PATH,
        "resample",
        image_path,
        output_path,
        "--jsonfile",
        "--size",
        str(downsample_size),
        "--format",
        ext,
        "--outname",
        output_file_name_noext,
    ]

    if page is not None:
        call.extend(["--page", str(page)])

    try:
        _ = tk_subproc.check_output(call, stderr=subprocess.STDOUT)
        # ^ the unused returned variable is warning/error messages.
        output_json_wpath = os.path.join(output_path, f"{output_file_name_noext}.json")
        output_json = json.load(open(output_json_wpath))
    except subprocess.CalledProcessError as ex:
        martian.exit(f"Could not generate downsampled image: \n{ex.output.decode()}")
    return output_json


def get_max_image_dimension(img_path) -> int:
    """Get the maximum dimension of an image."""
    image_info = call_tiffer_info(img_path)
    pages = image_info[TIFFER_INFO_PAGE_KEY]
    return max(max(page[TIFFER_INFO_HEIGHT_KEY], page[TIFFER_INFO_WIDTH_KEY]) for page in pages)


def call_tiffer_normalize(image_path, full_output_path, crop_box: BoundingBox | None = None):
    """Call tiffer normalize."""
    full_output_path = os.path.abspath(full_output_path)
    output_directory_path = os.path.dirname(full_output_path)
    image_file_name = os.path.basename(image_path)
    tiffer_output_file_name = os.path.join(output_directory_path, image_file_name)

    call = [
        _TIFFER_PATH,
        "normalize",
        image_path,
        output_directory_path,
    ]

    if crop_box:
        call.extend([f"--bounding-box={','.join(str(x) for x in crop_box.get_array())}"])

    print(f"command run: {' '.join(call)}")

    try:
        _ = tk_subproc.check_output(call, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as ex:
        martian.exit(f"Could not generate normalize image: \n{ex.output.decode()}")

    cr_io.hardlink_with_fallback(tiffer_output_file_name, full_output_path)
