#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Functions and classes for loading and manipulating images."""

from __future__ import annotations

import base64
import io
import os
import tempfile

from PIL import Image, ImageOps


def _base64_encode_image(filename, fmt="jpeg"):
    """Opens a file using PIL and returns it encoded as a base64 string.

    :param filename:
    :return:
    """
    with open(filename, "rb") as image_file:
        image_show = image_file.read()
    encoded_string = base64.b64encode(image_show).decode("utf-8")
    return f"data:image/{fmt};base64," + encoded_string


def base64_encode_image(fname):
    """base64 encode image."""
    _, fmt = os.path.splitext(fname)
    fmt = fmt.removeprefix(".")
    if fmt == "jpg":
        fmt = "jpeg"
    return _base64_encode_image(fname, fmt=fmt)


class WebImage:
    """A class for working with and investigating simple image files."""

    def __init__(self, filename: str, cropbox=None, markersize=None):
        """Create a new instance that can base64 encode the image and give coordinates.

        :param filename: Image file name
        :param cropbox: optional [ x0, y0, x1, y1 ] just held as an attribute,
        :               defaults to whole image
        :param markersize: optional marker size for plotly for plotting a capture area spot
        """
        self._base64 = _base64_encode_image(filename)
        self.filename = filename
        with Image.open(filename) as img:
            self.width, self.height = img.size
        self.cropbox = cropbox if cropbox is not None else [0, 0, self.width - 1, self.height - 1]
        self.markersize = markersize

    @property
    def base64_encoded_str(self):
        """:return: String for a web summary,i.e. "data:image/jpg;base64,..."."""
        return self._base64

    def base64_encoded_grayscale_image(self) -> str:
        """Get encoded string of grascaled image for web summary."""
        with Image.open(self.filename) as img:
            with io.BytesIO() as img_bytes:
                ImageOps.grayscale(img).save(img_bytes, format="PNG")
                image_show = img_bytes.getvalue()
        encoded_string = base64.b64encode(image_show).decode("utf-8")
        return f"data:image/png;base64,{encoded_string}"

    def resize_and_encode_image(self, new_width=None, new_height=None):
        """:param new_width: New image height.

        :param new_height: New image width
        :return: A base64 encoded string with the new image in it.
        """
        # TODO: We want to be able to encode this without saving to a file.
        if not new_width and not new_height:
            raise ValueError("Width and/or height must be set when resizing image.")
        elif not new_width:
            new_width = self.width * new_height // self.height
        elif not new_height:
            new_height = self.height * new_width // self.width

        _, fname = os.path.split(self.filename)
        tmp_file = os.path.join(tempfile.mkdtemp(), "tmp_" + fname)
        with Image.open(self.filename) as img:
            img2 = img.resize((new_width, new_height), Image.Resampling.LANCZOS)
        img2.save(tmp_file)
        img2.close()
        return WebImage(tmp_file)
