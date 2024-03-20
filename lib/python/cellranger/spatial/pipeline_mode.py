#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Define the pipeline mode for spatial imaging subpipeline."""


from __future__ import annotations

from enum import Enum
from typing import NamedTuple


class Product(str, Enum):
    VISIUM = "Visium"
    CYT = "CytAssist"
    VISIUM_HD_NOCYT_PD = "Visium-HD no CytAssist image"


class SlideType(str, Enum):
    VISIUM = "Visium-Slide"
    XL = "XL-Slide"
    VISIUM_HD = "Visium-HD-Slide"


class PipelineMode(NamedTuple):
    """Pipeline mode used for the spatial imaging subpipeline."""

    product: Product
    slide: SlideType

    def validate(self):
        """Validate the pipeline mode is valid."""
        # TODO: (dongyao) once we finalize the utility of the PipelineMode,
        # add check to prevent illegal product-slide combination

        try:
            _ = Product(self.product)
        except Exception as err:
            raise ValueError(
                f"invalid product '{self.product!s}' of type {type(self.product)}"
            ) from err

        try:
            _ = SlideType(self.slide)
        except Exception as err:
            raise ValueError(
                f"invalid slide type '{self.slide!s}' of type {type(self.slide)}"
            ) from err

    def is_visium_hd_with_fiducials(self) -> bool:
        return self.slide == SlideType.VISIUM_HD and self.product == Product.CYT

    def is_cytassist(self) -> bool:
        return self.product == Product.CYT

    def is_visium_hd(self) -> bool:
        return self.slide == SlideType.VISIUM_HD
