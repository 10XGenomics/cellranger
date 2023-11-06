#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Read environment variables."""

from __future__ import annotations

import os


def product():
    """Return the value of the environment variable TENX_PRODUCT if it exists.

    or "tenxranger" otherwise.
    """
    return os.getenv("TENX_PRODUCT", "tenxranger")
