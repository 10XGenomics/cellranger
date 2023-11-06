#!/usr/bin/env python3
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Conditionally include stage code."""

try:
    from .internal import main
except ImportError:
    from .stub import main

__all__ = ["main"]
