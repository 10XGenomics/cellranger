#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Definiton of the Metric class."""

from __future__ import annotations

REFERENCE_PATH = "reference_path"


class Metrics:
    def update(self, other):
        for k, v in other.__dict__.items():
            if v is not None:
                setattr(self, k, v)


class BarcodeFilterResults(Metrics):
    def __init__(self, default_value: int = 0):
        self.filtered_bcs = default_value
        self.filtered_bcs_lb = default_value
        self.filtered_bcs_ub = default_value
        self.filtered_bcs_var = default_value
        self.filtered_bcs_cv = float(default_value)
        self.filtered_bcs_cutoff = default_value

    @staticmethod
    def init_with_constant_call(n_bcs: int) -> BarcodeFilterResults:
        res = BarcodeFilterResults()
        res.filtered_bcs = n_bcs
        res.filtered_bcs_lb = n_bcs
        res.filtered_bcs_ub = n_bcs
        res.filtered_bcs_var = 0
        res.filtered_bcs_cv = 0
        res.filtered_bcs_cutoff = 0
        return res

    def to_dict_with_prefix(
        self, i: int, sample: str | None, method: str
    ) -> dict[str, int | float]:
        sample_prefix = "_" + sample if sample else ""
        return {
            "gem_group_%d%s_%s_%s" % (i, sample_prefix, key, method): value
            for key, value in self.__dict__.items()
        }
