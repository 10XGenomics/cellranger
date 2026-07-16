#!/usr/bin/env python
#
# Copyright (c) 2026 10x Genomics, Inc. All rights reserved.
#
"""Compute segmentation plots for Visium HD data."""


__MRO__ = """
stage STRUCTIFY_WEBSUMMARY_INPUTS(
    in  WebSummaryCellTypeInputs websummary_inputs,
    in  map<json>                cell_type_spatial_plot,
    out WebSummaryCellTypeInputs websummary_inputs_with_spatial_charts,
    src py                       "stages/cas_cell_typing/structify_websummary_inputs",
)
"""


def main(
    args,
    outs,
):
    if args.cell_type_spatial_plot:
        # Deep copy the input to avoid modifying the original
        websummary_inputs = args.websummary_inputs

        cell_type_spatial_plot = args.cell_type_spatial_plot

        # Get the bundles from the websummary_inputs
        bundles = websummary_inputs.get("cell_type_websummary_bundles", {})

        # Update each model family's spatial plot
        for model_family, file_path in cell_type_spatial_plot.items():
            if model_family in bundles:
                # Update the cell_type_spatial_plot field directly
                bundles[model_family]["cell_type_spatial_plot"] = file_path

        outs.websummary_inputs_with_spatial_charts = websummary_inputs
    else:
        outs.websummary_inputs_with_spatial_charts = args.websummary_inputs
