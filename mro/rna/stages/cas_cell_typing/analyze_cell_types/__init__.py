#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Analyze cell types generated from CALL_CELL_TYPES."""

import json
import os

import h5py
import martian

import cellranger.cell_typing.broad_tenx.cas_metrics as cas_metrics
import cellranger.cell_typing.common_cell_typing as common_cell_typing
import cellranger.feature.utils as feature_utils
import tenkit.safe_json as tk_safe_json
from cellranger.cell_typing.common_cell_typing import (
    COARSE_CELL_TYPES_KEY,
    FINE_CELL_TYPES_KEY,
)

__MRO__ = """
stage ANALYZE_CELL_TYPES(
    in  csv    cell_types,
    in  path   analysis,
    in  csv    tissue_positions,
    in  png    tissue_lowres_image,
    in  json   scale_factors,
    in  float  cas_frac_returned_bcs     "Fraction of barcodes that went into cell annotation that were returned",
    in  string cas_model_used,
    in  string pipestance_type,
    in  bool   cas_success,
    in  json   metadata,
    in  bool   is_pd,
    out json   cell_annotation_analysis_metrics,
    out json   cell_type_bar_chart,
    out json   spatial_cell_types_chart,
    src py     "stages/cas_cell_typing/analyze_cell_types",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    if not args.cell_types:
        martian.clear(outs)
        return
    cas_mets = {}
    metadata_dict = {}
    if args.analysis and "/clustering/_gene_expression_graphclust/clusters" in h5py.File(
        os.path.join(args.analysis, "analysis.h5"), "r"
    ):
        if args.metadata:
            with open(args.metadata) as f:
                metadata_dict = json.load(f)
        if cas_metrics.analysis_has_clustering(args.analysis):
            cas_mets = cas_metrics.get_cas_cluster_purity(
                cell_types=args.cell_types, analysis=args.analysis
            )
        # Make some diagnostic plots
        ## Stacked barchart of cluster purity
        if args.is_pd:
            cell_type_bar_chart = cas_metrics.cell_type_bar_chart(
                cell_types=args.cell_types, analysis=args.analysis, stacked=True
            )
            feature_utils.write_json_from_dict(cell_type_bar_chart, outs.cell_type_bar_chart)
        else:
            outs.cell_type_bar_chart = None
    else:
        outs.cell_type_bar_chart = None

    if args.cas_model_used and args.pipestance_type is not None and args.cas_success is not None:
        prefix = args.cas_model_used.lower().replace(" ", "_")
        cas_coarse_fine_mets = cas_metrics.get_cell_type_summary_mets(cell_types=args.cell_types)
        if args.cas_frac_returned_bcs:
            cas_mets["cas_frac_returned_bcs"] = args.cas_frac_returned_bcs
        cas_mets["cell_annotation_model"] = args.cas_model_used
        cas_mets["pipestance_type"] = args.pipestance_type
        cas_mets["cas_success"] = args.cas_success
        cas_mets["tree_version_used"] = metadata_dict.get("tree_version_used")
        cas_mets["display_map_version_used"] = metadata_dict.get("display_map_version_used")
        cas_mets[f"{prefix}_num_coarse_cell_categories"] = cas_coarse_fine_mets.get(
            "num_coarse_cell_categories"
        )
        cas_mets[f"{prefix}_num_fine_cell_categories"] = cas_coarse_fine_mets.get(
            "num_fine_cell_categories"
        )

        cas_mets[f"{prefix}_top_coarse_cell_types"] = common_cell_typing.get_top_cell_types(
            args.cell_types, COARSE_CELL_TYPES_KEY, 10
        )
        cas_mets[f"{prefix}_top_fine_cell_types"] = common_cell_typing.get_top_cell_types(
            args.cell_types, FINE_CELL_TYPES_KEY, 10
        )
        cas_mets[f"{prefix}_low_umi_barcodes"] = common_cell_typing.get_cell_type_counts(
            args.cell_types, FINE_CELL_TYPES_KEY, "Low UMI Barcode"
        )

        with open(martian.make_path(outs.cell_annotation_analysis_metrics), "w") as f:
            tk_safe_json.dump_numpy(
                cas_mets,
                f,
                indent=4,
                sort_keys=False,
                separators=(",", ": "),
            )
    else:
        outs.cell_annotation_analysis_metrics = None

    ## Spatial plot of cell types
    if (
        args.tissue_positions
        and args.scale_factors
        and args.tissue_lowres_image
        and args.cell_types
    ):
        spatial_plot = cas_metrics.spatial_cell_types_plot(
            tissue_positions_path=args.tissue_positions,
            tissue_lowres_image=args.tissue_lowres_image,
            cell_types=args.cell_types,
            analysis=args.analysis,
            scale_factors=args.scale_factors,
        )
        feature_utils.write_json_from_dict(spatial_plot, outs.spatial_cell_types_chart)
    else:
        outs.spatial_cell_types_chart = None
