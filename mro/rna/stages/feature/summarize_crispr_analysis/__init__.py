#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import martian

import cellranger.cr_io as cr_io
import cellranger.feature.crispr.protospacer_calling as protospacer_calling

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage SUMMARIZE_CRISPR_ANALYSIS(
    in  csv  feature_reference,
    in  csv  protospacer_calls_summary,
    in  csv  protospacer_calls_per_cell,
    in  json cells_per_protospacer,
    in  csv  protospacer_umi_thresholds_csv,
    in  json protospacer_umi_thresholds_json,
    in  csv  perturbation_efficiencies_by_feature,
    in  csv  perturbations_efficiencies_by_target,
    in  path perturbation_effects_by_feature,
    in  path perturbation_effects_by_target,
    out path crispr_analysis,
    src py   "stages/feature/summarize_crispr_analysis",
) using (
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_CRISPR_ANALYSIS stage inputs."""

    feature_reference: None | mro_filetypes.csv
    protospacer_calls_summary: None | mro_filetypes.csv
    protospacer_calls_per_cell: None | mro_filetypes.csv
    cells_per_protospacer: None | mro_filetypes.json
    protospacer_umi_thresholds_csv: None | mro_filetypes.csv
    protospacer_umi_thresholds_json: None | mro_filetypes.json
    perturbation_efficiencies_by_feature: None | mro_filetypes.csv
    perturbations_efficiencies_by_target: None | mro_filetypes.csv
    perturbation_effects_by_feature: None | mro_filetypes.path
    perturbation_effects_by_target: None | mro_filetypes.path


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_CRISPR_ANALYSIS stage outputs."""

    crispr_analysis: mro_filetypes.path


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_CRISPR_ANALYSIS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """SUMMARIZE_CRISPR_ANALYSIS chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    list_of_files = [
        args.feature_reference,
        args.protospacer_calls_summary,
        args.protospacer_calls_per_cell,
        args.cells_per_protospacer,
        args.protospacer_umi_thresholds_csv,
        args.protospacer_umi_thresholds_json,
        args.perturbation_efficiencies_by_feature,
        args.perturbations_efficiencies_by_target,
    ]

    os.makedirs(outs.crispr_analysis, exist_ok=True)

    for file_path, file_name in zip(list_of_files, protospacer_calling.CRISPR_ANALYSIS_FILE_NAMES):
        if file_path is None:
            continue
        cr_io.hardlink_with_fallback(file_path, os.path.join(outs.crispr_analysis, file_name))

    if args.perturbation_effects_by_feature and os.path.isdir(args.perturbation_effects_by_feature):
        perturbation_effects_by_feature_dir = os.path.join(
            outs.crispr_analysis, "perturbation_effects_by_feature"
        )
        os.makedirs(perturbation_effects_by_feature_dir, exist_ok=True)
        cr_io.hardlink_with_fallback(
            args.perturbation_effects_by_feature,
            perturbation_effects_by_feature_dir,
        )

    if args.perturbation_effects_by_target and os.path.isdir(args.perturbation_effects_by_target):
        perturbation_effects_by_target_dir = os.path.join(
            outs.crispr_analysis, "perturbation_effects_by_target"
        )
        os.makedirs(perturbation_effects_by_target_dir, exist_ok=True)
        cr_io.hardlink_with_fallback(
            args.perturbation_effects_by_target,
            perturbation_effects_by_target_dir,
        )
