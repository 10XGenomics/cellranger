#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#
import cellranger.io as cr_io
import os
import itertools
import cellranger.feature.crispr.protospacer_calling as protospacer_calling

__MRO__ = """
stage SUMMARIZE_CRISPR_ANALYSIS(
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
    mem_gb = 4,
)
"""
def main(args, outs):
    list_of_files = [args.protospacer_calls_summary, args.protospacer_calls_per_cell, args.cells_per_protospacer,
                     args.protospacer_umi_thresholds_csv, args.protospacer_umi_thresholds_json,
                     args.perturbation_efficiencies_by_feature, args.perturbations_efficiencies_by_target]

    cr_io.makedirs(outs.crispr_analysis, allow_existing=True)

    for (file_path, file_name) in itertools.izip(list_of_files, protospacer_calling.CRISPR_ANALYSIS_FILE_NAMES):
        if file_path is None:
            continue
        cr_io.copy(file_path, os.path.join(outs.crispr_analysis, file_name))

    if os.path.isdir(args.perturbation_effects_by_feature):
        perturbation_effects_by_feature_dir = os.path.join(outs.crispr_analysis, 'perturbation_effects_by_feature')
        cr_io.makedirs(perturbation_effects_by_feature_dir, allow_existing=True)
        cr_io.copytree(args.perturbation_effects_by_feature, perturbation_effects_by_feature_dir, allow_existing=True)

    if os.path.isdir(args.perturbation_effects_by_target):
        perturbation_effects_by_target_dir = os.path.join(outs.crispr_analysis, 'perturbation_effects_by_target')
        cr_io.makedirs(perturbation_effects_by_target_dir, allow_existing=True)
        cr_io.copytree(args.perturbation_effects_by_target, perturbation_effects_by_target_dir, allow_existing=True)
