#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#
import cellranger.feature.crispr.protospacer_calling as protospacer_calling
import json
import os
import cellranger.feature.utils as feature_utils
import cellranger.rna.library as rna_library
import cellranger.matrix as cr_matrix
import cellranger.feature.constants as feature_constants

__MRO__ = """
stage CALL_PROTOSPACERS(
    in  csv    filtered_barcodes,
    in  h5     molecule_info,
    in  h5     filtered_feature_counts_matrix,
    in  string feature_type,
    in  csv    feature_reference,
    in  json   counter_metrics_json,
    out csv    protospacer_calls_summary,
    out csv    protospacer_calls_per_cell,
    out json   protospacer_call_metrics_json,
    out json   cells_per_protospacer,
    out json   protospacer_umi_thresholds_json,
    out csv    protospacer_umi_thresholds_csv,
    src py     "stages/feature/call_protospacers",
) using (
    mem_gb = 6,
)
"""

def set_empty(outs):
    outs.protospacer_calls_per_cell = None
    outs.protospacer_calls_summary = None
    outs.cells_per_protospacer = None
    outs.protospacer_umi_thresholds_json = None
    outs.protospacer_umi_thresholds_csv = None
    outs.protospacer_call_metrics_json = None

def main(args, outs):

    if not(os.path.isfile(args.molecule_info) and os.path.isfile(args.filtered_feature_counts_matrix)):
        set_empty(outs)
        return

    with open(args.counter_metrics_json) as f:
        protospacer_call_metrics = json.load(f)
    report_prefix = feature_constants.PREFIX_FROM_FEATURE_TYPE.get(args.feature_type, 'FEATURE') + '_'

    filtered_feature_counts_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_feature_counts_matrix)
    filtered_guide_counts_matrix = filtered_feature_counts_matrix.select_features_by_type(rna_library.CRISPR_LIBRARY_TYPE)
    num_gex_cbs = len(filtered_feature_counts_matrix.bcs)

    if feature_utils.check_if_none_or_empty(filtered_guide_counts_matrix):
        set_empty(outs)
        return

    feature_defs = filtered_guide_counts_matrix.feature_ref.feature_defs
    feature_map = {feature_def.id:feature_def.tags.get('sequence') for feature_def in feature_defs}

    """Protospacer calling"""
    (perturbation_calls_table, presence_calls,
        cells_with_ps, ps_calls_summary, umi_thresholds)  = protospacer_calling.get_ps_calls_and_summary(filtered_guide_counts_matrix,
                                                                                        feature_map,)
    protospacer_call_metrics.update(protospacer_calling.get_protospacer_call_metrics(ps_calls_summary, num_gex_cbs, report_prefix))

    perturbation_calls_table.to_csv(outs.protospacer_calls_per_cell)
    ps_calls_summary.to_csv(outs.protospacer_calls_summary)
    feature_utils.write_json_from_dict(cells_with_ps, outs.cells_per_protospacer)
    feature_utils.write_json_from_dict(umi_thresholds, outs.protospacer_umi_thresholds_json)
    feature_utils.write_csv_from_dict(umi_thresholds, outs.protospacer_umi_thresholds_csv, "Protospacer, UMI threshold\n")
    feature_utils.write_json_from_dict(protospacer_call_metrics, outs.protospacer_call_metrics_json)
