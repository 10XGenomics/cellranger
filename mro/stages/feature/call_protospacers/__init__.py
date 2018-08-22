#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#
import cellranger.feature.crispr.analysis as crispr_analysis
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
    out json   umi_thresholds_json,
    out csv    umi_thresholds_csv,
    src py     "stages/feature/call_protospacers",
) using (
    mem_gb = 6,
)
"""

def main(args, outs):
    if not(os.path.isfile(args.molecule_info) and os.path.isfile(args.filtered_feature_counts_matrix)):
        (outs.protospacer_calls_summary, outs.protospacer_calls_per_cell,
            outs.protospacer_call_metrics_json) = (None, None, None, None)
        return
    protospacer_call_metrics = feature_utils.translate_feature_metrics_json_keys(args.counter_metrics_json, args.feature_type)
    report_prefix = feature_constants.PREFIX_FROM_FEATURE_TYPE.get(args.feature_type, 'FEATURE') + '_'

    filtered_feature_counts_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_feature_counts_matrix)
    filtered_guide_counts_matrix = filtered_feature_counts_matrix.select_features_by_type(rna_library.CRISPR_LIBRARY_TYPE)
    num_gex_cbs = len(filtered_feature_counts_matrix.bcs)

    if filtered_guide_counts_matrix is None:
        return

    feature_defs = filtered_guide_counts_matrix.feature_ref.feature_defs
    feature_map = {feature_def.id:feature_def.tags.get('sequence') for feature_def in feature_defs}

    """Protospacer calling"""
    (perturbation_calls_table, presence_calls,
        cells_with_ps, ps_calls_summary, umi_thresholds)  = crispr_analysis.get_ps_calls_and_summary(filtered_guide_counts_matrix,
                                                                                        feature_map,)
    protospacer_call_metrics.update(crispr_analysis.get_protospacer_call_metrics(ps_calls_summary, num_gex_cbs, report_prefix))

    perturbation_calls_table.to_csv(outs.protospacer_calls_per_cell)
    ps_calls_summary.to_csv(outs.protospacer_calls_summary)
    feature_utils.write_json_from_dict(cells_with_ps, outs.cells_per_protospacer)
    feature_utils.write_json_from_dict(umi_thresholds, outs.umi_thresholds_json)
    feature_utils.write_csv_from_dict(umi_thresholds, outs.umi_thresholds_csv, "Protospacer, UMI threshold\n")
    feature_utils.write_json_from_dict(protospacer_call_metrics, outs.protospacer_call_metrics_json)








