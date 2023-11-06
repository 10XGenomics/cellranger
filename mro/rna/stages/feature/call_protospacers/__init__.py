#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#
"""Assign protospacers to cells i.e.identify which cells express.

which protospacers above background
"""

import json

import numpy as np
from six import ensure_str

import cellranger.feature.utils as feature_utils
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.feature.feature_assigner import GuideAssigner
from cellranger.pandas_utils import sanitize_dataframe

__MRO__ = """
stage CALL_PROTOSPACERS(
    in  csv    filtered_barcodes,
    in  h5     filtered_feature_counts_matrix,
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
    """Set outs to empty."""
    outs.protospacer_calls_per_cell = None
    outs.protospacer_calls_summary = None
    outs.cells_per_protospacer = None
    outs.protospacer_umi_thresholds_json = None
    outs.protospacer_umi_thresholds_csv = None
    outs.protospacer_call_metrics_json = None


def write_csv_from_dict(input_dict, out_file_name, header=None):
    with open(out_file_name, "w") as f:
        if header is not None:
            f.write(header)
        for key, value in input_dict.items():
            if isinstance(key, bytes):
                key = ensure_str(key)
            elif not isinstance(key, str):
                key = str(key)
            f.write(key)
            f.write(",")
            if isinstance(value, bytes):
                value = ensure_str(value)
            elif not isinstance(value, str):
                value = str(value)
            f.write(value)
            f.write("\n")


def main(args, outs):
    input_files = [args.filtered_feature_counts_matrix, args.filtered_barcodes]
    all_inputs_present = feature_utils.all_files_present(input_files)

    no_cells_found = len(feature_utils.get_gex_cell_list(args.filtered_barcodes)) == 0

    if not (all_inputs_present) or no_cells_found:
        set_empty(outs)
        return

    with open(args.counter_metrics_json) as f:
        protospacer_call_metrics = json.load(f)

    filtered_feature_counts_matrix = cr_matrix.CountMatrix.load_h5_file(
        args.filtered_feature_counts_matrix
    )
    filtered_guide_counts_matrix = filtered_feature_counts_matrix.select_features_by_type(
        rna_library.CRISPR_LIBRARY_TYPE
    )

    if feature_utils.check_if_none_or_empty(filtered_guide_counts_matrix):
        set_empty(outs)
        return

    # Protospacer calling
    guide_assigner = GuideAssigner(
        matrix=filtered_feature_counts_matrix,
        feature_type=rna_library.CRISPR_LIBRARY_TYPE,
    )

    guide_assigner.assignments = guide_assigner.get_feature_assignments()
    guide_assigner.features_per_cell_table = guide_assigner.get_features_per_cell_table()
    guide_assigner.assignment_metadata = guide_assigner.compute_assignment_metadata()
    guide_assigner.feature_calls_summary = guide_assigner.get_feature_calls_summary()
    guide_assigner.cells_per_feature = guide_assigner.get_cells_per_feature()

    guide_assignments_matrix = guide_assigner.create_guide_assignments_matrix()
    features_per_cell_table = guide_assignments_matrix.get_features_per_cell_table()
    cells_per_feature = guide_assignments_matrix.get_cells_per_feature()
    feature_calls_summary = guide_assignments_matrix.get_feature_calls_summary()
    # temporary assertions to make sure these outputs do not change
    # when we switch to using feature assignments matrix
    assert cells_per_feature == guide_assigner.cells_per_feature
    for bc, row in guide_assigner.features_per_cell_table.iterrows():
        assert (row == features_per_cell_table.loc[bc]).all()
    for bc, row in guide_assigner.feature_calls_summary.iloc[4:, :].iterrows():
        old_values = np.array(row, dtype="float")
        new_values = np.array(feature_calls_summary.loc[ensure_str(bc)], dtype="float")
        assert np.allclose(
            old_values, new_values, atol=0.05  # new_values is rounded to 1 decimal point
        )

    protospacer_call_metrics.update(guide_assigner.get_feature_assignment_metrics())

    sanitize_dataframe(guide_assigner.features_per_cell_table, inplace=True).to_csv(
        outs.protospacer_calls_per_cell
    )
    sanitize_dataframe(guide_assigner.feature_calls_summary, inplace=True).to_csv(
        outs.protospacer_calls_summary
    )
    feature_utils.write_json_from_dict(guide_assigner.cells_per_feature, outs.cells_per_protospacer)
    feature_utils.write_json_from_dict(
        guide_assigner.assignment_metadata.umi_thresholds, outs.protospacer_umi_thresholds_json
    )

    write_csv_from_dict(
        guide_assigner.assignment_metadata.umi_thresholds,
        outs.protospacer_umi_thresholds_csv,
        "Protospacer,UMI threshold\n",
    )

    feature_utils.write_json_from_dict(protospacer_call_metrics, outs.protospacer_call_metrics_json)
