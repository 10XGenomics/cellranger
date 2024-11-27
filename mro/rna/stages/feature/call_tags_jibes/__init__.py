#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved
#
"""Assigns tags to cells."""

import pickle

import martian
import numpy as np
import pandas as pd

import cellranger.feature.multiplexing.metrics_utils as multi_metrics_utils
import cellranger.feature.utils as feature_utils
import cellranger.feature_ref as feature_ref
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.feature import jibes_tag_assigner as cr_jta
from cellranger.feature.feature_assigner import (
    CellEstimateCantConvergeException,
    calculate_fat_tail_frac,
)
from cellranger.feature.feature_assignments import CMO_CATEGORY_NAMES, CellsPerFeature

# pylint: disable=wrong-import-position, too-many-locals
from cellranger.pandas_utils import sanitize_dataframe

__MRO__ = """
stage CALL_TAGS_JIBES(
    in  csv    marginal_tag_calls_per_cell,
    in  csv    marginal_tag_frequencies,
    in  csv    filtered_barcodes,
    in  h5     filtered_feature_counts_matrix,
    in  h5     molecule_info,
    in  string throughput,
    in  string library_type,
    in  string multiplexing_method,
    in  float  min_assignment_confidence,
    out json   jibes_parameters,
    out csv    jibes_model_summary,
    out json   jibes_summary_data,
    out csv    assignment_confidence_table,
    out csv    tag_calls_summary,
    out csv    tag_calls_per_cell,
    out json   tag_call_metrics,
    out json   cells_per_tag,
    out json   tag_umi_thresholds_json,
    out csv    tag_umi_thresholds_csv,
    out pickle tag_assigner_pickle,
    out json   non_singlet_barcodes,
    src py     "stages/feature/call_tags_jibes",
) split (
) using (
    volatile = strict,
)
"""

np.random.seed(0)
MIN_MEM_GB = 6


def _add_empty_tag_assignments(cells_per_feature_object, all_tags_list):
    """Ensure that tags without assignments have an entry in the.

    cells_per_feature object with value []
    """
    tags_with_assignments = set(cells_per_feature_object.keys())
    tags_without_assignments = set(all_tags_list).difference(tags_with_assignments)

    for tag in tags_without_assignments:
        cells_per_feature_object[tag] = []

    return cells_per_feature_object


def set_empty(outs, all_tags):
    """Exit with empty functional outs."""
    outs.jibes_parameters = None
    outs.jibes_model_summary = None
    outs.jibes_summary_data = None
    outs.assignment_confidence_table = None
    outs.tag_calls_per_cell = None
    outs.tag_calls_summary = None

    if all_tags == []:
        outs.cells_per_tag = None
    else:
        empty_assignments = CellsPerFeature()
        empty_assignments = _add_empty_tag_assignments(empty_assignments, all_tags)
        empty_assignments.save_to_file(outs.cells_per_tag)
    outs.tag_umi_thresholds_json = None
    outs.tag_umi_thresholds_csv = None
    outs.tag_call_metrics = None
    outs.tag_assigner_pickle = None
    outs.non_singlet_barcodes = None


def split(args):
    return {"chunks": [], "join": {"__mem_gb": 16, "__threads": 4}}


def _get_input_files(args):
    return [
        args.marginal_tag_calls_per_cell,
        args.filtered_feature_counts_matrix,
        args.filtered_barcodes,
        args.molecule_info,
    ]


def join(args, outs, chunk_defs, chunk_outs):  # pylint: disable=too-many-locals,too-many-statements
    matrix_present = feature_utils.all_files_present([args.filtered_feature_counts_matrix])
    if not matrix_present:
        set_empty(outs=outs, all_tags=[])
        return

    filtered_feature_counts_matrix = cr_matrix.CountMatrix.load_h5_file(
        args.filtered_feature_counts_matrix
    )
    # in case of vdj cells need to subset barcodes in the filtered_matrix
    filtered_feature_counts_matrix = filtered_feature_counts_matrix.select_barcodes_by_seq(
        feature_utils.get_gex_cell_list(args.filtered_barcodes)
    )

    # exactly one of library_type or multiplexing_method is required to be set
    assert (args.library_type is None) ^ (args.multiplexing_method is None)
    if args.multiplexing_method is not None:
        multiplexing_method = rna_library.BarcodeMultiplexingType(args.multiplexing_method)
        assert multiplexing_method.is_cell_multiplexed(), "Unsupported multiplexing method!"
        library_type = multiplexing_method.multiplexing_library_type()
        if multiplexing_method.type == rna_library.CellLevel.Hashtag:
            filtered_tag_counts_matrix = (
                filtered_feature_counts_matrix.select_features_by_type_and_tag(
                    library_type, feature_ref.HASHTAG_TAG
                )
            )
        else:
            filtered_tag_counts_matrix = filtered_feature_counts_matrix.select_features_by_type(
                library_type
            )
    else:
        filtered_tag_counts_matrix = filtered_feature_counts_matrix.select_features_by_type(
            args.library_type
        )

    if filtered_tag_counts_matrix is None:
        set_empty(outs=outs, all_tags=[])
        return

    all_tags = [x.id for x in filtered_tag_counts_matrix.feature_ref.feature_defs]
    if feature_utils.check_if_none_or_empty(filtered_tag_counts_matrix):
        set_empty(outs=outs, all_tags=[])
        return

    input_files = _get_input_files(args)
    all_inputs_present = feature_utils.all_files_present(input_files)
    no_cells_found = len(feature_utils.get_gex_cell_list(args.filtered_barcodes)) == 0

    if not all_inputs_present or no_cells_found:
        set_empty(outs=outs, all_tags=[])
        return

    tag_call_metrics = {}

    # Tag calling
    try:
        tag_assigner = cr_jta.run_assignment_stage(
            filtered_tag_counts_matrix,
            args.marginal_tag_calls_per_cell,
            library_type,
            args.throughput,
            confidence=args.min_assignment_confidence,
        )
    except CellEstimateCantConvergeException as conv_err:
        martian.exit(f"{conv_err}")
    except cr_jta.JibesError as err:
        # this is a PD code path
        if library_type == rna_library.ANTIGEN_LIBRARY_TYPE:
            set_empty(outs=outs, all_tags=[])
            return
        else:
            martian.exit(f"Deplex Error: {err}")

    cr_jta.make_parameter_table_rows(tag_assigner.fitter, outs.jibes_parameters)

    sanitize_dataframe(tag_assigner.jibes_assignments).to_csv(outs.assignment_confidence_table)
    tag_call_metrics.update(tag_assigner.fitter.get_snr_dictionary())

    # TODO: Avoid pickle here, for now make tag_assigner.fitter a much smaller object
    tag_assigner.prepare_for_pickle()
    with open(outs.tag_assigner_pickle, "wb") as f:
        pickle.dump(tag_assigner, f)

    if tag_assigner.assignments is None or len(tag_assigner.assignments) == 0:
        set_empty(outs=outs, all_tags=all_tags)
        return
    tag_assigner.assignment_metadata = tag_assigner.compute_assignment_metadata()
    tag_assigner.feature_calls_summary = tag_assigner.get_feature_calls_summary()
    tag_assigner.cells_per_feature = tag_assigner.get_cells_per_feature()

    tag_assignments_matrix = tag_assigner.create_tag_assignments_matrix()
    features_per_cell_table = tag_assignments_matrix.get_features_per_cell_table()
    cells_per_feature = tag_assignments_matrix.get_cells_per_feature()
    feature_calls_summary = tag_assignments_matrix.get_feature_calls_summary(
        cat_names=CMO_CATEGORY_NAMES
    )

    # write out non-tag assignments so they can be used later while making plots
    # without having to use the tag assigner pickle.
    non_singlet_barcodes = CellsPerFeature()
    for assignment, barcode_indices in tag_assigner.non_singlet_barcodes.items():
        non_singlet_barcodes[assignment.encode()] = filtered_feature_counts_matrix.ints_to_bcs(
            barcode_indices
        )
    non_singlet_barcodes.save_to_file(outs.non_singlet_barcodes)

    # Save memory
    del tag_assigner.sub_matrix
    del filtered_feature_counts_matrix
    tag_call_metrics.update(tag_assigner.get_feature_assignment_metrics())

    sanitize_dataframe(features_per_cell_table, inplace=True).to_csv(outs.tag_calls_per_cell)
    sanitize_dataframe(feature_calls_summary, inplace=True).to_csv(outs.tag_calls_summary)

    # Ensure a [] entry in cells_per_feature for tags without assignments
    cells_per_feature = _add_empty_tag_assignments(cells_per_feature, all_tags)
    cells_per_feature.save_to_file(outs.cells_per_tag)
    sanitize_dataframe(tag_assigner.assignment_metadata.umi_thresholds, inplace=True).to_csv(
        outs.tag_umi_thresholds_csv, index_label="tag"
    )
    tag_assigner.assignment_metadata.umi_thresholds.to_json(
        path_or_buf=outs.tag_umi_thresholds_json, orient="index"
    )

    # calculate fraction of reads from contaminant tags
    report_prefix = rna_library.get_library_type_metric_prefix(
        rna_library.MULTIPLEXING_LIBRARY_TYPE
    )
    name_fr_contaminants = report_prefix + "observed_frac_contaminant_tag"
    tag_call_metrics[name_fr_contaminants] = multi_metrics_utils.get_frac_contaminant_tags(
        tag_assigner.assignments, args.molecule_info
    )

    # calculate fat_tail_frac
    if args.marginal_tag_frequencies:
        marginal_freqs_df = pd.read_csv(args.marginal_tag_frequencies, index_col=0)
        tag_call_metrics.update(
            calculate_fat_tail_frac(
                freqs_df=marginal_freqs_df,
                report_prefix=report_prefix,
            )
        )

    # add a "_jibes" suffix for the JIBES metrics
    tag_call_metrics = multi_metrics_utils.add_jibes_suffix(tag_call_metrics)

    # fix tag_umi_thresholds_csv for bytes
    umi_thresholds = sanitize_dataframe(
        tag_assigner.assignment_metadata.umi_thresholds, inplace=True
    )
    umi_thresholds.to_csv(outs.tag_umi_thresholds_csv, index_label="tag")
    tag_assigner.assignment_metadata.umi_thresholds.to_json(
        path_or_buf=outs.tag_umi_thresholds_json, orient="index"
    )
    feature_utils.write_json_from_dict(tag_call_metrics, outs.tag_call_metrics)

    return
