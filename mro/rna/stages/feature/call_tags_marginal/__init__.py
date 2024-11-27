#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved
#
"""Assigns tags to cells."""

import martian
import numpy as np

import cellranger.feature.utils as feature_utils
import cellranger.feature_ref as feature_ref
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.feature.feature_assigner import CellEstimateCantConvergeException, TagAssigner
from cellranger.feature.throughputs import G19_N_GEMS
from cellranger.pandas_utils import sanitize_dataframe

__MRO__ = """
stage CALL_TAGS_MARGINAL(
    in  csv    filtered_barcodes,
    in  h5     filtered_feature_counts_matrix,
    in  string throughput,
    in  string multiplexing_method,
    in  string library_type,
    out csv    marginal_tag_calls_per_cell,
    out csv    marginal_tag_frequencies,
    out json   tag_contaminant_info,
    src py     "stages/feature/call_tags_marginal",
) split (
)
"""

np.random.seed(0)
MIN_MEM_GB = 6


def set_empty(outs):
    """Exit with empty functional outs."""
    outs.marginal_tag_calls_per_cell = None
    outs.tag_contaminant_info = None
    outs.marginal_tag_frequencies = None


def split(args):
    return {"chunks": [], "join": {"__mem_gb": 16}}


def _get_input_files(args):
    return [
        args.filtered_feature_counts_matrix,
        args.filtered_barcodes,
    ]


def join(args, outs, chunk_defs, chunk_outs):
    input_files = _get_input_files(args)
    all_inputs_present = feature_utils.all_files_present(input_files)
    no_cells_found = len(feature_utils.get_gex_cell_list(args.filtered_barcodes)) == 0

    if not all_inputs_present or no_cells_found:
        set_empty(outs)
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

    if feature_utils.check_if_none_or_empty(filtered_tag_counts_matrix):
        set_empty(outs)
        return

    # Tag calling
    n_gems = G19_N_GEMS[args.throughput]
    tag_assigner = TagAssigner(
        matrix=filtered_tag_counts_matrix,
        feature_type=library_type,
        n_gems=n_gems,
    )
    tag_assigner.assignments = tag_assigner.get_feature_assignments()
    if tag_assigner.assignments is None or len(tag_assigner.assignments) == 0:
        set_empty(outs)
        return
    try:
        tag_assigner.assignment_metadata = tag_assigner.compute_assignment_metadata()
    except CellEstimateCantConvergeException as ex:
        martian.exit(f"{ex}")

    tag_assigner.feature_calls_summary = tag_assigner.get_feature_calls_summary()
    tag_assigner.cells_per_feature = tag_assigner.get_cells_per_feature()

    tag_assignments_matrix = tag_assigner.create_tag_assignments_matrix()
    features_per_cell_table = tag_assignments_matrix.get_features_per_cell_table()

    # Save memory
    del tag_assigner.sub_matrix
    del filtered_feature_counts_matrix
    sanitize_dataframe(features_per_cell_table, inplace=True).to_csv(
        outs.marginal_tag_calls_per_cell
    )

    tag_assigner.assignment_metadata.umi_thresholds.T.to_json(outs.tag_contaminant_info)

    # Diagnostic counts per k-let for microfluidics (marginal caller only)
    tag_assignment_freqs_df = tag_assigner.assignment_metadata.freq_counts
    sanitize_dataframe(tag_assignment_freqs_df, inplace=True).to_csv(
        outs.marginal_tag_frequencies, index_label="n_Tags"
    )
    return
