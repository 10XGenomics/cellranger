#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved
#
"""Assign antibodies to cells."""

import cellranger.feature.utils as feature_utils
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.fast_utils import MultiGraph

# from cellranger.feature.feature_assigner import AntibodyAssigner
from cellranger.websummary.histograms import make_antibody_histograms
from cellranger.websummary.treemaps import MIN_ANTIBODY_UMI, make_antibody_treemap_plot

# pylint: disable=invalid-name

__MRO__ = """
stage CALL_ANTIBODIES(
    in  h5   filtered_feature_counts_matrix,
    in  bool is_antibody,
    in  bool is_spatial,
    in  json multi_graph,
    in  string sample_id,
    out json antibody_histograms_json,
    out json antibody_treemap_json,
    src py   "stages/feature/call_antibodies",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    """TODO CELLRANGER-7261: Revise the memory reservation."""
    # Use number of cells and non-zero items in count matrix to predict memory usage
    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_feature_counts_matrix)
    num_cells = matrix_dims[1]
    nnz = matrix_dims[2]
    # This stage used to run a GMM model for single samples, and had a regression estimate then.  Now, it just makes
    # histograms and most of the memory is consumed by loading the matrix which can be meaningful for AGGR samples.
    # We'll use two methods of estimating, one based on the matrix the other based on the regression,
    # and pick the highest to avoid OOM errors.
    mem_gb_matrix = cr_matrix.CountMatrix.get_mem_gb_from_matrix_dim(num_cells, nnz)
    mem_gb_regression = 1.5 + 7e-09 * nnz + 9e-06 * num_cells
    mem_gb = max(mem_gb_matrix, mem_gb_regression)
    vmem_gb = mem_gb + 3.0

    return {
        "chunks": [],
        "join": {"__mem_gb": mem_gb, "__threads": 1, "__vmem_gb": vmem_gb},
    }


def main(args, _outs):
    pass


def join(args, outs, chunk_defs, chunk_outs):
    """This stage is the entry to the _ANTIBODY_ANALYZER pipeline. It calls the feature barcode.

    processing APIs such as AntibodyOrAntigenAssigner, AntibodyOrAntigenReporter, and AntibodyAssignmentsMatrix
    to assign certain antibodies to cells, as well as generate plots, metrics, and summary outputs
    related to these assignments.
    Inputs:
        filtered_feature_counts_matrix: path to the filtered count matrix h5 file
    Outputs:
        antibody_histograms_json: histograms for all antibodies above a certain threshold
    """
    if not args.filtered_feature_counts_matrix:
        outs.antibody_histograms_json = None
        outs.antibody_treemap_json = None
        return

    filtered_feature_counts_matrix = cr_matrix.CountMatrix.load_h5_file(
        args.filtered_feature_counts_matrix
    )
    # Determine whether the input is antigen or antibody
    feature_type = (
        rna_library.ANTIBODY_LIBRARY_TYPE if args.is_antibody else rna_library.ANTIGEN_LIBRARY_TYPE
    )

    filtered_ab_counts_matrix = filtered_feature_counts_matrix.select_features_by_type(feature_type)
    del filtered_feature_counts_matrix

    # TODO: Maybe someday use these results in CS
    # ### Run antibody assignments and generate various summary outputs
    # ab_assigner = AntibodyAssigner(
    #     matrix=filtered_feature_counts_matrix, feature_type=feature_type
    # )
    # ab_assigner.assignments = ab_assigner.get_feature_assignments()
    # ab_assigner.assignment_metadata = ab_assigner.compute_assignment_metadata()

    # Read in the multi graph

    is_hashtag_multiplexed = False
    if args.multi_graph:
        multi_graph = MultiGraph.from_path(args.multi_graph)
        is_hashtag_multiplexed = multi_graph.is_hashtag_multiplexed()

    hashtags = multi_graph.sample_tag_ids().get(args.sample_id) if is_hashtag_multiplexed else None

    antibody_histograms = make_antibody_histograms(
        filtered_ab_counts_matrix, is_spatial=args.is_spatial, hashtags=hashtags
    )
    antibody_treemap = make_antibody_treemap_plot(
        filtered_ab_counts_matrix, feature_type, MIN_ANTIBODY_UMI, args.is_spatial
    )
    if antibody_histograms is None:
        outs.antibody_histograms_json = None
    else:
        feature_utils.write_json_from_dict(antibody_histograms, outs.antibody_histograms_json)
    if antibody_treemap is None:
        outs.antibody_treemap_json = None
    else:
        feature_utils.write_json_from_dict(antibody_treemap, outs.antibody_treemap_json)
