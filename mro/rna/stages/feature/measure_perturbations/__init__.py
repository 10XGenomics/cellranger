#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#

import sys

import pandas as pd
from six import ensure_str

import cellranger.feature.crispr.measure_perturbations as measure_perturbations

pd.set_option("compute.use_numexpr", False)
import cellranger.analysis.diffexp as cr_diffexp
import cellranger.feature.utils as feature_utils
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library

SUMMARY_FILE_NAME = "transcriptome_analysis"

__MRO__ = """
stage MEASURE_PERTURBATIONS(
    in  csv  protospacer_calls_per_cell,
    in  h5   filtered_feature_counts_matrix,
    in  csv  feature_reference,
    in  bool by_feature,
    in  bool ignore_multiples,
    out csv  perturbation_efficiencies,
    out path perturbation_effects_path,
    src py   "stages/feature/measure_perturbations",
) split (
)
"""


def split(args):
    # Note that the complex data flow here makes it very difficult to
    # precisely determine the memory usage as a function of the inputs.
    mem_gib = 1 + cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(
        args.filtered_feature_counts_matrix, scale=10
    )

    # Sometimes this stage uses a lot of vmem. It has proven difficult to reproduce and measure.
    vmem_gib = 6 + 2 * mem_gib
    return {"chunks": [], "join": {"__mem_gb": mem_gib, "__threads": 5, "__vmem_gb": vmem_gib}}


def join(args, outs, _chunk_defs, _chunk_outs):
    list_file_paths = [
        args.protospacer_calls_per_cell,
        args.filtered_feature_counts_matrix,
        args.feature_reference,
    ]
    if not (feature_utils.all_files_present(list_file_paths)):
        outs.perturbation_efficiencies = None
        return

    feature_count_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_feature_counts_matrix)
    gex_count_matrix = feature_count_matrix.select_features_by_type(
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE
    )
    feature_ref_table = pd.read_csv(ensure_str(args.feature_reference), na_filter=False)

    protospacers_per_cell = pd.read_csv(
        ensure_str(args.protospacer_calls_per_cell), index_col=0, na_filter=False
    )
    if "target_gene_id" not in feature_ref_table.columns.tolist():
        sys.stderr.write(
            "Feature ref does not specify target gene IDs; that is a requirement for measuring perturbation efficiencies"
        )
        outs.perturbation_efficiencies = None
        return

    if "Non-Targeting" not in list(feature_ref_table["target_gene_id"].values):
        sys.stderr.write(
            "Non-Targeting guides required as controls for differential expression calculations"
        )
        outs.perturbation_efficiencies = None
        return

    (
        log2_fold_change,
        log2_fold_change_ci,
        num_cells_per_perturbation,
        results_per_perturbation,
        results_all_perturbations,
    ) = measure_perturbations.get_perturbation_efficiency(
        feature_ref_table,
        protospacers_per_cell,
        feature_count_matrix,
        args.by_feature,
        args.ignore_multiples,
    )

    # results_all_perturbations is an OrderedDict. The call to save_differential_expression_csv below assumes this when it assigns cluster_names from the keys of the dict

    if (log2_fold_change is None) or (log2_fold_change_ci is None):
        outs.perturbation_efficiencies = None
        return
    perturbation_efficiency_summary = (
        measure_perturbations.construct_perturbation_efficiency_summary(
            log2_fold_change,
            log2_fold_change_ci,
            num_cells_per_perturbation,
            args.by_feature,
        )
    )
    perturbation_efficiency_summary.to_csv(outs.perturbation_efficiencies, index=False)
    measure_perturbations.save_top_perturbed_genes(
        outs.perturbation_effects_path, results_per_perturbation
    )
    cr_diffexp.save_differential_expression_csv(
        None,
        results_all_perturbations,
        gex_count_matrix,
        outs.perturbation_effects_path,
        cluster_names=list(results_per_perturbation.keys()),
        file_name=SUMMARY_FILE_NAME,
    )
    # this call assumes that results_all_perturbations is an OrderedDict, hence can get ordered names from keys()
