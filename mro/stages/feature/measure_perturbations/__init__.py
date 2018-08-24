#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved
#
import cellranger.feature.crispr.measure_perturbations as measure_perturbations
import pandas as pd
pd.set_option("compute.use_numexpr", False)
import cellranger.matrix as cr_matrix
import sys
import cellranger.feature.utils as feature_utils

__MRO__ = """
stage MEASURE_PERTURBATIONS_PD(
    in  csv protospacer_calls_per_cell,
    in  h5  filtered_gene_bc_matrices_h5,
    in  csv feature_reference,
    in  bool by_feature,
    in  bool ignore_multiples,
    out csv perturbation_efficiencies,
    out path transcriptome_analysis_csv,
    src py  "stages/feature_pd/measure_perturbations_pd",
) using (
    mem_gb = 6,
)
"""

def main(args, outs):
    list_file_paths = [args.protospacer_calls_per_cell, args.filtered_feature_counts_matrix, args.feature_reference]
    if not(feature_utils.all_files_present(list_file_paths)):
        outs.perturbation_efficiencies = None
        return

    feature_count_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_feature_counts_matrix)
    feature_ref_table = pd.read_csv(args.feature_reference)

    protospacers_per_cell = pd.read_csv(args.protospacer_calls_per_cell, index_col = 0, na_filter = False)

    if "target_gene_id" not in feature_ref_table.columns.tolist():
        sys.stderr.write("Feature ref does not specify target gene IDs; that is a requirement for measuring perturbation efficiencies")
        outs.perturbation_efficiencies = None
        return

    if "Non-Targeting" not in list(feature_ref_table['target_gene_id'].values):
        sys.stderr.write("Non-Targeting guides required as controls for differential expression calculations")
        outs.perturbation_efficiencies = None
        return

    (log2_fold_change, log2_fold_change_ci,
        num_cells_per_perturbation,
        results_all_perturbations, results_per_perturbation) = measure_perturbations.get_perturbation_efficiency(
                                                                            feature_ref_table,
                                                                            protospacers_per_cell,
                                                                            feature_count_matrix,
                                                                            args.by_feature,
                                                                            args.ignore_multiples,
                                                                            )

    if (log2_fold_change is None) or (log2_fold_change_ci is None):
        outs.perturbation_efficiencies = None
        return
    summary_df = measure_perturbations.construct_df(log2_fold_change, log2_fold_change_ci,
                                                    num_cells_per_perturbation, args.by_feature)
    summary_df.to_csv(outs.perturbation_efficiencies, index=False)
    measure_perturbations.save_transcriptome_analysis_csv(outs.transcriptome_analysis_csv, results_per_perturbation)

