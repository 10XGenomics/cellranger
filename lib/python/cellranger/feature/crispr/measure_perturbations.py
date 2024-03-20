# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
from __future__ import annotations

import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

from cellranger.analysis.analysis_types import DifferentialExpression
from cellranger.constants import FILTER_LIST

pd.set_option("compute.use_numexpr", False)
import collections
import os
import sys

import cellranger.analysis.diffexp as cr_diffexp
import cellranger.rna.library as rna_library
import tenkit.stats as tk_stats

NUM_BOOTSTRAPS = 500  # number of bootstrap draws to do for calculating
# empirical confidence intervals for perturbation efficiencies
CI_LOWER_BOUND = 5.0  # CI lower bound (ie percentile value) for perturbation efficiencies
CI_UPPER_BOUND = 95.0  # CI upper bound (ie percentile value) for perturbation efficiencies
# for which we can't or won't compute perturbation efficiencies
CONTROL_LIST = ["Non-Targeting"]  # Target IDs used for specifying control perturbations
MIN_NUMBER_CELLS_PER_PERTURBATION = 10  # Minimum number of cells a perturbation has to have before we compute differential expression for it
MIN_COUNTS_PERTURBATION = 5
MIN_COUNTS_CONTROL = 5

UMI_NUM_TRIES = 10  # Number of initial points to try for GMM-fitting
UMI_MIX_INIT_SD = 0.25  # Intial standard deviation for GMM components

PERTURBATION_EFFICIENCY_SUMMARY_COLUMNS = [
    "Perturbation",
    "target_string",
    "Log2 Fold Change",
    "p Value",
    "Log2 Fold Change Lower Bound",
    "Log2 Fold Change Upper Bound",
    "Cells with Perturbation",
    "Mean UMI Count Among Cells with Perturbation",
    "Cells with Non-Targeting Guides",
    "Mean UMI Count Among Cells with Non-Targeting Guides",
]

TOP_GENES_SUMMARY_MAP = collections.OrderedDict(
    [
        ("Gene Name", "Gene Name"),
        ("Gene ID", "Gene ID"),
        ("log2_fold_change", "Log2 Fold Change"),
        ("adjusted_p_value", "Adjusted p-value"),
    ]
)
NUM_TOP_GENES = 10


def construct_perturbation_efficiency_summary(
    f_change,
    f_change_ci,
    num_cells_per_perturbation,
    by_feature,
    summary_columns=PERTURBATION_EFFICIENCY_SUMMARY_COLUMNS,
):
    if (f_change is None) or (f_change_ci is None):
        return None

    if by_feature:
        summary_columns[1] = "Target Guide"
    else:
        summary_columns[1] = "Target Gene"

    this_df = pd.DataFrame(columns=summary_columns)
    counter = 0
    control_num_cells = num_cells_per_perturbation["Non-Targeting"]
    for key in sorted(f_change.keys()):
        this_key_results = f_change.get(key)
        this_key_ci = f_change_ci.get(key)
        if this_key_results is None:
            continue

        this_num_cells = num_cells_per_perturbation[key]

        for ps, results in this_key_results.items():
            lower_bound = this_key_ci.get(ps)[0]
            upper_bound = this_key_ci.get(ps)[1]
            this_df.loc[counter] = (
                key,
                ps,
                results[0],
                results[1],
                lower_bound,
                upper_bound,
                this_num_cells,
                tk_stats.robust_divide(results[2], this_num_cells),
                control_num_cells,
                tk_stats.robust_divide(results[3], control_num_cells),
            )
            counter += 1
    this_df.sort_values(by=["Log2 Fold Change"], ascending=True, inplace=True)
    return this_df


def save_top_perturbed_genes(
    base_dir,
    results_per_perturbation,
    column_map=TOP_GENES_SUMMARY_MAP,
    num_genes_to_keep=NUM_TOP_GENES,
):
    if not results_per_perturbation:
        return
    os.makedirs(base_dir, exist_ok=True)
    fn = os.path.join(base_dir, "top_perturbed_genes.csv")

    list_df_results = []
    summary_df_columns = []
    for perturbation in results_per_perturbation:
        this_results = sanitize_perturbation_results(results_per_perturbation.get(perturbation))

        if this_results is None:
            continue

        this_results = this_results[list(column_map)[:num_genes_to_keep]]
        this_results.reset_index(drop=True, inplace=True)
        list_df_results.append(this_results)
        summary_df_columns += [f"Perturbation: {perturbation}, {s}" for s in column_map.values()]

    summary_df = pd.concat(list_df_results, ignore_index=True, axis=1)
    summary_df.columns = summary_df_columns
    summary_df.to_csv(fn, index=False)


def sanitize_perturbation_results(res_table):
    if res_table is None or res_table.empty:
        return None

    res_table = res_table[res_table["sum_b"] > 0]
    # at least 1 count amongst all the control cells
    res_table = res_table[
        (res_table["sum_a"] >= MIN_COUNTS_PERTURBATION) | (res_table["sum_b"] >= MIN_COUNTS_CONTROL)
    ]
    # at least the minimum number of counts in either category (perturbation or control)

    res_table["abs_log2_fold_change"] = np.abs(res_table["log2_fold_change"])
    res_table.sort_values(
        by=["abs_log2_fold_change", "adjusted_p_value", "Gene Name"],
        ascending=[False, True, True],
        inplace=True,
    )
    # sort by abs log2 fold change, adjusted_p_value, Gene Name, in that order
    return res_table


def _add_bcs_without_ps_calls(bc_targets_dict, bcs):
    str_bcs = [ensure_str(bc) for bc in bcs]
    bcs_without_ps_calls = list(set(str_bcs).difference(set(bc_targets_dict.keys())))
    for bc in bcs_without_ps_calls:
        bc_targets_dict[bc] = {
            "target_name": "None",
            "num_features": -1,
            "target_id": "None",
            "feature_call": "None",
        }

    return bc_targets_dict


def _get_bc_targets_dict(
    target_info,
    protospacers_per_cell,
    ignore_multiples=False,
    control_list=CONTROL_LIST,
    filter_list=FILTER_LIST,
):
    """Get the dict of barcode pairs.

    The values are dicts with information about the targets of the features assgined to the barcode.

    Returns:
        a dict of bc:{} pairs. All barcodes which have been assigned at least 1 protospacer are keys in this dict.

    Note:
        barcodes without any protospacers will not be present in this dict.
    """
    bc_targets_dict = {}

    for this_row in protospacers_per_cell.itertuples():
        feature_call = this_row.feature_call

        if feature_call in target_info:
            # single feature
            this_target_id = target_info.get(feature_call).get("target_id")
            this_target_name = target_info.get(feature_call).get("target_name")
        else:
            this_target_id = this_target_name = feature_call

        if this_row.num_features > 1:
            # multiple features
            if ignore_multiples:
                this_target_id = this_target_name = "Ignore"
            else:
                this_features = feature_call.split("|")
                this_target_ids = [target_info.get(x).get("target_id") for x in this_features]
                this_target_names = [target_info.get(x).get("target_name") for x in this_features]

                if set(this_target_ids) == set(control_list):
                    this_target_id = this_target_name = "Non-Targeting"
                    # each feature is a non-targeting guide, and so the cell is a control cell
                else:
                    this_target_ids = list(set(this_target_ids).difference(set(filter_list)))
                    this_target_names = list(set(this_target_names).difference(set(filter_list)))
                    if this_target_ids:
                        this_target_id = "|".join(this_target_ids)
                        this_target_name = "|".join(this_target_names)
                    else:
                        this_target_id = this_target_name = "Ignore"

        bc_targets_dict[this_row.Index] = {
            "num_features": this_row.num_features,
            "feature_call": feature_call,
            "target_id": this_target_id,
            "target_name": this_target_name,
        }
    return bc_targets_dict


def _get_target_name_from_id(target_id_name_map, target_id, sep="|"):
    if sep not in target_id:
        return target_id_name_map.get(target_id, target_id)
    return "|".join([target_id_name_map.get(x, x) for x in target_id.split(sep)])


def _should_filter(
    perturbation_name,
    feature_ref_table,
    target_info,
    filter_list=FILTER_LIST,
    by_feature=False,
):
    target_id_name_map = _get_target_id_name_map(feature_ref_table, target_info, by_feature)
    target_tuple = _get_target_id_from_name(perturbation_name, target_id_name_map, target_info)
    return all(x in filter_list for x in target_tuple[1])


def _get_target_id_from_name(this_perturbation_name, target_id_name_map, target_info):
    if "|" not in this_perturbation_name:
        return ([this_perturbation_name], [target_id_name_map[this_perturbation_name]])

    p_names = this_perturbation_name.split("|")

    return (p_names, [target_id_name_map[p_name] for p_name in p_names])


def get_perturbation_efficiency(
    feature_ref_table,
    protospacers_per_cell,
    feature_count_matrix,
    by_feature=False,
    ignore_multiples=False,
    filter_list=FILTER_LIST,
    num_bootstraps=NUM_BOOTSTRAPS,
    ci_lower=CI_LOWER_BOUND,
    ci_upper=CI_UPPER_BOUND,
    sseq_params=None,
):
    """Calculates log2 fold change and empirical confidence intervals for log2 fold change for target genes.

    Args:
        feature_ref_table (pd.DataFrame): obtained by reading the feature reference file
        protospacers_per_cell (pd.DataFrame): protospacer calls per cell. Output of protospacer calling stage
        feature_count_matrix (CountMatrix): Feature Barcode Matrix
        by_feature (bool): if True, cells are grouped by the combination of protospacers
                           present in them, rather than the gene targets of those protospacers.
                           Right now, the pipeline measures perturbations both by feature
                           (by_feature=True) or by target (by_feature=False) by calling
                           the MEASURE_PERTURBATIONS_PD stage twice using the "as" keyword in
                           Martian 3.0
        ignore_multiples (bool): if True, cells with multiple protospacers are ignored from the analysis. Default is False
        filter_list (list): list of target id values to be filtered out of differential expression calculations
        num_bootstraps (int): number of bootstrap draws for computing empirical confideence interval for log2 fold change
        ci_lower (float):  percentile value for calculating the lower bound of the emprirical confidence interval for log2 fold change
        ci_upper (float): percentile value for calculating the upper bound of the emprirical confidence interval for log2 fold change
        sseq_params (dict): global parameters for differential expression calculations. If None (default),
                    computed by cr_diffexp.compute_sseq_params

    Returns:
        log2_fold_change: keys are perturbations, grouped by features or by targets.
                            Values are a dictionary where the keys are individual features or targets,
                            and values are a tuple of floats, each with 4 entries.
                            First entry is the log2_fold_change, second entry is the p-value.=,
                            third entry is the sum of UMI counts for that target gene across all cells in the condition,
                            fourth entry is the sum of UMI counts for that target gene across all cells in the control.
                            Eg.:: {'COA3': {'COA3': (-1.4, 0.02, 5, 20)},
                                'COX18|COX6C': {'COX18': (-0.87, 1.0, 3, 40), ...} ...}
        fold_change_CI: similar nested dict, except the tuple indicates the confidence interval.
                        First entry is the lower percentile bound and second is the
                        upper percentile bound.
                        Eg. {'COA3': {'COA3': (-2.12, -0.93)},
                                'COX18|COX6C': {'COX18': (-1.16, -0.57),...}... }
    """
    matrix = feature_count_matrix.select_features_by_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
    target_info = _get_target_info(feature_ref_table)
    target_id_name_map = _get_target_id_name_map(feature_ref_table, target_info, by_feature)
    (cluster_per_bc, perturbation_keys) = _get_ps_clusters(
        target_info, protospacers_per_cell, matrix.bcs, by_feature, ignore_multiples
    )

    num_cells_per_perturbation = {
        name: sum(cluster_per_bc == cluster) for (cluster, name) in perturbation_keys.items()
    }

    diff_exp_results = _analyze_transcriptome(
        matrix,
        target_id_name_map,
        cluster_per_bc,
        perturbation_keys,
        feature_ref_table,
        by_feature,
        target_info,
        num_bootstraps,
        ci_lower,
        ci_upper,
        filter_list,
    )

    if diff_exp_results is None:
        return (None, None, None, None, None)

    results_per_perturbation = diff_exp_results["results_per_perturbation"]
    results_all_perturbations = diff_exp_results["results_all_perturbations"]
    log2_fold_change_CIs = diff_exp_results["log2_fold_change_CIs"]
    log2_fold_change = diff_exp_results["log2_fold_change"]

    return (
        log2_fold_change,
        log2_fold_change_CIs,
        num_cells_per_perturbation,
        results_per_perturbation,
        results_all_perturbations,
    )


def _analyze_transcriptome(
    matrix,
    target_id_name_map,
    target_calls,
    perturbation_keys,
    feature_ref_table,
    by_feature,
    target_info,
    num_bootstraps,
    ci_lower,
    ci_upper,
    filter_list=FILTER_LIST,
):
    """Compute differential expression for each perturbation vs non-targeting control.

    Args:
        matrix (CountMatrix): gene expression data, sub-selected from all feature expression data in the CountMatrix
        target_calls (np.array(int)): integer perturbation target labels.
                                        Numbering starts at 1
        perturbation_keys (dict): (cluster_number:perturbation_name) pairs
        feature_ref_table (pd.DataFrame): obtained by reading the feature reference file
        by_feature (bool): if True, cells are grouped by features (rather than by target)
        target_info (dict): Nested dict: {feature1: {'target_id': value1, 'target_name': value2}, ...}
        filter_list (list): list of target ids to be filtered out of the analysis

    Returns:
        dict: {'results_all_perturbations': named tuple containing a dataframe of differential expression outputs,
               'results_per_perturbation': dict of perturbation_name:diffexp_dataframe pairs,
               }
    """
    results_per_perturbation = collections.OrderedDict()
    # MEASURE_PERTURBATIONS assumes this dict to be ordered
    filter_cluster_indices = [x for x in perturbation_keys if perturbation_keys[x] in filter_list]
    n_clusters = len(perturbation_keys)
    n_effective_clusters = n_clusters - len(filter_cluster_indices)

    # Create a numpy array with 3*k columns, where k is the number of perturbations
    # k = n_clusters - 2 since we don't need to report results for ignore and non-targeting
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((matrix.features_dim, 3 * n_effective_clusters))

    nt_indices = [x for x in perturbation_keys if perturbation_keys[x] == "Non-Targeting"]

    if (nt_indices is None) or (nt_indices == []):
        return None
    nt_index = nt_indices[0]

    in_control_cluster = target_calls == nt_index
    feature_defs = matrix.feature_ref.feature_defs
    gene_ids = [ensure_str(feature_def.id) for feature_def in feature_defs]
    gene_names = [ensure_str(feature_def.name) for feature_def in feature_defs]

    log2_fold_change_CIs = {}
    log2_fold_change = {}
    cluster_counter = 1
    column_counter = 0
    for cluster in perturbation_keys:
        perturbation_name = perturbation_keys.get(cluster)
        if (cluster in filter_cluster_indices) or _should_filter(
            perturbation_name, feature_ref_table, target_info, filter_list, by_feature
        ):
            continue

        in_cluster = target_calls == cluster
        group_a = np.flatnonzero(in_cluster)

        if len(group_a) < MIN_NUMBER_CELLS_PER_PERTURBATION:
            continue

        group_b = np.flatnonzero(in_control_cluster)
        sys.stdout.flush()

        both_conditions = np.concatenate([group_a, group_b])
        local_matrix = matrix.select_barcodes(both_conditions)

        (
            local_sseq_params,
            new_group_a,
            new_group_b,
            matrix_groups,
        ) = cr_diffexp.get_local_sseq_params(matrix.m, group_a, group_b)

        de_result = cr_diffexp.sseq_differential_expression(
            matrix_groups, new_group_a, new_group_b, local_sseq_params
        )
        de_result["Gene ID"] = gene_ids
        de_result["Gene Name"] = gene_names
        results_per_perturbation[perturbation_name] = de_result

        all_de_results[:, 0 + 3 * (cluster_counter - 1)] = de_result["sum_a"] / len(group_a)
        all_de_results[:, 1 + 3 * (cluster_counter - 1)] = de_result["log2_fold_change"]
        all_de_results[:, 2 + 3 * (cluster_counter - 1)] = de_result["adjusted_p_value"]
        column_counter += 3

        (this_log2_fc, this_log2_cis) = _get_log2_fold_change(
            perturbation_name,
            de_result,
            target_id_name_map,
            target_info,
            local_matrix,
            new_group_a,
            new_group_b,
            local_sseq_params,
        )

        log2_fold_change[perturbation_name] = this_log2_fc
        log2_fold_change_CIs[perturbation_name] = this_log2_cis

        cluster_counter += 1

    all_de_results = all_de_results[:, 0:column_counter]

    return {
        "results_per_perturbation": results_per_perturbation,
        "results_all_perturbations": DifferentialExpression(all_de_results),
        "log2_fold_change_CIs": log2_fold_change_CIs,
        "log2_fold_change": log2_fold_change,
    }


def _get_log2_fold_change(
    perturbation_name,
    this_results,
    target_id_name_map,
    target_info,
    matrix,
    group_a,
    group_b,
    local_params,
    filter_list=FILTER_LIST,
):
    (this_names, this_ids) = _get_target_id_from_name(
        perturbation_name, target_id_name_map, target_info
    )

    log2_fc = {}
    log2_cis = {}

    for name, target in zip(this_names, this_ids):
        if target in filter_list:
            continue

        deg_result = _get_ko_per_target(this_results, target)
        if deg_result is not None:
            log2_fc[name] = deg_result
            log2_cis[name] = _get_fold_change_cis(matrix, target, group_a, group_b, local_params)

    return (log2_fc, log2_cis)


def _get_ko_per_target(results, target):
    deg_result = results.loc[results["Gene ID"] == target]
    if deg_result.empty:
        return None
    return (
        deg_result["log2_fold_change"].values[0],
        deg_result["p_value"].values[0],
        deg_result["sum_a"].values[0],
        deg_result["sum_b"].values[0],
    )


def _get_fold_change_cis(
    matrix,
    target,
    cond_a,
    cond_b,
    computed_params,
    num_bootstraps=NUM_BOOTSTRAPS,
    ci_lower=CI_LOWER_BOUND,
    ci_upper=CI_UPPER_BOUND,
):
    # filter the matrix to select only the target gene
    this_matrix = matrix.select_features_by_ids([ensure_binary(target)])

    x = this_matrix.m

    log2_fold_change_vals = list(range(num_bootstraps))
    for i in range(num_bootstraps):
        this_cond_a = np.random.choice(cond_a, size=len(cond_a), replace=True)
        this_cond_b = np.random.choice(cond_b, size=len(cond_b), replace=True)
        x_a = x[:, this_cond_a]
        x_b = x[:, this_cond_b]

        # Size factors
        size_factor_a = np.sum(computed_params["size_factors"][cond_a])
        size_factor_b = np.sum(computed_params["size_factors"][cond_b])

        gene_sums_a = np.squeeze(np.asarray(x_a.sum(axis=1)))
        gene_sums_b = np.squeeze(np.asarray(x_b.sum(axis=1)))

        log2_fold_change_vals[i] = np.log2((1 + gene_sums_a) / (1 + size_factor_a)) - np.log2(
            (1 + gene_sums_b) / (1 + size_factor_b)
        )

    return (
        np.percentile(log2_fold_change_vals, ci_lower),
        np.percentile(log2_fold_change_vals, ci_upper),
    )


def _get_target_id_name_map(feature_ref_table, target_info, by_feature=False):
    """Returns a dict from either gene name or protospacer to gene id.

    Returns:
        a dict of target_gene_name:target_gene_id pairs if by_feature = False
        and a dict of protospacer:target_gene_id pairs if by_feature = True
    """
    if not (by_feature):
        return {target_info[x]["target_name"]: target_info[x]["target_id"] for x in target_info}
    return {x: target_info[x]["target_id"] for x in target_info}


def _get_target_info(feature_ref_table):
    """Returns a nested dict.

    {feature1: {'target_id': value1, 'target_name': value2}, ...}
    """
    target_info = {}
    for this_row in feature_ref_table.itertuples():
        this_target_id = getattr(this_row, "target_gene_id")
        this_target_name = getattr(this_row, "target_gene_name", this_target_id)
        target_info[this_row.id] = {
            "target_id": this_target_id,
            "target_name": this_target_name,
        }

    return target_info


def _get_ps_clusters(
    target_info, protospacers_per_cell, bcs, by_feature=True, ignore_multiples=True
):
    """Returns a tuple (target_calls, perturbation_keys).

    Args:
        target_calls (np.array(int)): identifies the perturbation assigned to
            each cell in the gene-barcode matrix
        perturbation_keys (dict): (cluster_number:perturbation_name) pairs
    """
    bc_targets_dict = _get_bc_targets_dict(target_info, protospacers_per_cell, ignore_multiples)
    bc_targets_dict = _add_bcs_without_ps_calls(bc_targets_dict, bcs)
    if not (by_feature):
        return _get_ps_clusters_by_target(
            target_info, bc_targets_dict, protospacers_per_cell, bcs, ignore_multiples
        )
    return _get_ps_clusters_by_feature(
        target_info, bc_targets_dict, protospacers_per_cell, bcs, ignore_multiples
    )


def _get_ps_clusters_by_target(
    target_info, bc_targets_dict, protospacers_per_cell, bcs, ignore_multiples=True
):
    target_id_name_map = {v.get("target_id"): v.get("target_name") for v in target_info.values()}
    target_id_per_bc = [bc_targets_dict.get(ensure_str(bc)).get("target_id") for bc in bcs]
    unique_targets_list = list(set(target_id_per_bc))
    unique_targets_list.sort()

    target_to_int = {el: i + 1 for (i, el) in enumerate(unique_targets_list)}

    perturbation_from_cluster = {
        b: _get_target_name_from_id(target_id_name_map, a, sep="|")
        for (a, b) in target_to_int.items()
    }

    calls_vector = np.asarray([target_to_int.get(x) for x in target_id_per_bc])

    return (calls_vector, perturbation_from_cluster)


def _get_ps_clusters_by_feature(
    target_info, bc_targets_dict, protospacers_per_cell, bcs, ignore_multiples=True
):
    feature_per_bc = [_get_feature_from_bc(bc_targets_dict, ensure_str(bc)) for bc in bcs]
    unique_feature_calls = list(set(feature_per_bc))
    unique_feature_calls.sort()

    feature_to_int = {el: i + 1 for (i, el) in enumerate(unique_feature_calls)}
    perturbation_from_cluster = {b: a for (a, b) in feature_to_int.items()}
    calls_vector = np.asarray([feature_to_int.get(x) for x in feature_per_bc])

    return (calls_vector, perturbation_from_cluster)


def _get_feature_from_bc(bc_targets_dict, bc, filter_list=FILTER_LIST):
    if bc_targets_dict.get(bc).get("target_id") not in filter_list:
        return bc_targets_dict.get(bc).get("feature_call")

    if bc_targets_dict.get(bc).get("target_id") in ["None"]:
        return "Ignore"

    return bc_targets_dict.get(bc).get("target_id")
