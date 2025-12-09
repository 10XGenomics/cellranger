# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
# Keeping the annotation import for now until PEP 563 is mandatory in a future Python version
# Check: https://peps.python.org/pep-0563/
# Check: https://github.com/astral-sh/ruff/issues/7214
from __future__ import annotations

import csv
import sys
from collections import OrderedDict
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from cellranger.analysis.analysis_types import DifferentialExpression
from cellranger.analysis.diffexp import (
    get_local_sseq_params_from_presubset_matrix,
    sseq_differential_expression,
)
from cellranger.constants import FILTER_LIST
from cellranger.feature.feature_assignments import CELL_BARCODE, FEATURE_CALL, NUM_FEATURES
from cellranger.rna.feature_ref import TARGET_GENE_ID, TARGET_GENE_NAME
from tenkit.stats import robust_divide

if TYPE_CHECKING:
    from typing import Any

    from cellranger.feature_ref import FeatureReference

    fc_t = tuple[float, float, int, int]
    ci_t = tuple[np.float64, np.float64]
    np_1d_array_int64 = np.ndarray[tuple[int], np.dtype[np.int64]]
    np_1d_array_intp = np.ndarray[tuple[int], np.dtype[np.intp]]
    from cellranger.matrix import CountMatrix

pd.set_option("compute.use_numexpr", False)

NUM_THREADS = 4


@dataclass
class TargetInfo:
    target_id: str
    gene_id: str
    gene_name: str


@dataclass(frozen=True)
class ProtospacerCall:
    feature_call: str
    gene_id: str


@dataclass
class PerturbationResult:
    perturbation_name: str
    target_or_gene_name: str
    log2_fc: float
    p_val: float
    lower_bound: np.float64
    upper_bound: np.float64
    num_cells_with_perb: int
    mean_umi_count_perb: float
    num_cells_with_control: int
    mean_umi_count_control: float


PERT_EFFI_SUMM_COLS_BY_FEAT = (
    "Perturbation",
    "Target Guide",
    "Log2 Fold Change",
    "p Value",
    "Log2 Fold Change Lower Bound",
    "Log2 Fold Change Upper Bound",
    "Cells with Perturbation",
    "Mean UMI Count Among Cells with Perturbation",
    "Cells with Non-Targeting Guides",
    "Mean UMI Count Among Cells with Non-Targeting Guides",
)
PERT_EFFI_SUMM_COLS_BY_TRGT = (
    PERT_EFFI_SUMM_COLS_BY_FEAT[:1] + ("Target Gene",) + PERT_EFFI_SUMM_COLS_BY_FEAT[2:]
)
NUM_BOOTSTRAPS = 500  # number of bootstrap draws to do for calculating
# empirical confidence intervals for perturbation efficiencies
MIN_BOOTSTRAPS = 50  # minimum number of bootstrap iterations before convergence check
CONVERGENCE_WINDOW = 20  # window size for convergence check
CONVERGENCE_TOLERANCE = 0.01  # tolerance for CI convergence
CI_LOWER_BOUND = 5.0  # CI lower bound (ie percentile value) for perturbation efficiencies
CI_UPPER_BOUND = 95.0  # CI upper bound (ie percentile value) for perturbation efficiencies
# for which we can't or won't compute perturbation efficiencies
NON_TARGETING = "Non-Targeting"  # Target ID used to specify control perturbations
MIN_NUMBER_CELLS_PER_PERTURBATION = 10  # Minimum number of cells a perturbation has to have before we compute differential expression for it
MIN_COUNTS_PERTURBATION = 5
MIN_COUNTS_CONTROL = 5

UMI_NUM_TRIES = 10  # Number of initial points to try for GMM-fitting
UMI_MIX_INIT_SD = 0.25  # Intial standard deviation for GMM components


TOP_GENES_SUMMARY_MAP = OrderedDict(
    [
        ("Gene Name", "Gene Name"),
        ("Gene ID", "Gene ID"),
        ("log2_fold_change", "Log2 Fold Change"),
        ("adjusted_p_value", "Adjusted p-value"),
    ]
)


def read_and_validate_feature_ref(
    feature_reference: FeatureReference,
) -> None | dict[str, TargetInfo]:
    """Returns a dict of target_id to TargetInfo."""
    # This method used to take the feature_reference CSV as input
    if TARGET_GENE_ID not in feature_reference.all_tag_keys:
        sys.stderr.write(
            "feature_reference CSV does not have target_gene_id column which "
            + "is a requirement for measuring perturbation efficiencies"
        )
        return None
    assert (
        len(feature_reference.feature_defs) > 0
    ), "This feature_reference does not have any feature defs"
    # Not all features in the feature def will have a `target_gene_id`, only those relevant to CRISPR
    # will, subselect on those.  Also check if this feature reference has both "target_gene_id" and "target_gene_name",
    # if only the former, use it as an alias for the later
    first_target = next(
        (x for x in feature_reference.feature_defs if TARGET_GENE_ID in x.tags), None
    )
    name_column = TARGET_GENE_NAME
    if first_target is None:
        sys.stderr.write(
            "Did not find any feature reference definition with the 'target_gene_id' value."
        )
        return None
    elif TARGET_GENE_NAME not in first_target.tags:
        name_column = TARGET_GENE_ID

    target_info = {
        x.id.decode(): TargetInfo(
            target_id=x.id.decode(), gene_id=x.tags[TARGET_GENE_ID], gene_name=x.tags[name_column]
        )
        for x in feature_reference.feature_defs
        if TARGET_GENE_ID in x.tags
    }

    if not any(v.gene_id == NON_TARGETING for v in target_info.values()):
        sys.stderr.write(
            "Non-Targeting guides required as controls for differential expression calculations"
        )
        return None
    return target_info


def save_perturbation_efficiency_summary(
    outpath: str,
    fold_change_per_perturbation: dict[tuple[str, int], dict[str, PerturbationResult]],
    by_feature: bool,
) -> None:
    with open(outpath, "w", newline="") as outfile:
        writer = csv.writer(outfile, dialect="unix", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(PERT_EFFI_SUMM_COLS_BY_FEAT if by_feature else PERT_EFFI_SUMM_COLS_BY_TRGT)
        rows = sorted(
            (
                (
                    fc.perturbation_name,
                    fc.target_or_gene_name,
                    fc.log2_fc,
                    fc.p_val,
                    fc.lower_bound,
                    fc.upper_bound,
                    fc.num_cells_with_perb,
                    fc.mean_umi_count_perb,
                    fc.num_cells_with_control,
                    fc.mean_umi_count_control,
                )
                for fold_changes in fold_change_per_perturbation.values()
                for fc in fold_changes.values()
            ),
            key=lambda x: x[2],
        )
        writer.writerows(rows)


def sanitize_perturbation_results(
    res_table: pd.DataFrame,
) -> pd.DataFrame:
    # at least 1 count amongst all the control cells
    res_table = res_table[res_table["sum_b"] > 0]
    # at least the minimum number of counts in either category (perturbation or control)
    res_table = res_table[
        (res_table["sum_a"] >= MIN_COUNTS_PERTURBATION) | (res_table["sum_b"] >= MIN_COUNTS_CONTROL)
    ]
    res_table["abs_log2_fold_change"] = np.abs(res_table["log2_fold_change"])
    # sort by abs log2 fold change, adjusted_p_value, Gene Name, in that order
    res_table.sort_values(
        by=["abs_log2_fold_change", "adjusted_p_value", "Gene Name"],
        ascending=[False, True, True],
        inplace=True,
    )
    res_table = res_table[list(TOP_GENES_SUMMARY_MAP)]
    res_table.reset_index(drop=True, inplace=True)
    return res_table


def _get_bc_targets_dict(
    target_info: dict[str, TargetInfo],
    protospacers_per_cell_path: str,
) -> dict[str, ProtospacerCall]:
    """Get the dict of barcode pairs.

    The values are dicts with information about the targets of the features assigned to the barcode.

    Returns:
        a dict of bc:{} pairs. All barcodes which have been assigned at least 1 protospacer are keys in this dict.

    Note:
        barcodes without any protospacers will not be present in this dict.
    """
    bc_targets_dict: dict[str, ProtospacerCall] = {}

    with open(protospacers_per_cell_path) as csv_file:
        csv_reader = csv.reader(csv_file)

        # Read the header
        col_name_to_idx = {n: i for i, n in enumerate(next(csv_reader))}
        feature_call_idx = col_name_to_idx[FEATURE_CALL]
        num_features_idx = col_name_to_idx[NUM_FEATURES]
        cell_barcode_idx = col_name_to_idx[CELL_BARCODE]
        for row in csv_reader:
            feature_call = row[feature_call_idx]
            num_features = int(row[num_features_idx])
            cell_barcode = row[cell_barcode_idx]

            if feature_call in target_info:
                # single feature
                gene_id = target_info[feature_call].gene_id
            else:
                gene_id: str = feature_call

            if num_features > 1:
                # multiple features
                this_features = feature_call.split("|")
                gene_ids = [target_info[x].gene_id for x in this_features]
                if all(gene_id == NON_TARGETING for gene_id in gene_ids):
                    # each feature is a non-targeting guide, and so the cell is a control cell
                    gene_id = NON_TARGETING
                else:
                    gene_ids = sorted(set(gene_ids).difference(FILTER_LIST))
                    if gene_ids:
                        gene_id = "|".join(gene_ids)
                    else:
                        gene_id = "Ignore"

            bc_targets_dict[cell_barcode] = ProtospacerCall(
                feature_call=feature_call,
                gene_id=gene_id,
            )
    return bc_targets_dict


def _should_filter(
    perturbation_name: str,
    target_id_name_map: dict[str, str],
):
    target_tuple = _get_target_id_from_name(perturbation_name, target_id_name_map)
    return all(x in FILTER_LIST for x in target_tuple[1])


def _get_target_id_from_name(
    this_perturbation_name: str, target_id_name_map: dict[str, str]
) -> tuple[list[str], list[str]]:
    if "|" not in this_perturbation_name:
        return ([this_perturbation_name], [target_id_name_map[this_perturbation_name]])

    p_names = this_perturbation_name.split("|")

    return (p_names, [target_id_name_map[p_name] for p_name in p_names])


def get_perturbation_efficiency(
    which_perturbations: set[str],
    target_info: dict[str, TargetInfo],
    protospacers_per_cell_path: str,
    gex_matrix: CountMatrix,
    by_feature: bool,
) -> (
    tuple[
        list[str],
        DifferentialExpression,
        dict[tuple[str, int], dict[str, PerturbationResult]],
    ]
    | None
):
    """Calculates log2 fold change and empirical confidence intervals for log2 fold change for target genes.

    Args:
        targets (list[TargetInfo]): list of TargetInfo dataclass objects
        protospacers_per_cell (str): path to the protospacer calls per cell CSV file
        gex_matrix (CountMatrix): Feature Barcode Matrix w/ only GEX data
        by_feature (bool): if True, cells are grouped by the combination of protospacers
                           present in them, rather than the gene targets of those protospacers.

    Returns:
        perturbation_names (list[str]): List of perturbation_names
        results_all_perturbations (DifferentialExpression): All DE results
        fold_change_per_perturbation (dict[tuple[str, int], dict[str, FoldChange]]): ((perturbation_name, perturbation_idx), dict[target/feature name, FoldChange])
    """
    (target_calls, perturbation_keys) = _get_ps_clusters(
        target_info,
        protospacers_per_cell_path,
        [bc.decode() for bc in gex_matrix.bcs],
        by_feature,
    )

    target_to_gene_id = {
        target.target_id if by_feature else target.gene_name: target.gene_id
        for target in target_info.values()
    }
    perturbation_names: list[str] = []

    # Create a numpy array with 3*k columns, where k is the number of perturbations
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((gex_matrix.features_dim, 3 * len(which_perturbations)))

    nt_indices = [x for x in perturbation_keys if perturbation_keys[x] == NON_TARGETING]

    if len(nt_indices) == 0:
        return None
    nt_index = nt_indices[0]
    in_control_cluster = target_calls == nt_index
    control_num_cells = sum(in_control_cluster)

    # Early exit if insufficient control cells
    if control_num_cells < MIN_NUMBER_CELLS_PER_PERTURBATION:
        return None

    feature_defs = gex_matrix.feature_ref.feature_defs
    gene_ids = [feature_def.id.decode() for feature_def in feature_defs]
    gene_names = [feature_def.name for feature_def in feature_defs]

    fold_change_per_perturbation: dict[tuple[str, int], dict[str, PerturbationResult]] = {}
    cluster_counter = 1
    column_counter = 0

    # Pre-compute control group for reuse
    group_b = np.flatnonzero(in_control_cluster)
    lgb = len(group_b)
    for cluster, perturbation_name in perturbation_keys.items():
        if perturbation_name not in which_perturbations or _should_filter(
            perturbation_name, target_to_gene_id
        ):
            continue

        in_cluster = target_calls == cluster
        group_a = np.flatnonzero(in_cluster)

        if len(group_a) < MIN_NUMBER_CELLS_PER_PERTURBATION:
            continue

        sys.stdout.flush()

        both_conditions = np.concatenate([group_a, group_b])
        lga = len(group_a)
        del group_a
        local_matrix = gex_matrix.select_barcodes(both_conditions)

        (
            local_sseq_params,
            new_group_a,
            new_group_b,
            matrix_groups,
        ) = get_local_sseq_params_from_presubset_matrix(local_matrix.m, lga, lgb)

        de_result = sseq_differential_expression(
            matrix_groups, new_group_a, new_group_b, local_sseq_params, threads=NUM_THREADS
        )
        assert de_result is not None
        de_result["Gene ID"] = gene_ids
        de_result["Gene Name"] = gene_names

        all_de_results[:, 0 + 3 * (cluster_counter - 1)] = de_result["sum_a"] / lga
        all_de_results[:, 1 + 3 * (cluster_counter - 1)] = de_result["log2_fold_change"]
        all_de_results[:, 2 + 3 * (cluster_counter - 1)] = de_result["adjusted_p_value"]
        column_counter += 3

        num_cells = sum(target_calls == cluster)

        fold_change_per_perturbation[perturbation_name, cluster] = _get_log2_fold_change(
            perturbation_name,
            num_cells,
            control_num_cells,
            de_result,
            target_to_gene_id,
            local_matrix,
            new_group_a,
            new_group_b,
            local_sseq_params,
        )
        del de_result
        perturbation_names.append(perturbation_name)
        cluster_counter += 1

    all_de_results = all_de_results[:, 0:column_counter]
    results_all_perturbations = DifferentialExpression(all_de_results)

    return (
        perturbation_names,
        results_all_perturbations,
        fold_change_per_perturbation,
    )


def _get_log2_fold_change(
    perturbation_name: str,
    num_cells: int,
    control_num_cells: int,
    results: pd.DataFrame,
    target_to_gene_id: dict[str, str],
    matrix: CountMatrix,
    group_a: np_1d_array_intp,
    group_b: np_1d_array_intp,
    local_params: dict[str, Any],
) -> dict[str, PerturbationResult]:
    (this_names, this_ids) = _get_target_id_from_name(
        perturbation_name,
        target_to_gene_id,
    )

    perturbation_results: dict[str, PerturbationResult] = {}

    for name, target in zip(this_names, this_ids):
        if target in FILTER_LIST:
            continue
        deg_result = results.loc[results["Gene ID"] == target]
        if deg_result.empty:
            continue

        lower_bound, upper_bound = _get_fold_change_cis(
            matrix,
            target,
            group_a,
            group_b,
            local_params,
        )
        log2_fold_change = deg_result["log2_fold_change"].values[0]
        p_value = deg_result["p_value"].values[0]
        sum_a = deg_result["sum_a"].values[0]
        sum_b = deg_result["sum_b"].values[0]
        perturbation_results[name] = PerturbationResult(
            perturbation_name=perturbation_name,
            target_or_gene_name=name,
            log2_fc=log2_fold_change,
            p_val=p_value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            num_cells_with_perb=num_cells,
            mean_umi_count_perb=robust_divide(sum_a, num_cells),
            num_cells_with_control=control_num_cells,
            mean_umi_count_control=robust_divide(sum_b, control_num_cells),
        )

    return perturbation_results


def _get_fold_change_cis(
    matrix: CountMatrix,
    target: str,
    cond_a: np_1d_array_intp,
    cond_b: np_1d_array_intp,
    computed_params: dict[str, Any],
) -> tuple[np.float64, np.float64]:
    np.random.seed(0)

    # filter the matrix to select only the target gene
    this_matrix = matrix.select_features_by_ids([target.encode()])

    x = this_matrix.m

    # Use numpy array instead of list for better performance
    log2_fold_change_vals = np.zeros(NUM_BOOTSTRAPS)

    # Early convergence detection
    for i in range(NUM_BOOTSTRAPS):
        this_cond_a = np.random.choice(cond_a, size=len(cond_a), replace=True)
        this_cond_b = np.random.choice(cond_b, size=len(cond_b), replace=True)

        gene_sums_a = np.sum(x[:, this_cond_a], axis=1)
        gene_sums_b = np.sum(x[:, this_cond_b], axis=1)
        size_factor_a = np.sum(computed_params["size_factors"][this_cond_a])
        size_factor_b = np.sum(computed_params["size_factors"][this_cond_b])

        log2_fold_change_vals[i] = np.log2((1 + gene_sums_a) / (1 + size_factor_a)) - np.log2(
            (1 + gene_sums_b) / (1 + size_factor_b)
        )

        # Check for convergence after minimum iterations
        if i >= MIN_BOOTSTRAPS and i % CONVERGENCE_WINDOW == 0:
            # Calculate current confidence interval
            current_vals = log2_fold_change_vals[: i + 1]
            lower_current = np.percentile(current_vals, CI_LOWER_BOUND)
            upper_current = np.percentile(current_vals, CI_UPPER_BOUND)

            # Calculate previous window confidence interval
            if i >= MIN_BOOTSTRAPS + CONVERGENCE_WINDOW:
                prev_vals = log2_fold_change_vals[: i + 1 - CONVERGENCE_WINDOW]
                lower_prev = np.percentile(prev_vals, CI_LOWER_BOUND)
                upper_prev = np.percentile(prev_vals, CI_UPPER_BOUND)

                # Check convergence
                if (
                    abs(lower_current - lower_prev) < CONVERGENCE_TOLERANCE
                    and abs(upper_current - upper_prev) < CONVERGENCE_TOLERANCE
                ):
                    # Converged, use current values
                    log2_fold_change_vals = log2_fold_change_vals[: i + 1]
                    break

    return (
        np.percentile(log2_fold_change_vals, CI_LOWER_BOUND),
        np.percentile(log2_fold_change_vals, CI_UPPER_BOUND),
    )


def _get_ps_clusters(
    target_info: dict[str, TargetInfo],
    protospacers_per_cell_path: str,
    barcodes: list[str],
    by_feature: bool = True,
) -> tuple[np_1d_array_int64, dict[int, str]]:
    """Returns a tuple (target_calls, perturbation_keys).

    Args:
        target_calls (np.array(int)): identifies the perturbation assigned to
            each cell in the gene-barcode matrix
        perturbation_keys (dict): (cluster_number:perturbation_name) pairs

    Returns:
        target_calls (np_1d_array_int64): integer perturbation target labels starting from 1
        perturbation_keys (dict): dict[cluster_number, perturbation_name]
    """
    bc_targets_dict = _get_bc_targets_dict(target_info, protospacers_per_cell_path)

    # Optimize: Use set operations for faster lookups
    missing_barcodes = set(barcodes) - set(bc_targets_dict.keys())
    default_call = ProtospacerCall(feature_call="None", gene_id="None")
    for bc in missing_barcodes:
        bc_targets_dict[bc] = default_call

    if by_feature:
        return _get_ps_clusters_by_feature(bc_targets_dict, barcodes)
    else:
        return _get_ps_clusters_by_target(target_info, bc_targets_dict, barcodes)


def _get_ps_clusters_by_target(
    target_info: dict[str, TargetInfo],
    bc_targets_dict: dict[str, ProtospacerCall],
    barcodes: list[str],
) -> tuple[np_1d_array_int64, dict[int, str]]:
    gene_id_to_gene_name = {v.gene_id: v.gene_name for v in target_info.values()}

    def get_target_name(gene_id: str) -> str:
        sep = "|"
        if sep not in gene_id:
            return gene_id_to_gene_name.get(gene_id, gene_id)
        return "|".join([gene_id_to_gene_name.get(gid, gid) for gid in gene_id.split(sep)])

    gene_ids = [bc_targets_dict[bc].gene_id for bc in barcodes]
    unique_gene_ids = sorted(set(gene_ids))

    gene_id_to_idx = {gene_id: idx for (idx, gene_id) in enumerate(unique_gene_ids, start=1)}
    target_calls = np.asarray([gene_id_to_idx[x] for x in gene_ids], dtype=np.int64)

    perturbation_keys = {idx: get_target_name(gid) for (gid, idx) in gene_id_to_idx.items()}

    return (target_calls, perturbation_keys)


def _get_ps_clusters_by_feature(
    bc_targets_dict: dict[str, ProtospacerCall],
    barcodes: list[str],
) -> tuple[np_1d_array_int64, dict[int, str]]:
    def _get_feature_from_pc(protospacer_call: ProtospacerCall) -> str:
        if protospacer_call.gene_id not in FILTER_LIST:
            return protospacer_call.feature_call

        if protospacer_call.gene_id == "None":
            return "Ignore"

        return protospacer_call.gene_id

    features = [_get_feature_from_pc(bc_targets_dict[bc]) for bc in barcodes]
    unique_features = sorted(set(features))

    feature_to_idx = {feat: idx for (idx, feat) in enumerate(unique_features, start=1)}
    target_calls = np.asarray([feature_to_idx[feat] for feat in features], dtype=np.int64)

    perturbation_keys = {b: a for a, b in feature_to_idx.items()}

    return (target_calls, perturbation_keys)
