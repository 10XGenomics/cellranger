# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
# Keeping the annotation import for now until PEP 563 is mandatory in a future Python version
# Check: https://peps.python.org/pep-0563/
# Check: https://github.com/astral-sh/ruff/issues/7214
from __future__ import annotations

import csv
import sys
from collections import Counter, OrderedDict
from dataclasses import dataclass
from typing import TYPE_CHECKING, NamedTuple

import h5py as h5
import numpy as np
import pandas as pd

from cellranger.analysis.analysis_types import DifferentialExpression
from cellranger.analysis.diffexp import SSEQ_ZETA_QUANTILE
from cellranger.constants import FILTER_LIST
from cellranger.fast_utils import nb_de_from_stats_o3
from cellranger.feature.feature_assignments import CELL_BARCODE, FEATURE_CALL, NUM_FEATURES
from cellranger.rna.feature_ref import TARGET_GENE_ID, TARGET_GENE_NAME
from tenkit.stats import robust_divide

if TYPE_CHECKING:
    from cellranger.feature_ref import FeatureReference

    fc_t = tuple[float, float, int, int]
    ci_t = tuple[np.float64, np.float64]
    np_1d_array_int64 = np.ndarray[tuple[int], np.dtype[np.int64]]
    np_1d_array_intp = np.ndarray[tuple[int], np.dtype[np.intp]]
    np_1d_array_f64 = np.ndarray[tuple[int], np.dtype[np.float64]]
    np_1d_array_u64 = np.ndarray[tuple[int], np.dtype[np.uint64]]
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


class _PerturbationDEArrays(NamedTuple):
    """Per-gene sSeq DE result arrays for one perturbation (from nb_de_from_stats_o3).

    Bundled because the four arrays are all gene-indexed outputs of the same DE call and travel
    together into the fold-change assembly.
    """

    log2_fold_change: np_1d_array_f64
    p_values: np_1d_array_f64
    # sums_in/sums_out are the per-gene group count sums from the Rust DE kernel
    # (`DiffExpResult.sums_{in,out}: Array1<u64>`), so they arrive as uint64, not float.
    sums_in: np_1d_array_u64
    sums_out: np_1d_array_u64


@dataclass
class PerturbationStats:
    """Per-(gene x group) sSeq sufficient statistics emitted by SPLIT_PERTURBATION_MATRIX.

    Each row of the 2-D arrays is one group (the shared control in ``control_stats``; the chunk's
    perturbations in a ``chunk_stats`` file). These feed ``nb_de_from_stats_o3`` in place of the
    full per-chunk matrix; the streaming reducer that produced them is bit-identical to the in-memory
    reduction, so end-to-end the feature_sums/log2FC/bootstrap-CIs match the matrix path bit-for-bit
    and only the DE p-values drift (~1e-11, within the rtol=1e-6 gate).
    """

    n_genes: int
    names: list[str]
    feature_sums: np.ndarray  # [n_groups, n_genes] uint64
    sum_w: np.ndarray  # [n_groups, n_genes] float64
    sum_wsq: np.ndarray  # [n_groups, n_genes] float64
    sum_ct: np.ndarray  # [n_groups] float64
    sum_inv: np.ndarray  # [n_groups] float64
    n: np.ndarray  # [n_groups] int64
    count_totals: list[np.ndarray]  # per group, per-cell count_total (float64)
    name_to_idx: dict[str, int]


def read_perturbation_stats(path: str) -> PerturbationStats:
    """Load a sufficient-statistics H5 (``control_stats`` or a ``chunk_stats`` file)."""
    with h5.File(path, "r") as f:
        g = f["stats"]
        feature_sums = g["feature_sums"][:]
        sum_w = g["sum_w"][:]
        sum_wsq = g["sum_wsq"][:]
        sum_ct = g["sum_ct"][:]
        sum_inv = g["sum_inv"][:]
        n = g["n"][:]
        offsets = g["count_totals_offsets"][:]
        flat_count_totals = g["count_totals"][:]
        raw_names = g["perturbation_names"][:]
        n_genes = int(g.attrs["n_genes"])
    names = [name.decode() if isinstance(name, bytes) else str(name) for name in raw_names]
    count_totals = [flat_count_totals[offsets[i] : offsets[i + 1]] for i in range(len(names))]
    return PerturbationStats(
        n_genes=n_genes,
        names=names,
        feature_sums=feature_sums,
        sum_w=sum_w,
        sum_wsq=sum_wsq,
        sum_ct=sum_ct,
        sum_inv=sum_inv,
        n=n,
        count_totals=count_totals,
        name_to_idx={name: i for i, name in enumerate(names)},
    )


def read_bootstrap_count_total(path: str) -> np.ndarray:
    """Per-cell count_total (sum over all GEX genes) stored alongside a SPLIT bootstrap matrix.

    Written in the matrix's cell (column) storage order, so concatenating the control matrix's and a
    chunk matrix's arrays gives the combined per-cell size-factor numerators in the same order as the
    hstacked bootstrap matrix.
    """
    with h5.File(path, "r") as f:
        return f["matrix"]["count_total"][:]


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
    keep_barcodes: set[str],
) -> dict[str, ProtospacerCall]:
    """Get the dict of barcode pairs.

    The values are dicts with information about the targets of the features assigned to the barcode.

    Args:
        target_info: dict of target_id to TargetInfo.
        protospacers_per_cell_path: path to the protospacer calls per cell CSV.
        keep_barcodes: only parse rows whose barcode is in this set (the barcodes the caller
            actually needs). On the pre-split MEASURE chunk this is the chunk's ~15k cells, not
            all ~790k, so the parse is byte-identical but far cheaper — every downstream consumer
            only ever reads ``bc_targets_dict[bc]`` for barcodes the caller passes in.

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
            cell_barcode = row[cell_barcode_idx]
            if cell_barcode not in keep_barcodes:
                continue
            feature_call = row[feature_call_idx]
            num_features = int(row[num_features_idx])

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


def get_target(perturbation: str, target_info: dict[str, TargetInfo]) -> str:
    """Return the target gene names of the specified pipe-separated perturbation.

    Deduplicate by gene ID and return the gene names.
    Match the behaviour of `_get_ps_clusters_by_target`.
    """
    gene_id_names = sorted(
        set((target_info[x].gene_id, target_info[x].gene_name) for x in perturbation.split("|"))
    )
    if all(gene_id == NON_TARGETING for gene_id, _ in gene_id_names):
        return NON_TARGETING

    target = "|".join(
        gene_name for gene_id, gene_name in gene_id_names if gene_id not in FILTER_LIST
    )
    return target if target else "Ignore"


def get_perturbations(
    protospacer_calls_per_cell: str,
    target_info: dict[str, TargetInfo],
    by_feature: bool,
) -> list[str]:
    """Return a sorted list of perturbations from the protospacer calls per cell file.

    Translate protospacer IDs to target gene names when by_feature is false.
    """
    with open(protospacer_calls_per_cell) as f:
        feature_calls = (x["feature_call"] for x in csv.DictReader(f))
        perturbation_counts = (
            Counter(feature_calls)
            if by_feature
            else Counter(get_target(x, target_info) for x in feature_calls)
        )
        return sorted(
            perturbation
            for perturbation, cells in perturbation_counts.items()
            if cells >= MIN_NUMBER_CELLS_PER_PERTURBATION
        )


def get_perturbation_efficiency(
    which_perturbations: set[str],
    target_info: dict[str, TargetInfo],
    protospacers_per_cell_path: str,
    combined_bcs: list[bytes],
    gex_feature_ref: FeatureReference,
    control_stats: PerturbationStats,
    chunk_stats: PerturbationStats,
    bootstrap_matrix: CountMatrix,
    combined_col_sums: np_1d_array_f64,
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

    The full per-chunk matrix is gone: the sSeq DE is computed from per-(gene x group) sufficient
    statistics (``control_stats`` shared, ``chunk_stats`` this chunk's perturbations) via
    ``nb_de_from_stats_o3``, and the empirical-CI bootstrap resamples a single target gene out of the
    small ``bootstrap_matrix`` (target-gene rows x ``[control ++ perturbed]`` cells, in storage order
    — the same cell order the old ``hstack([control, slice])`` produced). ``combined_col_sums`` is the
    per-cell count_total over all GEX genes in that same cell order (the bootstrap size-factor
    numerator). Because the streaming reducer is bit-identical to the matrix reduction and the
    bootstrap cell order is preserved, both the DE and the CIs match the matrix path.

    Args:
        which_perturbations: perturbation names this chunk handles.
        target_info: dict of target_id to TargetInfo.
        protospacers_per_cell_path: path to the protospacer calls per cell CSV file.
        combined_bcs: barcodes of the bootstrap matrix cells (control then perturbed, storage order).
        gex_feature_ref: GEX feature reference (for gene_id->row and the result-array gene axis).
        control_stats: shared non-targeting control sufficient statistics.
        chunk_stats: this chunk's perturbations' sufficient statistics (keyed by name).
        bootstrap_matrix: target-gene rows x combined cells, for the empirical-CI bootstrap.
        combined_col_sums: per-cell count_total over all GEX genes, in combined-cell order.
        by_feature: if True, group by protospacer combination rather than gene target.

    Returns:
        perturbation_names (list[str]): List of perturbation_names
        results_all_perturbations (DifferentialExpression): All DE results
        fold_change_per_perturbation (dict[tuple[str, int], dict[str, FoldChange]]): ((perturbation_name, perturbation_idx), dict[target/feature name, FoldChange])
    """
    (target_calls, perturbation_keys) = _get_ps_clusters(
        target_info,
        protospacers_per_cell_path,
        [bc.decode() for bc in combined_bcs],
        by_feature,
    )

    target_to_gene_id = {
        target.target_id if by_feature else target.gene_name: target.gene_id
        for target in target_info.values()
    }
    perturbation_names: list[str] = []

    n_genes = len(gex_feature_ref.feature_defs)
    # Create a numpy array with 3*k columns, where k is the number of perturbations
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((n_genes, 3 * len(which_perturbations)))

    nt_indices = [x for x in perturbation_keys if perturbation_keys[x] == NON_TARGETING]

    if len(nt_indices) == 0:
        return None
    nt_index = nt_indices[0]
    in_control_cluster = target_calls == nt_index
    # Vectorized count instead of Python builtin sum() iterating an n_cells-long bool array.
    control_num_cells = int(in_control_cluster.sum())

    # Early exit if insufficient control cells
    if control_num_cells < MIN_NUMBER_CELLS_PER_PERTURBATION:
        return None

    # gene_id -> matrix row. Keep the FIRST occurrence to match the previous
    # results.loc[results["Gene ID"] == target].values[0] (first-match) semantics, even if a
    # feature reference ever carried a duplicate gene id.
    gene_id_to_row: dict[str, int] = {}
    for row_idx, feature_def in enumerate(gex_feature_ref.feature_defs):
        gene_id_to_row.setdefault(feature_def.id.decode(), row_idx)

    fold_change_per_perturbation: dict[tuple[str, int], dict[str, PerturbationResult]] = {}

    # Control group (Non-Targeting cells) is shared by every perturbation.
    group_b = np.flatnonzero(in_control_cluster)

    # Collect the perturbations that pass the which/filter/min-cell gates, in cluster order so the
    # DE result columns match the matrix path. The control group is shared, so its statistics are
    # reduced once upstream (in SPLIT) instead of once per perturbation.
    valid_perturbations: list[tuple[int, str, np_1d_array_intp]] = []
    for cluster, perturbation_name in perturbation_keys.items():
        if perturbation_name not in which_perturbations or _should_filter(
            perturbation_name, target_to_gene_id
        ):
            continue
        group_a = np.flatnonzero(target_calls == cluster)
        if len(group_a) < MIN_NUMBER_CELLS_PER_PERTURBATION:
            continue
        valid_perturbations.append((cluster, perturbation_name, group_a))

    if not valid_perturbations:
        return (
            perturbation_names,
            DifferentialExpression(all_de_results[:, 0:0]),
            fold_change_per_perturbation,
        )

    sys.stdout.flush()

    # One batched sSeq DE call for all valid perturbations from the sufficient statistics, sharing
    # the (once-reduced) control group. Perturbation stats are gathered by name in valid order.
    ctrl = control_stats.name_to_idx[NON_TARGETING]
    pidx = [chunk_stats.name_to_idx[name] for (_, name, _) in valid_perturbations]
    batched_results = nb_de_from_stats_o3(
        n_genes,
        control_stats.feature_sums[ctrl],
        control_stats.sum_w[ctrl],
        control_stats.sum_wsq[ctrl],
        float(control_stats.sum_ct[ctrl]),
        float(control_stats.sum_inv[ctrl]),
        int(control_stats.n[ctrl]),
        control_stats.count_totals[ctrl],
        chunk_stats.feature_sums[pidx],
        chunk_stats.sum_w[pidx],
        chunk_stats.sum_wsq[pidx],
        chunk_stats.sum_ct[pidx].tolist(),
        chunk_stats.sum_inv[pidx].tolist(),
        chunk_stats.n[pidx].tolist(),
        [chunk_stats.count_totals[i] for i in pidx],
        SSEQ_ZETA_QUANTILE,
        big_count=900,
        threads=NUM_THREADS,
    )

    for cluster_counter, ((cluster, perturbation_name, group_a), de_dict) in enumerate(
        zip(valid_perturbations, batched_results), start=1
    ):
        lga = len(group_a)
        de_arrays = _PerturbationDEArrays(
            log2_fold_change=np.asarray(de_dict["log2_fold_change"]),
            p_values=np.asarray(de_dict["p_values"]),
            sums_in=np.asarray(de_dict["sums_in"]),
            sums_out=np.asarray(de_dict["sums_out"]),
        )
        adjusted_p_values = np.asarray(de_dict["adjusted_p_values"])

        col = 3 * (cluster_counter - 1)
        all_de_results[:, col] = de_arrays.sums_in / lga
        all_de_results[:, col + 1] = de_arrays.log2_fold_change
        all_de_results[:, col + 2] = adjusted_p_values

        # The bootstrap CIs slice a single target gene out of bootstrap_matrix and index
        # group_a/group_b directly. group_a/group_b index combined (control ++ perturbed) cells in
        # storage order, and combined_col_sums[group_x] is each cell's all-GEX count_total in the
        # same order as the matrix path, so the bootstrap draws are byte-identical.
        size_factors = combined_col_sums / de_dict["median_ab"]
        fold_change_per_perturbation[perturbation_name, cluster] = _get_log2_fold_change(
            perturbation_name,
            lga,
            control_num_cells,
            target_to_gene_id,
            bootstrap_matrix,
            group_a,
            group_b,
            size_factors,
            gene_id_to_row,
            de_arrays,
        )
        perturbation_names.append(perturbation_name)

    all_de_results = all_de_results[:, 0 : 3 * len(perturbation_names)]
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
    target_to_gene_id: dict[str, str],
    bootstrap_matrix: CountMatrix,
    group_a: np_1d_array_intp,
    group_b: np_1d_array_intp,
    size_factors: np_1d_array_f64,
    gene_id_to_row: dict[str, int],
    de_arrays: _PerturbationDEArrays,
) -> dict[str, PerturbationResult]:
    (this_names, this_ids) = _get_target_id_from_name(
        perturbation_name,
        target_to_gene_id,
    )

    perturbation_results: dict[str, PerturbationResult] = {}

    for name, target in zip(this_names, this_ids):
        if target in FILTER_LIST:
            continue
        # Direct gene_id -> row lookup into the Rust DE result arrays replaces the previous
        # per-perturbation pandas DataFrame + .loc filter; gene_id_to_row keeps the first match,
        # so a missing target (row is None) is exactly the old empty-DataFrame skip.
        row = gene_id_to_row.get(target)
        if row is None:
            continue

        lower_bound, upper_bound = _get_fold_change_cis(
            bootstrap_matrix,
            target,
            group_a,
            group_b,
            size_factors,
        )
        perturbation_results[name] = PerturbationResult(
            perturbation_name=perturbation_name,
            target_or_gene_name=name,
            log2_fc=de_arrays.log2_fold_change[row],
            p_val=de_arrays.p_values[row],
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            num_cells_with_perb=num_cells,
            mean_umi_count_perb=robust_divide(de_arrays.sums_in[row], num_cells),
            num_cells_with_control=control_num_cells,
            mean_umi_count_control=robust_divide(de_arrays.sums_out[row], control_num_cells),
        )

    return perturbation_results


def _get_fold_change_cis(
    bootstrap_matrix: CountMatrix,
    target: str,
    cond_a: np_1d_array_intp,
    cond_b: np_1d_array_intp,
    size_factors: np_1d_array_f64,
) -> tuple[np.float64, np.float64]:
    np.random.seed(0)

    # filter the matrix to select only the target gene
    this_matrix = bootstrap_matrix.select_features_by_ids([target.encode()])

    # Densify the single target gene's count vector once and pre-gather the per-group
    # count + size-factor arrays. The previous code re-sliced a sparse matrix
    # (`x[:, this_cond_a]`) every bootstrap iteration, which dominates runtime on the
    # 1M-cell case where the control group spans hundreds of thousands of cells.
    # `cond_a`/`cond_b` index the full GEX matrix's cells; `size_factors` is the per-cell
    # size-factor array over all cells (count_total / median_ab). Counts are non-negative
    # integers, so a gathered subset sum is exact regardless of order, and the size-factor
    # sums are over the same float values in the same order as the original per-perturbation
    # path, so the bootstrap CIs are unchanged.
    counts = this_matrix.m.toarray().ravel()
    counts_a = counts[cond_a]
    counts_b = counts[cond_b]
    size_factors_a = size_factors[cond_a]
    size_factors_b = size_factors[cond_b]
    len_a = len(cond_a)
    len_b = len(cond_b)

    # Use numpy array instead of list for better performance
    log2_fold_change_vals = np.zeros(NUM_BOOTSTRAPS)

    # Early convergence detection
    for i in range(NUM_BOOTSTRAPS):
        # `np.random.choice(cond_x, size=len_x, replace=True)` is exactly
        # `cond_x.take(np.random.randint(0, len_x, len_x))`, so drawing the index arrays
        # directly keeps the RNG stream (seeded above) bit-for-bit identical while letting
        # us index the small pre-gathered group arrays instead of the full matrix.
        idx_a = np.random.randint(0, len_a, size=len_a)
        idx_b = np.random.randint(0, len_b, size=len_b)

        gene_sums_a = counts_a[idx_a].sum().item()
        gene_sums_b = counts_b[idx_b].sum().item()
        size_factor_a = size_factors_a[idx_a].sum()
        size_factor_b = size_factors_b[idx_b].sum()

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
    keep_barcodes = set(barcodes)
    bc_targets_dict = _get_bc_targets_dict(target_info, protospacers_per_cell_path, keep_barcodes)

    # Optimize: Use set operations for faster lookups
    missing_barcodes = keep_barcodes - set(bc_targets_dict.keys())
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
