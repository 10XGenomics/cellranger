#!/usr/bin/env python3
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import csv
import json
import os
import sys
from functools import reduce
from typing import TYPE_CHECKING, TypedDict

import numpy as np
import scipy.stats
from six import ensure_binary, ensure_str

import cellranger.analysis.constants as analysis_constants
import cellranger.library_constants as lib_constants
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats

if TYPE_CHECKING:
    from cellranger.matrix import CountMatrix

# Default number of bootstrap samples to use
NUM_MULTIPLET_BOOTSTRAP_SAMPLES = 1000


def compute_count_purity(
    counts0: np.ndarray,
    counts1: np.ndarray,
    classifications: np.ndarray[int, np.dtype[np.bytes_]] | None = None,
):
    """Compute fraction of counts in putative single-cell GEMs.

    originating from the non-cell transcriptome
    """
    if classifications is None:
        classifications = classify_gems(counts0, counts1)
    frac0 = counts0.astype(float) / (counts0 + counts1).astype(float)
    purity0 = frac0[classifications == analysis_constants.GEM_CLASS_GENOME0]
    purity1 = 1 - frac0[classifications == analysis_constants.GEM_CLASS_GENOME1]

    # Compute number of purity outliers.
    # Note: This is not used for any important metrics, just for diagnostics.
    threshold0, threshold1 = 1.0, 1.0
    fit_purity0 = purity0[np.logical_and(purity0 > 0, purity0 < 1)]
    fit_purity1 = purity1[np.logical_and(purity1 > 0, purity1 < 1)]
    if len(fit_purity0) > 1 and len(fit_purity1) > 1:
        try:
            alpha0, beta0, _, _ = scipy.stats.beta.fit(fit_purity0, floc=0, fscale=1)
            alpha1, beta1, _, _ = scipy.stats.beta.fit(fit_purity1, floc=0, fscale=1)
            threshold0 = scipy.stats.beta.ppf(
                analysis_constants.COUNT_PURITY_OUTLIER_PROB_THRESHOLD, alpha0, beta0
            )
            threshold1 = scipy.stats.beta.ppf(
                analysis_constants.COUNT_PURITY_OUTLIER_PROB_THRESHOLD, alpha1, beta1
            )
        except scipy.stats._continuous_distns.FitSolverError as e:
            print(e, file=sys.stderr)
            threshold0, threshold1 = 1.0, 1.0
        except scipy.stats._continuous_distns.FitDataError as e:
            print(e, file=sys.stderr)
            threshold0, threshold1 = 1.0, 1.0

    outlier0 = np.logical_and(
        classifications == analysis_constants.GEM_CLASS_GENOME0, frac0 < threshold0
    )
    outlier1 = np.logical_and(
        classifications == analysis_constants.GEM_CLASS_GENOME1, (1 - frac0) < threshold1
    )
    n_outlier0: int = sum(outlier0)
    n_outlier1: int = sum(outlier1)
    frac_outlier0 = tk_stats.robust_divide(n_outlier0, len(purity0))
    frac_outlier1 = tk_stats.robust_divide(n_outlier1, len(purity1))
    is_outlier: np.ndarray[int, np.dtype[np.int_]] = np.logical_or(outlier0, outlier1).astype(int)

    # Let the UMI count purity be the total fraction of counts that don't belong,
    #   for all barcodes classified as non-multiplets.
    gems0 = classifications == analysis_constants.GEM_CLASS_GENOME0
    mean_purity0 = tk_stats.robust_divide(
        counts0[gems0].sum(), (counts0[gems0] + counts1[gems0]).sum()
    )

    gems1 = classifications == analysis_constants.GEM_CLASS_GENOME1
    mean_purity1 = tk_stats.robust_divide(
        counts1[gems1].sum(), (counts0[gems1] + counts1[gems1]).sum()
    )

    single_cell = (classifications == analysis_constants.GEM_CLASS_GENOME0) | (
        classifications == analysis_constants.GEM_CLASS_GENOME1
    )
    mean_overall_purity = tk_stats.robust_divide(
        np.maximum(counts0[single_cell], counts1[single_cell]).sum(),
        (counts0 + counts1)[single_cell].sum(),
    )

    return (
        mean_purity0,
        mean_purity1,
        mean_overall_purity,
        n_outlier0,
        n_outlier1,
        frac_outlier0,
        frac_outlier1,
        is_outlier,
        classifications,
    )


def infer_multiplets_from_observed(n_obs_multiplets, n_cells0, n_cells1):
    """Given a number of observed multiplets and cell counts for two transcriptomes,.

    infer the total number of multiplets (observed + unobserved)
    """
    if n_cells0 == 0 or n_cells1 == 0:
        return 0

    # Prior probability of a doublet given counts for each cell type (ignore N_cells > 2)
    p_obs_multiplet = (
        2
        * (float(n_cells0) / float(n_cells0 + n_cells1))
        * (float(n_cells1) / float(n_cells0 + n_cells1))
    )

    # Analytical MOM/MLE of binomial N given p, k
    mle = float(n_obs_multiplets) / p_obs_multiplet
    # In some (artificial) datasets, the mle can be higher than the total number of cells
    # observed.  This occurs when n_obs_multiplets > n_cells0|1. The right way to fix that would be to
    # do inference in a full model that didn't fix some parameters.  In practice, multigenomes
    # are a rare analysis and most data isn't artificial, so we are implementing
    # a small hack instead.
    return min(mle, float(n_obs_multiplets + n_cells0 + n_cells1))


def classify_gems(
    counts0: np.ndarray[int, np.dtype[np.int_]], counts1: np.ndarray[int, np.dtype[np.int_]]
) -> np.ndarray[int, np.dtype[np.bytes_]]:
    """Classify counts by inferred number of distinct transcriptomes present in each GEM (1 or 2).

    Report analysis_constants.GEM_CLASS_GENOME0 for a single cell w/ transcriptome 0,
    report analysis_constants.GEM_CLASS_GENOME1 for a single cell w/ transcriptome 1,
    report analysis_constants.GEM_CLASS_MULTIPLET for multiple transcriptomes.
    """
    # Assumes that most of the GEMs are single-cell; model counts independently
    thresh0, thresh1 = [analysis_constants.DEFAULT_MULTIPLET_THRESHOLD] * 2
    if sum(counts0 > counts1) >= 1 and sum(counts1 > counts0) >= 1:
        thresh0 = np.percentile(
            counts0[counts0 > counts1], analysis_constants.MULTIPLET_PROB_THRESHOLD * 100.0
        )
        thresh1 = np.percentile(
            counts1[counts1 > counts0], analysis_constants.MULTIPLET_PROB_THRESHOLD * 100.0
        )

    # If input is a pure species instead of a mixture then modeling counts independently
    # can result in an absurdly low threshold for the missing species causing FP labels that
    # show up as an inflated number of Multiplets.
    thresholds = sorted([thresh0, thresh1])
    fold_change = thresholds[1] / thresholds[0]
    if (thresholds[0] < 50) and (fold_change > 25):
        thresh0 = thresh1 = np.percentile(
            counts0 + counts1, analysis_constants.MULTIPLET_PROB_THRESHOLD * 100.0
        )
    doublet = np.logical_and(counts0 >= thresh0, counts1 >= thresh1)
    dtype = np.dtype(("S", max(len(cls) for cls in analysis_constants.GEM_CLASSES)))
    result = np.where(
        doublet, analysis_constants.GEM_CLASS_MULTIPLET, analysis_constants.GEM_CLASS_GENOME0
    ).astype(dtype)
    result[
        np.logical_and(
            np.logical_not(result == analysis_constants.GEM_CLASS_MULTIPLET), counts1 > counts0
        )
    ] = analysis_constants.GEM_CLASS_GENOME1

    return result


class MultiGenomeAnalysisResult(TypedDict, total=False):
    barcode: list[bytes]
    call: list[bytes]
    count0: list[int]
    count1: list[int]
    genome0: str
    genome1: str
    purity_outlier: list[int]


class MultiGenomeAnalysis:
    """Analysis of matrices when >1 genome is present."""

    def __init__(self, filtered_matrix: CountMatrix | None = None):
        self.filtered_matrix = filtered_matrix

        self.summary = {}
        self.result: MultiGenomeAnalysisResult = {}
        self.top_two_txomes = None
        self.suffix = ""
        self.n_gems = None

    def is_zero_matrix(self):
        if self.filtered_matrix is None:
            return True
        if self.filtered_matrix.bcs_dim == 0 or self.filtered_matrix.features_dim == 0:
            return True
        return False

    def _infer_multiplets(
        self,
        counts0: np.ndarray[int, np.dtype[np.int_]],
        counts1: np.ndarray[int, np.dtype[np.int_]],
        bootstraps: int = NUM_MULTIPLET_BOOTSTRAP_SAMPLES,
        **kwargs,
    ):
        """An overridable method to determine the number of multiplets, should return.

        Args:
            counts0:
            counts1:
            bootstraps: Number of bootstrap iterations to run

        Returns:
            int: observed multiplets
            np.ndarray: the number of inferred multiplets,
                either as a single value, or as a list/vector of values if bootstrapping is used.
            np.ndarray: strings with elements of either GEM_CLASS_GENOME0,
                GEM_CLASS_GENOME1 or GEM_CLASS_MULTIPLET
        """
        classifications = classify_gems(counts0, counts1)
        n_multiplet_obs: int = sum(classifications == analysis_constants.GEM_CLASS_MULTIPLET)
        assert bootstraps > 0
        assert len(counts0) == len(counts1)
        # Fix random seed
        np.random.seed(0)

        n_multiplet_boot: np.ndarray[int, np.dtype[np.float64]] = np.zeros(bootstraps)
        for i in range(bootstraps):
            boot_idx = np.random.choice(len(counts0), len(counts0))
            counts0_boot = counts0[boot_idx]
            counts1_boot = counts1[boot_idx]
            gem_cls_boot = classify_gems(counts0_boot, counts1_boot)
            n_obs_multiplet_boot = sum(gem_cls_boot == analysis_constants.GEM_CLASS_MULTIPLET)
            n_cells0_boot = sum(gem_cls_boot == analysis_constants.GEM_CLASS_GENOME0)
            n_cells1_boot = sum(gem_cls_boot == analysis_constants.GEM_CLASS_GENOME1)
            n_multiplet_boot[i] = infer_multiplets_from_observed(
                n_obs_multiplet_boot, n_cells0_boot, n_cells1_boot
            )
        return n_multiplet_obs, n_multiplet_boot, classifications

    def run_all(self):
        d = {}

        # Compute N_cells > 1 between top two genomes
        assert self.filtered_matrix is not None
        genomes = self.filtered_matrix.get_genomes()
        genome_mats = [self.filtered_matrix.select_features_by_genome(g) for g in genomes]

        txome_counts = [mat.m.sum() for mat in genome_mats]
        top_txome_idx = sorted(np.argsort(txome_counts)[::-1][0:2])

        top_two_txomes = [genomes[i] for i in top_txome_idx]
        self.top_two_txomes = top_two_txomes

        top_txome_cell_bc_seqs = [genome_mats[i].bcs for i in top_txome_idx]
        use_barcodes: list[bytes] = sorted(
            reduce(lambda a, x: a | set(x), top_txome_cell_bc_seqs, set())
        )
        n_barcodes = len(use_barcodes)

        # Don't compute multiplet / purity metrics if no cells detected
        if n_barcodes == 0:
            return

        # Sum the filtered matrices by barcode and stack the two genomes together
        top_two_filt_mats = [genome_mats[i] for i in top_txome_idx]
        top_txome_reads_per_bc: np.ndarray[int, np.dtype[np.int_]] = np.vstack(
            tuple(
                m.select_barcodes_by_seq(use_barcodes).get_counts_per_bc()
                for m in top_two_filt_mats
            )
        )

        n_multiplet_obs, n_multiplet_boot, gem_class_call = self._infer_multiplets(
            top_txome_reads_per_bc[0, :], top_txome_reads_per_bc[1, :], n_gems=self.n_gems
        )
        d["filtered_bcs_observed_all"] = n_barcodes
        d["filtered_bcs_observed_multiplets"] = int(round(n_multiplet_obs))
        d["filtered_bcs_inferred_multiplets"] = int(round(n_multiplet_boot.mean()))
        multiplet_rate = tk_stats.robust_divide(n_multiplet_boot.mean(), len(use_barcodes))
        d["filtered_bcs_inferred_multiplet_rate"] = multiplet_rate
        d["filtered_bcs_inferred_normalized_multiplet_rate"] = 1000 * tk_stats.robust_divide(
            multiplet_rate, len(use_barcodes)
        )
        if n_multiplet_boot.size > 1:
            d["filtered_bcs_inferred_multiplet_rate_lb"] = tk_stats.robust_divide(
                np.percentile(n_multiplet_boot, 2.5), len(use_barcodes)
            )
            d["filtered_bcs_inferred_multiplet_rate_ub"] = tk_stats.robust_divide(
                np.percentile(n_multiplet_boot, 97.5), len(use_barcodes)
            )

        (
            purity0,
            purity1,
            overall_purity,
            n_purity_outlier0,
            n_purity_outlier1,
            frac_purity_outlier0,
            frac_purity_outlier1,
            is_outlier,
            _,
        ) = compute_count_purity(
            top_txome_reads_per_bc[0, :], top_txome_reads_per_bc[1, :], gem_class_call
        )
        d[f"{top_two_txomes[0]}_filtered_bcs_mean_count_purity"] = purity0
        d[f"{top_two_txomes[1]}_filtered_bcs_mean_count_purity"] = purity1
        d[f"{lib_constants.MULTI_REFS_PREFIX}_filtered_bcs_mean_count_purity"] = overall_purity
        d[f"{top_two_txomes[0]}_filtered_bcs_purity_outliers"] = n_purity_outlier0
        d[f"{top_two_txomes[1]}_filtered_bcs_purity_outliers"] = n_purity_outlier1
        d[f"{top_two_txomes[0]}_filtered_bcs_frac_purity_outlier"] = frac_purity_outlier0
        d[f"{top_two_txomes[1]}_filtered_bcs_frac_purity_outlier"] = frac_purity_outlier1
        d[f"{lib_constants.MULTI_REFS_PREFIX}_filtered_bcs_frac_purity_outlier"] = (
            frac_purity_outlier0 + frac_purity_outlier1
        )
        self.result = {
            "barcode": use_barcodes,
            "call": gem_class_call.tolist(),
            "count0": top_txome_reads_per_bc[0, :].tolist(),
            "count1": top_txome_reads_per_bc[1, :].tolist(),
            "genome0": top_two_txomes[0],
            "genome1": top_two_txomes[1],
            "purity_outlier": is_outlier.tolist(),
        }
        self.summary = d
        self._add_suffix_to_metrics()

    def _add_suffix_to_metrics(self):
        """Method to update the result and suffix to add a suffix if need be."""
        if self.suffix == "":
            return

        sum_keys = list(self.summary.keys())
        for key in sum_keys:
            self.summary[key + self.suffix] = self.summary.pop(key)
        res_keys = list(self.result.keys())
        for key in res_keys:
            self.result[key + self.suffix] = self.result.pop(key)

    def _get_wo_suffix(self, column: str):
        """Get a value from the result with the suffix removed."""
        if self.suffix != "":
            column = column + self.suffix
        return self.result[column]

    def save_gem_class_csv(self, base_dir: str):
        csv_file_path = os.path.join(base_dir, "gem_classification.csv")
        os.makedirs(os.path.dirname(csv_file_path), exist_ok=True)
        with open(csv_file_path, "w") as f:
            writer = csv.writer(f, lineterminator=os.linesep)
            writer.writerow(
                ["barcode", self._get_wo_suffix("genome0"), self._get_wo_suffix("genome1"), "call"]
            )
            genome0 = ensure_str(analysis_constants.GEM_CLASS_GENOME0)
            genome1 = ensure_str(analysis_constants.GEM_CLASS_GENOME1)
            for i in range(len(self._get_wo_suffix("barcode"))):
                call = ensure_str(self._get_wo_suffix("call")[i])

                call = call.replace(
                    genome0,
                    self._get_wo_suffix("genome0"),
                )
                call = call.replace(
                    genome1,
                    self._get_wo_suffix("genome1"),
                )
                writer.writerow(
                    [
                        ensure_str(self._get_wo_suffix("barcode")[i]),
                        self._get_wo_suffix("count0")[i],
                        self._get_wo_suffix("count1")[i],
                        call,
                    ]
                )

    def save_gem_class_json(self, base_dir: str):
        json_file_path = MultiGenomeAnalysis.json_path(base_dir)
        os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
        with open(json_file_path, "w") as f:
            tk_safe_json.dump_numpy(self.result, f, indent=4, sort_keys=True)

    def save_summary_json(self, filename):
        with open(filename, "w") as f:
            tk_safe_json.dump_numpy(self.summary, f, indent=4, sort_keys=True)

    @staticmethod
    def load_default_format(base_dir: str):
        json_file_path = MultiGenomeAnalysis.json_path(base_dir)
        if os.path.exists(json_file_path):
            return MultiGenomeAnalysis.load_json(json_file_path)
        else:
            return None

    @staticmethod
    def load_json(filename):
        analysis = MultiGenomeAnalysis()
        with open(filename) as f:
            analysis.result = json.load(f)
        if "call" in analysis.result:
            # Convert calls to bytes.
            # Most of these will be one of just a few strings, so avoid
            # duplicating the str/byte conversions and also duplicating in
            # memory.
            # We'd use six.moves.intern=sys.intern here except this is bytes.
            bytes_intern = {}
            for call in analysis.result["call"]:
                if call not in bytes_intern:
                    bytes_intern[call] = ensure_binary(call)
            analysis.result["call"] = [bytes_intern[call] for call in analysis.result["call"]]
        return analysis

    @staticmethod
    def json_path(base_path: str):
        return os.path.join(base_path, "analysis.json")
