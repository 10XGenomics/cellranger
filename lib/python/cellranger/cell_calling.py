#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Functions for calling cell-associated barcodes."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import numpy.ma as ma

import cellranger.sgt as cr_sgt
import cellranger.stats as cr_stats
from cellranger.analysis.diffexp import adjust_pvalue_bh
from cellranger.chemistry import (
    CHEMISTRY_DESCRIPTION_FIELD,
    CHEMISTRY_SC3P_LT,
    HT_CHEMISTRIES,
    SC3P_V3_CHEMISTRIES,
    SC3P_V4_CHEMISTRIES,
    SC5P_CHEMISTRIES,
    SC5P_V3_CHEMISTRIES,
)

# Minimum number of UMIs for a barcode to be called as a cell
MIN_GLOBAL_UMIS = 0

# Maximum percentage of mitochondrial reads allowed for a barcode to be called a cell
MAX_MITO_PCT = 100.0

# Minimum number of UMIS per barcode to consider after the initial cell calling
MIN_UMIS = 500

# Default number of background simulations to make
NUM_SIMS = 100000

# Minimum number of UMIS per barcode to consider after the initial cell calling for targeted GEX
TARGETED_CC_MIN_UMIS_ADDITIONAL_CELLS = 10

# Minimum number of target gene UMIs per barcode required to be called a cell for targeted GEX
TARGETED_CC_MIN_UMIS_FROM_TARGET_GENES = 1


def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
    """Estimate a gene expression profile by Simple Good Turing.

    Args:
      raw_mat (sparse matrix): Sparse matrix of all counts
      barcode_indices (np.array(int)): Barcode indices to use
      nz_feat (np.array(int)): Indices of features that are non-zero at least once

    Returns:
      profile (np.array(float)): Estimated probabilities of length len(nz_feat).
    """
    # Initial profile estimate
    prof_mat = matrix[:, barcode_indices]

    profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
    zero_feat = np.flatnonzero(profile == 0)

    # Simple Good Turing estimate
    p_smoothed, p0 = cr_sgt.sgt_proportions(profile[np.flatnonzero(profile)])

    n0 = len(zero_feat)
    if n0 == 0:
        # Renormalize in absence of 0 class
        p_smoothed = p_smoothed / p_smoothed.sum()
        p0_i = -1.0  # unused
    else:
        # Distribute p0 equally among the zero elements.
        p0_i = p0 / n0

    profile_p = np.repeat(p0_i, len(nz_feat))
    profile_p[np.flatnonzero(profile)] = p_smoothed

    assert np.isclose(profile_p.sum(), 1.0)
    return profile_p


# Construct a background expression profile from barcodes with <= T UMIs
def est_background_profile_sgt(matrix, use_bcs):
    """Estimate a gene expression profile on a given subset of barcodes.

    Use Good-Turing to smooth the estimated profile.

    Args:
      matrix (scipy.sparse.csc_matrix): Sparse matrix of all counts
      use_bcs (np.array(int)): Indices of barcodes to use (col indices into matrix)

    Returns:
      profile (use_features, np.array(float)): Estimated probabilities of length use_features.
    """
    # Use features that are nonzero anywhere in the data
    use_feats = np.flatnonzero(np.asarray(matrix.sum(1)))

    # Estimate background profile
    bg_profile_p = estimate_profile_sgt(matrix, use_bcs, use_feats)

    return (use_feats, bg_profile_p)


class NonAmbientBarcodeResult(NamedTuple):
    eval_bcs: np.ndarray  # Candidate barcode indices (n)
    log_likelihood: np.ndarray  # Ambient log likelihoods (n)
    pvalues: np.ndarray  # pvalues (n)
    pvalues_adj: np.ndarray  # B-H adjusted pvalues (n)
    is_nonambient: np.ndarray  # Boolean nonambient calls (n)
    emptydrops_minimum_umis: int  # Min UMI threshold for empty drops (1)


def get_empty_drops_fdr(chemistry_description: str) -> float:
    """Gets the maximum adjusted p-value to call a barcode as non-ambient."""
    # The chips used with V4 have roughly double the GEMs as the older V3 chips
    chemistries = SC3P_V4_CHEMISTRIES + SC3P_V3_CHEMISTRIES + SC5P_CHEMISTRIES + HT_CHEMISTRIES
    chem_names = [chem[CHEMISTRY_DESCRIPTION_FIELD] for chem in chemistries]
    return 0.001 if chemistry_description in chem_names else 0.01


def get_empty_drops_range(chemistry_description: str, num_probe_bcs: int | None) -> tuple[int, int]:
    """Gets the range of values to use for empty drops background given a chemistry description.

    Args:
        chemistry_description: A string describing the chemistry
        num_probe_bcs: The number of probe or OCM multiplexing barcodes

    Returns:
        (lower_range, upper_range)
    """
    if chemistry_description == CHEMISTRY_SC3P_LT[CHEMISTRY_DESCRIPTION_FIELD]:
        n_partitions = 9000
    elif chemistry_description in [
        chem[CHEMISTRY_DESCRIPTION_FIELD]
        for chem in HT_CHEMISTRIES + SC3P_V4_CHEMISTRIES + SC5P_V3_CHEMISTRIES
    ]:
        if num_probe_bcs is None:
            n_partitions = 160000
        else:
            # OCM
            n_partitions = 40000 * num_probe_bcs
    elif num_probe_bcs is not None:
        # Flex
        n_partitions = 90000 * num_probe_bcs
    else:
        n_partitions = 90000
    return (n_partitions // 2, n_partitions)


def find_nonambient_barcodes(
    matrix,
    orig_cell_bcs,
    chemistry_description,
    num_probe_bcs,
    *,
    emptydrops_minimum_umis=MIN_UMIS,
    num_sims=NUM_SIMS,
    method="multinomial",
) -> NonAmbientBarcodeResult | None:
    """Call barcodes as being sufficiently distinct from the ambient profile.

    Args:
      matrix (CountMatrix): Full expression matrix.
      orig_cell_bcs (iterable of str): Strings of initially-called cell barcodes.
      chemistry_description: Change ambient RNA estimation for LT chemistry
      num_probe_bcs: The number of probe or OCM multiplexing barcodes
      emptydrops_minimum_umis: Minimum UMI threshold
      num_sims: Number of simulations

    Returns:
      NonAmbientBarcodeResult
    """
    assert method in ["dirichlet", "multinomial"]

    # Estimate an ambient RNA profile
    umis_per_bc = matrix.get_counts_per_bc()
    bc_order = np.argsort(umis_per_bc)

    low, high = get_empty_drops_range(chemistry_description, num_probe_bcs)
    max_adj_pvalue = get_empty_drops_fdr(chemistry_description)

    # Take what we expect to be the barcodes associated w/ empty partitions.
    print(f"Range empty barcodes: {low} - {high}")
    empty_bcs = bc_order[::-1][low:high]
    empty_bcs.sort()

    # Require non-zero barcodes
    nz_bcs = np.flatnonzero(umis_per_bc)
    nz_bcs.sort()

    ambient_bcs = np.intersect1d(empty_bcs, nz_bcs, assume_unique=True)

    if len(ambient_bcs) > 0:
        try:
            eval_features, ambient_profile_p = est_background_profile_sgt(matrix.m, ambient_bcs)
            if method == "dirichlet":
                alpha = cr_stats.estimate_dirichlet_overdispersion(
                    matrix.m, ambient_bcs, ambient_profile_p
                )

        except cr_sgt.SimpleGoodTuringError as e:
            print(str(e))
            return None
    else:
        eval_features = np.zeros(0, dtype=np.int64)
        ambient_profile_p = np.zeros(0)

    # Choose candidate cell barcodes
    orig_cell_bc_set = set(orig_cell_bcs)
    orig_cells = np.flatnonzero(
        np.fromiter(
            (bc in orig_cell_bc_set for bc in matrix.bcs), count=len(matrix.bcs), dtype=bool
        )
    )

    # No good incoming cell calls
    if orig_cells.sum() == 0:
        return None

    # Look at non-cell barcodes above a minimum UMI count
    eval_bcs = np.ma.array(np.arange(matrix.bcs_dim))
    eval_bcs[orig_cells] = ma.masked

    max_background_umis = np.max(umis_per_bc[empty_bcs], initial=0)
    emptydrops_minimum_umis = max(emptydrops_minimum_umis, 1 + max_background_umis)
    print(f"Max background UMIs: {max_background_umis}")

    eval_bcs[umis_per_bc < emptydrops_minimum_umis] = ma.masked
    n_unmasked_bcs = len(eval_bcs) - eval_bcs.mask.sum()

    eval_bcs = np.argsort(ma.masked_array(umis_per_bc, mask=eval_bcs.mask))[0:n_unmasked_bcs]
    # SORT the barcodes. This is a critical step; eval_bcs is a list of integers, and these indices are used
    # to get the counts via matrix.select_features_by_genome(genome).select_barcodes(eval_bcs).get_counts_per_bc(),
    # which sorts the input, and also to get the string barcode sequences via np.array(matrix.bcs)[eval_bcs],
    # which doesn't. Without sorting here, the matching between the sequences and their counts from both
    # genomes is essentially random, thus the species assignments will be wrong.
    eval_bcs.sort()

    if len(eval_bcs) == 0:
        return None

    assert not np.any(np.isin(eval_bcs, orig_cells))
    assert not np.any(np.isin(eval_bcs, empty_bcs))
    print(f"Number of candidate bcs: {len(eval_bcs)}")
    print(f"Range candidate bc umis: {umis_per_bc[eval_bcs].min()}, {umis_per_bc[eval_bcs].max()}")
    print(f"Number of empty bcs: {len(empty_bcs)}")
    print(f"Number of original cell calls: {len(orig_cells)}")

    eval_mat = matrix.m[eval_features, :][:, eval_bcs]

    if len(ambient_profile_p) == 0:
        obs_loglk = np.repeat(np.nan, len(eval_bcs))
        pvalues = np.repeat(1, len(eval_bcs))
        sim_loglk = np.repeat(np.nan, len(eval_bcs))
        return None

    # Compute observed log-likelihood of barcodes being generated from ambient RNA and
    # simulate log-likelihoods
    if method == "dirichlet":
        obs_loglk = cr_stats.eval_dirichlet_multinomial_loglikelihoods(
            eval_mat, alpha * ambient_profile_p
        )
        distinct_ns, sim_loglk = cr_stats.simulate_dirichlet_multinomial_loglikelihoods(
            alpha * ambient_profile_p, umis_per_bc[eval_bcs], num_sims=num_sims
        )
    else:
        obs_loglk = cr_stats.eval_multinomial_loglikelihoods(eval_mat, np.log(ambient_profile_p))
        distinct_ns, sim_loglk = cr_stats.simulate_multinomial_loglikelihoods(
            ambient_profile_p, umis_per_bc[eval_bcs], num_sims=num_sims
        )

    # Compute p-values
    pvalues = cr_stats.compute_ambient_pvalues(
        umis_per_bc[eval_bcs], obs_loglk, distinct_ns, sim_loglk
    )

    pvalues_adj = adjust_pvalue_bh(pvalues)

    print(f"Max adjusted P-value: {max_adj_pvalue}")
    print(f"Min observed P-value: {min(pvalues_adj)}")
    is_nonambient = pvalues_adj <= max_adj_pvalue

    print(f"Non-ambient bcs identified by empty drops: {sum(is_nonambient)}")

    return NonAmbientBarcodeResult(
        eval_bcs=eval_bcs,
        log_likelihood=obs_loglk,
        pvalues=pvalues,
        pvalues_adj=pvalues_adj,
        is_nonambient=is_nonambient,
        emptydrops_minimum_umis=emptydrops_minimum_umis,
    )
