#!/usr/bin/env python3
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Helper statistical methods for cell calling."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.special import gammaln

if TYPE_CHECKING:

    from cellranger.matrix import CountMatrixView

# Pre-seed the RNG on package import
RNG = np.random.default_rng(seed=42)


def _incremental_counts_from_sample(
    sample_draws: np.ndarray[tuple[int], np.dtype[np.intp]],
) -> np.ndarray[tuple[int], np.dtype[np.int_]]:
    """Produces an incremental count vector from a categorical sample vector.

    Takes a set of drawn samples from a multinomial or similar distribution
    distribution and produces a vector of the same size, with each element the
    the number of times the sample at that position has been seen
    indices.  I.e., it would turn the following vector of sample draws:
    [1, 1, 2, 1, 0, 2, 1, 3] into
    [1, 2, 1, 3, 1, 2, 4, 1]

    Args:
      sample_draws (np.ndarray(int)): A vector of sample indices.

    Returns:
      inc_counts (np.ndarray(int)): The incremental counts of those sample indices.

    """
    # Vectorized implementation: for each position i, the rank of sample_draws[i]
    # among sample_draws[0..i] (1-indexed). Stable argsort groups equal values in
    # original positional order; within each group, positions take ranks 1, 2, ...
    # Replaces a Python dict + loop. Byte-identical output to the prior impl
    # (verified by smoke test pvalue equality on CC-230 / CC-510).
    n = sample_draws.shape[0]
    if n == 0:
        return np.empty(0, dtype=sample_draws.dtype)
    order = np.argsort(sample_draws, kind="stable")
    sorted_vals = sample_draws[order]
    arange_n = np.arange(n)
    is_new_group = np.empty(n, dtype=bool)
    is_new_group[0] = True
    np.not_equal(sorted_vals[1:], sorted_vals[:-1], out=is_new_group[1:])
    # group_start[i] = position of the first element in i's value group (in sorted order).
    group_start = np.maximum.accumulate(np.where(is_new_group, arange_n, -1))
    rank_within = arange_n - group_start + 1
    inc_counts = np.empty(n, dtype=sample_draws.dtype)
    inc_counts[order] = rank_within
    return inc_counts


def eval_multinomial_loglikelihoods(
    matrix: CountMatrixView,
    logp: np.ndarray[tuple[int], np.dtype[np.float64]],
    total_umis: np.ndarray[tuple[int], np.dtype[np.int_]] | None = None,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Computes the multinomial log-likelihood for many barcodes.

    Multinomial log-likelihood for a single barcode where the count of UMIs for
    feature i is x_i is:
    l = log(gamma(sum_i(x_i) + 1)) - sum_i(log(gamma(x_i + 1)) + sum_i(log(p_i) * x_i)

    Because the input matrix is sparse and zero counts do not contribute to the
    likelihood, we can rapidly remove zero elements from the computation.

    Args:
      matrix (CountMatrixView): Matrix of UMI counts (feature x barcode)
      logp (np.ndarray(float)): The natural log of the multinomial probability
        vector across features
      total_umis (np.ndarray(int), optional): The total number of UMIs per barcode.  Saves computation
        if it can be precomputed.

    Returns:
      loglk (np.ndarray(float)): Log-likelihood for each barcode
    """
    loglk: np.ndarray[tuple[int], np.dtype[np.float64]] = np.zeros(matrix.bcs_dim, dtype=float)

    if total_umis is None:
        total_umis = matrix.sum(axis=0)
    assert total_umis is not None
    consts = gammaln(total_umis + 1)
    for bc_idx, (bc_feature_cnts, bc_feature_idxs) in enumerate(
        matrix.get_per_barcode_counts_and_features()
    ):
        short_logp = logp.take(bc_feature_idxs, axis=None, mode="clip")
        loglk[bc_idx] = (
            consts[bc_idx]
            - gammaln(bc_feature_cnts + 1).sum()
            + (bc_feature_cnts * short_logp).sum()
        )
    return loglk


def _eval_multinomial_loglikelihood_cumulative(
    sample_draws: np.ndarray[tuple[int], np.dtype[np.intp]],
    logp: np.ndarray[tuple[int], np.dtype[np.float64]],
    log_nvals: np.ndarray[tuple[int], np.dtype[np.float64]] | None = None,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Computes the cumulative multinomial log-likelihood for a vector.

    Given a vector of sample draws from a multinomial distribution, computes the
    log-likelihood of all of the draws up to each draw.  This is done incrementally,
    by first pre-computing the rank R of each draw as the number of times that
    feature was seen when it was drawn.  Then the incremental log-likelihood from
    that draw is

    dl = log(i) - log(R_i) + log(p_i)

    Args:
      sample_draws (np.ndarray(int)): A vector of sample draws
      logp (np.ndarray(float)): The log-probability of sampling each feature
      log_nvals: Optional precomputed `np.log(np.arange(1, n + 1))` of length at
        least `len(sample_draws)`. `compute_multinomial_ambient_pvalues` hoists
        this constant outside its `num_sims` simulation loop to skip ~num_sims
        recomputations.

    Returns:
      loglk (np.ndarray(float): The cumulative log-likelihood up to each draw
    """
    n = sample_draws.shape[0]
    marginal_counts = _incremental_counts_from_sample(sample_draws)
    if log_nvals is None or log_nvals.shape[0] < n:
        log_nvals = np.log(np.arange(1, n + 1))
    else:
        log_nvals = log_nvals[:n]
    loglk = log_nvals - np.log(marginal_counts) + logp[sample_draws]
    return np.cumsum(loglk)


def _draw_multinomial_sample(
    num_draws: int,
    p_cumulative: np.ndarray[tuple[int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.intp]]:
    """Returns a fixed size array of multinomial draws from the input probabilities.

    Given a number of samples to draw and a probability vector for all features
    to sample from, produces an array of feature indices of the desired size,
    drawn according to a multinomial distribution.  Because it is significantly
    faster to generate random numbers on the unit interval, we use that and then
    find the feature index using the cumulative probabilities of the probability
    vector.

    Args:
      num_draws (int): The number of samples to draw from the distribution
      p_cumulative (np.ndarray(float)): The cumulative probability of all possible
        features.  The length of this array should be the number of features.
        It should be sorted from smallest to largest, and the last entry should
        be 1.

    Returns:
      sample_draws (np.ndarray(intp)): A random set of draws from a multinomial
        distribution.
    """
    rng_nums = RNG.random(size=num_draws)
    result = np.searchsorted(p_cumulative, rng_nums)
    return result


def compute_multinomial_ambient_pvalues(
    ambient_profile_p: np.ndarray[tuple[int], np.dtype[np.float64]],
    umis_per_bc: np.ndarray[tuple[int], np.dtype[np.int_]],
    obs_loglk: np.ndarray[tuple[int], np.dtype[np.float64]],
    num_sims: int,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Compute p-values for observed multinomial log-likelihoods.

    Args:
      ambient_profile_p (nd.array(float)): Probability of observing each feature.
      umis_per_bc (nd.array(int)): UMI counts per barcode
      obs_loglk (nd.array(float)): Observed log-likelihoods of each barcode deriving from an ambient profile
      num_sims (int): Number of simulations to perform.

    Returns:
      pvalues (nd.array(float)): p-values
    """
    distinct_ns = np.flatnonzero(np.bincount(umis_per_bc))
    max_n = int(max(distinct_ns))
    p_cumulative = np.cumsum(ambient_profile_p)
    logp = np.log(ambient_profile_p)
    sim_n_idx = np.searchsorted(distinct_ns, umis_per_bc)
    num_barcodes = len(umis_per_bc)
    num_lower_loglk = np.zeros(num_barcodes, dtype=np.int32)
    # `log_nvals` is constant across the simulation loop; hoist it so
    # `_eval_multinomial_loglikelihood_cumulative` doesn't recompute it `num_sims`
    # times (memray flagged this as a top transient allocator).
    log_nvals = np.log(np.arange(1, max_n + 1))
    for _ in range(num_sims):
        draw = _draw_multinomial_sample(max_n, p_cumulative)
        # Note we can use the distinct N values as indices to extract the log-likelihood
        # at each N value.
        sim_loglk = _eval_multinomial_loglikelihood_cumulative(draw, logp, log_nvals)[
            distinct_ns - 1
        ]
        num_lower_loglk += (sim_loglk[sim_n_idx] < obs_loglk).astype(int)

    pvalues = (1 + num_lower_loglk) / (1 + num_sims)
    assert len(pvalues) == num_barcodes
    return pvalues
