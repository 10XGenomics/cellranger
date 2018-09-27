#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

""" Functions for calling cell-associated barcodes """

from collections import namedtuple
import numpy as np
import numpy.ma as ma
from cellranger.analysis.diffexp import adjust_pvalue_bh
import cellranger.stats as cr_stats
import cellranger.sgt as cr_sgt

# Number of additional barcodes to consider after the initial cell calling
N_CANDIDATE_BARCODES=20000

# Number of partitions (max number of barcodes to consider for ambient estimation)
N_PARTITIONS=90000

# Drop this top fraction of the barcodes when estimating ambient.
MAX_OCCUPIED_PARTITIONS_FRAC = 0.5

# Minimum number of UMIS per barcode to consider after the initial cell calling
MIN_UMIS = 100

# Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling
MIN_UMI_FRAC_OF_MEDIAN = 0.01

# Maximum adjusted p-value to call a barcode as non-ambient
MAX_ADJ_PVALUE = 0.01

def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
    """ Estimate a gene expression profile by Simple Good Turing.
    Args:
      raw_mat (sparse matrix): Sparse matrix of all counts
      barcode_indices (np.array(int)): Barcode indices to use
      nz_feat (np.array(int)): Indices of features that are non-zero at least once
    Returns:
      profile (np.array(float)): Estimated probabilities of length len(nz_feat).
    """
    # Initial profile estimate
    prof_mat = matrix[:,barcode_indices]

    profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
    zero_feat = np.flatnonzero(profile == 0)

    # Simple Good Turing estimate
    p_smoothed, p0 = cr_sgt.sgt_proportions(profile[np.flatnonzero(profile)])

    # Distribute p0 equally among the zero elements.
    p0_i = p0/len(zero_feat)

    profile_p = np.repeat(p0_i, len(nz_feat))
    profile_p[np.flatnonzero(profile)] = p_smoothed

    assert np.isclose(profile_p.sum(), 1.0)
    return profile_p


# Construct a background expression profile from barcodes with <= T UMIs
def est_background_profile_sgt(matrix, use_bcs):
    """ Estimate a gene expression profile on a given subset of barcodes.
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

NonAmbientBarcodeResult = namedtuple('NonAmbientBarcodeResult',
                                     ['eval_bcs',      # Candidate barcode indices (n)
                                      'log_likelihood',# Ambient log likelihoods (n)
                                      'pvalues',       # pvalues (n)
                                      'pvalues_adj',   # B-H adjusted pvalues (n)
                                      'is_nonambient', # Boolean nonambient calls (n)
                                      ])

def find_nonambient_barcodes(matrix, orig_cell_bcs,
                             min_umi_frac_of_median=MIN_UMI_FRAC_OF_MEDIAN,
                             max_adj_pvalue=MAX_ADJ_PVALUE,):
    """ Call barcodes as being sufficiently distinct from the ambient profile

    Args:
      matrix (CountMatrix): Full expression matrix.
      orig_cell_bcs (iterable of str): Strings of initially-called cell barcodes.
    Returns:
    TBD
    """
    # Estimate an ambient RNA profile
    umis_per_bc = matrix.get_counts_per_bc()
    bc_order = np.argsort(umis_per_bc)

    # Take what we expect to be the barcodes associated w/ empty partitions.
    empty_bcs = bc_order[::-1][(N_PARTITIONS/2):N_PARTITIONS]
    empty_bcs.sort()

    # Require non-zero barcodes
    nz_bcs = np.flatnonzero(umis_per_bc)
    nz_bcs.sort()

    use_bcs = np.intersect1d(empty_bcs, nz_bcs, assume_unique=True)

    if len(use_bcs) > 0:
        try:
            eval_features, ambient_profile_p = est_background_profile_sgt(matrix.m, use_bcs)
        except cr_sgt.SimpleGoodTuringError as e:
            print str(e)
            return None
    else:
        eval_features = np.zeros(0, dtype=int)
        ambient_profile_p = np.zeros(0)

    # Choose candidate cell barcodes
    orig_cell_bc_set = set(orig_cell_bcs)
    orig_cells = np.flatnonzero(np.fromiter((bc in orig_cell_bc_set for bc in matrix.bcs),
                                            count=len(matrix.bcs), dtype=bool))

    # No good incoming cell calls
    if orig_cells.sum() == 0:
        return None

    # Look at non-cell barcodes above a minimum UMI count
    eval_bcs = np.ma.array(np.arange(matrix.bcs_dim))
    eval_bcs[orig_cells] = ma.masked

    median_initial_umis = np.median(umis_per_bc[orig_cells])
    min_umis = int(max(MIN_UMIS, round(np.ceil(median_initial_umis * min_umi_frac_of_median))))
    print('Median UMIs of initial cell calls: {}'.format(median_initial_umis))
    print('Min UMIs: {}'.format(min_umis))

    eval_bcs[umis_per_bc < min_umis] = ma.masked
    n_unmasked_bcs = len(eval_bcs) - eval_bcs.mask.sum()

    # Take the top N_CANDIDATE_BARCODES by UMI count, of barcodes that pass the above criteria
    eval_bcs = np.argsort(ma.masked_array(umis_per_bc, mask=eval_bcs.mask))[0:n_unmasked_bcs][-N_CANDIDATE_BARCODES:]

    if len(eval_bcs) == 0:
        return None

    assert not np.any(np.isin(eval_bcs, orig_cells))
    print('Number of candidate bcs: {}'.format(len(eval_bcs)))
    print('Range candidate bc umis: {}, {}'.format(umis_per_bc[eval_bcs].min(), umis_per_bc[eval_bcs].max()))

    eval_mat = matrix.m[eval_features, :][:, eval_bcs]

    if len(ambient_profile_p) == 0:
        obs_loglk = np.repeat(np.nan, len(eval_bcs))
        pvalues = np.repeat(1, len(eval_bcs))
        sim_loglk = np.repeat(np.nan, len(eval_bcs))
        return None

    # Compute observed log-likelihood of barcodes being generated from ambient RNA
    obs_loglk = cr_stats.eval_multinomial_loglikelihoods(eval_mat, ambient_profile_p)

    # Simulate log likelihoods
    distinct_ns, sim_loglk = cr_stats.simulate_multinomial_loglikelihoods(ambient_profile_p, umis_per_bc[eval_bcs], num_sims=10000, verbose=True)

    # Compute p-values
    pvalues = cr_stats.compute_ambient_pvalues(umis_per_bc[eval_bcs], obs_loglk, distinct_ns, sim_loglk)

    pvalues_adj = adjust_pvalue_bh(pvalues)

    is_nonambient = pvalues_adj <= max_adj_pvalue

    return NonAmbientBarcodeResult(
        eval_bcs=eval_bcs,
        log_likelihood=obs_loglk,
        pvalues=pvalues,
        pvalues_adj=pvalues_adj,
        is_nonambient=is_nonambient,
    )
