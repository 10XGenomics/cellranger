#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Code to fit two component mixture model."""
import math

import numpy as np
import scipy
from sklearn import mixture as sk_mix

from cellranger.targeted.utils import (
    BOTH_SPHERICAL,
    BOTH_TIED,
    MIN_OFFTARGET_GENES,
    OFFTARGETS_ONLY,
    SPHERICAL,
    TIED,
    ClassStats,
    RpuFitParams,
    fit_offtarget_gaussian,
    get_class_stats,
)


# pylint: disable=singleton-comparison, too-many-locals
def fit_2comp_gmm(labels, values, covariance_type=SPHERICAL):
    """Fits a 2-component gaussian model to the on- and off-target enrichments.

    Assumes first component should be the larger/higher (usually targeted) one.
    """
    MIN_SD = 1e-2
    # the number of SDs covered by the IQR for the normal distribution
    NORM_IQR_NUM_SDS = 1.35
    RPU_NUM_TRIES = 100

    n_targeted_genes = (labels == True).sum()
    n_offtarget_genes = (labels == False).sum()
    total_genes = n_targeted_genes + n_offtarget_genes

    # arbitrarily set targeted genes to be component 1
    init_weights = [
        n_targeted_genes / float(total_genes),
        n_offtarget_genes / float(total_genes),
    ]
    init_means = [np.median(values[labels == True]), np.median(values[labels == False])]
    init_sds = [
        max(scipy.stats.iqr(values[labels == True]) / NORM_IQR_NUM_SDS, MIN_SD),
        max(scipy.stats.iqr(values[labels == False]) / NORM_IQR_NUM_SDS, MIN_SD),
    ]

    reshaped_vals = np.reshape(values, (-1, 1))

    max_lk = float("-inf")
    best_gmm = None
    best_means = [np.nan] * 2
    best_sds = [np.nan] * 2
    best_weights = [np.nan] * 2

    paired_means_to_try = [init_means]
    vals_to_try = np.linspace(0, np.nanmax(values), RPU_NUM_TRIES + 1)
    for init_mu_tgt in vals_to_try:
        for init_mu_offtgt in vals_to_try:
            if init_mu_tgt < init_mu_offtgt:
                continue
            paired_means_to_try.append([init_mu_tgt, init_mu_offtgt])

    for paired_means in paired_means_to_try:
        gmm = sk_mix.GaussianMixture(
            n_components=2,
            covariance_type=covariance_type,
            weights_init=init_weights,
            means_init=np.reshape(paired_means, (-1, 1)),
            precisions_init=(
                [1.0 / math.pow(init_sds[0], 2), 1.0 / math.pow(init_sds[1], 2)]
                if covariance_type == SPHERICAL
                else np.reshape([1.0 / math.pow(init_sds[1], 2)], (-1, 1))
            ),
            random_state=0,
        )
        gmm.fit(reshaped_vals)
        lk = np.sum(gmm.score_samples(reshaped_vals))
        if lk > max_lk:
            best_gmm = gmm
            max_lk = lk
            best_means = np.reshape(gmm.means_, (2)).tolist()
            best_sds = (
                np.sqrt(np.reciprocal(gmm.precisions_)).tolist()
                if covariance_type == SPHERICAL
                else [np.squeeze(gmm.precisions_).tolist()] * 2
            )
            best_weights = gmm.weights_.tolist()

    rpu_posterior = best_gmm.predict_proba(reshaped_vals)
    high_rpu_component = np.argmax(best_gmm.means_)
    low_rpu_component = np.argmin(best_gmm.means_)
    in_high_rpu_component = rpu_posterior[:, high_rpu_component] > 0.5
    in_low_rpu_component = rpu_posterior[:, low_rpu_component] >= 0.5

    if np.sum(in_high_rpu_component) > 0 and np.sum(in_low_rpu_component) > 0:
        log_rpu_threshold = np.min(values[in_high_rpu_component])
    else:
        log_rpu_threshold = np.nan
    params = RpuFitParams(
        log_rpu_threshold=log_rpu_threshold,
        mu_high=best_means[high_rpu_component],
        mu_low=best_means[low_rpu_component],
        sd_high=best_sds[high_rpu_component],
        sd_low=best_sds[low_rpu_component],
        alpha_high=best_weights[high_rpu_component],
        alpha_low=best_weights[low_rpu_component],
    )
    return params


# pylint: disable=singleton-comparison
def fit_enrichments(tgt_label, values, method, min_genes=MIN_OFFTARGET_GENES, FPR=0.001):
    """Wrapper for calculation of enrichments bases on RPU distributions.

    Bail if there aren't enough genes.
    Only includes genes considered quantifiable by the time we get to this function.
    Can function both in "offtargets_only" mode (estimate off-target distr, only
    compute enrichments for targeted genes) or in the "both" (2-comp GMM) mode.
    FPR only applies to the "offtargets_only" method.
    """
    n_targeted_genes = np.sum(tgt_label)
    n_offtarget_genes = tgt_label.shape[0] - n_targeted_genes

    # default values to return in case undeterminable for misc reasons
    params = RpuFitParams(
        log_rpu_threshold=np.nan,
        mu_high=np.nan,
        mu_low=np.nan,
        sd_high=np.nan,
        sd_low=np.nan,
        alpha_high=np.nan,
        alpha_low=np.nan,
    )
    class_stats = ClassStats(
        n_targeted_enriched=np.nan,
        n_targeted_not_enriched=np.nan,
        n_offtgt_enriched=np.nan,
        n_offtgt_not_enriched=np.nan,
    )

    if n_targeted_genes == 0 and n_offtarget_genes == 0:
        # nothing to do
        pass
    elif n_offtarget_genes == 0:
        # unlikely fortuitous case where no offtgts so completely enriched, maybe n_offtgt should be slightly > 0
        params = RpuFitParams(
            log_rpu_threshold=0,
            mu_high=np.nan,
            mu_low=np.nan,
            sd_high=np.nan,
            sd_low=np.nan,
            alpha_high=np.nan,
            alpha_low=np.nan,
        )
        log_rpu_threshold = params.log_rpu_threshold
        enrichments = values >= log_rpu_threshold
        class_stats = get_class_stats(tgt_label, enrichments, method=method)
    elif method == OFFTARGETS_ONLY and n_offtarget_genes < min_genes:
        # both here and below -- too few points to do anything useful
        pass
    elif method in [BOTH_SPHERICAL, BOTH_TIED] and (
        max(n_offtarget_genes, n_targeted_genes) < min_genes
        or min(n_targeted_genes, n_offtarget_genes) == 0
    ):
        pass
    else:
        # actually compute things
        if method == OFFTARGETS_ONLY:
            params = fit_offtarget_gaussian(values[~tgt_label], FPR=FPR)
        elif method == BOTH_SPHERICAL:
            params = fit_2comp_gmm(tgt_label, values, covariance_type=SPHERICAL)
        elif method == BOTH_TIED:
            params = fit_2comp_gmm(tgt_label, values, covariance_type=TIED)
        else:
            raise RuntimeError(f"Unknown enrichment calc method {method}")
        log_rpu_threshold = params.log_rpu_threshold
        enrichments = values >= log_rpu_threshold
        class_stats = get_class_stats(tgt_label, enrichments, method=method)
    return params, class_stats
