#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Utilities for dealing with targeted and features."""

# Do not add new things to this module.
# Instead, either find or create a module with a name that better describes
# the functionality implemented by the methods or classes you want to add.

import math
from collections import namedtuple

import numpy as np
import scipy
import sklearn.mixture as sk_mix

from cellranger.constants import OFF_TARGET_SUBSAMPLE, ON_TARGET_SUBSAMPLE

# pylint: disable=invalid-name

######################################################################
# Constants (labels, column names, colors) used for targeted WS
######################################################################

# Targeting column name constants in feature_metrics table
TARGETING_COLNAME = "is_targeted"
READS_COLNAME = "num_reads"
UMIS_COLNAME = "num_umis"
READS_IN_CELLS_COLNAME = "num_reads_cells"
UMIS_IN_CELLS_COLNAME = "num_umis_cells"
LOG_RPU_CELLS_COLNAME = "mean_reads_per_umi_cells_log10"
RPU_COLNAME = "mean_reads_per_umi"
ENRICHMENT_COLNAME = "enriched_rpu"
FEATURE_ID_COLNAME = "feature_id"
FEATURE_NAME_COLNAME = "feature_name"
TARGET_FEATURE_METRICS_COLS = [
    FEATURE_ID_COLNAME,
    FEATURE_NAME_COLNAME,
    TARGETING_COLNAME,
    READS_IN_CELLS_COLNAME,
    UMIS_IN_CELLS_COLNAME,
    LOG_RPU_CELLS_COLNAME,
    RPU_COLNAME,
    ENRICHMENT_COLNAME,
]

# Targeting color scheme for websummary
ON_TARGET_COLOR = "#0071D9"
OFF_TARGET_COLOR = "darkgrey"

# Labels used in targeted WS (possible metrics suffixes and corresponding label)
ON_TARGET_METRIC_SUFFIX = "_on_target"
OFF_TARGET_METRIC_SUFFIX = "_off_target"
WsLabel = namedtuple("ws_label", ["subsample_suffix", "metric_suffix", "label", "color"])
TARGETED_WS_LABEL = WsLabel(
    subsample_suffix=ON_TARGET_SUBSAMPLE,
    metric_suffix=ON_TARGET_METRIC_SUFFIX,
    label="Targeted",
    color=ON_TARGET_COLOR,
)
OFFTARGET_WS_LABEL = WsLabel(
    subsample_suffix=OFF_TARGET_SUBSAMPLE,
    metric_suffix=OFF_TARGET_METRIC_SUFFIX,
    label="Non-Targeted",
    color=OFF_TARGET_COLOR,
)
TARGETING_WS_LABELS = [OFFTARGET_WS_LABEL, TARGETED_WS_LABEL]

######################################################################
# Functions and constants related to enrichment calculations
######################################################################

EST_PREFILTER_RPU_COLNAME = "rpu_estimate_pre_umi_filtering"
POSTFILTER_RPU_COLNAME = "rpu_estimate_post_umi_filtering"
RPU_IS_CONF_EST_COLNAME = "is_rpu_estimate_confident"
# reasonable threshold for fitting a normal
MIN_OFFTARGET_GENES = 30
# min on-target rpu threshold under which to consider disabling rpu stats
MIN_RPU_THRESHOLD = 3

RpuFitParams = namedtuple(
    "RpuFitParams",
    ["log_rpu_threshold", "mu_high", "mu_low", "sd_high", "sd_low", "alpha_high", "alpha_low"],
)
ClassStats = namedtuple(
    "ClassStats",
    [
        "n_targeted_enriched",
        "n_targeted_not_enriched",
        "n_offtgt_enriched",
        "n_offtgt_not_enriched",
    ],
)

# keywords for scikit GMM wrapper fn
SPHERICAL = "spherical"
TIED = "tied"
# enrichment method calculations allowed
OFFTARGETS_ONLY = "offtargets_only"
BOTH_TIED = "both_tied"
BOTH_SPHERICAL = "both_spherical"
ENRICH_METHODS = [OFFTARGETS_ONLY, BOTH_TIED, BOTH_SPHERICAL]


# pylint: disable=singleton-comparison
def get_class_stats(labels, calls, method):
    """Get number of targeted and off-target genes classified as enriched and non-enriched."""
    n_targeted_enriched = ((labels == True) & (calls == True)).sum()
    n_targeted_not_enriched = ((labels == True) & (calls == False)).sum()
    if method in [BOTH_TIED, BOTH_SPHERICAL]:
        n_offtgt_enriched = ((labels == False) & (calls == True)).sum()
        n_offtgt_not_enriched = ((labels == False) & (calls == False)).sum()
    else:
        n_offtgt_enriched, n_offtgt_not_enriched = np.nan, np.nan
    class_stats = ClassStats(
        n_targeted_enriched=n_targeted_enriched,
        n_targeted_not_enriched=n_targeted_not_enriched,
        n_offtgt_enriched=n_offtgt_enriched,
        n_offtgt_not_enriched=n_offtgt_not_enriched,
    )
    return class_stats


# pylint: disable=singleton-comparison, too-many-locals
def fit_2comp_GMM(labels, values, covariance_type=SPHERICAL):
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


def fit_offtarget_gaussian(values, FPR=0.001):
    """Estimates the mean and sd of the provided values (log reads per umi).

    using the median and quartiles. Sets the threshold at an FPR of 0.001 for
    calling genes as significantly different from this distribution.
    """
    MIN_SD = 1e-2
    # the number of SDs covered by the 10th-90th percentile for the normal distribution
    NORM_10_90_NUM_SDS = 2.56

    mu = np.median(values)
    sd = max(MIN_SD, (np.percentile(values, 90) - np.percentile(values, 10)) / NORM_10_90_NUM_SDS)

    log_rpu_threshold = scipy.stats.norm.ppf(1.0 - FPR, loc=mu, scale=sd)

    params = RpuFitParams(
        log_rpu_threshold=log_rpu_threshold,
        mu_high=np.nan,
        mu_low=mu,
        sd_high=np.nan,
        sd_low=sd,
        alpha_high=np.nan,
        alpha_low=np.nan,
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
            params = fit_2comp_GMM(tgt_label, values, covariance_type=SPHERICAL)
        elif method == BOTH_TIED:
            params = fit_2comp_GMM(tgt_label, values, covariance_type=TIED)
        else:
            raise RuntimeError(f"Unknown enrichment calc method {method}")
        log_rpu_threshold = params.log_rpu_threshold
        enrichments = values >= log_rpu_threshold
        class_stats = get_class_stats(tgt_label, enrichments, method=method)
    return params, class_stats


def reformat_targeted_label(targeting_group):
    """Reformat internal labels to nicer labels for websummary."""
    if targeting_group.lower() in ["ontarget", "on-target"]:
        return TARGETED_WS_LABEL.label
    elif targeting_group.lower() in ["offtarget", "off-target"]:
        return OFFTARGET_WS_LABEL.label
    return targeting_group.title()
