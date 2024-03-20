#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Utilities for dealing with targeted and features."""

# Do not add new things to this module.
# Instead, either find or create a module with a name that better describes
# the functionality implemented by the methods or classes you want to add.

from collections import namedtuple

import numpy as np
import scipy

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


def reformat_targeted_label(targeting_group):
    """Reformat internal labels to nicer labels for websummary."""
    if targeting_group.lower() in ["ontarget", "on-target"]:
        return TARGETED_WS_LABEL.label
    elif targeting_group.lower() in ["offtarget", "off-target"]:
        return OFFTARGET_WS_LABEL.label
    return targeting_group.title()
