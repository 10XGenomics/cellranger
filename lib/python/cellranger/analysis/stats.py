#!/usr/bin/env python
import numpy as np
import scipy
from sklearn.utils import sparsefuncs

def normalize_by_umi(matrix):
    counts_per_bc = matrix.get_counts_per_bc()
    median_counts_per_bc = max(1.0, np.median(counts_per_bc))
    scaling_factors = median_counts_per_bc / counts_per_bc

    # Normalize each barcode's total count by median total count
    m = matrix.m.copy().astype(np.float64)
    sparsefuncs.inplace_column_scale(m, scaling_factors)

    return m

def normalize_by_idf(matrix):
    numbcs_per_feature = matrix.get_numbcs_per_feature()
    scaling_factors_row = np.log(matrix.bcs_dim + 1) - np.log(1 + numbcs_per_feature) 

    m = matrix.m.copy().astype(np.float64)
    sparsefuncs.inplace_row_scale(m, scaling_factors_row)

    return m

def summarize_columns(matrix):
    ''' Calculate mean and variance of each column, in a sparsity-preserving way.'''
    mu = matrix.mean(axis=0).A

    # sparse variance = E(col^2) - E(col)^2
    mu2 = matrix.multiply(matrix).mean(axis=0).A
    var = mu2 - mu**2

    return mu, var

def get_normalized_dispersion(mat_mean, mat_var, nbins=20):
    """ Calculates the normalized dispersion.  The dispersion is calculated for each feature
        and then normalized to see how its dispersion compares to samples that had a
        similar mean value.
    """
    # See equation in https://academic.oup.com/nar/article/40/10/4288/2411520
    # If a negative binomial is parameterized with mean m, and variance = m + d * m^2
    # then this d = dispersion as calculated below
    mat_disp = (mat_var - mat_mean) / np.square(mat_mean)

    quantiles = np.percentile(mat_mean, np.arange(0, 100, 100 / nbins))
    quantiles = np.append(quantiles, mat_mean.max())

    # merge bins with no difference in value
    quantiles = np.unique(quantiles)

    if len(quantiles) <= 1:
        # pathological case: the means are all identical. just return raw dispersion.
        return mat_disp

    # calc median dispersion per bin
    (disp_meds, _, disp_bins) = scipy.stats.binned_statistic(mat_mean, mat_disp, statistic='median', bins=quantiles)

    # calc median absolute deviation of dispersion per bin
    disp_meds_arr = disp_meds[disp_bins-1] # 0th bin is empty since our quantiles start from 0
    disp_abs_dev = abs(mat_disp - disp_meds_arr)
    (disp_mads, _, disp_bins) = scipy.stats.binned_statistic(mat_mean, disp_abs_dev, statistic='median', bins=quantiles)

    # calculate normalized dispersion
    disp_mads_arr = disp_mads[disp_bins-1]
    disp_norm = (mat_disp - disp_meds_arr) / disp_mads_arr
    return disp_norm
