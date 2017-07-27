#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utils for math operations
#

import math
import numpy as np
from scipy.stats import norm

def generalized_iqr(x, left_cutoff=0.25, right_cutoff=0.25):
    """Returns the width of the distribution given by x, cutting
    off the left and right tails at the specified cutoff values.
    """
    if len(x) == 0:
        return float('NaN')

    sorted_vals = sorted(x)
    left_index = int(len(x)*left_cutoff)
    right_index = int(len(x)*(1 - right_cutoff))
    iqr = sorted_vals[right_index] - sorted_vals[left_index]
    return iqr

def entropy(ps):
    """Calculates the entropy (log 2) of the distribution given by p
    """
    entropy = 0.0
    for p in ps:
        if not(p == 0):
            entropy -= p*math.log(p, 2)
    return entropy

def kl_divergence(ps, qs):
    """ Calculates the relative entropy of distribution of ps to qs
    """
    rel_entropy = 0.0
    for (p,q) in zip(ps, qs):
        if p > 0:
            rel_entropy += p*math.log(p/q, 2)
    return rel_entropy

def N50(lengths):
    """ Calculates the N50 for a given set of lengths
    """
    return NX(lengths, 0.5)

def NX(lengths, fraction):
    """ Calculates the N50 for a given set of lengths
    """
    lengths = sorted(lengths, reverse=True)
    tot_length = sum(lengths)
    cum_length = 0
    for length in lengths:
        cum_length += length
        if cum_length >= fraction*tot_length:
            return length

def mean_var_from_counts(val_counts):
    total_count = 0.0
    total_sum = 0.0
    total_sum_squares = 0.0
    for val, count in val_counts.iteritems():
        total_count += count
        total_sum += count*val
        total_sum_squares += count*val*val

    if total_count == 0.0:
        return float('NaN'), float('NaN')

    mean_val = total_sum/total_count
    var_val = total_sum_squares/total_count - mean_val*mean_val

    return (mean_val, var_val)

def robust_divide(a,b):
    """Handles 0 division and conversion to floats automatically
    """
    a = float(a)
    b = float(b)
    if b == 0:
        return float('NaN')
    else:
        return a/b

def robust_percentile(arr, p):
    if len(arr) == 0:
        return float('NaN')
    return np.percentile(np.array(arr), p)

def log_1minus(x):
    """Computes log(1 - x). More accurate than doing np.log(1-x)."""
    return np.log1p(-x)

def log_prob_correct_from_qual(q):
    """Computes the probability of no error given a phred quality."""
    return np.log1p(- 10**(-0.1 * q))

def log_prob_wrong_from_qual(q):
    """Computes the probability of error given a phred quality."""
    return (-0.1 * q) / np.log10(np.exp(1))

def qual_from_prob_correct(p):
    """Computes phred quality given probability of no error."""
    qual = np.round(-10.0 * np.log10(1-p))
    if np.isnan(qual):
        return 0
    elif qual > 60:
        return 60
    else:
        return int(qual)

def logaddexp(arr):
    """Computes log(exp(arr[0]) + exp(arr[1]) + ...). """
    assert(len(arr) >= 2)
    res = np.logaddexp(arr[0], arr[1])
    for i in arr[2:]:
        res = np.logaddexp(res, i)
    return res

def norm_std_from_iqr(mu, iqr):
    """Computes the standard deviation of a normal distribution with mean mu and IQR irq."""
    # Assuming normal distribution (so symmetric), Q1 = m - (iqr / 2)
    Q1 = mu - (iqr / 2.0)
    # Assuming a normal distribution with mean m, std s, and first quartile Q1,
    # we have (Q1 - m)/s = ppf(0.25)
    return (Q1 - mu) / norm.ppf(0.25)

def numpy_logical_and_list(list_of_logicals):
    assert(len(list_of_logicals) >= 2)
    output = list_of_logicals[0]
    for i in range(1,len(list_of_logicals)):
        output = np.logical_and(output, list_of_logicals[i])
    return output
