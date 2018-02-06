#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import collections
import numpy as np
import os
import pandas as pd
from scipy.misc import logsumexp
from scipy.special import gammaln
import scipy.stats
from sklearn.utils import sparsefuncs
import sys
import tenkit.stats as tk_stats
import cellranger.analysis.io as cr_io
import cellranger.analysis.clustering as cr_clustering
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

SSEQ_ZETA_QUANTILE = 0.995

DIFFERENTIAL_EXPRESSION = collections.namedtuple('DIFFERENTIAL_EXPRESSION', ['data'])

def estimate_size_factors(x):
    """ Estimate size factors (related to cell RNA content and GEM-to-GEM technical variance)
    Args:
      x - Sparse matrix (csc) of counts (gene x cell)
    Returns:
      Array of floats, one per cell.
    """
    counts_per_cell = np.squeeze(np.asarray(x.sum(axis=0)))
    size_factors = counts_per_cell.astype(np.float64) / np.median(counts_per_cell)
    return size_factors

def compute_sseq_params(x, zeta_quantile=SSEQ_ZETA_QUANTILE):
    """ Compute global parameters for the sSeq differential expression method.
    The key parameters are the shrunken gene-wise dispersions.

    This method was published in:
    Yu D, et al. (2013) Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq experiments with small sample size.
    Bioinformatics. 29: 1275-1282. doi: 10.1093/bioinformatics/btt143
    Args:
      x - Sparse matrix (csc) of counts (gene x cell)
      zeta_quantile (float) - Quantile of method-of-moments dispersion estimates to
                              use as the shrinkage target zeta.
    Returns:
      A dictionary containing the sSeq parameters and some diagnostic info.
    """
    # Number of cells
    N = x.shape[1]

    # Number of genes
    G = x.shape[0]

    # Estimate size factors and normalize the matrix for quick mean/var calcs
    size_factors = estimate_size_factors(x)
    # Cast to float to prevent truncation of 1 -> 0 for size factors < 1
    x_norm = x.copy().astype(np.float64)
    sparsefuncs.inplace_column_scale(x_norm, 1.0/size_factors)

    ## Estimate genewise mean, variance, and dispersion by the method of moments
    ## assuming that each gene follows a negative-binomial distribution.
    mean_g = np.squeeze(np.asarray(x_norm.mean(axis=1)))
    # V[X] = E[X^2] - E[X]^2
    mean_sq_g = np.squeeze(np.asarray(x_norm.multiply(x_norm).mean(axis=1)))
    var_g = mean_sq_g - np.square(mean_g)

    # Method of moments estimate of gene-wise dispersion (phi)
    # Only use genes with non-zero variance in the following estimation
    use_g = var_g > 0
    phi_mm_g = np.zeros(G)
    phi_mm_g[use_g] = np.maximum(0, (float(N)*var_g[use_g] - mean_g[use_g]*np.sum(1.0/size_factors)) / \
                                 (np.square(mean_g[use_g])*np.sum(1.0/size_factors)))

    ## Estimate the optimal global target dispersion (zeta_hat).
    ## The true optimal zeta is that which minimizes the MSE vs the true dispersions.
    ## The genewise dispersions will be "shrunk" towards our estimate of zeta.

    # Use a high quantile of the MoM dispersion as our shrinkage target
    # per the rule of thumb in Yu, et al.
    zeta_hat = np.nanpercentile(phi_mm_g[use_g], 100.0 * zeta_quantile)

    ## Compute delta, the optimal shrinkage towards zeta_hat
    ## This defines a linear function that shrinks the MoM dispersion estimates
    mean_phi_mm_g = np.mean(phi_mm_g[use_g])
    delta = (np.sum(np.square(phi_mm_g[use_g] - mean_phi_mm_g))/float(G-1)) / \
            (np.sum(np.square(phi_mm_g[use_g] - zeta_hat))/float(G-2))


    ## Compute the shrunken dispersion estimates
    ## Interpolate between the MoM estimates and zeta_hat by delta
    phi_g = np.full(G, np.nan)
    if np.any(phi_mm_g[use_g] > 0):
        phi_g[use_g] = (1-delta)*phi_mm_g[use_g] + delta*zeta_hat
    else:
        phi_g[use_g] = 0.0

    return {
        'N': N,
        'G': G,
        'size_factors': size_factors,
        'mean_g': mean_g,
        'var_g': var_g,
        'use_g': use_g,
        'phi_mm_g': phi_mm_g,
        'eval_zeta': None,
        'eval_asd': None,
        'asd_slope': None,
        'zeta_hat': zeta_hat,
        'delta': delta,
        'phi_g': phi_g,
    }

def neg_bin_log_pmf(k, mu, phi):
    """ Log(PMF) of negative binomial distribution with mean mu and dispersion phi,
    conveniently parameterized.
    Args:
      k (int) - NB random variable
      mu (float) - mean
      phi (float) - dispersion
    Returns:
      The log of the pmf at k. """
    r = 1.0/phi
    return gammaln(r+k) - (gammaln(r) + gammaln(k+1)) + \
        k * np.log(mu/(r+mu)) + \
        r * np.log(r/(r+mu))

def nb_exact_test(x_a, x_b, size_factor_a, size_factor_b, mu, phi):
    """ Compute p-value for a pairwise exact test using the negative binomial.
    Args:
      x_a (int) - Total count for a single gene in group A
      x_b (int) - Total count for a single gene in group B
      size_factor_a (float) - Sum of size factors for group A
      size_factor_b (float) - Sum of size factors for group B
      mu (float) - Common mean count for this gene
      phi (float) - Common dispersion for this gene
    Returns:
      p-value (float); the probability that a random pair of counts under the null hypothesis is more extreme
      than the observed pair of counts. """
    size_factor_a = float(size_factor_a)
    size_factor_b = float(size_factor_b)
    mu = float(mu)
    phi = float(phi)

    if (x_a + x_b) == 0:
        return 1.0
    if phi == 0:
        return 1.0
    if size_factor_a == 0 or size_factor_b == 0:
        return 1.0

    all_x_a = np.arange(0, 1+x_a+x_b, 1)
    all_x_b = np.arange(x_a+x_b, -1, -1)

    def log_prob(x, size_factor):
        return neg_bin_log_pmf(x, size_factor*mu, phi/size_factor)

    log_p_obs = log_prob(x_a, size_factor_a) + log_prob(x_b, size_factor_b)
    log_p_all = log_prob(all_x_a, size_factor_a) + log_prob(all_x_b, size_factor_b)

    more_extreme = log_p_all <= log_p_obs
    if np.sum(more_extreme) == 0:
        return 0.0

    return np.exp(logsumexp(log_p_all[log_p_all <= log_p_obs]) - logsumexp(log_p_all))

def nb_asymptotic_test(x_a, x_b, size_factor_a, size_factor_b, mu, phi):
    """ Compute p-value for a pairwise exact test using a fast beta approximation
    to the conditional joint distribution of (x_a, x_b).
    Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. Biostatistics, 9, 321-332
    "It is based a method-of-moments gamma approximation to the negative binomial distribution."
      - Personal communication w/ author
    Adapted from implementation in the "edgeR" package:
      https://github.com/Bioconductor-mirror/edgeR/blob/1ab290c9585335cf99bb41f50cfce2ce4d40f907/R/exactTestBetaApprox.R
    This function is vectorized. It always returns a vector even if the inputs are scalar.
    Args:
      x_a (int/np.array) - Total count for a single gene in group A
      x_b (int/np.array) - Total count for a single gene in group B
      size_factor_a (float) - Sum of size factors for group A
      size_factor_b (float) - Sum of size factors for group B
      mu (float/np.array) - Common mean count for this gene
      phi (float/np.array) - Common dispersion for this gene
    Returns:
      p-value (np.array); the probability that a random pair of counts under the null hypothesis is more extreme
      than the observed pair of counts. """
    x_a = np.array(x_a, ndmin=1, copy=False)
    x_b = np.array(x_b, ndmin=1, copy=False)
    mu = np.array(mu, ndmin=1, copy=False)
    phi = np.array(phi, ndmin=1, copy=False)

    alpha = size_factor_a * mu / (1 + phi*mu)
    beta = (size_factor_b/size_factor_a) * alpha

    total = x_a + x_b

    median = scipy.stats.beta.median(alpha, beta)
    left = ((x_a + 0.5) / total) < median
    right = np.logical_not(left)

    p = np.empty(len(x_a))
    p[left]  = 2 *    scipy.stats.beta.cdf((x_a[left]  + 0.5)/total[left],
                                           alpha[left], beta[left])
    p[right] = 2 * (1-scipy.stats.beta.cdf((x_a[right] - 0.5)/total[right],
                                           alpha[right], beta[right]))
    return p

def adjust_pvalue_bh(p):
    """ Multiple testing correction of p-values using the Benjamini-Hochberg procedure """
    descending = np.argsort(p)[::-1]
    # q = p * N / k where p = p-value, N = # tests, k = p-value rank
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

    # Return to original order
    return q[np.argsort(descending)]

def sseq_differential_expression(x, cond_a, cond_b, sseq_params, big_count=900):
    """ Run sSeq pairwise differential expression test.
      Args:
        x - Sparse matrix (csc) of counts (gene x cell)
        cond_a (np.array(int)): Indices of cells in group A
        cond_b (np.array(int)): Indices of cells in group B
        sseq_params (dict): Precomputed global parameters
        big_count (int): Use asymptotic approximation if both counts > this
      Returns:
        A pd.DataFrame with DE results for group A relative to group B """
    x_a = x[:,cond_a]
    x_b = x[:,cond_b]

    # Number of genes
    G = x.shape[0]

    # Size factors
    size_factor_a = np.sum(sseq_params['size_factors'][cond_a])
    size_factor_b = np.sum(sseq_params['size_factors'][cond_b])

    # Compute p-value for each gene
    p_values = np.ones(G)

    gene_sums_a = np.squeeze(np.asarray(x_a.sum(axis=1)))
    gene_sums_b = np.squeeze(np.asarray(x_b.sum(axis=1)))

    big = tk_stats.numpy_logical_and_list([sseq_params['use_g'], gene_sums_a > big_count, gene_sums_b > big_count])
    small = np.logical_and(sseq_params['use_g'], np.logical_not(big))

    sys.stderr.write("Computing %d exact tests and %d asymptotic tests.\n" %
                     (np.sum(small), np.sum(big)))

    # Compute exact test for small-count genes
    for i in np.flatnonzero(small):
        p_values[i] = nb_exact_test(gene_sums_a[i], gene_sums_b[i],
                                                         size_factor_a, size_factor_b,
                                                         sseq_params['mean_g'][i],
                                                         sseq_params['phi_g'][i])
    # Compute asymptotic approximation for big-count genes
    p_values[big] = nb_asymptotic_test(gene_sums_a[big],
                                                            gene_sums_b[big],
                                                            size_factor_a, size_factor_b,
                                                            sseq_params['mean_g'][big],
                                                            sseq_params['phi_g'][big])
    # Adjust p-values for multiple testing correction
    # Only adjust the genes that were actually tested
    adj_p_values = p_values.copy()
    adj_p_values[sseq_params['use_g']] = adjust_pvalue_bh(p_values[sseq_params['use_g']])

    de_result = pd.DataFrame({
        'tested': sseq_params['use_g'],
        'sum_a': gene_sums_a,
        'sum_b': gene_sums_b,
        'common_mean': sseq_params['mean_g'],
        'common_dispersion': sseq_params['phi_g'],
        'norm_mean_a': gene_sums_a / size_factor_a,
        'norm_mean_b': gene_sums_b / size_factor_b,
        'p_value': p_values,
        'adjusted_p_value': adj_p_values,
        # Introduce a pseudocount into log2(fold_change)
        'log2_fold_change': np.log2((1+gene_sums_a)/(1+size_factor_a)) - \
        np.log2((1+gene_sums_b)/(1+size_factor_b))
    })

    return de_result

def run_differential_expression(matrix, clusters, sseq_params=None):
    """ Compute differential expression for each cluster vs all other cells
        Args: matrix      - GeneBCMatrix  :  gene expression data
              clusters    - np.array(int) :  1-based cluster labels
              sseq_params - dict          :  params from compute_sseq_params """

    n_clusters = np.max(clusters)

    if sseq_params is None:
        print "Computing params..."
        sys.stdout.flush()
        sseq_params = compute_sseq_params(matrix.m)

    # Create a numpy array with 3*K columns;
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((matrix.genes_dim, 3*n_clusters))

    for cluster in xrange(1, 1+n_clusters):
        in_cluster = clusters == cluster
        group_a = np.flatnonzero(in_cluster)
        group_b = np.flatnonzero(np.logical_not(in_cluster))
        print 'Computing DE for cluster %d...' % cluster
        sys.stdout.flush()

        de_result = sseq_differential_expression(
            matrix.m, group_a, group_b, sseq_params)
        all_de_results[:, 0+3*(cluster-1)] = de_result['norm_mean_a']
        all_de_results[:, 1+3*(cluster-1)] = de_result['log2_fold_change']
        all_de_results[:, 2+3*(cluster-1)] = de_result['adjusted_p_value']

    return DIFFERENTIAL_EXPRESSION(all_de_results)

def save_differential_expression_csv(clustering_key, de, matrix, base_dir):
    out_dir = os.path.join(base_dir, clustering_key)
    cr_utils.makedirs(out_dir, allow_existing=True)

    diff_expression_fn = os.path.join(out_dir, 'differential_expression.csv')
    diff_expression_header = ['Gene ID', 'Gene Name']

    n_clusters = de.data.shape[1] / 3
    for i in xrange(n_clusters):
        diff_expression_header += [
            'Cluster %d Mean UMI Counts' % (i+1),
            'Cluster %d Log2 fold change' % (i+1),
            'Cluster %d Adjusted p value' % (i+1),
            ]

    diff_expression_prefixes = [(gene.id, gene.name) for gene in matrix.genes]
    cr_io.save_matrix_csv(diff_expression_fn,
                          de.data,
                          diff_expression_header,
                          diff_expression_prefixes)

def save_differential_expression_h5(f, clustering_key, de):
    group = f.create_group(f.root, cr_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP)

    cr_io.save_h5(f, group, clustering_key, de)

    cr_clustering.create_legacy_kmeans_nodes(f,
                                             cr_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP,
                                             cr_constants.ANALYSIS_H5_KMEANS_DIFFERENTIAL_EXPRESSION_GROUP,
                                             DIFFERENTIAL_EXPRESSION,
                                             clustering_key)
