#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import array
import numpy as np
import scipy.stats as sp_stats
import sys
import cellranger.constants as cr_constants
import tenkit.constants as tk_constants
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import collections
import sklearn.mixture as sk_mix

def to_col_vec(a):
    """ Convert a 1-d array to a column vector """
    return np.reshape(a, (len(a), 1))

def create_gmm(weights, means, sd):
    """ Create a 2-component GMM with tied variance and given initialization
        This uses the sklearn 0.17.1 interface - it changes in 0.18.x """
    gmm = sk_mix.GMM(n_components=2,
                     covariance_type='tied',
                     init_params='',
                     params='wmc')
    gmm.weights_ = np.array(weights)
    gmm.means_ = np.reshape(means, (len(means), 1))
    gmm._set_covars(np.reshape(sd, (1,1)))
    return gmm

def multistart_gmm(data, weights, means_list, sd):
    """ Sweep over the given initial mean vectors
        and return the result with the highest log-likelihood """
    best_gmm = None
    max_loglk = float('-inf')
    for means in means_list:
        gmm = create_gmm(weights=weights, means=means, sd=sd)
        gmm.fit(data)

        # sklearn 0.17 return type
        loglk = np.sum(gmm.score_samples(data)[0])

        if loglk > max_loglk:
            best_gmm = gmm
            max_loglk = loglk

    return best_gmm

# Inverse Simpson Index, or the effective diversity of power 2
def effective_diversity(counts):
    numerator = np.sum(counts)**2
    denominator = np.sum(v**2 for v in counts)
    effective_diversity = tk_stats.robust_divide(float(numerator), float(denominator))
    return effective_diversity

def compute_percentile_from_distribution(counter, percentile):
    """ Takes a Counter object (or value:frequency dict) and computes a single percentile.
    Uses Type 7 interpolation from:
      Hyndman, R.J.; Fan, Y. (1996). "Sample Quantiles in Statistical Packages".
    """
    assert 0 <= percentile <= 100

    n = np.sum(counter.values())
    h = (n - 1) * (percentile / 100.0)
    lower_value = None

    cum_sum = 0
    for value, freq in sorted(counter.items()):
        cum_sum += freq
        if cum_sum > np.floor(h) and lower_value is None:
            lower_value = value
        if cum_sum > np.ceil(h):
            return lower_value + (h - np.floor(h)) * (value - lower_value)

# Test for compute_percentile_from_distribution()
# def test_percentile(x, p):
#    c = Counter()
#    for xi in x:
#        c[xi] += 1
#    my_res = np.array([compute_percentile_from_distribution(c, p_i) for p_i in p], dtype=float)
#    numpy_res = np.percentile(x, p)
#    print np.sum(np.abs(numpy_res - my_res))

def compute_iqr_from_distribution(counter):
    p25 = compute_percentile_from_distribution(counter, 25)
    p75 = compute_percentile_from_distribution(counter, 75)
    return p75 - p25

def compute_median_from_distribution(counter):
    return compute_percentile_from_distribution(counter, 50)


def correct_bc_error(bc_confidence_threshold, seq, qual, wl_dist):
    '''Attempt to correct an incorrect BC sequence by computing
    the probability that a Hamming distance=1 BC generated
    the observed sequence, accounting for the prior distribution
    of the whitelist barcodes (wl_dist), and the QV of the base
    that must have been incorrect'''

    # QV values
    qvs = np.fromstring(qual, dtype=np.byte) - tk_constants.ILLUMINA_QUAL_OFFSET

    # Char array of read
    a = array.array('c', seq)

    # Likelihood of candidates
    wl_cand = []
    likelihoods = []

    # Enumerate Hamming distance 1 sequences - if a sequence
    # is on the whitelist, compute it's likelihood.
    for pos in range(len(a)):
        existing = a[pos]
        for c in tk_seq.NUCS:
            if c == existing:
                continue
            a[pos] = c
            test_str = a.tostring()

            # prior probability of this BC
            p_bc = wl_dist.get(test_str)
            if p_bc is not None:
                # probability of the base error
                edit_qv = min(33.0, float(qvs[pos]))
                p_edit = 10.0**(-edit_qv / 10.0)
                wl_cand.append(test_str)
                likelihoods.append(p_bc * p_edit)

        a[pos] = existing

    posterior = np.array(likelihoods)
    posterior /= posterior.sum()
    if len(posterior) > 0:
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)]

    return None

def determine_max_filtered_bcs(total_diversity, recovered_cells):
    """ Determine the max # of cellular barcodes to consider """
    return float(recovered_cells) * cr_constants.FILTER_BARCODES_MAX_RECOVERED_CELLS_MULTIPLE

def init_barcode_filter_result():
    return {
        'filtered_bcs': 0,
        'filtered_bcs_lb': 0,
        'filtered_bcs_ub': 0,
        'max_filtered_bcs': 0,
        'filtered_bcs_var': 0,
        'filtered_bcs_cv': 0,
    }

def summarize_bootstrapped_top_n(top_n_boot):
    top_n_bcs_mean = np.mean(top_n_boot)
    top_n_bcs_sd = np.std(top_n_boot)
    top_n_bcs_var = np.var(top_n_boot)
    result = {}
    result['filtered_bcs_var'] = top_n_bcs_var
    result['filtered_bcs_cv'] = tk_stats.robust_divide(top_n_bcs_sd, top_n_bcs_mean)
    result['filtered_bcs_lb'] = round(sp_stats.norm.ppf(0.025, top_n_bcs_mean, top_n_bcs_sd))
    result['filtered_bcs_ub'] = round(sp_stats.norm.ppf(0.975, top_n_bcs_mean, top_n_bcs_sd))
    result['filtered_bcs'] = int(round(top_n_bcs_mean))
    return result

def find_within_ordmag(x, baseline_idx):
    x_ascending = np.sort(x)
    baseline = x_ascending[-baseline_idx]
    cutoff = max(1, round(0.1 * baseline))
    # Return the index corresponding to the cutoff in descending order
    return len(x) - np.searchsorted(x_ascending, cutoff)

def filter_cellular_barcodes_ordmag(bc_counts, recovered_cells, total_diversity):
    """ Simply take all barcodes that are within an order of magnitude of a top barcode
        that likely represents a cell
    """
    if recovered_cells is None:
        recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    metrics = init_barcode_filter_result()
    max_filtered_bcs = determine_max_filtered_bcs(total_diversity, recovered_cells)
    metrics['max_filtered_bcs'] = max_filtered_bcs

    nonzero_bc_counts = bc_counts[bc_counts > 0]
    if len(nonzero_bc_counts) == 0:
        msg = "WARNING: All barcodes do not have enough reads for ordmag, allowing no bcs through"
        return [], metrics, msg

    baseline_bc_idx = int(round(float(recovered_cells) * (1 - cr_constants.ORDMAG_RECOVERED_CELLS_QUANTILE)))
    baseline_bc_idx = min(baseline_bc_idx, len(nonzero_bc_counts) - 1)
    assert baseline_bc_idx < max_filtered_bcs

    # Bootstrap sampling; run algo with many random samples of the data
    top_n_boot = np.array([
        find_within_ordmag(np.random.choice(nonzero_bc_counts, len(nonzero_bc_counts)), baseline_bc_idx)
        for i in xrange(cr_constants.ORDMAG_NUM_BOOTSTRAP_SAMPLES)
    ])

    metrics.update(summarize_bootstrapped_top_n(top_n_boot))

    # Get the filtered barcodes
    top_n = metrics['filtered_bcs']
    top_bc_idx = np.sort(np.argsort(bc_counts)[::-1][0:top_n])
    return top_bc_idx, metrics, None

def filter_cellular_barcodes_fixed_cutoff(bc_counts, cutoff):
    nonzero_bcs = len(bc_counts[bc_counts > 0])
    top_n = min(cutoff, nonzero_bcs)
    top_bc_idx = np.sort(np.argsort(bc_counts)[::-1][0:top_n])
    metrics = {
        'filtered_bcs': top_n,
        'filtered_bcs_lb': top_n,
        'filtered_bcs_ub': top_n,
        'max_filtered_bcs': 0,
        'filtered_bcs_var': 0,
        'filtered_bcs_cv': 0,
    }
    return top_bc_idx, metrics, None

def filter_cellular_barcodes_manual(matrix, cell_barcodes):
    """ Take take all barcodes that were given as cell barcodes """
    barcodes = list(set(matrix.bcs) & set(cell_barcodes))

    metrics = {
        'filtered_bcs': len(barcodes),
        'filtered_bcs_lb': len(barcodes),
        'filtered_bcs_ub': len(barcodes),
        'max_filtered_bcs': 0,
        'filtered_bcs_var': 0,
        'filtered_bcs_cv': 0,
    }

    return barcodes, metrics, None


def merge_filtered_metrics(filtered_metrics):
    result = {
        'filtered_bcs': 0,
        'filtered_bcs_lb': 0,
        'filtered_bcs_ub': 0,
        'max_filtered_bcs': 0,
        'filtered_bcs_var': 0,
        'filtered_bcs_cv': 0,
    }
    for i, fm in enumerate(filtered_metrics):
        # Add per-gem group metrics
        result.update({'gem_group_%d_%s' % (i + 1, key): value for key, value in fm.iteritems()})

        # Compute metrics over all gem groups
        result['filtered_bcs'] += fm['filtered_bcs']
        result['filtered_bcs_lb'] += fm['filtered_bcs_lb']
        result['filtered_bcs_ub'] += fm['filtered_bcs_ub']
        result['max_filtered_bcs'] += fm['max_filtered_bcs']
        result['filtered_bcs_var'] += fm['filtered_bcs_var']

    # Estimate CV based on sum of variances and means
    result['filtered_bcs_cv'] = tk_stats.robust_divide(
        np.sqrt(result['filtered_bcs_var']), fm['filtered_bcs'])

    return result

def correct_umis(dupe_keys):
    corrected_dupe_keys = collections.defaultdict(dict)
    for dupe_key, umis in dupe_keys.iteritems():
        for umi in umis:
            new_umi = correct_umi(umi, umis)
            if not (new_umi == umi):
                corrected_dupe_keys[dupe_key][umi] = new_umi

    return corrected_dupe_keys

def correct_umi(seq, counts):
    corrected_seq = seq
    count = counts.get(seq, 0)

    a = array.array('c', seq)
    for pos in xrange(len(a)):
        existing = a[pos]
        for c in tk_seq.NUCS:
            if c == existing:
                continue
            a[pos] = c
            test_str = a.tostring()

            value = counts.get(test_str, 0)
            if value > count or (value == count and corrected_seq < test_str):
                corrected_seq = test_str
                count = value

        a[pos] = existing
    return corrected_seq

def est_background_profile_bottom(matrix, bottom_frac):
    """Construct a background expression profile from the barcodes that make up the bottom b% of the data
    Args:
      matrix (scipy.sparse.csc_matrix): Feature x Barcode matrix
      bottom_frac (float): Use barcodes making up the bottom x fraction of the counts (0-1)
    Returns:
      (nz_feat (ndarray(int)), profile_p (ndarray(float)): Indices of nonzero features and background profile
    """
    assert bottom_frac >= 0 and bottom_frac <= 1
    umis_per_bc = np.ravel(np.asarray(matrix.sum(0)))
    barcode_order = np.argsort(umis_per_bc)

    cum_frac = np.cumsum(umis_per_bc[barcode_order]) / float(umis_per_bc.sum())
    max_bg_idx = np.searchsorted(cum_frac, bottom_frac, side='left')
    bg_mat = matrix[:, barcode_order[0:max_bg_idx]]

    nz_feat = np.flatnonzero(np.asarray(bg_mat.sum(1)))
    bg_profile = np.ravel(bg_mat[nz_feat, :].sum(axis=1))
    bg_profile_p = bg_profile / float(np.sum(bg_profile))
    assert np.isclose(bg_profile_p.sum(), 1)

    return (nz_feat, bg_profile_p)

def eval_multinomial_loglikelihoods(matrix, profile_p, max_mem_gb=0.1):
    """Compute the multinomial log PMF for many barcodes
    Args:
      matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
      profile_p (np.ndarray(float)): Multinomial probability vector
      max_mem_gb (float): Try to bound memory usage.
    Returns:
      log_likelihoods (np.ndarray(float)): Log-likelihood for each barcode
    """
    gb_per_bc = float(matrix.shape[0] * matrix.dtype.itemsize) / (1024**3)
    bcs_per_chunk = max(1, int(round(max_mem_gb/gb_per_bc)))
    num_bcs = matrix.shape[1]

    loglk = np.zeros(num_bcs)

    for chunk_start in xrange(0, num_bcs, bcs_per_chunk):
        chunk = slice(chunk_start, chunk_start+bcs_per_chunk)
        matrix_chunk = matrix[:,chunk].transpose().toarray()
        n = matrix_chunk.sum(1)
        loglk[chunk] = sp_stats.multinomial.logpmf(matrix_chunk, n, p=profile_p)
    return loglk

def simulate_multinomial_loglikelihoods(profile_p, umis_per_bc,
                                        num_sims=1000, jump=1000,
                                        n_sample_feature_block=1000000, verbose=False):
    """Simulate draws from a multinomial distribution for various values of N.

       Uses the approximation from Lun et al. ( https://www.biorxiv.org/content/biorxiv/early/2018/04/04/234872.full.pdf )

    Args:
      profile_p (np.ndarray(float)): Probability of observing each feature.
      umis_per_bc (np.ndarray(int)): UMI counts per barcode (multinomial N).
      num_sims (int): Number of simulations per distinct N value.
      jump (int): Vectorize the sampling if the gap between two distinct Ns exceeds this.
      n_sample_feature_block (int): Vectorize this many feature samplings at a time.
    Returns:
      (distinct_ns (np.ndarray(int)), log_likelihoods (np.ndarray(float)):
      distinct_ns is an array containing the distinct N values that were simulated.
      log_likelihoods is a len(distinct_ns) x num_sims matrix containing the
        simulated log likelihoods.
    """
    distinct_n = np.flatnonzero(np.bincount(umis_per_bc))

    loglk = np.zeros((len(distinct_n), num_sims), dtype=float)
    num_all_n = np.max(distinct_n) - np.min(distinct_n)
    if verbose:
        print 'Number of distinct N supplied: %d' % len(distinct_n)
        print 'Range of N: %d' % num_all_n
        print 'Number of features: %d' % len(profile_p)

    sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
    k = 0

    log_profile_p = np.log(profile_p)

    for sim_idx in xrange(num_sims):
        if verbose and sim_idx % 100 == 99:
            sys.stdout.write('.')
            sys.stdout.flush()
        curr_counts = np.ravel(sp_stats.multinomial.rvs(distinct_n[0], profile_p, size=1))

        curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[0], p=profile_p)

        loglk[0, sim_idx] = curr_loglk

        for i in xrange(1, len(distinct_n)):
            step = distinct_n[i] - distinct_n[i-1]
            if step >= jump:
                # Instead of iterating for each n, sample the intermediate ns all at once
                curr_counts += np.ravel(sp_stats.multinomial.rvs(step, profile_p, size=1))
                curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[i], p=profile_p)
                assert not np.isnan(curr_loglk)
            else:
                # Iteratively sample between the two distinct values of n
                for n in xrange(distinct_n[i-1]+1, distinct_n[i]+1):
                    j = sampled_features[k]
                    k += 1
                    if k >= n_sample_feature_block:
                        # Amortize this operation
                        sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
                        k = 0
                    curr_counts[j] += 1
                    curr_loglk += log_profile_p[j] + np.log(float(n)/curr_counts[j])

            loglk[i, sim_idx] = curr_loglk

    if verbose:
        sys.stdout.write('\n')

    return distinct_n, loglk

def compute_ambient_pvalues(umis_per_bc, obs_loglk, sim_n, sim_loglk):
    """Compute p-values for observed multinomial log-likelihoods
    Args:
      umis_per_bc (nd.array(int)): UMI counts per barcode
      obs_loglk (nd.array(float)): Observed log-likelihoods of each barcode deriving from an ambient profile
      sim_n (nd.array(int)): Multinomial N for simulated log-likelihoods
      sim_loglk (nd.array(float)): Simulated log-likelihoods of shape (len(sim_n), num_simulations)
    Returns:
      pvalues (nd.array(float)): p-values
    """
    assert len(umis_per_bc) == len(obs_loglk)
    assert sim_loglk.shape[0] == len(sim_n)

    # Find the index of the simulated N for each barcode
    sim_n_idx = np.searchsorted(sim_n, umis_per_bc)
    num_sims = sim_loglk.shape[1]

    num_barcodes = len(umis_per_bc)

    pvalues = np.zeros(num_barcodes)

    for i in xrange(num_barcodes):
        num_lower_loglk = np.sum(sim_loglk[sim_n_idx[i],:] < obs_loglk[i])
        pvalues[i] = float(1 + num_lower_loglk) / (1 + num_sims)
    return pvalues
