#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import array
import numpy as np
import scipy.stats
import cellranger.constants as cr_constants
import tenkit.constants as tk_constants
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats

# Inverse Simpson Index, or the effective diversity of power 2
def effective_diversity(counts):
    numerator = np.sum(counts)**2
    denominator = np.sum(v**2 for v in counts)
    effective_diversity = tk_stats.robust_divide(float(numerator), float(denominator))
    return effective_diversity

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
                p_edit = 10.0**(-edit_qv/10.0)
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

def compute_percentile_from_distribution(counter, percentile):
    """ Takes a Counter object (or value:frequency dict) and computes a single percentile.
    Uses Type 7 interpolation from:
      Hyndman, R.J.; Fan, Y. (1996). "Sample Quantiles in Statistical Packages".
    """
    assert 0 <= percentile <= 100

    n = np.sum(counter.values())
    h = (n-1)*(percentile/100.0)
    lower_value = None

    cum_sum = 0
    for value, freq in sorted(counter.items()):
        cum_sum += freq
        if cum_sum > np.floor(h) and lower_value is None:
            lower_value = value
        if cum_sum > np.ceil(h):
            return lower_value + (h-np.floor(h)) * (value-lower_value)

# Test for compute_percentile_from_distribution()
#def test_percentile(x, p):
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

# barcode filtering methods

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

def find_within_ordmag(x, baseline_idx):
    x_ascending = np.sort(x)
    baseline = x_ascending[-baseline_idx]
    cutoff = max(1, round(0.1*baseline))
    # Return the index corresponding to the cutoff in descending order
    return len(x) - np.searchsorted(x_ascending, cutoff)

def summarize_bootstrapped_top_n(top_n_boot):
    top_n_bcs_mean = np.mean(top_n_boot)
    top_n_bcs_sd = np.std(top_n_boot)
    top_n_bcs_var = np.var(top_n_boot)
    result = {}
    result['filtered_bcs_var'] = top_n_bcs_var
    result['filtered_bcs_cv'] = tk_stats.robust_divide(top_n_bcs_sd, top_n_bcs_mean)
    result['filtered_bcs_lb'] = round(scipy.stats.norm.ppf(0.025, top_n_bcs_mean, top_n_bcs_sd))
    result['filtered_bcs_ub'] = round(scipy.stats.norm.ppf(0.975, top_n_bcs_mean, top_n_bcs_sd))
    result['filtered_bcs'] = int(round(top_n_bcs_mean))
    return result

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

    baseline_bc_idx = int(round(float(recovered_cells) * (1-cr_constants.ORDMAG_RECOVERED_CELLS_QUANTILE)))
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
