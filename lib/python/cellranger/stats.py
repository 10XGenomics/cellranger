#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

import numpy as np
from scipy import optimize
from scipy.special import betaln, gammaln, loggamma

import tenkit.stats as tk_stats

# Pre-seed the RNG on package import
RNG = np.random.default_rng(seed=42)


def effective_diversity(counts):
    """Inverse Simpson Index, or the effective diversity of power 2."""
    numerator = np.sum(counts) ** 2
    denominator = np.sum(v**2 for v in counts)
    return tk_stats.robust_divide(float(numerator), float(denominator))


def incremental_counts_from_sample(sample_draws):
    """Produces an incremental count vector from a categorial sample vector.

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
    idxs, c = np.unique(sample_draws, return_counts=True)
    mask = np.isin(sample_draws, idxs[c > 1])
    inc_counts = np.ones_like(sample_draws)
    counts = {}
    for i, j in zip(np.arange(len(sample_draws))[mask], sample_draws[mask]):
        if j not in counts:
            counts[j] = 2
        else:
            inc_counts[i] = counts[j]
            counts[j] += 1
    return inc_counts


def collapse_draws_to_counts(sample_draws, num_features):
    """Produces a count of items in the sample draws.

    Takes an array of samples drawn from a distribution and returns a new
    array of counts of how many times each feature was seen in the draws.

    Args:
      sample_draws (np.ndarray(int)): A vector of sample indices.

    Returns:
      counts (np.ndarray(int): A vector of counts for each feature index.

    """
    return np.bincount(sample_draws, minlength=num_features)


def eval_multinomial_loglikelihoods(matrix, logp, n=None):
    """Computes the multinomial log-likelihood for many barcodes.

    Multinomial log-likelihood for a single barcode where the count of UMIs for
    feature i is x_i is:
    l = log(gamma(sum_i(x_i) + 1)) - sum_i(log(gamma(x_i + 1)) + sum_i(log(p_i) * x_i)

    Because the input matrix is sparse and zero counts do not contribute to the
    likelihood, we can rapidly remove zero elements from the computation.

    Args:
      matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
      logp (np.ndarray(float)): The natural log of the multinomial probability
        vector across features
      n (int, optional): The total number of UMIs per barcode.  Saves computation
        if it can be precomputed.

    Returns:
      loglk (np.ndarray(float)): Log-likelihood for each barcode
    """
    num_bcs = matrix.shape[1]
    loglk = np.zeros(num_bcs, dtype=float)
    if n is None:
        n = np.asarray(matrix.sum(axis=0))[0]

    consts = gammaln(n + 1)
    for i in range(num_bcs):
        idx_start, idx_end = matrix.indptr[i], matrix.indptr[i + 1]
        idxs = matrix.indices[idx_start:idx_end]
        row = matrix.data[idx_start:idx_end]
        short_logp = logp.take(idxs, axis=None, mode="clip")
        loglk[i] = consts[i] - gammaln(row + 1).sum() + (row * short_logp).sum()
    return loglk


def eval_multinomial_loglikelihood_cumulative(sample_draws, logp):
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

    Returns:
      loglk (np.ndarray(float): The cumulative log-likelihood up to each draw
    """
    marginal_counts = incremental_counts_from_sample(sample_draws)
    nvals = np.arange(1, len(sample_draws) + 1)
    loglk = np.log(nvals) - np.log(marginal_counts) + logp[sample_draws]
    return np.cumsum(loglk)


def eval_dirichlet_multinomial_loglikelihoods(matrix, alpha, n=None):
    """Computes the Dirichlet-multinomial log-likelihood for many barcodes.

    Dirichlet-Multinomial log-likelihood for a single barcode where the count
    of UMIs for feature i is x_i is:
    l = log(sum_i(x_i)) + log(beta(sum_i(a_i), sum_i(x_i))) -
        sum_i(log(x_i)) - sum_i(log(beta(a_i, x_i)))

    Because the input matrix is sparse and zero counts do not contribute to the
    likelihood, we can rapidly remove zero elements from the computation.

    Args:
      matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
      alpha (np.ndarray(float)): The vector of Dirichlet parameters for each feature
      n (int, optional): The total number of UMIs per barcode.  Saves computation
        if it can be precomputed.

    Returns:
      loglk (np.ndarray(float)): Log-likelihood for each barcode
    """
    num_bcs = matrix.shape[1]
    loglk = np.zeros(num_bcs)
    if n is None:
        n = np.asarray(matrix.sum(axis=0))[0]
    consts = np.log(n) + betaln(np.sum(alpha), n)
    for bc_index in range(matrix.shape[1]):
        idx_start, idx_end = matrix.indptr[bc_index], matrix.indptr[bc_index + 1]
        idxs = matrix.indices[idx_start:idx_end]
        row = matrix.data[idx_start:idx_end]
        short_alpha = alpha.take(idxs, axis=None, mode="clip")
        loglk[bc_index] = consts[bc_index] - np.log(row).sum() - betaln(short_alpha, row).sum()
    return loglk


def eval_dirichlet_multinomial_loglikelihood_cumulative(sample_draws, alpha):
    """Computes the cumulative Dirichlet-multinomial log-likelihood for a vector.

    Given a vector of sample draws from a Dirichlet-multinomial distribution, computes the
    log-likelihood of all of the draws up to each draw.  This is done incrementally,
    by first pre-computing the rank R of each draw as the number of times that
    feature was seen when it was drawn.  Then the incremental log-likelihood from
    that draw is

    dl = log(i) - log(R_i) - log(i + sum_i(a_i) - 1) + log(R_i + a_i - 1)

    Args:
      sample_draws (np.ndarray(int)): A vector of sample draws
      alpha (np.ndarray(float)): The Dirichlet parameter for each feature

    Returns:
      loglk (np.ndarray(float): The cumulative log-likelihood up to each draw
    """
    marginal_counts = incremental_counts_from_sample(sample_draws)
    nvals = np.arange(1, len(sample_draws) + 1)
    alpha_0 = np.sum(alpha)
    loglk = (
        np.log(nvals)
        - np.log(marginal_counts)
        - np.log(nvals + alpha_0 - 1)
        + np.log((marginal_counts - 1) + alpha[sample_draws])
    )
    return np.cumsum(loglk)


def estimate_dirichlet_overdispersion(matrix, ambient_bcs, p):
    """Estimates the best-fit overdispersion parameter for data.

    Uses a Dirichlet-multinomial to maximize the log-likelihood of the ambient
    barcode signal, given an input probability per feature that is scaled by
    a fixed overdispersion parameter.

    Args:
      matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
      ambient_bcs (np.array): Array of barcode indexes to use for the ambient background
      p (np.ndarray(float)): The estimated multinomial probability vector across features

    Returns:
      alpha (float): Best-fit overdispersion for the data
    """
    matrix = matrix[:, ambient_bcs]

    umis_per_bc = np.asarray(matrix.sum(axis=0)).flatten()

    def ambient_loglk(alpha, matrix, p, umis_per_bc):
        # We return the negative so function minimization gets the max likelihood
        return np.sum(-eval_dirichlet_multinomial_loglikelihoods(matrix, alpha * p, n=umis_per_bc))

    bounds = [0.001, 10000]
    result = optimize.minimize_scalar(
        ambient_loglk,
        bounds=bounds,
        args=(matrix, p, umis_per_bc),
    )
    if not result.success:
        raise ValueError(f"Could not find valid alpha: {result.message}")
    else:
        print(
            f"Alpha = {result.x}: maximizes the likelihood of the ambient barcodes. Search bounds: {bounds}"
        )
    return result.x


def draw_multinomial_sample(num_draws, p_cumulative):
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
      sample_draws (np.ndarray(float)): A random set of draws from a multinomial
        distribution.
    """
    rng_nums = RNG.random(size=num_draws)
    return np.searchsorted(p_cumulative, rng_nums)


def draw_dirichlet_multinomial_sample(num_draws, alpha):
    """Returns an array of Dirichlet-multinomial draws from the input parameters.

    Given a number of samples to draw and a vector of alpha parameters for all
    features to sample from, produces an array of feature indices of the desired
    size, drawn according to a Dirichlet-multinomial distribution.  Note that
    Dirichlet-multinomial draws are *not* independent and identically-distributed,
    so that all draws for a sample must be done with a single call or the
    sample variance will be under-stated.

    Args:
      num_draws (int): The number of samples to draw from the distribution
      alpha (np.ndarray(float)): The Dirichlet parameters for each feature.
      Normalized to a sum of 1, these are the mean probabilities of each feature
      over many draws of many samples.  The length of this array should be the
      number of features.

    Returns:
      sample_draws (np.ndarray(float)): A random set of draws from a Dirichlet-
        multinomial distribution.
    """
    probs = RNG.dirichlet(alpha)
    return draw_multinomial_sample(num_draws, np.cumsum(probs))


def simulate_multinomial_loglikelihoods(p, umis_per_bc, num_sims):
    """Simulate draws from a multinomial distribution for many values of N.

    Note that the samples within each simulation are not independent; we generate
    N draws for the largest value of N and use subsamples of the largest simulation
    for smaller values of N.

    Args:
      p (np.ndarray(float)): Probability of observing each feature.
      umis_per_bc (np.ndarray(int)): UMI counts per barcode.
      num_sims (int): Number of simulations to perform.

    Returns:
      distinct_ns (np.ndarray(int)): an array containing the distinct N values
          that were simulated.
      log_likelihoods (np.ndarray(float)): a len(distinct_ns) x num_sims matrix
          containing the simulated log likelihoods.
    """
    distinct_n = np.flatnonzero(np.bincount(umis_per_bc))
    max_n = max(distinct_n)
    loglk = np.zeros((len(distinct_n), num_sims), dtype=float)
    p_cumulative = np.cumsum(p)
    logp = np.log(p)

    for i in range(num_sims):
        draw = draw_multinomial_sample(max_n, p_cumulative)
        # Note we can use the distinct N values as indices to extract the log-likelihood
        # at each N value.
        loglk[:, i] = eval_multinomial_loglikelihood_cumulative(draw, logp)[distinct_n - 1]
    return distinct_n, loglk


def simulate_dirichlet_multinomial_loglikelihoods(alpha, umis_per_bc, num_sims):
    """Simulate draws from a Dirichlet-multinomial distribution for many values of N.

    Note that the samples within each simulation are not independent; we generate
    N draws for the largest value of N and use subsamples of the largest simulation
    for smaller values of N.

    Args:
      alpha (np.ndarray(float)): Dirichlet parameter for each feature.
      umis_per_bc (np.ndarray(int)): UMI counts per barcode.
      num_sims (int): Number of simulations to perform.

    Returns:
      distinct_ns (np.ndarray(int)): an array containing the distinct N values
          that were simulated.
      log_likelihoods (np.ndarray(float)): a len(distinct_ns) x num_sims matrix
          containing the simulated log likelihoods.
    """
    distinct_n = np.flatnonzero(np.bincount(umis_per_bc))
    max_n = max(distinct_n)
    loglk = np.zeros((len(distinct_n), num_sims), dtype=float)

    for i in range(num_sims):
        draw = draw_dirichlet_multinomial_sample(max_n, alpha)
        # Note we can use the distinct N values as indices to extract the log-likelihood
        # at each N value.
        loglk[:, i] = eval_dirichlet_multinomial_loglikelihood_cumulative(draw, alpha)[
            distinct_n - 1
        ]
    return distinct_n, loglk


def compute_ambient_pvalues(umis_per_bc, obs_loglk, sim_n, sim_loglk):
    """Compute p-values for observed multinomial log-likelihoods.

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

    for i in range(num_barcodes):
        num_lower_loglk = np.sum(sim_loglk[sim_n_idx[i], :] < obs_loglk[i])
        pvalues[i] = float(1 + num_lower_loglk) / (1 + num_sims)
    return pvalues


class Curve:
    """Curve information for plotting.

    Attributes:
        x     (list of int): Curve x-axis (number of cells)
        y     (list of float): Curve y-axis (expected value of number of
              unique clonotypes)
        y_std (list of float): Standard deviation of number of unique
              clonotypes in a curve
        y_ciu (list of float): Upper bound of 95% confidence interval
              for number of unique clonotypes in a curve
        y_cil (list of float): Lower bound of 95% confidence interval
              for number of unique clonotypes in a curve
    """

    def __init__(self):
        """Initialize an empty curve."""
        self.x: list[int] = []
        self.y: list[float] = []
        self.y_std: list[float] | list[None] = []
        self.y_ciu: list[float] | list[None] = []
        self.y_cil: list[float] | list[None] = []

    def is_empty(self):
        """Check if curve is empty or calculated.

        Returns:
            bool
        """
        if len(self.y) == 0 or len(self.x) == 0:
            return True
        return False

    def is_consistent(self):
        """Check if calculated curve is consistent (length match).

        Returns:
            bool
        """
        len_list = [len(i) for i in [self.x, self.y, self.y_cil, self.y_ciu]]
        if len_list.count(len_list[0]) != len(len_list):
            return False
        return True


class Diversity:
    """Represents diversity with rarefaction and extrapolation curves.

    Attributes:
        sorted_hist         (list of tuple): A histogram sorted by abundance
                                            (most abundant has lowest index)
        freq_counts         (dict of int: int): Frequency counts table
        N                   (int): total number of counts in histogram
        rarefaction_curve   (Curve): Rarefaction curve information
        extrapolation_curve (Curve): Extrapolation curve information
    """

    def __init__(self, hist):
        self.sorted_hist = sorted(hist, reverse=True)
        self.N: int = self._calc_n()
        self.freq_counts: dict[int, int] = self._get_freq_counts()
        # Rarefaction curve
        self.rarefaction_curve: Curve = Curve()
        # Extrapolation curve
        self.extrapolation_curve: Curve = Curve()

    def is_diversity_curve_possible(self):
        if self.N == 0 or self.f_0_chao1() in [0, -1]:
            return False
        return True

    def _get_freq_counts(self) -> dict[int, int]:
        """Get frequency counts (f_k in Colwell et al 2012).

        Returns:
            dict
        """
        ret_dict: dict[int, int] = {}
        for abund in self.sorted_hist:
            if abund not in ret_dict:
                ret_dict[abund] = 1
            else:
                ret_dict[abund] += 1
        return ret_dict

    def _calc_n(self) -> int:
        """Calculates number of samples.

        (length of input histogram)

        Returns:
            int
        """
        return sum(self.sorted_hist)

    def _alpha(self, n, k):
        """Calculates alpha parameter in Colwell et al 2012.

        Args:
            n (int): rarefaction point n
            k (int): Frequency

        Returns:
            int
        """
        if k > self.N - n:
            return 0
        log_alpha = (
            loggamma(self.N - k + 1)
            + loggamma(self.N - n + 1)
            - loggamma(self.N + 1)
            - loggamma(self.N - k - n + 1)
        )
        return np.exp(np.real(log_alpha))

    # Minimum variance unbiased estimator (MVUE) for both hypergeometric
    # and multinomial models
    # eq (4) in Colwell et al 2012
    # eq (5) in Colwell et al 2012
    def _rarefaction(self, n):
        """Calculates rarefaction and standard deviation of rarefaction for a single point n.

        Args:
            n (int): rarefaction point n

        Returns:
            float, float
        """
        # number of observed clonotypes (S_obs in paper)
        num_clon = float(sum(self.freq_counts.values()))
        sum_helper_exp_val = 0
        for k, count in self.freq_counts.items():
            sum_helper_exp_val += self._alpha(n, k) * count
        exp_val = num_clon - sum_helper_exp_val
        sum_helper_std_dev = 0
        for k, count in self.freq_counts.items():
            sum_helper_std_dev += (1 - self._alpha(n, k)) ** 2 * count - float(
                exp_val**2
            ) / self.assemblage_size_estimate()
        return exp_val, np.sqrt(sum_helper_std_dev)

    def calc_rarefaction_curve(self, num_steps: int = 40):
        """Calculates rarefaction curve for num_steps between 1 and N.

        Args:
            num_steps (int): Number of steps between 1 and N
        """
        step_size = int(self.N // (num_steps - 1))
        if step_size == 0:
            step_size = 1
        rc_x = range(1, self.N, step_size)
        if len(rc_x) < num_steps:
            rc_x = range(1, self.N + step_size, step_size)

        placeholder = 0.0
        rc_y_exp = [placeholder] * len(rc_x)
        rc_y_std = [placeholder] * len(rc_x)
        rc_y_ciu = [placeholder] * len(rc_x)
        rc_y_cil = [placeholder] * len(rc_x)
        for idx, n in enumerate(rc_x):
            rc_y_exp[idx], rc_y_std[idx] = self._rarefaction(n)
            rc_y_cil[idx] = rc_y_exp[idx] - 1.96 * rc_y_std[idx]
            rc_y_ciu[idx] = rc_y_exp[idx] + 1.96 * rc_y_std[idx]
        self.rarefaction_curve.x = list(rc_x)
        self.rarefaction_curve.y = rc_y_exp
        self.rarefaction_curve.y_std = rc_y_std
        self.rarefaction_curve.y_cil = rc_y_cil
        self.rarefaction_curve.y_ciu = rc_y_ciu

    def calc_rarefaction_curve_plotly(
        self, origin: str, color, num_steps: int = 40, stdev: bool = True
    ):
        """Create a plotly curve for rarefaction.

        Args:
            origin (str): Origin string (see vdj inputs)
            color: color of the curve
            num_steps (int): number of extrapolation steps
            stdev (bool): flag for calculation of standard deviation (and confidence intarvals)

        Returns:
            [dict]
        """
        self.calc_rarefaction_curve(num_steps)
        lines = []
        lines.append(
            {
                "x": self.rarefaction_curve.x,
                "y": self.rarefaction_curve.y,
                "type": "scatter",
                "name": origin,
                "mode": "lines",
                "line_color": color,
            }
        )
        if stdev:
            lines.append(
                {
                    "x": self.rarefaction_curve.x + self.rarefaction_curve.x[::-1],
                    "y": self.rarefaction_curve.y_ciu + self.rarefaction_curve.y_cil[::-1],
                    "type": "scatter",
                    "name": origin + " (Stdev)",
                    "fill": "toself",
                    "line_color": "rgba(255,255,255,0)",
                    "fillcolor": color,
                    "opacity": 0.2,
                }
            )
        return lines

    # eq (9) in Colwell et al 2012
    # not using the approximation used in eq (9)
    def _extrapolation(self, n_plus_m):
        """Calculates extrapolation for a single point (N + m).

        standard deviation is not implemented.

        Args:
            n_plus_m (int): Distance of extrapolation point to N
                (how many more samples collected)

        Returns:
            float: Expected value
            None: Standard deviation (not implemented)
        """
        # number of observed clonotypes (s_obs in paper)
        s_obs = float(sum(self.freq_counts.values()))
        f_0 = float(self.f_0_chao1())
        f_1 = float(self.freq_counts[1])

        N = float(self.N)  # n in paper
        m = float(n_plus_m - self.N)

        brackets = 1.0 - (1.0 - f_1 / N / f_0) ** m

        exp_val = s_obs + f_0 * brackets
        # Standard deviation for extrapolation is not implemented
        std_dev = None

        return exp_val, std_dev

    def calc_extrapolation_curve(self, num_steps: int = 40, max_n: int = 2000):
        """Calculates extrapolation curve for num_steps between N and max_n.

        Args:
            num_steps (int): Number of steps between N and max_n
            max_n (int): The limit to which the function extrapolates

        Returns:
            None
        """
        step_size = int((max_n - self.N) / (num_steps - 1))
        ec_x = range(self.N, max_n, step_size)
        if len(ec_x) < num_steps:
            ec_x = range(self.N, max_n + step_size, step_size)
        ec_y_exp = [0.0] * len(ec_x)
        ec_y_std = [None] * len(ec_x)
        ec_y_ciu = [None] * len(ec_x)
        ec_y_cil = [None] * len(ec_x)
        for idx, n_plus_m in enumerate(ec_x):
            ec_y_exp[idx], ec_y_std[idx] = self._extrapolation(n_plus_m)
            if ec_y_std[idx] is not None:
                ec_y_cil[idx] = ec_y_exp[idx] - 1.96 * ec_y_std[idx]
                ec_y_ciu[idx] = ec_y_exp[idx] + 1.96 * ec_y_std[idx]
        self.extrapolation_curve.x = list(ec_x)
        self.extrapolation_curve.y = ec_y_exp
        self.extrapolation_curve.y_std = ec_y_std
        self.extrapolation_curve.y_cil = ec_y_cil
        self.extrapolation_curve.y_ciu = ec_y_ciu

    def f_0_chao1(self):
        """Estimate f0 (number of clonotypes that we haven't seen yet) using chao1 estimator.

        Based on: Colwell et al 2012

        Returns:
            float
        """
        if 1 in self.freq_counts and 2 in self.freq_counts:
            return float(self.freq_counts[1] ** 2) / float(2 * self.freq_counts[2])
        elif 1 in self.freq_counts:
            return float(self.freq_counts[1] * (self.freq_counts[1] - 1)) / 2.0
        else:
            # print "Error in f_0_chao1" # TODO log message
            return -1.0

    def assemblage_size_estimate(self):
        """Estimate assemblage size (number of total unique clonotypes) using chao1 estimator.

        Based on: Colwell et al 2012

        Returns:
            float
        """
        return self.N + self.f_0_chao1()

    def create_both_curves(
        self, outfile, rc_steps: int = 40, ec_steps: int = 40, max_n: int = 2000
    ):
        """Calculate rarefaction and extrapolation diversity curves.

        Save them in csv format in outfile.

        Args:
            outfile (str): Path to output file
            rc_steps (int): Rarefaction curve steps between 1 and N
            ec_steps (int): Extrapolation steps between N and max_n
            max_n (int): Maximum number of cells extrapolating to
        """
        # Check if extrapolation_curve is calculated. If not, call the appropriate func to calc
        if self.extrapolation_curve.is_empty():
            self.calc_extrapolation_curve(ec_steps, max_n)
        # Check if rarefactoion_curve is calculated. If not, call the appropriate func to calc
        if self.rarefaction_curve.is_empty():
            self.calc_rarefaction_curve(rc_steps)
        if ~self.rarefaction_curve.is_consistent():
            return None  # TODO add appropriate logging message
        if ~self.extrapolation_curve.is_consistent():
            return None  # TODO add appropriate logging message
        header = "x,y,y_cil,y_ciu,type\n"
        with open(outfile, "w") as outhandle:
            outhandle.write(header)
            for idx, x in enumerate(self.rarefaction_curve.x):
                numbers = [
                    x,
                    self.rarefaction_curve.y[idx],
                    self.rarefaction_curve.y_cil[idx],
                    self.rarefaction_curve.y_ciu[idx],
                ]
                line = ",".join([str(n) for n in numbers] + ["rarefaction"]) + "\n"
                outhandle.write(line)
            for idx, x in enumerate(self.extrapolation_curve.x):
                numbers = [
                    x,
                    self.extrapolation_curve.y[idx],
                    self.extrapolation_curve.y_cil[idx],
                    self.extrapolation_curve.y_ciu[idx],
                ]
                line = ",".join([str(n) for n in numbers] + ["extrapolation"]) + "\n"
                outhandle.write(line)
        return 1
