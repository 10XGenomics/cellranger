# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
"""Python implementation of classes also implemented in Rust.

 This was the original implementation
before these classes were migrated to PyO3.  They are kept around to enable easier
prototyping of model changes and because they are occasionally passed as pickles.
"""
from __future__ import annotations

import numpy as np
import statsmodels.api as sm
from scipy.stats import norm
from six import ensure_str

import cellranger.feature.feature_assigner as cr_fa
import cellranger.feature.throughputs
from cellranger.analysis.combinatorics import generate_all_multiplets, multinomial_comb
from cellranger.analysis.jibes_data import JibesData
from tenkit.stats import robust_divide

DEFAULT_BLANK_PROB = 0.04
_MAX_K_LETS_TO_CONSIDER = 3


class JibesModelPy:
    """A class that stores parameters and can simulate from a JIBES model."""

    # pylint: disable=too-few-public-methods
    def __init__(
        self,
        backgrounds: list[float],
        foregrounds,
        std_devs,
        frequencies=None,
        blank_prob=DEFAULT_BLANK_PROB,
        n_gems=cellranger.feature.throughputs.N_G,
    ):
        num_tags = len(backgrounds)
        assert len(backgrounds) == num_tags
        assert len(foregrounds) == num_tags
        assert len(std_devs) == num_tags
        if frequencies:
            assert len(frequencies) == num_tags
            self.frequencies = frequencies
        else:
            self.frequencies = np.repeat(robust_divide(1, num_tags), num_tags)
        self.blank_freq = blank_prob
        self.num_tags = num_tags
        self.background: np.ndarray[int, np.dtype[np.float64]] = np.array(backgrounds)
        self.foreground: np.ndarray[int, np.dtype[np.float64]] = np.array(foregrounds)
        self.B = np.vstack((self.background, np.diag(self.foreground)))
        self.std_devs: np.ndarray[int, np.dtype[np.float64]] = np.array(std_devs)
        self.n_gems = n_gems
        if np.sum(np.isnan(self.B)) != 0:
            raise ValueError("NaN Detected in Coefficient Matrix")
        if np.sum(np.isnan(self.std_devs)):
            raise ValueError("NaN Detected in Std. Devs.")

    def simulate(self, N):
        """Simulate a set of observations given an input number of "observed" cells.

        :param N: How many cells would be recovered across all k-lets?
        :return: A tuple with matrix Y, giving the observed counts, and a matrix X, with the
                 simulated latent states
        """
        # TODO: Blank state not yet implemented
        cnts = [int(x) for x in cr_fa.get_multiplet_counts(N, self.n_gems)][: self.num_tags]
        total = np.sum(cnts)
        x = np.zeros((total, self.num_tags + 1))
        x[:, 0] = np.ones(total)
        cur_row = 0
        for index, cnt in enumerate(cnts):
            to_sample = (index + 1) * cnt
            samps = np.random.choice(self.num_tags, to_sample)
            start_row = cur_row
            for i in range(index + 1):
                cur_row = start_row
                for j in range(cnt):
                    x[cur_row, 1 + samps[j + i * j]] += 1
                    cur_row += 1
        means = np.matmul(x, self.B)
        Y = np.zeros((total, self.num_tags))
        for k in range(self.num_tags):
            Y[:, k] = np.random.normal(means[:, k], self.std_devs[k])
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        assert self.num_tags < len(alphabet)
        data = JibesData(Y, list(alphabet[: self.num_tags]), None)
        return (data, x)

    def foreground_as_vector(self):
        """Return the diagonal matrix as a vector of the diagnoals instead."""
        return np.diag(self.foreground)


class JibesEMPy:  # pylint: disable=too-many-instance-attributes
    """Class to perform EM fitting of JIBES model."""

    def __init__(
        self,
        data: JibesData,
        model: JibesModelPy,
        optimize_pop_freqs=False,
        max_k_lets=_MAX_K_LETS_TO_CONSIDER,
        n_gems=cellranger.feature.throughputs.N_G,
    ):
        """Fit a Joint Inference by Exploiting Stoichiometry.

        :param data:
        :param model:
        :param optimize_pop_freqs: Optimize tag populations (otherwise all equal)
        :param max_k_lets: Maximum k-lets to consider
        :param n_gems: Number of GEMs used
        """
        assert isinstance(model, JibesModelPy)
        assert isinstance(data, JibesData)
        self.optimize_pop_freqs = optimize_pop_freqs
        n = data.counts.shape[0]
        self.estimated_cells = cr_fa.calculate_expected_total_cells(n, n_gems=n_gems)
        self.max_k_let_setting = max_k_lets
        self.n_gems = n_gems
        exp_cnts = cr_fa.get_multiplet_counts(self.estimated_cells, n_gems=n_gems)
        # Only integrate up to a multiplet number
        # with expected count >= 1, as multiplets of higher integer values lead to
        # combinatorial explosions
        max_multiplets = np.max(np.nonzero(exp_cnts)) + 1
        max_multiplets = max(max_multiplets, 2)
        if max_multiplets > self.max_k_let_setting:
            # For this version of the algorithm, to save memory we'll have a "catch-all" category
            # that has all the states, and limit the set considered to _MAX_K_LETS_TO_CONSIDER,
            # from scanning data, it doesn't appear the doublets/triplets of one cell type are
            # that different from the singlets
            # TODO: Verify this assumption and add tests
            print(
                f"Limiting Latent States to K-lets of {self.max_k_let_setting} with a blank state"
            )
            self.k_let_limited = True
            max_multiplets = self.max_k_let_setting
        else:
            self.k_let_limited = False
        latent_states = generate_all_multiplets(model.num_tags, max_multiplets, self.k_let_limited)
        self.max_modeled_k_let = max(model.num_tags, max_multiplets)
        self.latent_states = np.array(latent_states)
        print("Total Latent States = " + str(len(self.latent_states)))
        print(f"Observed Barcodes = {n}")
        print(f"Inferred Cells = {self.estimated_cells}")
        self.data = data
        self.model = model
        ones = np.ones((self.z, 1))
        self.X = np.hstack((ones, self.latent_states))
        # pylint: disable=invalid-name
        self.data_full, self.X_full = self._reshape_for_regression()
        self.posterior = None
        self.LL = -np.inf
        self.converged = False
        self.iterations = 0

    @property
    def col_names(self):
        """Returns the names of the vector."""
        return self.data.column_names

    @property
    def n(self):
        """Number of observed barcodes."""
        return self.data.counts.shape[0]

    @property
    def z(self):
        """Number of latent states in the model."""
        return self.latent_states.shape[0]

    @property
    def k(self):
        """Number of observations per data point."""
        return self.model.num_tags

    def _calculate_latent_state_weights(self):
        """Calculate the probabilities for each type of singlet or multiplet.

        Based on our poisson expectation.
        """
        # cnts here starts at 1 cell, the 0 group is not included
        cnts = cr_fa.get_multiplet_counts_unrounded(self.estimated_cells, n_gems=self.n_gems)[
            : self.max_modeled_k_let
        ]
        p_k_let = cnts / np.sum(cnts)
        if self.k_let_limited:
            # We mush all the probability for the higher states into the last state that we want to account for all that data
            # Note no +1 here as it's 1 based array for multiple
            p_k_let[-1] = np.sum(p_k_let[self.max_k_let_setting :])
        x = self.X[:, 1:]
        klet = x.sum(axis=1).astype(np.int32)
        state = p_k_let[klet[1:] - 1]
        state = np.log(state)
        pis = np.log(self.model.frequencies)
        z = x.shape[0]
        comb_probs = np.zeros(z)
        comb_probs[0] = np.log(self.model.blank_freq)
        log_not_blank = np.log(1.0 - self.model.blank_freq)
        for zi in range(1, z):  # pylint: disable=invalid-name
            relevant_ps = np.nonzero(x[zi, :])
            cnts = x[zi, relevant_ps].flatten()
            ps = pis[relevant_ps]  # pylint: disable=invalid-name
            p = np.sum(cnts * ps)
            # TODO: Investigate precision issues when multiplier is large
            multiplier = float(multinomial_comb(cnts))
            comb_probs[zi] = p + np.log(multiplier) + state[zi - 1] + log_not_blank
        return comb_probs

    def _calculate_posterior_by_state(self):
        # Critical bit for the E step
        mu = np.matmul(self.X, self.model.B)
        z_prior = self._calculate_latent_state_weights()
        ll_posterior = np.tile(z_prior, (self.n, 1))
        for i in range(self.n):
            ll_posterior[i, :] += norm.logpdf(self.data.counts[i, :], mu, self.model.std_devs).sum(
                axis=1
            )
        ll_max = np.max(ll_posterior, axis=1)
        # TODO: There must be a better way than transposing here
        ll_posterior_altered = ll_posterior.T - ll_max
        posterior = np.exp(ll_posterior_altered.T)
        marginal = posterior.sum(axis=1, keepdims=True)
        self.posterior = posterior / marginal
        self.LL = np.log(marginal).sum() + ll_max.sum()

    def _maximize_parameters(self):
        """Maximize the parameters by weighted regression."""
        # The M-step
        W = self.posterior.flatten()
        new_std_devs = np.zeros(self.k)
        new_foregrounds = np.zeros(self.k)
        new_backgrounds = np.zeros(self.k)
        for k in range(self.k):
            res_wls = sm.WLS(self.data_full[:, k], self.X_full[:, [0, k + 1]], weights=W).fit()
            # pylint: disable=no-member
            residuals = res_wls.fittedvalues - np.squeeze(np.asarray(self.data_full[:, k]))
            # TODO: Ensure this is never 0
            var = np.sum(W * np.power(residuals, 2.0))
            new_std_devs[k] = np.sqrt(var / self.n)
            new_backgrounds[k] = res_wls.params[0]
            new_foregrounds[k] = res_wls.params[1]
        # Update frequencies of each population
        if self.optimize_pop_freqs:
            wght_cnts = np.matmul(W.reshape((self.n * self.z, 1)), self.X_full[:, 1:])
            total = wght_cnts.sum(axis=0)
            freqs = total / np.sum(total)
        else:
            freqs = None
        self.model = JibesModelPy(
            new_backgrounds,
            new_foregrounds,
            new_std_devs,
            freqs,
            blank_prob=self.model.blank_freq,
            n_gems=self.n_gems,
        )

    def _reshape_for_regression(self):
        # Make a new data matrix by replicating each row K times
        to_cat = []
        for i in range(self.n):
            to_cat.append(np.repeat(i, self.z))
        replicated_indices = np.concatenate(to_cat)
        data_new = self.data.counts[replicated_indices, :]
        X_full = np.tile(self.X, (self.n, 1))  # pylint: disable=invalid-name
        return (data_new, X_full)

    def one_EM_step(self):  # pylint: disable=invalid-name
        """Perform one step of the EM algorithm.

        :returns The current LL of the model.
        """
        if self.posterior is None:
            self._calculate_posterior_by_state()
        print("LL = " + str(self.LL))
        self._maximize_parameters()
        # Calculate this after to ensure we are current with the parameters
        self._calculate_posterior_by_state()
        self.iterations += 1
        return self.LL

    # Note these default parameters are mirrored in rust code and should be updated at the same time
    def perform_EM(  # pylint: disable=invalid-name
        self, max_reps=50000, abs_tol=1e-2, rel_tol=1e-7
    ):
        """Run the EM algorithm.

        Runs until either max_reps is exceeded or the change in  the LL after
        a run is <= tol or 1 - new/old <= rel_tol
        :param max_reps: Maximum number of repetitions
        :param abs_tol: Absolute difference in LL which would terminate the algorithm
        :param rel_tol: Relative difference in LL which would terminate the algorithm
        :return: The LL of the model after fitting
        """
        last_ll = self.LL
        rep = 0
        while True:
            self.one_EM_step()
            rep += 1
            rel_change = 1.0 - self.LL / last_ll
            abs_change = self.LL - last_ll
            if rep > max_reps:
                print(f"Did not converge in {max_reps} reps")
                break
            if not np.isneginf(last_ll) and ((abs_change <= abs_tol) or (rel_change <= rel_tol)):
                print("EM algorithm has converged")
                self.converged = True
                break
            last_ll = self.LL
        return self.LL

    def get_snr_dictionary(self):
        """Get a dictionary of "snr_TAG_NAME": snr_value.

        Returns:
            dict: Dictionary intended for json output
        """
        if not self.converged:
            raise ArithmeticError("Model must converge before SNRs can be reported.")
        snrs = self.model.foreground / self.model.std_devs
        return {"snr_" + ensure_str(x): y for x, y in zip(self.data.column_names, snrs)}
