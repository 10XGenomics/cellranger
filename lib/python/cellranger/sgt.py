#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

"""Simple Good-Turing estimator.

Based on S implementation in::

  William A. Gale & Geoffrey Sampson (1995) Good-turing frequency estimation without tears,
  Journal of Quantitative Linguistics, 2:3, 217-237, DOI: 10.1080/09296179508590051
"""

from __future__ import annotations

import numpy as np
import scipy.stats as sp_stats


class SimpleGoodTuringError(Exception):
    pass


def _averaging_transform(r, nr):
    d = np.concatenate((np.ones(1, dtype=int), np.diff(r)))
    dr = np.concatenate((0.5 * (d[1:] + d[0:-1]), np.array((d[-1],), dtype=float)))
    return nr.astype(float) / dr


def _rstest(r, coef):
    return r * np.power(1 + 1.0 / r, 1 + coef)


def simple_good_turing(xr: np.ndarray, xnr: np.ndarray):
    """Make a Simple Good-Turing estimate of the frequencies.

    Args:
      xr (np.array(int)): Non-zero item frequencies
      xnr (np.array(int)): Non-zero frequencies of frequencies
    Returns:
      (rstar (np.array(float)), p0 (float)):
        rstar: The adjusted non-zero frequencies
        p0: The total probability of unobserved items
    """
    xr = xr.astype(float)
    xnr = xnr.astype(float)

    xN = np.sum(xr * xnr)

    # Get Linear Good-Turing estimate
    xnrz = _averaging_transform(xr, xnr)
    slope, _, _, _, _ = sp_stats.linregress(np.log(xr), np.log(xnrz))

    if slope < -1:
        print(f"The SGT slope is {slope}.")
    else:
        # The Simple Good-Turing (SGT) method is not applicable when its slope is larger than -1,
        # in which case use a slope of -1, which results in r* = r, and which is the MLE.
        # Use Turing estimates at low frequencies, and switch to MLE estimates for higher frequencies.
        slope = -1
        print("The SGT slope {slope} is > -1. The SGT estimator is not applicable to these data.")

    xrst = _rstest(xr, slope)
    xrstrel = xrst / xr

    # Get traditional Good-Turing estimate
    xrtry = xr == np.concatenate((xr[1:] - 1, np.zeros(1)))
    xrstarel = np.zeros(len(xr))
    xrstarel[xrtry] = (
        (xr[xrtry] + 1) / xr[xrtry] * np.concatenate((xnr[1:], np.zeros(1)))[xrtry] / xnr[xrtry]
    )

    # Determine when to switch from GT to LGT estimates
    tursd = np.ones(len(xr))
    for i in range(len(xr)):
        if xrtry[i]:
            tursd[i] = float(i + 2) / xnr[i] * np.sqrt(xnr[i + 1] * (1 + xnr[i + 1] / xnr[i]))

    xrstcmbrel = np.zeros(len(xr))
    useturing = True
    for r in range(len(xr)):
        if not useturing:
            xrstcmbrel[r] = xrstrel[r]
        elif np.abs(xrstrel[r] - xrstarel[r]) * (1 + r) / tursd[r] > 1.65:
            xrstcmbrel[r] = xrstarel[r]
        else:
            useturing = False
            xrstcmbrel[r] = xrstrel[r]

    # Renormalize the probabilities for observed objects
    sumpraw = np.sum(xrstcmbrel * xr * xnr / xN)

    xrstcmbrel = xrstcmbrel * (1 - xnr[0] / xN) / sumpraw
    p0 = xnr[0] / xN

    return (xr * xrstcmbrel, p0)


def sgt_proportions(frequencies: np.ndarray):
    """Use Simple Good-Turing estimate to adjust for unobserved items.

    Args:
      frequencies (np.array(int)): Nonzero frequencies of items

    Returns:
        pstar (np.array[float]): The adjusted non-zero proportions
        p0 (float): The total probability of unobserved items
    """
    if len(frequencies) == 0:
        raise ValueError("Input frequency vector is empty")
    if np.count_nonzero(frequencies) != len(frequencies):
        raise ValueError("Frequencies must be greater than zero")

    freqfreqs = np.bincount(frequencies)
    assert freqfreqs[0] == 0
    use_freqs = np.flatnonzero(freqfreqs)

    if len(use_freqs) < 10:
        raise SimpleGoodTuringError(
            "Too few non-zero frequency items (%d). Aborting SGT." % len(use_freqs)
        )

    rstar, p0 = simple_good_turing(use_freqs, freqfreqs[use_freqs])

    # rstar contains the smoothed frequencies.
    # Map each original frequency r to its smoothed rstar.
    rstar_dict = dict(zip(use_freqs, rstar))

    rstar_sum = np.sum(freqfreqs[use_freqs] * rstar)
    rstar_i = np.fromiter((rstar_dict[f] for f in frequencies), dtype=float, count=len(frequencies))
    pstar = (1 - p0) * (rstar_i / rstar_sum)

    assert np.isclose(p0 + np.sum(pstar), 1)
    return (pstar, p0)
