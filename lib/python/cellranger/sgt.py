#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

"""
Simple Good-Turing estimator.
Based on S implementation in
  William A. Gale & Geoffrey Sampson (1995) Good-turing frequency estimation without tears,
  Journal of Quantitative Linguistics, 2:3, 217-237, DOI: 10.1080/09296179508590051
"""

import numpy as np
import scipy.stats as sp_stats
import itertools

class SimpleGoodTuringError(Exception):
    pass

def _averaging_transform(r, nr):
    d = np.concatenate((np.ones(1, dtype=int), np.diff(r)))
    dr = np.concatenate((
        0.5 * (d[1:] + d[0:-1]),
        np.array((d[-1],), dtype=float),
        ))
    return nr.astype(float)/dr

def _rstest(r, coef):
    return r * np.power(1 + 1/r, 1 + coef)

def simple_good_turing(xr, xnr):
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

    xN = np.sum(xr*xnr)

    # Get Linear Good-Turing estimate
    xnrz = _averaging_transform(xr, xnr)
    slope, intercept, _, _, _ = sp_stats.linregress(np.log(xr), np.log(xnrz))

    if slope > -1:
        raise SimpleGoodTuringError("The log-log slope is > -1 (%d); the SGT estimator is not applicable to these data." % slope)

    xrst = _rstest(xr,slope)
    xrstrel = xrst/xr

    # Get traditional Good-Turing estimate
    xrtry = xr == np.concatenate((xr[1:]-1, np.zeros(1)))
    xrstarel = np.zeros(len(xr))
    xrstarel[xrtry] = (xr[xrtry]+1) / xr[xrtry] * \
                      np.concatenate((xnr[1:], np.zeros(1)))[xrtry] / xnr[xrtry]

    # Determine when to switch from GT to LGT estimates
    tursd = np.ones(len(xr))
    for i in xrange(len(xr)):
        if xrtry[i]:
            tursd[i] = float(i+2) / xnr[i] * np.sqrt(xnr[i+1] * (1 + xnr[i+1]/xnr[i]))

    xrstcmbrel = np.zeros(len(xr))
    useturing = True
    for r in xrange(len(xr)):
        if not useturing:
            xrstcmbrel[r]  = xrstrel[r]
        else:
            if np.abs(xrstrel[r]-xrstarel[r]) * (1+r)/tursd[r] > 1.65:
                xrstcmbrel[r] = xrstarel[r]
            else:
                useturing = False
                xrstcmbrel[r] = xrstrel[r]

    # Renormalize the probabilities for observed objects
    sumpraw = np.sum(xrstcmbrel * xr * xnr / xN)

    xrstcmbrel = xrstcmbrel * (1 - xnr[0] / xN) / sumpraw
    p0 = xnr[0]/xN

    return (xr * xrstcmbrel, p0)

def sgt_proportions(frequencies):
    """Use Simple Good-Turing estimate to adjust for unobserved items

    Args:
      frequencies (np.array(int)): Nonzero frequencies of items
    Returns:
      (pstar (np.array(float)), p0 (float)):
        pstar: The adjusted non-zero proportions
        p0: The total probability of unobserved items
    """
    if len(frequencies) == 0:
        raise ValueError("Input frequency vector is empty")
    if np.count_nonzero(frequencies) != len(frequencies):
        raise ValueError("Frequencies must be greater than zero")

    freqfreqs = np.bincount(frequencies)
    assert freqfreqs[0] == 0
    use_freqs = np.flatnonzero(freqfreqs)

    if len(use_freqs) < 10:
        raise SimpleGoodTuringError("Too few non-zero frequency items (%d). Aborting SGT." % len(use_freqs))


    rstar, p0 = simple_good_turing(use_freqs, freqfreqs[use_freqs])

    # rstar contains the smoothed frequencies.
    # Map each original frequency r to its smoothed rstar.
    rstar_dict = dict(itertools.izip(use_freqs, rstar))

    rstar_sum = np.sum(freqfreqs[use_freqs] * rstar)
    rstar_i = np.fromiter((rstar_dict[f] for f in frequencies),
                          dtype=float, count=len(frequencies))
    pstar = (1 - p0) * (rstar_i / rstar_sum)

    assert np.isclose(p0 + np.sum(pstar), 1)
    return (pstar, p0)


def test_prosody():
    data = (
        (1, 120),
        (2, 40),
        (3, 24),
        (4, 13),
        (5, 15),
        (6, 5),
        (7, 11),
        (8, 2),
        (9, 2),
        (10, 1),
        (12, 3),
        (14, 2),
        (15, 1),
        (16, 1),
        (17, 3),
        (19, 1),
        (20, 3),
        (21, 2),
        (23, 3),
        (24, 3),
        (25, 3),
        (26, 2),
        (27, 2),
        (28, 1),
        (31, 2),
        (32, 2),
        (33, 1),
        (34, 2),
        (36, 2),
        (41, 3),
        (43, 1),
        (45, 3),
        (46, 1),
        (47, 1),
        (50, 1),
        (71, 1),
        (84, 1),
        (101, 1),
        (105, 1),
        (121, 1),
        (124, 1),
        (146, 1),
        (162, 1),
        (193, 1),
        (199, 1),
        (224, 1),
        (226, 1),
        (254, 1),
        (257, 1),
        (339, 1),
        (421, 1),
        (456, 1),
        (481, 1),
        (483, 1),
        (1140, 1),
        (1256, 1),
        (1322, 1),
        (1530, 1),
        (2131, 1),
        (2395, 1),
        (6925, 1),
        (7846, 1),
    )

    # Computed using R 3.5.1 w/ the Gale S code
    expect_p0 = 0.003883244
    expect_rstar = np.array((
        0.7628079,
        1.706448,
        2.679796,
        3.663988,
        4.653366,
        5.645628,
        6.63966,
        7.634856,
        8.63086,
        9.627446,
        11.62182,
        13.61725,
        14.61524,
        15.61336,
        16.6116,
        18.60836,
        19.60685,
        20.6054,
        22.60264,
        23.60133,
        24.60005,
        25.5988,
        26.59759,
        27.59639,
        30.59294,
        31.59183,
        32.59073,
        33.58964,
        35.58751,
        40.58235,
        42.58035,
        44.57836,
        45.57738,
        46.57641,
        49.57351,
        70.55399,
        83.54229,
        100.5272,
        104.5237,
        120.5097,
        123.507,
        145.4879,
        161.474,
        192.4472,
        198.4421,
        223.4205,
        225.4188,
        253.3947,
        256.3922,
        338.3218,
        420.2514,
        455.2215,
        480.2,
        482.1983,
        1138.636,
        1254.537,
        1320.48,
        1528.302,
        2128.788,
        2392.562,
        6918.687,
        7838.899,
    ))

    xr = np.array([d[0] for d in data], dtype=int)
    xnr = np.array([d[1] for d in data], dtype=int)

    rstar, p0 = simple_good_turing(xr, xnr)

    assert np.abs(p0 - expect_p0) < 1e-9
    assert np.all(np.abs(rstar - expect_rstar) < 1e-3)
    assert np.all((np.abs(rstar - expect_rstar))/expect_rstar < 1e-4)
