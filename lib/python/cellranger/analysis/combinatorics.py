# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Some useful combinatorics functions."""

from __future__ import annotations

from collections.abc import Generator

import numpy as np
from scipy.special import comb


# pylint: disable=invalid-name
def multinomial_comb(cnts):
    """Number of N choose k1, k2, ...k_n combinations possible.

    from a vector of counts k_i summing to N
    """
    N = np.sum(cnts)
    numerator = np.math.factorial(N)
    denominator = np.array([np.math.factorial(x) for x in cnts])
    return numerator / np.prod(denominator)


def calc_all_possible_nonnegative_solutions(n, r):
    """Calculate number of solutions.

    Solutions to::

        x_1 + x_2 + x_3 + ... x_n = r
    """
    N = n + r - 1
    k = r
    return comb(N, k)


def calc_all_possible_positive_solutions(n, r):
    """Calculate number of integer solutions.

    Solutions to:

        x_1 + x_2 + x_3 + ... x_n = r

    with all x_i > 0

    Args:
        n: number of elements to choose from
        r: The total
    """
    return calc_all_possible_nonnegative_solutions(n, r - n)


def generate_all_multiplets(n_tags: int, max_multiplets: int, add_unit_vector_at_end: bool = False):
    """Effectively generate all solutions.

    Solutions to::

        x_1 + x_2 + x_3 + ... = k for k = 1..max_multiplets
    """
    solutions: list[list[int]] = []
    for k in range(max_multiplets + 1):
        cur_solutions = list(generate_solutions(n_tags, k))
        assert int(calc_all_possible_nonnegative_solutions(n_tags, k)) == len(cur_solutions)
        solutions += cur_solutions
    # only relevant to JIBES
    if add_unit_vector_at_end:
        solutions.append([1] * n_tags)
    return solutions


def generate_all_multiplets_as_array(
    n_tags: int, max_multiplets: int, add_unit_vector_at_end: bool = False
) -> np.ndarray[tuple[int, int], np.dtype[np.int32]]:
    """Helper method for pyo3.

    Effectively generate all solutions to::

        x_1 + x_2 + x_3 + ... = k for k = 1..max_multiplets
    """
    solutions = generate_all_multiplets(n_tags, max_multiplets, add_unit_vector_at_end)
    return np.array(solutions, dtype=np.int32)


def generate_solutions(
    elements: int, k: int, only_positive: bool = False
) -> Generator[list[int], None, None]:
    """A generator with all integer solutions to x_1 + x_2 + ... x_n = k.

    Args:
        elements: Number of variables that can be added
        k: total
        only_positive: Only return solutions where all x_i are > 0
    """
    assert elements > 0

    if elements == 1:
        yield [k]
    else:
        start = 0
        end = k + 1
        if only_positive:
            start = 1
            end = k
        for i in range(start, end):
            current = [i]
            for val in generate_solutions(elements - 1, k - i, only_positive):
                yield current + val
