# Copyright (c) 2020 10x Genomics, Inc. All rights reserved.

"""A wrapper to avoid importing tables directly.

We need to monkey patch it on the fly to avoid issues with it changing thread
usage under the hood. See CELLRANGER-2936
"""


from __future__ import annotations

import sys

import numexpr
import tables

# Hold if we've queried martian and the current threads to use
_CACHED_MARTIAN_STATUS = [False, numexpr.utils.detect_number_of_cores(), 0]


def _get_n_threads():
    """Memoized function to query martian for thread usage,.

    assuming it's loaded.
    """
    if not _CACHED_MARTIAN_STATUS[0]:
        _CACHED_MARTIAN_STATUS[2] += 1
        try:
            import martian  # pylint: disable=bad-option-value,import-outside-toplevel

            n_threads = martian.get_threads_allocation()
            _CACHED_MARTIAN_STATUS[1] = n_threads
            tables.parameters.MAX_NUMEXPR_THREADS = n_threads
            _CACHED_MARTIAN_STATUS[0] = True
            sys.stdout.write(f"Martian runtime detected. Try Num. = {_CACHED_MARTIAN_STATUS[2]}\n")
        except:  # pylint: disable=bare-except
            sys.stdout.write(
                f"Warning: No martian runtime detected. Try Num. = {_CACHED_MARTIAN_STATUS[2]}\n"
            )
    return _CACHED_MARTIAN_STATUS[1]


def file_constructor_wrapper(function):
    """A decorator around the initialization function.

    to make sure it doesn't set the threads above what we want.
    """

    def wrapper(*args, **kwargs):  # pylint: disable=missing-docstring
        n_threads = _get_n_threads()
        kwargs["MAX_NUMEXPR_THREADS"] = n_threads
        kwargs["MAX_BLOSC_THREADS"] = n_threads
        return function(*args, **kwargs)

    return wrapper


# in case this value is used anywhere else

# specifically avoid any alternate settings by wrapping the constructor
tables.File.__init__ = file_constructor_wrapper(tables.File.__init__)
