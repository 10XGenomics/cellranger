# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
"""Do not increase the use of this module."""

MEM_GB_PER_THREAD = 8


def get_thread_request_from_mem_gb(mem_gb: float | int) -> int:
    """Do not increase the use of this function. Users should instead use the mrp option --mempercore.

    Return 1, 2, or 4 threads.
    """
    est_threads = round(mem_gb / MEM_GB_PER_THREAD)
    return next((threads for threads in (1, 2, 4) if threads >= est_threads), 4)
