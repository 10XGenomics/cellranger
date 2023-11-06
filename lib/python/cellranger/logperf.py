#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

# Performance logging
import resource
import time


class LogPerf:
    """Print before/after maxrss and elapsed time for code blocks to stdout."""

    def __init__(self, note) -> None:
        self.note = note
        self.start: int = None

    def __enter__(self) -> None:
        self.start = time.monotonic_ns()
        print(
            self.note,
            f"self_postmaxrss_mb\t{resource.getrusage(resource.RUSAGE_SELF)[2] / 1e3:0.1f}",
            sep="\t",
        )
        print(
            self.note,
            f"children_postmaxrss_mb\t{resource.getrusage(resource.RUSAGE_CHILDREN)[2] / 1e3:0.1f}",
            sep="\t",
            flush=True,
        )

    def __exit__(self, e_type, e_value, e_trace):
        end = time.monotonic_ns()
        assert self.start is not None
        print(
            self.note,
            f"self_postmaxrss_mb\t{resource.getrusage(resource.RUSAGE_SELF)[2] / 1e3:0.1f}",
            sep="\t",
        )
        print(
            self.note,
            f"children_postmaxrss_mb\t{resource.getrusage(resource.RUSAGE_CHILDREN)[2] / 1e3:0.1f}",
            sep="\t",
        )
        print(
            self.note,
            f"elapsed_sec\t{(end - self.start)/1e9:.0f}",
            sep="\t",
            flush=True,
        )

    @staticmethod
    def mem() -> None:
        print("self_maxrss_mb\t%0.1f" % (resource.getrusage(resource.RUSAGE_SELF)[2] / 1e3))
        print(
            "children_maxrss_mb\t%0.1f" % (resource.getrusage(resource.RUSAGE_CHILDREN)[2] / 1e3),
            flush=True,
        )
