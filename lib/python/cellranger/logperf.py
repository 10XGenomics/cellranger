#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

# Performance logging

import time
import resource
import sys

class LogPerf:
    """ Print before/after maxrss and elapsed time for code blocks to stdout """
    def __init__(self, note):
        self.note = note

    def __enter__(self):
        self.start = time.time()
        print '%s\tself_premaxrss_mb\t%0.1f' % (self.note, resource.getrusage(resource.RUSAGE_SELF)[2]/1e3)
        print '%s\tchildren_premaxrss_mb\t%0.1f' % (self.note, resource.getrusage(resource.RUSAGE_SELF)[2]/1e3)
        sys.stdout.flush()

    def __exit__(self, e_type, e_value, e_trace):
        print '%s\tself_postmaxrss_mb\t%0.1f' % (self.note, resource.getrusage(resource.RUSAGE_SELF)[2]/1e3)
        print '%s\tchildren_postmaxrss_mb\t%0.1f' % (self.note, resource.getrusage(resource.RUSAGE_CHILDREN)[2]/1e3)
        print '%s\telapsed_sec\t%d' % (self.note, time.time() - self.start)
        sys.stdout.flush()

    @staticmethod
    def mem():
        print 'self_maxrss_mb\t%0.1f' % (resource.getrusage(resource.RUSAGE_SELF)[2]/1e3)
        print 'children_maxrss_mb\t%0.1f' % (resource.getrusage(resource.RUSAGE_CHILDREN)[2]/1e3)
        sys.stdout.flush()
