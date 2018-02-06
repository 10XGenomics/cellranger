#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Helper for functions for logging calls to subprocesses
#

import subprocess
import logging
import sys

# On linux, provide a method to set PDEATHSIG on child processes.
if sys.platform.startswith('linux'):
    import ctypes
    import ctypes.util
    from signal import SIGKILL

    _LIBC = ctypes.CDLL(ctypes.util.find_library('c'))

    _PR_SET_PDEATHSIG = ctypes.c_int(1)  # <sys/prctl.h>

    def child_preexec_set_pdeathsig():
        """When used as the preexec_fn argument for subprocess.Popen etc,
        causes the subprocess to recieve SIGKILL if the parent process
        terminates."""
        zero = ctypes.c_ulong(0)
        _LIBC.prctl(_PR_SET_PDEATHSIG, ctypes.c_ulong(SIGKILL),
                    zero, zero, zero)
else:
    child_preexec_set_pdeathsig = None  # pylint: disable=invalid-name


def check_call(*args, **kwargs):
    logging.log(logging.INFO, "subprocess check_call: %s" % " ".join(*args))
    if not 'preexec_fn' in kwargs:
        kwargs['preexec_fn'] = child_preexec_set_pdeathsig
    return subprocess.check_call(*args, **kwargs)

def check_output(*args, **kwargs):
    logging.log(logging.INFO, "subprocess check_output: %s" % " ".join(*args))
    if not 'preexec_fn' in kwargs:
        kwargs['preexec_fn'] = child_preexec_set_pdeathsig
    return subprocess.check_output(*args, **kwargs)

def call(*args, **kwargs):
    logging.log(logging.INFO, "subprocess call: %s" % " ".join(*args))
    if not 'preexec_fn' in kwargs:
        kwargs['preexec_fn'] = child_preexec_set_pdeathsig
    return subprocess.call(*args, **kwargs)

def Popen(*args, **kwargs):
    logging.log(logging.INFO, "subprocess Popen: %s" % " ".join(*args))
    if not 'preexec_fn' in kwargs:
        kwargs['preexec_fn'] = child_preexec_set_pdeathsig
    return subprocess.Popen(*args, **kwargs)
