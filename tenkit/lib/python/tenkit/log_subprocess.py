#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Helper for functions for logging calls to subprocesses
#

import subprocess
import logging

def check_call(*args, **kwargs):
    logging.log(logging.INFO, "subprocess check_call: %s" % " ".join(*args))
    return subprocess.check_call(*args, **kwargs)

def call(*args, **kwargs):
    logging.log(logging.INFO, "subprocess call: %s" % " ".join(*args))
    return subprocess.call(*args, **kwargs)
