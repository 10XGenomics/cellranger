#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import socket
from cellranger.preflight import PreflightException
import cellranger.vdj.constants as vdj_constants

def check_refdata(reference_path, denovo):
    hostname = socket.gethostname()

    if reference_path is None and not denovo:
        raise PreflightException("Must specify --reference unless --denovo is specified.")

    if reference_path is None:
        return

    print "Checking reference_path (%s) on %s..." % (reference_path, hostname)

    required_files = [
        vdj_constants.REFERENCE_FASTA_PATH,
    ]

    for filename in required_files:
        p = os.path.join(reference_path, filename)

        if not os.path.isfile(p):
            raise PreflightException("Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname)

def check_chain(chain):
    if chain not in vdj_constants.CHAIN_TYPE_SPECS:
        raise PreflightException("Must specify --chain as one of: " + ", ".join(vdj_constants.CHAIN_TYPE_SPECS) + ".")
