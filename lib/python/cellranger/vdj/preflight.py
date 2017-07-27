#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import socket
import cellranger.vdj.constants as vdj_constants

def check_refdata(reference_path):
    hostname = socket.gethostname()

    print "Checking reference_path (%s) on %s..." % (reference_path, hostname)

    required_files = [
        vdj_constants.REFERENCE_FASTA_PATH,
    ]
    for filename in required_files:
        p = os.path.join(reference_path, filename)
        if not os.path.isfile(p):
            return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname

    return True, None
