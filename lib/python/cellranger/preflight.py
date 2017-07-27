#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import socket
import subprocess
import cellranger.constants as cr_constants

def check_refdata(reference_path):
    hostname = socket.gethostname()

    print "Checking reference_path (%s) on %s..." % (reference_path, hostname)

    required_files = [
        cr_constants.REFERENCE_METADATA_FILE,
        cr_constants.REFERENCE_FASTA_PATH,
        cr_constants.REFERENCE_GENES_GTF_PATH,
        cr_constants.REFERENCE_GENES_INDEX_PATH,
    ]
    for filename in required_files:
        p = os.path.join(reference_path, filename)
        if not os.path.isfile(p):
            return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname

    fasta_path = os.path.join(reference_path, cr_constants.REFERENCE_FASTA_PATH)
    for filename in cr_constants.STAR_REQUIRED_FILES:
        p = os.path.join(reference_path, cr_constants.REFERENCE_STAR_PATH, filename)
        if not os.path.exists(p):
            return False, "Your reference doesn't appear to be indexed. Please run the mkreference tool"

        # catch if fasta has been modified - buffer accounts for potential changed modification times based on copy/tarball extraction order
        mtime_buffer = 60
        if os.path.getmtime(p) + mtime_buffer < os.path.getmtime(fasta_path):
             return False, "Timestamps on the STAR index files have changed. Please reindex your reference using the mkreference tool"

    return True, None

def record_package_versions():
    for package in cr_constants.PACKAGE_VERSION_CMDS:
        name = package['name']
        cmd = package['cmd']

        version = subprocess.check_output(cmd, shell=True)
        print '%s: %s' % (name, version)
