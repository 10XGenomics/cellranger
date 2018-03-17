#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import socket
import tenkit.log_subprocess as tk_subproc
import tenkit.preflight as tk_preflight
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants

class PreflightException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def check(result):
    """ Translate (ok,msg) style from tenkit into an exception """
    ok, msg = result
    if not ok:
        raise PreflightException(msg)

def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True

def check_sample_def(sample_def):
    hostname = socket.gethostname()

    check(tk_preflight.check_gem_groups(sample_def))

    print "Checking FASTQ folder..."
    for sample_def in sample_def:
        read_path = sample_def["read_path"]
        if not read_path.startswith('/'):
            raise PreflightException("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            raise PreflightException("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            raise PreflightException("On machine: %s, cellranger does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            raise PreflightException("Specified FASTQ folder is empty: " + read_path)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not is_int(lane):
                    raise PreflightException("Lanes must be a comma-separated list of numbers.")

        check(tk_preflight.check_sample_indices(sample_def))

def check_refdata(reference_path):
    hostname = socket.gethostname()

    if reference_path is None:
        raise PreflightException("Must specify a transcriptome reference path.")

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
            raise PreflightException("Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname)

    for filename in cr_constants.STAR_REQUIRED_FILES:
        p = os.path.join(reference_path, cr_constants.REFERENCE_STAR_PATH, filename)
        if not os.path.exists(p):
            raise PreflightException("Your reference doesn't appear to be indexed. Please run the mkreference tool")


def check_chemistry(name, custom_def, allowed_chems):
    check(cr_chem.check_chemistry_defs())
    check(cr_chem.check_chemistry_arg(name, allowed_chems))

    if name == cr_chem.CUSTOM_CHEMISTRY_NAME:
        check(cr_chem.check_chemistry_def(custom_def))

def check_environment():
    check(tk_preflight.check_open_fh())

def record_package_versions():
    for package in cr_constants.PACKAGE_VERSION_CMDS:
        name = package['name']
        cmd = package['cmd']

        version = tk_subproc.check_output(cmd, shell=True)
        print '%s: %s' % (name, version)

def check_read_length(x):
    if x < 1:
        raise PreflightException("Specified read length must be greater than or equal to 1.")
