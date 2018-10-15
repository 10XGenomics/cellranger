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
import cellranger.library_constants as lib_constants
import cellranger.rna.library
from cellranger.feature_ref import FeatureDefException
import cellranger.rna.feature_ref as rna_feature_ref
import itertools
import csv

ALLOWED_LIBRARY_TYPES = [
    lib_constants.GENE_EXPRESSION_LIBRARY_TYPE,
    cellranger.rna.library.CRISPR_LIBRARY_TYPE,
    cellranger.rna.library.ANTIBODY_LIBRARY_TYPE, 
]

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

def check_sample_def(sample_def, feature_ref=None, pipeline=None):
    hostname = socket.gethostname()

    check(tk_preflight.check_gem_groups(sample_def))

    # Check uniqueness of sample_def entries
    sd_entries = sorted([
        (sd.get("read_path"),
         sd.get("sample_names"),
         sd.get("sample_indices"),
         sd.get("lanes")) for sd in sample_def
    ])

    for i in range(len(sd_entries) - 1):
        if sd_entries[i] == sd_entries[i + 1]:
            msg = "Duplicated entry in the input FASTQ data. Please use a unique combination of fastq path and sample name."
            msg += "\nPath: %s" % sd_entries[i][0]
            msg += "\nNote in demux mode, a unique combination fastq path, sample indices, and lanes is required."
            raise PreflightException(msg)

    print "Checking FASTQ folder..."
    for sample_def in sample_def:
        read_path = sample_def["read_path"]
        if read_path.strip() == "":
            raise PreflightException("Empty fastq path specifed. Please specify an absolute path.")
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

        if pipeline == "count":
            options = ", ".join(("'%s'" % x for x in ALLOWED_LIBRARY_TYPES))
            library_type = sample_def.get("library_type", None)

            # Check for empty library_type
            if library_type == '':
                msg = ("library_type field may not be an empty string."
                    "\nThe 'library_type' field in the libraries csv"
                    " must be one of %s, or start with '%s'") % \
                    (options, cellranger.rna.library.CUSTOM_LIBRARY_TYPE_PREFIX)
                raise PreflightException(msg)

            # Check for a valid library_type
            if not (library_type is None or library_type in ALLOWED_LIBRARY_TYPES or \
            library_type.startswith(cellranger.rna.library.CUSTOM_LIBRARY_TYPE_PREFIX)):

                msg = ("Unknown library_type: '%s'."
                    "\nThe 'library_type' field in the libraries csv"
                    " must be one of %s, or start with '%s'") % \
                    (library_type, options, cellranger.rna.library.CUSTOM_LIBRARY_TYPE_PREFIX)
                raise PreflightException(msg)

            # Check that the library_type exists in the feature_ref
            if feature_ref is not None and \
            library_type is not None and \
            library_type != cr_constants.GENE_EXPRESSION_LIBRARY_TYPE:

                if not any(x.feature_type == library_type for x in feature_ref.feature_defs):
                    msg = "You declared a library with library_type = '%s', but there are no features declared with that feature_type in the feature reference." % library_type
                    msg += "\nCheck that the 'library_type' field in the libraries csv matches at least 1 entry in the 'feature_type' field in the feature reference csv"
                    raise PreflightException(msg)

        elif pipeline == "vdj":
            # library type can be missing, or VDJ
            library_type = sample_def.get("library_type", None)
            if library_type is not None and not (library_type == lib_constants.VDJ_LIBRARY_TYPE):
                msg = "You declared a library with library_type = '%s'. For the vdj pipeline, the library_type field in sample_def must be missing or '%s'" % (library_type, lib_constants.VDJ_LIBRARY_TYPE)
                raise PreflightException(msg)


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

def expand_libraries_csv(csv_path):

    if not os.path.isfile(csv_path):
        raise PreflightException("Could not find the libraries csv file %s" % csv_path)

    if not os.access(csv_path, os.R_OK):
        raise PreflightException("libraries csv is not readable, please check file permissions: %s" % csv_path)

    with open(csv_path, 'rU') as f:

        rows = itertools.ifilter(lambda x: not x.startswith('#'), f)
        reader = csv.DictReader(rows)
        required_cols = set(cr_constants.LIBRARIES_CSV_FIELDS)

        libraries = []
        row_num = 1

        for row in reader:
            print row.keys()

            if not set(row.keys()).issuperset(required_cols):
                raise PreflightException('Libraries csv file must contain the following comma-delimited fields: "%s".' % ', '.join(required_cols))

            for key in row.iterkeys():
                if key is None:
                    msg = "Invalid libraries CSV file: incorrrect number of columns on line %d" % row_num
                    raise PreflightException(msg)
                row[key] = row[key].strip()

            library = {
                "fastq_mode": "ILMN_BCL2FASTQ",
                "gem_group": None,
                "lanes": None,
                "read_path": row["fastqs"],
                "sample_names": [row["sample"]],
                "library_type": row["library_type"],
                "sample_indices": ["any"],
            }

            libraries.append(library)
            row_num += 1

    return libraries


def check_chemistry(name, custom_def, allowed_chems):
    check(cr_chem.check_chemistry_defs())
    check(cr_chem.check_chemistry_arg(name, allowed_chems))

    if name == cr_chem.CUSTOM_CHEMISTRY_NAME:
        check(cr_chem.check_chemistry_def(custom_def))

def check_read_lengths_vs_chemistry(name, allowed_chems, r1_length, r2_length):
    if name == cr_chem.CUSTOM_CHEMISTRY_NAME:
        return
    else:
        chem = cr_chem.get_chemistry(name)

    sets = [("R1", "--r1-length", r1_length), ("R2", "--r2-length", r2_length)]

    for (read_type, flag, user_value) in sets:
        if user_value is not None:
            read_def = None
            if cr_chem.get_rna_read_def(chem).read_type == read_type:
                read_def = cr_chem.get_rna_read_def(chem)
            elif cr_chem.get_rna_read2_def(chem).read_type == read_type:
                read_def = cr_chem.get_rna_read2_def(chem)

            min_rna_len = cr_constants.MINIMUM_TRIMMED_READ_LENGTH_PREFLIGHT
            # Check that alignable sequence after bc/umi removal and hard-trimming is long enough
            if read_def is not None and read_def.offset + min_rna_len > user_value:
                msg = "You selected %s %d. In the selected chemistry is '%s', the RNA read on %s starts at position %d, leaving %d bp of alignable sequence after trimming.  At least %d bp of alignable sequence are required." % \
                (flag, user_value, name, read_type, read_def.offset, max(0, user_value - read_def.offset), cr_constants.MINIMUM_TRIMMED_READ_LENGTH_PREFLIGHT)
                msg += "\n"
                msg += "Please check your %s setting" % flag
                return msg

def check_feature_ref(transcriptome_ref_path, feature_ref_path):

    if not os.path.isfile(feature_ref_path):
        raise PreflightException("Could not find the feature reference file %s" % feature_ref_path)

    if not os.access(feature_ref_path, os.R_OK):
        raise PreflightException("feature reference is not readable, please check file permissions: %s" % feature_ref_path)

    try:
        feature_ref = rna_feature_ref.from_transcriptome_and_csv(
            transcriptome_ref_path,
            feature_ref_path)

        rna_feature_ref.FeatureExtractor(feature_ref)
        return feature_ref

    except FeatureDefException as e:
        raise PreflightException(str(e))

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
