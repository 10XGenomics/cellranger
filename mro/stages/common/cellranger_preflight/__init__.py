#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import martian
import os
import socket
import tenkit.preflight as tk_preflight
import cellranger.chemistry as cr_chem
import cellranger.preflight as cr_preflight
import cellranger.vdj.preflight as vdj_preflight

__MRO__ = """
stage CELLRANGER_PREFLIGHT_LOCAL(
    in map[]  sample_def,
    in string chemistry_name,
    in map    custom_chemistry_def,
    in path   reference_path,
    in path   vdj_reference_path,
    in bool   check_executables,
    in int    recovered_cells,
    in int    force_cells,
    src py     "stages/common/cellranger_preflight",
)
"""

def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True

def main(args, outs):
    hostname = socket.gethostname()

    print "Checking sample info..."
    ok, msg = tk_preflight.check_gem_groups(args.sample_def)
    if not ok:
        martian.exit(msg)

    print "Checking FASTQ folder..."
    for sample_def in args.sample_def:
        read_path = sample_def["read_path"]
        if not read_path.startswith('/'):
            martian.exit("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            martian.exit("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            martian.exit("On machine: %s, cellranger does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            martian.exit("Specified FASTQ folder is empty: " + read_path)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not is_int(lane):
                    martian.exit("Lanes must be a comma-separated list of numbers.")

        ok, msg = tk_preflight.check_sample_indices(sample_def)
        if not ok:
            martian.exit(msg)

    if args.reference_path is None and args.vdj_reference_path is None:
        martian.exit("Must specify either reference_path or vdj_reference_path.")

    print "Checking transcriptome..."
    if args.reference_path is not None:
        ok, msg = cr_preflight.check_refdata(args.reference_path)
        if not ok:
            martian.exit(msg)

    if args.vdj_reference_path is not None:
        ok, msg = vdj_preflight.check_refdata(args.vdj_reference_path)
        if not ok:
            martian.exit(msg)

    print "Checking chemistry..."
    ok, msg = cr_chem.check_chemistry_defs()
    if not ok:
        martian.exit(msg)

    ok, msg = cr_chem.check_chemistry_arg(args.chemistry_name)
    if not ok:
        martian.exit(msg)

    if args.chemistry_name == cr_chem.CUSTOM_CHEMISTRY_NAME:
        ok, msg = cr_chem.check_chemistry_def(args.custom_chemistry_def)
        if not ok:
            martian.exit(msg)

    # Open file handles limit - per CELLRANGER-824, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking at args.check_executables
    if args.check_executables:
        print "Checking system environment..."
        ok, msg = tk_preflight.check_open_fh()
        if not ok:
            martian.exit(msg)

    print "Checking optional arguments..."
    if args.recovered_cells is not None and args.force_cells is not None:
        martian.exit("Cannot specify both --force-cells and --expect-cells (or --cells) in the same run.")

    cr_preflight.record_package_versions()
