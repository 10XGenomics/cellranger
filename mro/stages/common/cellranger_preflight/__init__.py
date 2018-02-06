#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import martian
import cellranger.preflight as cr_preflight

__MRO__ = """
stage CELLRANGER_PREFLIGHT_LOCAL(
    in map[]  sample_def,
    in string chemistry_name,
    in map    custom_chemistry_def,
    in path   reference_path,
    in bool   check_executables,
    in int    recovered_cells,
    in int    force_cells,
    src py     "stages/common/cellranger_preflight",
)
"""

def run_preflight_checks(args):
    print "Checking sample info..."
    cr_preflight.check_sample_def(args.sample_def)

    print "Checking reference..."
    cr_preflight.check_refdata(args.reference_path)

    print "Checking chemistry..."
    cr_preflight.check_chemistry(args.chemistry_name, args.custom_chemistry_def,
                                 args.allowed_chems)

    if args.r1_length is not None:
        print "Checking read 1 length..."
        cr_preflight.check_read_length(args.r1_length)
    if args.r2_length is not None:
        print "Checking read 2 length..."
        cr_preflight.check_read_length(args.r2_length)

    # Open file handles limit - per CELLRANGER-824, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking at args.check_executables
    if args.check_executables:
        print "Checking system environment..."
        cr_preflight.check_environment()

    print "Checking optional arguments..."
    if args.recovered_cells is not None and args.force_cells is not None:
        raise cr_preflight.PreflightException("Cannot specify both --force-cells and --expect-cells in the same run.")

def main(args, outs):
    try:
        run_preflight_checks(args)
    except cr_preflight.PreflightException as e:
        martian.exit(e.msg)

    cr_preflight.record_package_versions()
