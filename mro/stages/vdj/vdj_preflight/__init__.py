#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import martian
import cellranger.preflight as cr_preflight
import cellranger.vdj.preflight as vdj_preflight

__MRO__ = """
stage VDJ_PREFLIGHT(
    in  map[]  sample_def,
    in  string chemistry_name,
    in  map    custom_chemistry_def,
    in  path   vdj_reference_path,
    in  bool   denovo,
    in  bool   check_executables,
    in  int    force_cells,
    in  string chain_type,
    src py     "stages/vdj/vdj_preflight",
)
"""

def run_preflight_checks(args):
    print "Checking sample info..."
    cr_preflight.check_sample_def(args.sample_def)

    print "Checking reference..."
    vdj_preflight.check_refdata(args.vdj_reference_path, args.denovo)

    print "Checking chain..."
    vdj_preflight.check_chain(args.chain_type)

    print "Checking chemistry..."
    cr_preflight.check_chemistry(args.chemistry_name, args.custom_chemistry_def, None)

    # Open file handles limit - per CELLRANGER-824, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking at args.check_executables
    if args.check_executables:
        print "Checking system environment..."
        cr_preflight.check_environment()

def main(args, outs):
    try:
        run_preflight_checks(args)
    except cr_preflight.PreflightException as e:
        martian.exit(e.msg)

    cr_preflight.record_package_versions()
