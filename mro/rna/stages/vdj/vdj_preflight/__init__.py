#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import martian

import cellranger.constants as cr_constants
import cellranger.preflight as cr_preflight
import cellranger.vdj.preflight as vdj_preflight

__MRO__ = """
stage VDJ_PREFLIGHT(
    in  map[]  sample_def,
    in  path   vdj_reference_path,
    in  bool   denovo,
    in  bool   full_check,
    in  path   inner_enrichment_primers,
    in  string chain_type,
    src py     "stages/vdj/vdj_preflight",
)
"""


def run_preflight_checks(args):
    print("Checking sample info...")
    cr_preflight.check_sample_def(args.sample_def, pipeline=cr_constants.PIPELINE_VDJ)

    print("Checking reference...")
    vdj_preflight.check_refdata(args.vdj_reference_path, args.denovo)

    print("Checking inner enrichment primers...")
    vdj_preflight.check_inner_enrichment_primers(
        args.inner_enrichment_primers, args.vdj_reference_path
    )

    print("Checking chain...")
    vdj_preflight.check_chain(args.chain_type)

    # Open file handles limit - per CELLRANGER-824, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking testing `not args.full_check`
    if args.full_check:
        print("Checking system environment...")
        cr_preflight.check_environment()


def main(args, outs):
    try:
        run_preflight_checks(args)
    except cr_preflight.PreflightException as e:
        martian.exit(e.msg)

    cr_preflight.record_package_versions()
