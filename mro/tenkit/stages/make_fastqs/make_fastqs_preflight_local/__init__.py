#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
import socket

import martian
from six import ensure_str

import tenkit.bcl as tk_bcl
import tenkit.preflight as tk_preflight
import tenkit.samplesheet as tk_sheet

__MRO__ = """
stage MAKE_FASTQS_PREFLIGHT_LOCAL(
    in path     run_path,
    in bool     check_executables,
    in int[]    lanes,
    in map[]    specs,
    in string   bcl2fastq2_args,
    in string   bc_read_type,
    in int      bc_start_index,
    in int      bc_length,
    in string   si_read_type,
    in string   umi_read_type,
    in int      umi_start_index,
    in int      umi_length,
    in string   project_name,
    in string   bases_mask,
    in string   all_mkfastq_args,
    src py      "stages/make_fastqs/make_fastqs_preflight_local",
)
"""


def check_spec(spec):
    # rule: if spec doesn't contain csv, it must contain lane/sample/index, with optional project.
    if not spec.get("csv"):
        # allow None for lanes, downstream will default to all
        if not spec.get("sample") or not spec.get("indices"):
            martian.exit(
                "Samplesheet spec without CSV must include lanes, sample and indices keys."
            )


def check_specs(args):
    hostname = socket.gethostname()
    specs = args.specs

    if not specs:
        martian.exit("Cannot create samplesheet with empty specs.")
    for spec in specs:
        check_spec(spec)

    if len(specs) > 1 and any(spec.get("csv") is not None for spec in specs):
        martian.exit("Cannot combine specs for CSV plus additional entries")

    # check for samplesheet
    csv_specs = [spec for spec in specs if spec.get("csv")]
    if csv_specs:
        csv_spec = csv_specs[0]
        csv_path = csv_spec["csv"]
        tk_preflight.check_file("samplesheet", csv_path, hostname)
        is_iem = tk_sheet.file_is_iem_samplesheet(csv_path)
        is_csv = tk_sheet.file_is_simple_samplesheet(csv_path)
        if not (is_iem or is_csv):
            martian.exit(f"Formatting error in sample sheet: {csv_path}")


def check_read_params(args, runinfo):
    read_info = tk_bcl.load_run_info(runinfo)[0]
    read_info_by_read_type = {r["read_name"]: r for r in read_info}

    # verify barcode
    if args.bc_read_type is None:
        martian.exit("Barcode read must be specified.")
    if args.bc_read_type not in read_info_by_read_type:
        martian.exit(f"Barcode read not found in run folder: {args.bc_read_type}")

    if args.bc_start_index is not None and args.bc_length is not None:
        if (
            args.bc_start_index + args.bc_length
            > read_info_by_read_type[args.bc_read_type]["read_length"]
        ):
            martian.exit(
                "Barcode out of bounds (%s:%d-%d)"
                % (args.bc_read_type, args.bc_start_index, args.bc_start_index + args.bc_length)
            )

    # if sample index reads not generated, must specify lanes to demux
    if args.si_read_type not in read_info_by_read_type:
        if not args.lanes or len(args.lanes) == 0:
            martian.exit("Lanes must be specified if no sample index reads were generated")

    # if UMI present, do bounds check
    if args.umi_read_type is not None and args.umi_read_type not in read_info_by_read_type:
        martian.exit(f"UMI read type not found in run folder: {args.umi_read_type}")
    if args.umi_start_index is not None and args.umi_length is not None:
        if (
            args.umi_start_index + args.umi_length
            > read_info_by_read_type[args.umi_read_type]["read_length"]
        ):
            martian.exit(
                "UMI out of bounds (%s:%d-%d)"
                % (args.umi_read_type, args.umi_start_index, args.umi_start_index + args.umi_length)
            )


def emit_info(args):
    martian.log_info("args:")
    martian.log_info(args.all_mkfastq_args)
    martian.log_info("runParameters.xml parameters:")
    tk_bcl.RTAVersionInformation(args.run_path).log_version_info()
    csv_specs = [spec for spec in args.specs if spec.get("csv")]
    if csv_specs:
        csv_spec = csv_specs[0]
        csv_path = csv_spec["csv"]
        if os.path.exists(csv_path):
            martian.log_info("samplesheet:")
            try:
                with open(csv_path) as sheet:
                    martian.log_info(f"\n{sheet.read()}")
            except PermissionError:
                martian.exit(f"Permission Denied: {csv_path}")

        else:
            martian.log_info(f"samplesheet not found: {csv_path}")
    else:
        martian.log_info("specs:")
        martian.log_info(args.specs)


def main(args, outs):
    hostname = socket.gethostname()

    print("Checking run folder...")
    tk_preflight.check_rta_complete(args.run_path)

    print("Checking RunInfo.xml...")
    runinfo = tk_preflight.check_runinfo_xml(args.run_path)

    print("Checking system environment...")
    ok, msg = tk_preflight.check_ld_library_path()
    if not ok:
        martian.exit(msg)

    if args.check_executables:
        print("Checking bcl2fastq...")
        rta_info = tk_bcl.RTAVersionInformation(args.run_path)
        rta_info.log_version_info()
        (major_ver, full_ver) = rta_info.check_bcl2fastq(hostname)
        martian.log_info(f"Running bcl2fastq mode: {major_ver}.  Version: {ensure_str(full_ver)}")

    if "--no-lane-splitting" in args.bcl2fastq2_args:
        martian.exit("The --no-lane-splitting option is not supported.")

    print("Emitting run information...")
    martian.log_info("-------mkfastq diagnostic start-------")
    emit_info(args)

    print("Checking read specification...")
    check_read_params(args, runinfo)
    martian.log_info("-------mkfastq diagnostic end-------")

    print("Checking samplesheet specs...")
    check_specs(args)

    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)
