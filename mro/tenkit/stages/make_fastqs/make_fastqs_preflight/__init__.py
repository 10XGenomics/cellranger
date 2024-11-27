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

__MRO__ = """
stage MAKE_FASTQS_PREFLIGHT(
    in path     run_path,
    in path     output_path,
    in path     interop_output_path,
    in bool     check_executables,
    in int      max_bcl2fastq_threads,
    src py      "stages/make_fastqs/make_fastqs_preflight",
)
"""


def main(args, outs):
    hostname = socket.gethostname()

    print("Checking run folder...")
    tk_preflight.check_rta_complete(args.run_path)

    print("Checking RunInfo.xml...")
    tk_preflight.check_runinfo_xml(args.run_path)

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

    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)

    if args.output_path is not None:
        tk_preflight.check_folder_or_create(
            "--output-dir", args.output_path, hostname, permission=os.W_OK | os.X_OK
        )

    if args.interop_output_path is not None:
        tk_preflight.check_folder_or_create(
            "--interop-dir", args.interop_output_path, hostname, permission=os.W_OK | os.X_OK
        )

    if args.max_bcl2fastq_threads < 1:
        msg = "Cannot run bcl2fastq with zero threads."
        martian.exit(msg)
