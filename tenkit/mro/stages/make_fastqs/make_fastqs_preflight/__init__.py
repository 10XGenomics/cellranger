#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import socket
import martian
import tenkit.preflight as tk_preflight
import tenkit.bcl as tk_bcl
import tenkit.constants as tk_constants
import os

__MRO__ = """
stage MAKE_FASTQS_PREFLIGHT(
    in path   run_path,
    in path   output_path,
    in path   interop_output_path,
    in string barcode_whitelist,
    in bool   check_executables,
    in int    max_bcl2fastq_threads,
    src py    "stages/make_fastqs/make_fastqs_preflight",
)
"""

def main(args, outs):
    hostname = socket.gethostname()

    print "Checking run folder..."
    tk_preflight.check_rta_complete(args.run_path)

    print "Checking RunInfo.xml..."
    tk_preflight.check_runinfo_xml(args.run_path)

    print "Checking system environment..."
    ok, msg = tk_preflight.check_ld_library_path()
    if not ok:
        martian.exit(msg)

    print "Checking barcode whitelist..."
    tk_preflight.check_barcode_whitelist(args.barcode_whitelist)

    if args.check_executables:
        print "Checking bcl2fastq..."
        (rta_version, rc_i2_read, bcl_params) = tk_bcl.get_rta_version(args.run_path)
        martian.log_info("RTA Version: %s" % rta_version)
        martian.log_info("BCL Params: %s" % str(bcl_params))

        (major_ver, full_ver) = tk_bcl.check_bcl2fastq(hostname, rta_version)
        martian.log_info("Running bcl2fastq mode: %s.  Version: %s" % (major_ver, full_ver))

    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)

    if args.output_path is not None:
        tk_preflight.check_folder_or_create("--output-dir", args.output_path, hostname, permission=os.W_OK|os.X_OK)

    if args.interop_output_path is not None:
        tk_preflight.check_folder_or_create("--interop-dir", args.interop_output_path, hostname, permission=os.W_OK|os.X_OK)

    if args.max_bcl2fastq_threads < 1:
        msg = "Cannot run bcl2fastq with zero threads."
        martian.exit(msg)

