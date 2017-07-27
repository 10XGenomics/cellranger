#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import socket
import subprocess
import xml.etree.ElementTree
import martian
import tenkit.constants as tk_constants
import tenkit.preflight as tk_preflight
import tenkit.bcl

__MRO__ = """
stage BCL_PROCESSOR_PREFLIGHT(
    in path run_path,
    in bool allow_no_barcodes,
    in bool check_executables,
    src py   "stages/preflight/bcl_processor",
)
"""

def main(args, outs):
    hostname = socket.gethostname()

    print "Checking run folder..."
    tk_preflight.check_rta_complete(args.run_path)

    print "Checking RunInfo.xml..."
    runinfo = tk_preflight.check_runinfo_xml(args.run_path)

    if not args.allow_no_barcodes:
        ok, msg = check_reads(runinfo)
        if not ok:
            martian.exit(msg)

    print "Checking system environment..."
    ok, msg = tk_preflight.check_ld_library_path()
    if not ok:
        martian.exit(msg)

    # Presence of SampleSheet.csv interferes with demux.
    # Ask customer to move it. Under older RTA, bcl2fastq looks for it
    # in Data/Intensities/BaseCalls while under newer RTA, it looks for it
    # at the top of the run folder.
    bc_dir = os.path.join(args.run_path, "Data", "Intensities", "BaseCalls")
    for ss_dir in [ args.run_path, bc_dir ]:
        ilmn_sample_sheet = os.path.join(ss_dir, "SampleSheet.csv")

        external = True
        try:
            import kitten
            external = False
        except ImportError:
            pass

        if external and os.path.exists(ilmn_sample_sheet):
            martian.exit("On machine: %s, SampleSheet.csv found in run folder that would interfere with demux:\n%s\nPlease move, rename, or delete the file and run demux again." % (hostname, ilmn_sample_sheet))

    if args.check_executables:
        print "Checking bcl2fastq..."
        # Determine the RTA version of the run and whether this instrument
        # requires i2 to RC'd
        (rta_version, rc_i2_read, bcl_params) = tenkit.bcl.get_rta_version(args.run_path)
        martian.log_info("RTA Version: %s" % rta_version)
        martian.log_info("BCL Params: %s" % str(bcl_params))


        # Determine the best available bcl2fastq version to use
        # Will call martian.exit() with an error message if there isn't
        # a compatible version available
        (major_ver, full_ver) = tenkit.bcl.check_bcl2fastq(hostname, rta_version)
        martian.log_info("Running bcl2fastq mode: %s.  Version: %s" % (major_ver, full_ver))

    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)

def check_bcl2fastq_v1(hostname):
    try:
        subprocess.check_call(["which", "configureBclToFastq.pl"])
    except subprocess.CalledProcessError:
        martian.exit("On machine: %s, bcl2fastq or configureBclToFastq.pl not found on PATH." % hostname)
    print "configureBclToFastq.pl version on %s: < 2.0" % hostname
    try:
        subprocess.check_call(["which", "perl"])
    except subprocess.CalledProcessError:
        martian.exit("On machine: %s, perl not found on PATH." % hostname)

def check_bcl2fastq_v2(hostname):
    try:
        output = subprocess.check_output(["bcl2fastq", "--version"])
    except subprocess.CalledProcessError:
        martian.exit("On machine: %s, bcl2fastq not found on PATH." % hostname)
    print "bcl2fastq version on %s: %s" % (hostname, output)

def check_bcl2fastq(hostname):
    try:
        subprocess.check_call(["which", "bcl2fastq"])
        check_bcl2fastq_v2(hostname)
    except subprocess.CalledProcessError:
        check_bcl2fastq_v1(hostname)

# TODO: this is different for Whole Genome than for Single Cell;
# figure out a more generic/permissive algorithm for preflight
# check (or parameterize the number of reads)
def check_reads(runinfo):
    tree = xml.etree.ElementTree.parse(runinfo)

    reads = tree.getroot().findall('./Run/Reads/Read')

    if len(reads) < tk_constants.REQUIRED_MIN_READS:
        return False, "Run does not match Chromium or GemCode read structure (2 template reads + 1 index read or 2 template reads + 2 index reads)."

    # Search for the GemCode barcode read if there are 4 reads
    if len(reads) == 4:
        for read in reads:
            if read.get('IsIndexedRead') == 'Y' and read.get('NumCycles') == str(tk_constants.DEMULTIPLEX_GEMCODE_BARCODE_LENGTH):
                return True, None

            # Allow a 1-base dummy read
            if read.get('IsIndexedRead') == 'Y' and read.get('NumCycles') == str(1):
                return True, None

        return False, "GemCode barcode read not found."

    return True, None
