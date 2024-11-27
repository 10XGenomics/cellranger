#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
import shlex
import shutil
import socket

import martian
from six import ensure_str

import tenkit.bcl as tk_bcl
import tenkit.log_subprocess as tk_proc

__MRO__ = """
stage BCL2FASTQ_WITH_SAMPLESHEET(
    in  path     run_path,
    in  path     output_path,
    in  path     interop_output_path,
    in  csv      samplesheet_path,
    in  string   bases_mask,
    in  string   bcl2fastq1_args,
    in  string   bcl2fastq2_args,
    in  string   bc_read_type,
    in  string   si_read_type,
    in  int      max_bcl2fastq_threads,
    in  bool     dual_indexed_flowcell,
    in  bool     filter_single_index,
    out path     fastq_path,
    out path     interop_path,
    out bool     rc_i2_read,
    out map      file_read_types_map,
    out string   bcl2fastq_version,
    out string   bcl2fastq_args,
    src py       "stages/make_fastqs/bcl2fastq_with_samplesheet"
)
split using (
)
"""


def split(args):
    # Illumina bcl2fastq2 guide recommends 32GB RAM
    return {"chunks": [{"__threads": args.max_bcl2fastq_threads, "__mem_gb": 32}]}


def main(args, outs):
    run_bcl2fastq(args, outs)


def remove_deprecated_args(argstr, major_ver: tk_bcl.Bcl2FastqVersion, full_ver: bytes):
    """Given the argstr to pass into bcl2fastq, parse out.

    any deprecated arguments and return the tokenized list.
    """
    arglist = shlex.split(argstr.strip())
    if major_ver == tk_bcl.Bcl2FastqVersion.V2:
        # TENKIT-58 -d/--demultiplex-threads no longer supported
        # bcl2fastq 2.19 says that --ignore-missing-controls is also deprecated
        if full_ver >= b"2.19":
            bad_param_args = ("-d", "--demultiplexing-threads")
            for arg in bad_param_args:
                if arg in arglist:
                    arg_idx = arglist.index(arg)
                    # get rid of arg and param
                    arglist = arglist[:arg_idx] + arglist[arg_idx + 2 :]
            bad_single_args = ("--ignore-missing-controls",)
            for arg in bad_single_args:
                if arg in arglist:
                    arg_idx = arglist.index(arg)
                    # get rid of arg
                    arglist = arglist[:arg_idx] + arglist[arg_idx + 1 :]
    return arglist


def join(args, outs, chunk_args, chunk_outs):
    chunk = chunk_outs[0]

    # if --output-dir and --interop-dir not overridden, move
    # into standard output location
    if not args.output_path:
        shutil.move(chunk.fastq_path, outs.fastq_path)
    else:
        outs.fastq_path = args.output_path

    if not args.interop_output_path:
        shutil.move(chunk.interop_path, outs.interop_path)
    else:
        outs.interop_path = args.interop_output_path

    outs.rc_i2_read = chunk.rc_i2_read
    outs.file_read_types_map = chunk.file_read_types_map
    outs.bcl2fastq_version = chunk.bcl2fastq_version
    outs.bcl2fastq_args = chunk.bcl2fastq_args


def run_bcl2fastq(args, outs):
    if args.output_path:
        outs.fastq_path = args.output_path

    output_dir = outs.fastq_path

    if args.interop_output_path:
        outs.interop_path = args.interop_output_path

    interop_dir = outs.interop_path

    martian.log_info(f"Running bcl2fastq on run: {args.run_path}")
    martian.log_info(f"FASTQ output dir: {output_dir}")

    run_info_xml = os.path.join(args.run_path, "RunInfo.xml")
    read_info = tk_bcl.load_run_info(run_info_xml)[0]
    if not args.bases_mask:
        use_bases_mask_val = tk_bcl.make_bases_mask_val(
            read_info,
            barcode_read=args.bc_read_type,
            sample_index_read=args.si_read_type,
            dual_indexed=args.dual_indexed_flowcell,
            ignore_dual_index=args.filter_single_index,
        )
    else:
        use_bases_mask_val = args.bases_mask

    outs.file_read_types_map = tk_bcl.get_bcl2fastq_read_type_map(
        read_info,
        barcode_read=args.bc_read_type,
        sample_index_read=args.si_read_type,
        dual_indexed=args.dual_indexed_flowcell,
        ignore_dual_index=args.filter_single_index,
    )

    # Determine the RTA version of the run and whether this instrument
    # requires i2 to be RC'd
    rta_info = tk_bcl.RTAVersionInformation(args.run_path)
    rta_info.log_version_info()
    outs.rc_i2_read = rta_info.rc_i2

    # Determine the best available bcl2fastq version to use
    # Will call martian.exit() with an error message if there isn't
    # a compatible version available
    hostname = socket.gethostname()
    (major_ver, full_ver) = rta_info.check_bcl2fastq(hostname)
    outs.bcl2fastq_version = full_ver

    martian.log_info(f"Using bcl2fastq version: {ensure_str(full_ver)}")
    martian.log_info(f"RC'ing i2 read: {rta_info.rc_i2!s}")

    # Restore the LD_LIBRARY_PATH set aside by sourceme.bash/shell10x.
    # Only do this for the environment in which BCL2FASTQ will run.
    new_environ = dict(os.environ)
    new_environ["LD_LIBRARY_PATH"] = os.environ["_TENX_LD_LIBRARY_PATH"]

    if major_ver == tk_bcl.Bcl2FastqVersion.V1:
        martian.exit(
            "bcl2fastq 1.8.4 is not currently supported. Please install bcl2fastq 2.17 or higher."
        )
        raise SystemExit()

    elif major_ver == tk_bcl.Bcl2FastqVersion.V2:
        if not os.path.exists(outs.interop_path):
            os.makedirs(outs.interop_path)
        if not os.path.exists(outs.fastq_path):
            os.makedirs(outs.fastq_path)

        # minimum-trimmed-read-length and mask-short-adapter-reads must be our call (SIs, UMIs)
        min_read_length = min(x["read_length"] for x in read_info)
        if min_read_length > 8:
            # ensure min is at sample-index, if extra base grabbed for QC purposes (I8n, for example)
            min_read_length = 8

        cmd = [
            "bcl2fastq",
            "--minimum-trimmed-read-length",
            str(min_read_length),
            "--mask-short-adapter-reads",
            str(min_read_length),
            "--create-fastq-for-index-reads",
            "--ignore-missing-positions",
            "--ignore-missing-filter",
            "--ignore-missing-bcls",
            #'-r', str(args.__threads), '-w', str(args.__threads),
            "--use-bases-mask=" + use_bases_mask_val,
            "-R",
            args.run_path,
            "--output-dir=" + output_dir,
            "--interop-dir=" + interop_dir,
            "--sample-sheet=" + args.samplesheet_path,
        ]
        cmd += remove_deprecated_args(args.bcl2fastq2_args, major_ver, full_ver)
        outs.bcl2fastq_args = " ".join(cmd)

        martian.log_info(f"Running bcl2fastq2: {outs.bcl2fastq_args}")

        try:
            ret = tk_proc.call(cmd, env=new_environ)
        except OSError:
            martian.throw(
                "bcl2fastq not found on PATH -- make sure you've added it to your environment"
            )

        if ret > 0:
            files_path = os.path.abspath(martian.make_path("_stderr"))
            enclosing_path = os.path.dirname(os.path.dirname(files_path))
            stderr_path = os.path.join(enclosing_path, b"_stderr")
            martian.exit(
                f"bcl2fastq exited with an error. You may have specified an invalid command-line option. See the full error here:\n{ensure_str(stderr_path)}"
            )
        elif ret < 0:
            # subprocess.call returns negative code (on UNIX): bcl2fastq killed by external signal
            martian.exit("bcl2fastq was killed with signal %d." % ret)
