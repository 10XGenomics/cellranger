#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# 1. Run bcl2fastq with the appropriate settings to get raw fastq files for
#    all the reads (if neccesary)
# 2. Demultiplex reads using the I1 reads, if present. Initially we will detect
#    common sample indicies by looking at the reads. In the future we will
#    accept a sample sheet
# 3. Put the output FASTQ files in a canonical location
#
import os
import re
import glob
import socket
import tenkit.log_subprocess
import tenkit.bcl as tk_bcl
import tenkit.lane as tk_lane
import shutil
import xml.etree.ElementTree as etree
import martian
from tenkit.constants import DEMULTIPLEX_BARCODE_LENGTH, DEMULTIPLEX_DEFAULT_SAMPLE_INDEX_LENGTH

__MRO__ = """
stage BARCODE_AWARE_BCL2FASTQ(
    in  path   run_path,
    in  int    num_threads,
    in  bool   split_by_tile,
    out path   raw_fastq_path,
    src py     "stages/bcl_processor/barcode_aware_bcl2fastq",
)
split using (
    in  string tile_suffix,
)
"""
def split(args):
    # Illumina bcl2fastq2 guide recommends 32GB RAM
    # empirical measurement on tile-divided NovaSeq S2 run = 8GB or so, 16 to be safe

    if args.split_by_tile:
        chunk_defs = [{'tile_suffix': char, '__threads': 6, '__mem_gb': 16} for char in "0123456789"]
    else:
        chunk_defs = [{'tile_suffix': '*', '__threads': 6, '__mem_gb': 32}]
    return {'chunks': chunk_defs}

def main(args, outs):
    process_raw_ilmn_data(args, outs)

def join(args, outs, chunk_args, chunk_outs):
    os.makedirs(outs.raw_fastq_path)
    for chunk_out in chunk_outs:
        for f in os.listdir(chunk_out.raw_fastq_path):
            in_file = os.path.join(chunk_out.raw_fastq_path, f)
            shutil.move(in_file, outs.raw_fastq_path)

def process_raw_ilmn_data(args, outs):
    """
    run_path must be the top-level Illumina run directory
    """
    input_dir =  os.path.join(args.run_path, "Data", "Intensities", "BaseCalls")
    output_dir = outs.raw_fastq_path

    martian.log_info("Running bcl2fastq on run: %s" % args.run_path)
    martian.log_info("FASTQ output dir: %s" % output_dir)

    if not os.path.exists(args.run_path):
        martian.throw("Run directory does not exist: %s" % args.run_path)

    run_info_xml = os.path.join(args.run_path, "RunInfo.xml")
    read_info, flowcell = tk_bcl.load_run_info(run_info_xml)
    use_bases_mask_val = tk_bcl.make_bases_mask_val(read_info)

    # Determine the RTA version of the run and whether this instrument
    # requires i2 to RC'd
    (rta_version, rc_i2_read, bcl_params) = tk_bcl.get_rta_version(args.run_path)
    martian.log_info("BCL folder RTA Version: %s" % rta_version)
    martian.log_info("BCL params: %s" % str(bcl_params))

    # Determine the best available bcl2fastq version to use
    # Will call martian.exit() with an error message if there isn't
    # a compatible version available
    hostname = socket.gethostname()
    (major_ver, full_ver) = tk_bcl.check_bcl2fastq(hostname, rta_version)

    martian.log_info("Using bcl2fastq version: %s" % full_ver)

    tile_split = args.tile_suffix != '*'

    try:
        # Internal use only. Move aside Illumina sample sheet so
        # bcl2fastq doesn't use it. For customers, there is a pre-flight
        # check to make sure there is no sample sheet in the places
        # bcl2fastq looks for it.
        import kitten

        # Older RTA put sheet into Data/Intensities/BaseCalls while
        # newer RTA put sheet at top of the BCL folder. Check both.
        for ss_dir in [ args.run_path, input_dir ]:
            ilmn_sample_sheet = os.path.join(ss_dir, "SampleSheet.csv")
            mv_sample_sheet = os.path.join(ss_dir, "IlluminaSampleSheet.csv")
            if os.path.exists(ilmn_sample_sheet):
                martian.log_info("Renaming the Illumina sample sheet")
                os.rename(ilmn_sample_sheet, mv_sample_sheet)
    except ImportError:
        pass

    # Restore the LD_LIBRARY_PATH set aside by sourceme.bash/shell10x.
    # Only do this for the environment in which BCL2FASTQ will run.
    new_environ = dict(os.environ)
    new_environ['LD_LIBRARY_PATH'] = os.environ['_TENX_LD_LIBRARY_PATH']

    if major_ver == tk_bcl.BCL2FASTQ_V1:
        if tile_split:
            martian.throw("Cannot support NovaSeq demux scheme on bcl2fastq v1.  Exiting.")

        # configure
        # write bigger fastq chunks to avoid blow-up of chunks
        cmd = [ "configureBclToFastq.pl",  "--fastq-cluster-count", "20000000",
            "--no-eamss", "--use-bases-mask=" + use_bases_mask_val,
            "--input-dir=" + input_dir, "--output-dir=" + output_dir ]

        martian.log_info("Running bcl2fastq setup command:")
        martian.log_info(" ".join(cmd))

        try:
            ret = tenkit.log_subprocess.call(cmd, env=new_environ)
        except OSError:
            martian.throw("configureBclToFastq.pl not found on path -- make sure you've added it to your environment")

        if ret != 0:
            martian.throw("configureBclToFastq.pl failed. Exiting.")

        # Run the actual makefiles
        makefile = os.path.join(output_dir, "Makefile")
        if not os.path.exists(makefile):
            martian.throw("BclToFastq Makefile not found where expected: %s" % makefile)

        martian.log_info("Running Makefile...")
        mk_cmd = [ "make", "-C", output_dir, "-j", str(args.num_threads)]
        martian.log_info(" ".join(mk_cmd))
        ret = tenkit.log_subprocess.call(mk_cmd, env=new_environ)

        if ret > 0:
            martian.throw("running the BclToFastq Makefile failed with code: %d. Exiting" % ret)
        elif ret < 0:
            martian.throw("Bcl2Fastq was killed with signal %d." % ret)

    elif major_ver == tk_bcl.BCL2FASTQ_V2:
        if tile_split:
            proj_output_dir = os.path.join(output_dir, "Tile%s" % args.tile_suffix, "Project_%s" % flowcell)
        else:
            proj_output_dir = os.path.join(output_dir, "Project_%s" % flowcell)

        fastq_output_dir = os.path.join(proj_output_dir, "fastq")
        interop_output_dir = os.path.join(proj_output_dir, "interop")

        if not os.path.exists(fastq_output_dir):
            os.makedirs(fastq_output_dir)

        if not os.path.exists(interop_output_dir):
            os.makedirs(interop_output_dir)

        min_read_length = min([x["read_length"] for x in read_info])

        if tile_split:
            flowcell_info = tk_lane.get_flowcell_layout(run_info_xml)
            if flowcell_info.tile_length is None:
                martian.throw("Cannot determine tile name length from RunInfo.xml")

            tiles_regex_prefix = "[0-9]"*(flowcell_info.tile_length-1)
            tiles_regex = "%s%s" % (tiles_regex_prefix, args.tile_suffix)
            cmd = [ "bcl2fastq" ,
                    "--minimum-trimmed-read-length", str(min_read_length),
                    # PIPELINES-1140 - required in bcl2fastq 2.17 to generate correct index read fastqs
                    "--mask-short-adapter-reads", str(min_read_length),
                    # LONGRANGER-121 - ignore missing bcl data
                    "--ignore-missing-bcls", "--ignore-missing-filter", "--ignore-missing-positions", "--ignore-missing-controls",
                    '-r', str(args.__threads), '-w', str(args.__threads),
                    # TENKIT-72 avoid CPU oversubscription
                    '-p', str(args.__threads),
                    "--use-bases-mask=" + use_bases_mask_val, "-R", args.run_path,
                    "--output-dir=" + fastq_output_dir,
                    "--interop-dir=" + interop_output_dir,
                    "--tiles=" + tiles_regex]
        else:
            cmd = ["bcl2fastq",
                   "--minimum-trimmed-read-length", str(min_read_length),
                   # PIPELINES-1140 - required in bcl2fastq 2.17 to generate correct index read fastqs
                   "--mask-short-adapter-reads", str(min_read_length),
                   # LONGRANGER-121 - ignore missing bcl data
                   "--ignore-missing-bcls", "--ignore-missing-filter", "--ignore-missing-positions",
                   "--ignore-missing-controls",
                   '-r', str(args.__threads), '-w', str(args.__threads),
                   # TENKIT-72 avoid CPU oversubscription
                   '-p', str(args.__threads),
                   "--use-bases-mask=" + use_bases_mask_val, "-R", args.run_path,
                   "--output-dir=" + fastq_output_dir,
                   "--interop-dir=" + interop_output_dir]

        martian.log_info("Running bcl2fastq 2: %s" % (" ".join(cmd)))

        try:
            ret = tenkit.log_subprocess.call(cmd, env=new_environ)
        except OSError:
            martian.throw("bcl2fastq not found on PATH -- make sure you've added it to your environment")

        if ret > 0:
            martian.exit("bcl2fastq failed. Exiting.")
        elif ret < 0:
            martian.exit("bcl2fastq was killed with signal %d." % ret)

   # Glob over all lanes - demultiplex handles whether to collapse them
    if tile_split:
        fastq_glob = os.path.join(output_dir, "Tile*", "Project_" + flowcell, "*", "*.fastq*")
    else:
        fastq_glob = os.path.join(output_dir, "Project_" + flowcell, "*", "*.fastq*")
    start_fastq_files = glob.glob(fastq_glob)

    # File renaming -- bcl2fastq names the reads R1, R2, R3, R4
    # Use our conventions to make them R1, I1, I2, R2, as the case may be.
    rename_fastq_files(read_info, start_fastq_files)


def rename_fastq_files(read_info, file_list):
    tmp_dirs = {}
    tmp_names = {}

    # Figure out new name for file, move it to temp directory with that name
    for f in file_list:
        path = os.path.dirname(f)
        basename = os.path.basename(f)

        if not (path in tmp_dirs.keys()):
            tmp_dir = os.path.join(path, "mv_tmp")
            os.mkdir(tmp_dir)
            tmp_dirs[path] = tmp_dir

        tmp_dir = tmp_dirs[path]

        for read_type in read_info:
            old_key = "_" + read_type["original_read_name"] + "_"
            new_key =  "_" + read_type["read_name"] + "_"

            if old_key in basename:
                new_name = re.sub(old_key, new_key, basename)
                break

        tmp_f = os.path.join(tmp_dir, new_name)
        new_f = os.path.join(path, new_name)
        tmp_names[tmp_f] = new_f

        martian.log_info("renaming: %s to %s" % (f, tmp_f))

        os.rename(f, tmp_f)

    # Now move all the temp file to their final location
    for (tmp_f, new_f) in tmp_names.items():
        martian.log_info("renaming: %s to %s" % (tmp_f, new_f))

        os.rename(tmp_f, new_f)

    # Remove the tmp directories
    for tmp_dir in tmp_dirs.values():
        os.rmdir(tmp_dir)
