#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import errno
import glob
import os
from collections import defaultdict

import martian
import pandas as pd
from six import ensure_binary, ensure_str

import tenkit.fasta as tk_fasta
import tenkit.log_subprocess as log_subprocess
import tenkit.samplesheet as tk_sheet

__MRO__ = """
stage MERGE_FASTQS_BY_LANE_SAMPLE(
    in  path     fastq_path,
    in  csv      samplesheet_path,
    in  bool     remove_undetermined_fastqs,
    in  bool     remove_split_fastqs,
    in  string   bcl2fastq_version,
    out bool     files_merged,
    out string[] merged_file_paths,
    src py       "stages/make_fastqs/merge_fastqs_by_lane_sample",
) split using (
    in string   project,
    in int      lane,
    in string   sample,
    in string[] input_files,
)"""


def split(args):
    chunk_dict = defaultdict(list)

    # find projects inside run dir
    dirnames = next(os.walk(args.fastq_path))[1]

    # get sample ID to original sample ID mapping
    iem_df = tk_sheet.file_get_iem_data_frame(args.samplesheet_path)
    rows = pd.DataFrame(iem_df, columns=["Sample_ID", "Original_Sample_ID", "Lane"]).to_numpy()

    # construct sample id structure (Sn)
    sample_count = 1
    sample_target_count = 1
    sample_internal_id_dict = {}
    sample_target_id_dict = {}
    internal_id_attr_dict = {}

    # An "S-value" here is defined as the unique numerical index assigned to a sample ID
    # by bcl2fastq when iterating over a sample sheet.  The first time a new Sample_ID is
    # encountered (in lane order-- see CSI-293), it is given a monotonically increasing value.
    # This S-value will be included in any FASTQs associated with that sample_ID, in the '_S#_'
    # part of the file name.
    #
    # Consider two monotonically increasing sequences, M and N, with indices i and j.
    # When 10x sample index rows are split into four, the S-values assigned to those
    # split sample IDs (if seen for the first time in a samplesheet) will be M[i]...M[i+3].
    # However, these all refer to the same sample, which should be given a single
    # S-value of N[j].
    #
    # Create that mapping here of i->j, which stores the desired S-value of the merged
    # output for each file with the S-value i.  For convenience, add to the map value
    # the original sample ID, to yield our final dict internal_id_attr_dict i->(j, Sample_ID[j])
    #
    # This will be used to determine the correct S-value and Sample_ID directory of any
    # output file.

    lane_row_dict = defaultdict(list)
    # determine lane ordering first for bcl2fastq 2.17 or lower;
    # bcl2fastq 2.18/2.19 use samplesheet appearance first
    for idx, row in enumerate(rows):
        # bcl2fastq 2.17 and earlier: sample number by lane appearance
        # or miseq: where there is no lane (len(row) <= 2)
        if args.bcl2fastq_version < "2.18" and len(row) > 2:
            _, _, lane = row
            lane_row_dict[lane].append(idx)
        else:
            # bcl2fastq 2.18 and later: sample number by samplesheet order
            # put all samples in the same "lane" in order to number
            # them by appearance in the sample sheet
            #
            # or if miseq, there is no lane
            lane_row_dict["all"].append(idx)

    bytes_intern = {}
    for lane, row_idxs in sorted(lane_row_dict.items()):
        for row_idx in row_idxs:
            if len(rows[row_idx]) > 2:
                demux_sample_id, original_sample_id, _ = rows[row_idx]
            else:
                demux_sample_id, original_sample_id = rows[row_idx]
            if original_sample_id not in sample_target_id_dict:
                sample_target_count_str = str(sample_target_count)
                sample_target_count_bytes = bytes_intern.get(sample_target_count_str, None)
                if sample_target_count_bytes is None:
                    sample_target_count_bytes = ensure_binary(sample_target_count_str)
                    bytes_intern[sample_target_count_str] = sample_target_count_bytes
                sample_target_id_dict[original_sample_id] = sample_target_count_bytes
                sample_target_count += 1

            if demux_sample_id not in sample_internal_id_dict:
                id_str = str(sample_count)
                id_bytes = bytes_intern.get(id_str, None)
                if id_bytes is None:
                    id_bytes = ensure_binary(id_str)
                    bytes_intern[id_str] = id_bytes
                sample_internal_id_dict[demux_sample_id] = id_bytes
                internal_id_attr_dict[id_bytes] = (
                    sample_target_id_dict[original_sample_id],
                    original_sample_id,
                )
                sample_count += 1

    # this should be fine even if --reports-dir and --stats-dir are elsewhere
    # theoretically, you could name your sample 'Reports' or 'Stats' and
    # move the --reports-dir and --stats-dir somewhere else but at that
    # point you're just being daft
    project_dirs = [dn for dn in dirnames if dn not in ("Reports", "Stats")]

    # this can be the case if no project is specified
    if not project_dirs:
        # add root folder
        project_dirs = ["."]

    for project in sorted(project_dirs):
        project_path = os.path.join(args.fastq_path, project)
        # split by detected sample (all types)
        r1_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(project_path, "R1", None, None)
        r2_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(project_path, "R2", None, None)
        r3_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(project_path, "R3", None, None)
        r4_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(project_path, "R4", None, None)
        i1_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(project_path, "I1", None, None)
        i2_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(project_path, "I2", None, None)

        # group files together by like sample
        all_files = r1_files + r2_files + r3_files + r4_files + i1_files + i2_files

        for path in sorted(all_files):
            file_spec = tk_fasta.IlmnFastqFile(path)

            assert isinstance(file_spec.prefix, bytes)
            # do not consider chunk if Undetermined
            if file_spec.prefix == b"Undetermined":
                continue

            assert isinstance(file_spec.s, bytes)
            # if file merged already, do not consider
            if file_spec.s.startswith(b"0"):
                continue

            # chunk by target output sample identifier for downstream processing
            target_internal_id, original_sample_id = internal_id_attr_dict[file_spec.s]

            chunk_dict[(project, file_spec.lane, original_sample_id, target_internal_id)].append(
                path
            )

    if len(chunk_dict) == 0:
        martian.exit(
            """
No FASTQs matched your sample sheet's lane, sample, and indexes.
This error is usually caused by one of the following:
1. You are demultiplexing a NovaSeq flowcell, but you are using an older
   version of bcl2fastq. Please upgrade to bcl2fastq 2.20+.
2. There may be a mistake in your samplesheet, please double check that
   you are specifying the right lanes, samples, and index sets.
3. Your sequencing data might be of such low quality that all the reads
   ended up in the "Undetermined" FASTQs. Please double check the quality
   of your sequence data.
"""
        )

    chunk_defs = [
        {
            "input_files": file_list,
            "project": tup[0],
            "lane": tup[1],
            "sample_id": tup[2],
            "output_snum": tup[3],
        }
        for tup, file_list in sorted(chunk_dict.items())
    ]

    # MARTIAN-396 FIXME: dummy chunk if no chunk_defs
    if not chunk_defs:
        dummy_def = {
            "input_files": None,
            "project": None,
            "lane": None,
            "sample_id": None,
            "output_snum": None,
        }
        chunk_defs.append(dummy_def)

    return {"chunks": chunk_defs}


def main(args, outs):
    if not args.sample_id:
        # dummy chunk, no merge necessary
        outs.files_merged = False
        outs.merged_file_paths = []
        return

    if not args.input_files:
        # no files to merge
        outs.files_merged = False
        outs.merged_file_paths = []
        return

    # group files by read
    file_infos = [tk_fasta.IlmnFastqFile(path) for path in args.input_files]
    paths_by_type = defaultdict(list)  # type: defaultdict[bytes, list[str]]
    for info in file_infos:
        paths_by_type[info.read].append(info.filename)

    # determine if files are spread out in different subfolders (e.g., project/sample_id)
    subfolders = {os.path.dirname(path) for path in args.input_files}
    if len(subfolders) > 1:
        # create unified output folder
        output_folder = os.path.join(
            os.path.dirname(os.path.dirname(args.input_files[0])), args.sample_id
        )

        # multiple threads creating same directory at once has defeated
        # the os.path.exists/makedirs check on my end; so just try/catch
        # for the path exists error
        try:
            os.makedirs(output_folder)
        except OSError as ex:
            errnum = ex.args[0]
            if errnum != errno.EEXIST:
                raise ex
    else:
        # just use the single folder all files live in
        output_folder = subfolders.pop()

    out_files = []
    for paths in paths_by_type.values():
        # if there is only one file, just do a move
        file_info = tk_fasta.IlmnFastqFile(paths[0])
        if len(paths) == 1:
            out_name = "%s_S0%s_L%03d_%s_%03d.fastq.gz" % (
                ensure_str(file_info.prefix),
                args.output_snum,
                file_info.lane,
                ensure_str(file_info.read),
                file_info.group,
            )
            out_path = os.path.join(output_folder, os.path.basename(out_name))
            # TODO: Don't be evil
            os.rename(paths[0], out_path)

        else:
            out_file = "%s_S0%s_L%03d_%s_%03d.fastq.gz" % (
                ensure_str(file_info.prefix),
                args.output_snum,
                file_info.lane,
                ensure_str(file_info.read),
                file_info.group,
            )
            out_path = os.path.join(output_folder, os.path.basename(out_file))
            with open(out_path, "wb") as out_file:
                log_subprocess.check_call(["cat"] + paths, stdout=out_file)
        out_files.append(out_path)

    # need something for non-blank chunk_outs (Martian foible)
    outs.files_merged = True
    outs.merged_file_paths = out_files


def join(args, outs, chunk_args, chunk_outs):
    if args.remove_split_fastqs:
        unique_folders = set()
        for chunk_arg, chunk_out in zip(chunk_args, chunk_outs):
            # if no files present or no files need merging, no folders necessary to remove
            if not chunk_out.files_merged:
                continue
            for input_file in chunk_arg.input_files:
                if os.path.exists(input_file):
                    martian.log_info(f"Removing split file: {input_file}")
                    os.remove(input_file)
            unique_folders.update(
                [os.path.dirname(input_file) for input_file in chunk_arg.input_files]
            )

        # remove empty split folders
        for folder in unique_folders:
            if not glob.glob(f"{folder}/*fastq*"):
                martian.log_info(f"Removing empty folder: {folder}")
                os.rmdir(folder)

    # rename merged files (which have a _S0n FASTQ name for distinction in flat directories)
    # to final name, with _Sn
    for chunk_out in chunk_outs:
        merged_files = chunk_out.merged_file_paths
        for merged_file in merged_files:
            file_info = tk_fasta.IlmnFastqFile(merged_file)
            out_name = "%s_S%s_L%03d_%s_%03d.fastq.gz" % (
                ensure_str(file_info.prefix),
                ensure_str(file_info.s[1:]),
                file_info.lane,
                ensure_str(file_info.read),
                file_info.group,
            )
            out_path = os.path.join(os.path.dirname(merged_file), os.path.basename(out_name))
            # TODO: Don't be evil
            os.rename(merged_file, out_path)

    if args.remove_undetermined_fastqs:
        undetermined_fastqs = glob.glob(f"{args.fastq_path}/Undetermined*fastq.gz")
        for undetermined_fastq in undetermined_fastqs:
            martian.log_info(f"Removing undetermined file: {undetermined_fastq}")
            os.remove(undetermined_fastq)
