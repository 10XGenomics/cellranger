#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import os
import glob
import martian
import errno
from tenkit import fasta as tk_fasta
from tenkit import samplesheet as tk_sheet
from tenkit import log_subprocess

__MRO__="""
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
    in string project,
    in int    lane,
    in string sample,
    in string input_files,
)"""

def split(args):
    chunk_dict = defaultdict(list)

    # find projects inside run dir
    dirnames = next(os.walk(args.fastq_path))[1]

    # get sample ID to original sample ID mapping
    iem_df = tk_sheet.file_get_iem_data_frame(args.samplesheet_path)
    rows = iem_df.as_matrix(['Sample_ID', 'Original_Sample_ID', 'Lane'])

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
    
    for lane, row_idxs in sorted(lane_row_dict.items()):
        for row_idx in row_idxs:
            if len(rows[row_idx]) > 2:
                demux_sample_id, original_sample_id, _ = rows[row_idx]
            else:
                demux_sample_id, original_sample_id = rows[row_idx]
            if original_sample_id not in sample_target_id_dict:
                sample_target_id_dict[original_sample_id] = str(sample_target_count)
                sample_target_count += 1

            if demux_sample_id not in sample_internal_id_dict:
                sample_internal_id_dict[demux_sample_id] = str(sample_count)
                internal_id_attr_dict[str(sample_count)] = (sample_target_id_dict[original_sample_id], original_sample_id)
                sample_count += 1

    # this should be fine even if --reports-dir and --stats-dir are elsewhere
    # theoretically, you could name your sample 'Reports' or 'Stats' and
    # move the --reports-dir and --stats-dir somewhere else but at that
    # point you're just being daft
    project_dirs = [dn for dn in dirnames if dn not in ('Reports', 'Stats')]

    # this can be the case if no project is specified
    if not project_dirs:
        # add root folder
        project_dirs = ['.']

    for project in sorted(project_dirs):
        project_path = os.path.join(args.fastq_path, project)
        # split by detected sample (all types)
        r1_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R1', None, None)
        r2_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R2', None, None)
        r3_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R3', None, None)
        r4_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R4', None, None)
        i1_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'I1', None, None)
        i2_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'I2', None, None)

        # group files together by like sample
        all_files = r1_files + r2_files + r3_files + r4_files + i1_files + i2_files

        for path in sorted(all_files):
            file_spec = tk_fasta.IlmnFastqFile(path)

            # do not consider chunk if Undetermined
            if file_spec.prefix == 'Undetermined':
                continue

            # if file merged already, do not consider
            if file_spec.s[0] == '0':
                continue

            # chunk by target output sample identifier for downstream processing
            target_internal_id, original_sample_id = internal_id_attr_dict[file_spec.s]

            chunk_dict[(project, file_spec.lane, original_sample_id, target_internal_id)].append(path)

    chunk_defs = [{'input_files': file_list,
                   'project': tup[0],
                   'lane': tup[1],
                   'sample_id': tup[2],
                   'output_snum': tup[3]} for tup, file_list in sorted(chunk_dict.items())]

    # MARTIAN-396 FIXME: dummy chunk if no chunk_defs
    if not chunk_defs:
        dummy_def = {'input_files': None, 'project': None, 'lane': None, 'sample_id': None, 'output_snum': None}
        chunk_defs.append(dummy_def)

    return { 'chunks': chunk_defs }


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
    paths_by_type = defaultdict(list)
    for info in file_infos:
        paths_by_type[info.read].append(info.filename)

    # determine if files are spread out in different subfolders (e.g., project/sample_id)
    subfolders = set([os.path.dirname(path) for path in args.input_files])
    if len(subfolders) > 1:
        # create unified output folder
        output_folder = os.path.join(os.path.dirname(os.path.dirname(args.input_files[0])),args.sample_id)

        # multiple threads creating same directory at once has defeated
        # the os.path.exists/makedirs check on my end; so just try/catch
        # for the path exists error
        try:
            os.makedirs(output_folder)
        except OSError, e:
            errnum, message = e.args
            if errnum != errno.EEXIST:
                raise e
    else:
        # just use the single folder all files live in
        output_folder = subfolders.pop()

    out_files = []
    for read_type, paths in paths_by_type.items():
        # if there is only one file, just do a move
        file_info = tk_fasta.IlmnFastqFile(paths[0])
        if len(paths) == 1:
            out_name = "%s_S0%s_L%s_%s_%s.fastq.gz" % (
                file_info.prefix,
                args.output_snum,
                str(file_info.lane).zfill(3),
                file_info.read,
                str(file_info.group).zfill(3)
            )
            out_path = os.path.join(output_folder, os.path.basename(out_name))
            os.rename(paths[0], out_path)

        else:
            out_file = "%s_S0%s_L%s_%s_%s.fastq.gz" % (
                file_info.prefix,
                args.output_snum,
                str(file_info.lane).zfill(3),
                file_info.read,
                str(file_info.group).zfill(3))
            out_path = os.path.join(output_folder, os.path.basename(out_file))
            subprocess_args = ["cat"] + paths + [">", out_path]
            log_subprocess.check_call(" ".join(subprocess_args), shell=True)
        out_files.append(out_path)

    # need something for non-blank chunk_outs (Martian foible)
    outs.files_merged = True
    outs.merged_file_paths = out_files


def join(args, outs, chunk_args, chunk_outs):
    if args.remove_split_fastqs:
        unique_folders = set([])
        for chunk_arg, chunk_out in zip(chunk_args, chunk_outs):
            # if no files present or no files need merging, no folders necessary to remove
            if not chunk_out.files_merged:
                continue
            for input_file in chunk_arg.input_files:
                if os.path.exists(input_file):
                    martian.log_info("Removing split file: %s" % input_file)
                    os.remove(input_file)
            unique_folders.update([os.path.dirname(input_file) for input_file in chunk_arg.input_files])

        # remove empty split folders
        for folder in unique_folders:
            if not glob.glob("%s/*fastq*" % folder):
                martian.log_info("Removing empty folder: %s" % folder)
                os.rmdir(folder)

    # rename merged files (which have a _S0n FASTQ name for distinction in flat directories)
    # to final name, with _Sn
    for chunk_out in chunk_outs:
        merged_files = chunk_out.merged_file_paths
        for merged_file in merged_files:
            file_info = tk_fasta.IlmnFastqFile(merged_file)
            out_name = "%s_S%s_L%s_%s_%s.fastq.gz" % (
                file_info.prefix,
                file_info.s[1:],
                str(file_info.lane).zfill(3),
                file_info.read,
                str(file_info.group).zfill(3)
            )
            out_path = os.path.join(os.path.dirname(merged_file), os.path.basename(out_name))
            os.rename(merged_file, out_path)

    if args.remove_undetermined_fastqs:
        undetermined_fastqs = glob.glob("%s/Undetermined*fastq.gz" % args.fastq_path)
        for undetermined_fastq in undetermined_fastqs:
            martian.log_info("Removing undetermined file: %s" % undetermined_fastq)
            os.remove(undetermined_fastq)

