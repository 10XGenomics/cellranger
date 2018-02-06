#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Determine assembly subsampling rate.

import numpy as np
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.report as vdj_report
import cellranger.vdj.umi_info as vdj_umi_info

__MRO__ = """
stage FILTER_UMIS(
    in  h5     umi_info,
    in  path   vdj_reference_path,
    in  int[]  gem_groups,
    in  int    target_n50,
    out pickle chunked_reporter,
    out json   summary,
    src py     "stages/vdj/filter_umis",
) split using (
    in  int    gem_group,
    in  int    start_row,
    in  int    end_row,
)
"""

def split(args):
    """ Chunk the UMI info HDF5 file by gem group """

    num_entries = vdj_umi_info.get_num_rows(args.umi_info)
    if num_entries > 1e9:
        print 'Warning: There are >1e9 entries in the umi_info - this could potentially cause an out-of-memory error.'

    # This will cause an OOM if there are >1.5e9 UMIs
    barcode_indices = vdj_umi_info.get_column(args.umi_info, 'barcode_idx')
    barcodes = vdj_umi_info.get_column(args.umi_info, 'barcodes')

    chunks = []

    start_row = 0
    prev_gem_group = None
    prev_barcode_idx = None

    for row, barcode_idx in enumerate(barcode_indices):
        if barcode_idx == prev_barcode_idx:
            continue

        _, gem_group = cr_utils.split_barcode_seq(barcodes[barcode_idx])

        if prev_gem_group is not None and gem_group != prev_gem_group:
            # Write complete chunk
            end_row = row
            mem_gb = max(cr_constants.MIN_MEM_GB,
                         2*int(np.ceil(vdj_umi_info.get_mem_gb(args.umi_info,
                                                               start_row=start_row,
                                                               end_row=end_row))))

            chunks.append({
                'gem_group': prev_gem_group,
                'start_row': start_row,
                'end_row': end_row,
                '__mem_gb': mem_gb,
            })

            start_row = end_row

        prev_gem_group = gem_group
        prev_barcode_idx = barcode_idx

    # Write final chunk
    end_row = vdj_umi_info.get_num_rows(args.umi_info)
    mem_gb = max(cr_constants.MIN_MEM_GB,
                 2*int(np.ceil(vdj_umi_info.get_mem_gb(args.umi_info,
                                                       start_row=start_row,
                                                       end_row=end_row))))

    # Handle case where umi info is empty by supplying a dummy gem group
    if prev_gem_group is None:
        prev_gem_group = args.gem_groups[0]

    chunks.append({
        'gem_group': prev_gem_group,
        'start_row': start_row,
        'end_row': end_row,
        '__mem_gb': mem_gb,
    })


    return {'chunks': chunks}

# Subsampling strategy:
# 1. Compute N50 read pairs per UMI. Compute sampling rate based on this N50 and
# the target reads/UMI.
# 2. Pass the subsampling rate to the assembler.
def main(args, outs):
    np.random.seed(0)

    unique_gem_groups = np.unique(args.gem_groups).tolist()

    reporter = vdj_report.VdjReporter(gem_groups=unique_gem_groups)

    # Load the umi info
    umi_info = vdj_umi_info.read_umi_info(args.umi_info, args.start_row, args.end_row)

    # Compute N50 read pairs per UMI for this gem group
    rppu_n50 = tk_stats.NX(umi_info['reads'], 0.5)
    if rppu_n50 is None:
        rppu_n50 = float('NaN')

    reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_n50',
                              cr_constants.MULTI_REFS_PREFIX, args.gem_group).set_value(rppu_n50)

    reporter.save(outs.chunked_reporter)


def join(args, outs, chunk_defs, chunk_outs):
    # Merge metrics summaries
    reporters = [chunk_out.chunked_reporter for chunk_out in chunk_outs]
    final_report = cr_report.merge_reporters(reporters)
    final_report.report_summary_json(outs.summary)
    outs.chunked_reporter = None
