#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Filter out noise UMIs.

import itertools
import numpy as np
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.report as vdj_report
import cellranger.vdj.stats as vdj_stats
import cellranger.vdj.umi_info as vdj_umi_info

__MRO__ = """
stage FILTER_UMIS(
    in  h5     umi_info,
    in  path   vdj_reference_path,
    in  float  intra_barcode_nx,
    in  int[]  gem_groups,
    in  int    target_n50,
    out map    min_readpairs_per_umi,
    out map    subsample_rate,
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

# Filtering strategy:
# 1. Remove very low support UMIs by applying a filter on reads/UMI for each barcode separately.
# 2. Across all remaining UMIs, get N50. Compute sampling rate based on this N50 and
# the target reads/UMI.
# 3. For each chain get the read/UMI distribution, subsample it according to the rate from above,
# and then get the read/UMI cutoff on the subsampled distribution
# as the low component of a 2-class K-means.
# 4. Get the min threshold across all chains (that represent N90 of reads).
# 5. Pass the subsampling rate and read/UMI cutoff from the previous step to the assembler.
# The cutoff is based on subsampled data and will be applied after the assembler does it's internal subsampling.
def main(args, outs):
    np.random.seed(0)

    unique_gem_groups = np.unique(args.gem_groups).tolist()

    reporter = vdj_report.VdjReporter(gem_groups=unique_gem_groups)

    # Load the umi info
    umi_info = vdj_umi_info.read_umi_info(args.umi_info, args.start_row, args.end_row)


    # Compute initial within-barcode thresholds
    # Assumes the fraction of noise-UMI reads is < some fraction (NX)
    barcode_nx = np.zeros(len(umi_info['barcodes']), dtype=int)

    # Assume grouped by barcode
    for bc, bc_reads in itertools.groupby(itertools.izip(umi_info['barcode_idx'],
                                                         umi_info['reads']),
                                          key=lambda x: x[0]):
        bc_reads_arr = np.fromiter((reads for bc,reads in bc_reads), umi_info['reads'].dtype)
        barcode_nx[bc] = tk_stats.NX(bc_reads_arr, args.intra_barcode_nx)


    # Filter out UMIs below the within-BC threshold (in-place)
    top_in_bc = umi_info['reads'] >= barcode_nx[umi_info['barcode_idx']]

    for col in vdj_umi_info.UMI_INFO_COLS.iterkeys():
        umi_info[col] = np.compress(top_in_bc, umi_info[col])


    # Compute N50 read pairs per UMI for this gem group
    # and use it to subsample to the target N50.
    rppu_n50 = tk_stats.NX(umi_info['reads'], 0.5)
    if rppu_n50 is None:
        rppu_n50 = float('NaN')

    reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_n50',
                              cr_constants.MULTI_REFS_PREFIX, args.gem_group).set_value(rppu_n50)
    if rppu_n50 == 0:
        subsample_rate = 1.0
    else:
        subsample_rate = min(1.0, tk_stats.robust_divide(args.target_n50, rppu_n50))

    reporter._get_metric_attr('vdj_assembly_subsample_rate',
                              args.gem_group).set_value(subsample_rate, 1.0)

    # Weighted average of subsample rates where weight = sum of readpairs on UMIs for each gem-group
    reporter._get_metric_attr('vdj_assembly_overall_subsample_rate').set_value(subsample_rate*sum(umi_info['reads']),
                                                                               sum(umi_info['reads']))

    # Find the global (per-chain) thresholds
    thresholds = {}
    chain_totals = {}

    # Sort the chains alphabetically for determinism in e.g. multi-library vs single-library
    #   runs.
    chain_tuples = list(enumerate(umi_info['chains']))
    sorted_chain_tuples = sorted(chain_tuples, key=lambda x: x[1])

    for chain_idx, chain in sorted_chain_tuples:
        chain_reads = umi_info['reads'][umi_info['chain_idx'] == chain_idx]
        chain_totals[chain] = chain_reads.sum()

        # Record the per-chain N50 read pairs per UMI (but don't use it)
        chain_n50 = tk_stats.NX(chain_reads, 0.5)
        if chain_n50 is None:
            chain_n50 = float('NaN')
        reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_n50',
                                  chain, args.gem_group).set_value(chain_n50)

        print "Computing per-chain threshold for %s" % chain

        thresholds[chain] = vdj_stats.compute_readpairs_per_umi_threshold(chain_reads, subsample_rate)

        print "  %d" % thresholds[chain]

        reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_threshold',
                                  chain, args.gem_group).set_value(thresholds[chain])

    # Take the min threshold among the chains that make up N90 of all reads
    chain_n90 = tk_stats.NX(chain_totals.values(), 0.9)

    use_chains = [chain for chain in thresholds.iterkeys() if chain_totals[chain] >= chain_n90]
    use_thresholds = [thresholds[c] for c in use_chains]

    print "Using thresholds from " + str(use_chains) + ": " + str(use_thresholds)

    # Handle case where no chains were detected
    if len(use_chains) == 0:
        threshold = 1
    else:
        threshold = min(use_thresholds)

    outs.min_readpairs_per_umi = {args.gem_group: int(threshold)}
    outs.subsample_rate = {args.gem_group: float(subsample_rate)}

    reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_threshold',
                              cr_constants.MULTI_REFS_PREFIX, args.gem_group).set_value(threshold)

    reporter.save(outs.chunked_reporter)


def join(args, outs, chunk_defs, chunk_outs):
    reporters = [chunk_out.chunked_reporter for chunk_out in chunk_outs]
    final_report = cr_report.merge_reporters(reporters)
    final_report.report_summary_json(outs.summary)
    outs.chunked_reporter = None

    outs.min_readpairs_per_umi = {}
    outs.subsample_rate = {}
    for chunk_def, chunk_out in zip(chunk_defs, chunk_outs):
        outs.min_readpairs_per_umi.update(chunk_out.min_readpairs_per_umi)
        outs.subsample_rate.update(chunk_out.subsample_rate)
