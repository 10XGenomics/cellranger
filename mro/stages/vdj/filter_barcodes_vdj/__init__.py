#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import itertools
import json
import numpy as np
import tenkit.safe_json as tk_safe_json
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.report as vdj_report
import cellranger.vdj.stats as vdj_stats
import cellranger.vdj.umi_info as vdj_umi_info

__MRO__ = """
stage FILTER_BARCODES_VDJ(
    in  json   barcode_counts,
    in  json   extract_reads_summary,
    in  int[]  gem_groups,
    in  string barcode_whitelist,
    in  int    recovered_cells,
    in  int    force_cells,
    in  h5     umi_info,
    in  int    min_readpairs_per_umi  "Loose threshold used in assembly",
    out json   cell_barcodes,
    out csv    barcode_support,
    out json   summary,
    out csv    barcode_umi_summary,
    src py     "stages/vdj/filter_barcodes_vdj",
) split using (
)
"""

DICT_BCS_PER_MEM_GB = 1750000

RPU_MIX_INIT_SD = 0.25
UMI_MIX_INIT_SD = 0.25

def split(args):
    # Need to store umi_info and a json with a dict containing 1 key per barcode
    umi_info_mem_gb = 2*int(np.ceil(vdj_umi_info.get_mem_gb(args.umi_info)))

    return {
        'chunks': [{
            '__mem_gb': int(np.ceil(max(cr_constants.MIN_MEM_GB, umi_info_mem_gb))),
        }]
    }


def write_barcode_umi_summary(umi_info_filename, reporter, filename, threshold, cell_barcode_set):
    """ Write a summary of UMI readpair-counts per (barcode, chain) tuple.
        Args: filename - output filename
              threshold (int) - min read pairs per UMI used in asm
              barcodes - set of barcode strings """

    # Load the umi info
    umi_info = vdj_umi_info.read_umi_info(umi_info_filename)
    chains = umi_info['chains']
    barcodes = umi_info['barcodes']

    sep = ','

    with open(filename, 'w') as writer:
        field_names = ["bc"]
        field_names += [chain + "_all_umis" for chain in reporter.vdj_genes] + \
                       [chain + "_good_umis" for chain in reporter.vdj_genes]
        writer.write(sep.join(field_names))
        writer.write("\n")

        # Assume sorted by barcode
        for bc_idx, umi_iter in itertools.groupby(itertools.izip(umi_info['barcode_idx'],
                                                                 umi_info['chain_idx'],
                                                                 umi_info['reads']),
                                                  key=lambda x: x[0]):
            bc = barcodes[bc_idx]
            if bc not in cell_barcode_set:
                continue

            # Count UMIs
            umis = list(umi_iter)
            chain_counts = defaultdict(int)
            good_chain_counts = defaultdict(int)
            for bc_idx, chain_idx, reads in umis:
                chain = chains[chain_idx]
                chain_counts[chain] += 1
                chain_counts[cr_constants.MULTI_REFS_PREFIX] += 1

                _, gem_group = cr_utils.split_barcode_seq(barcodes[bc_idx])

                if reads >= threshold:
                    good_chain_counts[chain] += 1
                    good_chain_counts[cr_constants.MULTI_REFS_PREFIX] += 1

            # Report barcode totals
            flds = {}
            flds["bc"] = bc

            num_good_umis = good_chain_counts[cr_constants.MULTI_REFS_PREFIX]
            reporter._get_metric_attr('vdj_recombinome_total_umis_per_cell_distribution').add(num_good_umis)
            reporter._get_metric_attr('vdj_recombinome_total_umis_per_cell_median').add(num_good_umis)

            # Report per-chain totals for this barcode
            for chain in reporter.vdj_genes:
                chain_all_umis = chain_counts[chain]
                chain_good_umis = good_chain_counts[chain]

                flds[chain + "_all_umis"] = chain_all_umis
                flds[chain + "_good_umis"] = chain_good_umis

                reporter._get_metric_attr('vdj_recombinome_umis_per_cell_distribution', chain).add(chain_good_umis)
                reporter._get_metric_attr('vdj_recombinome_umis_per_cell_median', chain).add(chain_good_umis)

            writer.write(sep.join([str(flds[name]) for name in field_names]))
            writer.write("\n")


def save_cell_barcodes_json(barcodes, filename):
    with open(filename, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(sorted(list(barcodes))), f, indent=4, sort_keys=True)

def call_cell_barcodes(umi_info_path, gem_group):
    """ Call cell barcodes by UMI support.
        Args: umi_info_path (str) - path to umi info h5
              gem_group (int) -  gem group
        Returns: (bc_support, cell_bcs, rt, ut)
                 where bc_support = dict of { barcode: umi_count },
                       cell_bcs = list(str) of cell barcodes)
                       rt = read pair per umi threshold used
                       ut = umi threshold """

    # Get umi info for this gem group only
    bc_idx = vdj_umi_info.get_column(umi_info_path, 'barcode_idx')
    bc_str = vdj_umi_info.get_column(umi_info_path, 'barcodes')
    bc_gg = np.array([int(cr_utils.split_barcode_seq(bc)[1]) for bc in bc_str])
    bc_in_gg = bc_gg == gem_group
    umi_in_gg = bc_in_gg[bc_idx]

    umi_read_pairs = vdj_umi_info.get_column(umi_info_path, 'reads')
    rpu_threshold, umi_threshold, bc_support, confidence = vdj_stats.call_vdj_cells(
        umi_barcode_idx=bc_idx[umi_in_gg],
        umi_read_pairs=umi_read_pairs[umi_in_gg],
        barcodes=bc_str,
        rpu_mix_init_sd=RPU_MIX_INIT_SD,
        umi_mix_init_sd=UMI_MIX_INIT_SD,
        verbosity=1,
    )

    cell_bcs = [bc for bc, umis in bc_support.iteritems() if umis >= umi_threshold ]

    return bc_support, cell_bcs, rpu_threshold, umi_threshold, confidence

def main(args, outs):
    np.random.seed(0)

    unique_gem_groups = np.unique(args.gem_groups).tolist()
    reporter = vdj_report.VdjReporter(gem_groups=unique_gem_groups)

    cell_barcodes = set()
    bc_support = defaultdict(int)

    # Load barcode whitelist
    barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)

    all_gem_groups = sorted(set(args.gem_groups))

    if args.recovered_cells:
        recovered_cells = args.recovered_cells
    else:
        recovered_cells = cr_constants.DEFAULT_TOP_BARCODE_CUTOFF * len(all_gem_groups)

    for gem_group in all_gem_groups:
        if barcode_whitelist is None:
            break

        # Load barcode raw read count distribution
        barcode_dist = cr_utils.load_barcode_dist(args.barcode_counts,
                                                  barcode_whitelist,
                                                  gem_group,
                                                  proportions=False)
        counts = np.array(barcode_dist.values())

        # Append gem group to barcode seqs
        barcodes = np.array([cr_utils.format_barcode_seq(seq, gem_group) for seq in barcode_dist.keys()])

        # Call cell barcodes
        gg_bc_support, gg_cell_bcs, rpu_threshold, umi_threshold, confidence = call_cell_barcodes(
            args.umi_info,
            int(gem_group))

        # Record the RPU and UMI thresholds
        reporter._get_metric_attr('vdj_filter_bcs_rpu_threshold',
                                  gem_group).set_value(rpu_threshold)
        reporter._get_metric_attr('vdj_filter_bcs_umi_threshold',
                                  gem_group).set_value(umi_threshold)
        reporter._get_metric_attr('vdj_filter_bcs_confidence',
                                  gem_group).set_value(confidence)

        if len(gg_bc_support) > 0:
            if args.force_cells is not None:
                sorted_bcs = map(lambda kv: kv[0],
                                 sorted(gg_bc_support.items(), key=lambda kv: kv[1], reverse=True))
                gg_cell_bcs = sorted_bcs[:min(len(sorted_bcs), args.force_cells)]

            # Update set of BCs called as cells
            cell_barcodes.update(set(gg_cell_bcs))

            # Sum BC support
            for bc, count in gg_bc_support.iteritems():
                bc_support[bc] += count

        # Load the extract_reads summary to get the total raw reads
        total_read_pairs = cr_utils.get_metric_from_json(args.extract_reads_summary, 'total_read_pairs')

        reporter.vdj_filter_barcodes_cb(cell_barcodes, barcodes, counts, total_read_pairs, recovered_cells)

    save_cell_barcodes_json(cell_barcodes, outs.cell_barcodes)

    with open(outs.barcode_support, 'w') as f:
        f.write('barcode,count\n')
        for k,v in bc_support.iteritems():
            f.write('%s,%d\n' % (k,v))

    write_barcode_umi_summary(args.umi_info,
                              reporter,
                              outs.barcode_umi_summary,
                              args.min_readpairs_per_umi,
                              cell_barcodes)

    reporter.report_summary_json(outs.summary)

def join(args, outs, chunk_defs, chunk_outs):
    cr_utils.copy(chunk_outs[0].cell_barcodes, outs.cell_barcodes)
    cr_utils.copy(chunk_outs[0].barcode_support, outs.barcode_support)
    cr_utils.copy(chunk_outs[0].summary, outs.summary)
    cr_utils.copy(chunk_outs[0].barcode_umi_summary, outs.barcode_umi_summary)
