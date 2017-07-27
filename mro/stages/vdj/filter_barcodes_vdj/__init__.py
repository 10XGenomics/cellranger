#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import csv
import itertools
import json
import numpy as np
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.report as vdj_report
import cellranger.vdj.umi_info as vdj_umi_info

__MRO__ = """
stage FILTER_BARCODES_VDJ(
    in  json   barcode_counts,
    in  tsv    contig_summary            "Assembler contig summary",
    in  tsv    umi_summary               "Assembler UMI summary",
    in  int    min_umis                  "Minimum passing UMIs in the top contig for a cell",
    in  json   extract_reads_summary,
    in  int[]  gem_groups,
    in  string barcode_whitelist,
    in  int    recovered_cells,
    in  int    force_cells,
    in  h5     umi_info,
    in  map    min_readpairs_per_umi     "Loose thresholds used in assembly",
    in  float  readpairs_per_umi_nx      "Estimate of max readpairs per umi",
    in  float  readpairs_per_umi_ratio   "Divide above estimate by this for strict-thresh",
    in  json   assemble_metrics_summary,
    out json   cell_barcodes,
    out csv    barcode_support,
    out json   summary,
    out csv    barcode_umi_summary,
    src py     "stages/vdj/filter_barcodes_vdj",
) split using (
)
"""
DICT_BCS_PER_MEM_GB = 1750000

def split(args):
    # Need to store umi_info and a json with a dict containing 1 key per barcode
    umi_info_mem_gb = 2*int(np.ceil(vdj_umi_info.get_mem_gb(args.umi_info)))

    bc_diversity = len(cr_utils.load_barcode_whitelist(args.barcode_whitelist))
    assemble_summary_mem_gb = tk_stats.robust_divide(bc_diversity, DICT_BCS_PER_MEM_GB)

    return {
        'chunks': [{
            '__mem_gb': int(np.ceil(max(cr_constants.MIN_MEM_GB, umi_info_mem_gb + assemble_summary_mem_gb))),
        }]
    }


def write_barcode_umi_summary(umi_info_filename, reporter, filename, thresholds, cell_barcode_set):
    """ Write a summary of UMI readpair-counts per (barcode, chain) tuple.
        Args: filename - output filename
              thresholds - good UMI threshold in readpairs (dict of str(gem_group) : int(threshold))
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
                threshold = thresholds[str(gem_group)]

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

def call_cell_barcodes(umi_summary_filename, gem_group, min_umis, threshold_nx, threshold_ratio):
    """ Call cell barcodes by contig/UMI read support.
        umi_summary_filename (str) - path to umi summary tsv generated by vdj_asm
        gem_group (int) -  gem group
        min_umis (int) - min passing UMIs on highest-passing-UMI-contig to call cell
        Returns: (d,b,t)
                 where d = dict of { barcode: best_contig_kth_umi_readpairs },
                            k = min_umis and
                            kth_umi_readpairs = 0 if best_contig has <k umis,
                       b = list(str) of cell barcodes)
                       t = read pair threshold used """


    with open(umi_summary_filename) as f:
        # First pass: compute threshold
        reader = csv.reader(f, delimiter='\t')

        hdr = next(reader)
        bc_col = hdr.index('barcode')
        umi_col = hdr.index('umi')
        reads_col = hdr.index('reads')
        thresh_col = hdr.index('min_umi_reads')
        good_col = hdr.index('good_umi')
        contigs_col = hdr.index('contigs')

        def use_umi(row):
            return (row[umi_col] != '') and \
                (row[contigs_col] != '') and \
                (row[good_col] == 'True')

        read_pairs = []
        assembly_rppu_threshold = 1

        bc_support = {}

        for row in reader:
            # Only take this gem group
            _, gg = cr_utils.split_barcode_seq(row[bc_col])
            if str(gg) != str(gem_group):
                continue

            # Initialize all barcodes
            bc_support[row[bc_col]] = 0

            if not use_umi(row):
                continue

            # Get the RPPU threshold that was used in assembly
            # The tsv reports reads per UMI, so divide by 2 for pairs.
            assembly_rppu_threshold = int(row[thresh_col])/2
            read_pairs.append(int(row[reads_col])/2)

        read_pairs = np.array(read_pairs, dtype=int)

        # Estimate the high end of the distribution
        if len(read_pairs) > 0:
            high_rppu = tk_stats.NX(read_pairs, threshold_nx)
        else:
            high_rppu = 1

        # Take UMIs within X of the high end, roughly corresponding the to highest mode
        # and therefore to molecules amplified from the first cycle.
        threshold = int(round(tk_stats.robust_divide(high_rppu, threshold_ratio)))

        # Don't drop below the looser threshold that was used in assembly.
        threshold = max(assembly_rppu_threshold, threshold)


        # Second pass: Call as cell BCs those with at least k UMIs
        # passing the strict threshold computed above.
        f.seek(0)
        reader = csv.reader(f, delimiter='\t')
        next(reader)

        cell_barcodes = []

        good_umi_iter = itertools.ifilter(use_umi, reader)
        bc_group_iter = itertools.groupby(good_umi_iter, key=lambda row: row[bc_col])

        for bc, rows in bc_group_iter:
            # Restrict to the current gem group
            bc_seq, gg = cr_utils.split_barcode_seq(bc)
            if str(gg) != str(gem_group):
                continue

            # Collect readpair support for all UMIs for all contigs
            contig_umis_readpairs = defaultdict(list)
            for row in rows:
                contig_umis_readpairs[row[contigs_col]].append(int(row[reads_col])/2)

            # Get the max (contig-kth-umi)
            best_kth_umi_readpairs = 0

            for contig, umi_readpairs in contig_umis_readpairs.iteritems():
                # Sort UMIs by readpairs, descending
                umi_readpairs = np.array(umi_readpairs, dtype=int)
                umi_readpairs[::-1].sort()

                # Get the kth UMI's readpair support or 0
                if len(umi_readpairs) >= min_umis:
                    kth_umi_readpairs = umi_readpairs[min_umis-1]
                else:
                    kth_umi_readpairs = 0

                best_kth_umi_readpairs = max(best_kth_umi_readpairs, kth_umi_readpairs)

            bc_support[bc] = best_kth_umi_readpairs

            if best_kth_umi_readpairs >= threshold:
                cell_barcodes.append(bc)

        return bc_support, cell_barcodes, threshold


def main(args, outs):
    unique_gem_groups = np.unique(args.gem_groups).tolist()
    reporter = vdj_report.VdjReporter(gem_groups=unique_gem_groups)

    cell_barcodes = set()
    bc_support = {}

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
        gg_bc_support, gg_cell_bcs, threshold = call_cell_barcodes(args.umi_summary, int(gem_group),
                                                                   args.min_umis,
                                                                   args.readpairs_per_umi_nx,
                                                                   args.readpairs_per_umi_ratio)

        # Record the threshold
        reporter._get_metric_attr('vdj_filtered_bc_contig_kth_umi_readpair_threshold',
                                  gem_group).set_value(threshold)

        if len(gg_bc_support) > 0:
            if args.force_cells is not None:
                sorted_bcs = map(lambda kv: kv[0],
                                 sorted(gg_bc_support.items(), key=lambda kv: kv[1], reverse=True))
                gg_cell_bcs = sorted_bcs[:min(len(sorted_bcs), args.force_cells)]

            cell_barcodes.update(set(gg_cell_bcs))
            bc_support.update(gg_bc_support)

        # Load the extract_reads summary to get the total raw reads
        total_read_pairs = cr_utils.get_metric_from_json(args.extract_reads_summary, 'total_read_pairs')

        # Load the assembly metrics summary to get the total assemblable reads
        assemblable_read_pairs_by_bc = cr_utils.get_metric_from_json(args.assemble_metrics_summary, 'assemblable_read_pairs_by_bc')
        assemblable_read_pairs = sum(assemblable_read_pairs_by_bc.get(bc, 0) for bc in cell_barcodes)

        reporter.vdj_filter_barcodes_cb(cell_barcodes, barcodes, counts, total_read_pairs, assemblable_read_pairs, recovered_cells)

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
