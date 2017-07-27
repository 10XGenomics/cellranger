#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import csv
import json
import numpy as np
import random
import martian
import tenkit.safe_json as tk_safe_json
import cellranger.chemistry as cr_chem
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import cellranger.stats as cr_stats
import cellranger.constants as cr_constants

__MRO__ = """
stage FILTER_BARCODES(
    in  string sample_id,
    in  h5     matrices_h5,
    in  json   raw_fastq_summary,
    in  json   attach_bcs_summary,
    in  int    recovered_cells,
    in  int    force_cells,
    in  h5     barcode_summary,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    in  map    chemistry_def,
    in  json   cell_barcodes          "Cell barcode override",
    out json   summary,
    out csv    filtered_barcodes,
    out h5     filtered_matrices_h5,
    out path   filtered_matrices_mex,
    src py     "stages/counter/filter_barcodes",
) split using (
)
"""

def split(args):
    mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups)
    chunks = [{
        '__mem_gb': mem_gb,
    }]
    return {'chunks': chunks}

def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]
    cr_utils.copy(chunk_out.summary, outs.summary)
    cr_utils.copy(chunk_out.filtered_matrices_h5, outs.filtered_matrices_h5)
    cr_utils.copy(chunk_out.filtered_barcodes, outs.filtered_barcodes)
    cr_utils.copytree(chunk_out.filtered_matrices_mex, outs.filtered_matrices_mex)

def main(args, outs):
    filtered_matrices = filter_barcodes(args, outs)
    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, args.gem_groups, cr_chem.get_description(args.chemistry_def))
    filtered_matrices.save_h5(outs.filtered_matrices_h5, extra_attrs=matrix_attrs)
    filtered_matrices.save_mex(outs.filtered_matrices_mex)

def filter_barcodes(args, outs):
    random.seed(0)
    np.random.seed(0)

    matrices = cr_matrix.GeneBCMatrices.load_h5(args.matrices_h5)

    summary = {}

    total_diversity = len(matrices.matrices.values()[-1].bcs)

    if args.cell_barcodes is not None:
        method_name = cr_constants.FILTER_BARCODES_MANUAL
    elif args.force_cells is not None:
        method_name = cr_constants.FILTER_BARCODES_FIXED_CUTOFF
    else:
        method_name = cr_constants.FILTER_BARCODES_ORDMAG

    summary['total_diversity'] = total_diversity
    summary['filter_barcodes_method'] = method_name

    # Initialize filtered matrices object
    filtered_matrices = cr_matrix.GeneBCMatrices(matrices.matrices.keys(),
                                                 [m.genes for m in matrices.matrices.values()],
                                                 [m.bcs for m in matrices.matrices.values()][0])

    # Get unique gem groups
    unique_gem_groups = sorted(list(set(args.gem_groups)))

    # Get per-gem group cell load
    if args.recovered_cells is not None:
        gg_recovered_cells = int(float(args.recovered_cells) / float(len(unique_gem_groups)))
    else:
        gg_recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    if args.force_cells is not None:
        gg_force_cells = int(float(args.force_cells) / float(len(unique_gem_groups)))

    filtered_metrics = []
    filtered_bcs = []

    # Track filtered barcodes for each genome
    bcs_per_genome = collections.defaultdict(list)

    # Filter each genome's matrix
    for genome, matrix in matrices.matrices.iteritems():
        filtered_metrics = []
        filtered_bcs = []

        # Filter each gem group individually
        for gem_group in unique_gem_groups:
            gg_matrix = matrix.select_barcodes_by_gem_group(gem_group)
            if method_name == cr_constants.FILTER_BARCODES_ORDMAG:
                gg_total_diversity = len(gg_matrix.bcs)
                gg_bc_counts = gg_matrix.get_reads_per_bc()
                gg_filtered_indices, gg_filtered_metrics, msg = cr_stats.filter_cellular_barcodes_ordmag(
                    gg_bc_counts, gg_recovered_cells, gg_total_diversity)
                gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)
            elif method_name == cr_constants.FILTER_BARCODES_MANUAL:
                with(open(args.cell_barcodes)) as f:
                    cell_barcodes = json.load(f)
                gg_filtered_bcs, gg_filtered_metrics, msg = cr_stats.filter_cellular_barcodes_manual(
                    gg_matrix, cell_barcodes)
            elif method_name == cr_constants.FILTER_BARCODES_FIXED_CUTOFF:
                gg_bc_counts = gg_matrix.get_reads_per_bc()
                gg_filtered_indices, gg_filtered_metrics, msg = cr_stats.filter_cellular_barcodes_fixed_cutoff(
                    gg_bc_counts, gg_force_cells)
                gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)
            else:
                martian.exit("Unsupported BC filtering method: %s" % method_name)

            if msg is not None:
                martian.log_info(msg)

            filtered_metrics.append(gg_filtered_metrics)
            filtered_bcs.extend(gg_filtered_bcs)

            bcs_per_genome[genome].extend(gg_filtered_bcs)

        # Merge metrics over all gem groups
        txome_summary = cr_stats.merge_filtered_metrics(filtered_metrics)

        # Append method name to metrics
        summary.update({
            ('%s_%s_%s' % (genome, key, method_name)): txome_summary[key] \
            for (key,_) in txome_summary.iteritems()})

        txome_filtered_matrix = matrix.select_barcodes_by_seq(filtered_bcs)
        filtered_matrices.matrices[genome] = txome_filtered_matrix
        summary['%s_filtered_bcs' % genome] = txome_summary['filtered_bcs']
        summary['%s_filtered_bcs_cv' % genome] = txome_summary['filtered_bcs_cv']

    # Re-compute various metrics on the filtered matrices
    matrix_summary = matrices.report(
        summary_json_paths=[args.raw_fastq_summary, args.attach_bcs_summary],
        barcode_summary_h5_path=args.barcode_summary,
        recovered_cells=args.recovered_cells,
        cell_bc_seqs=[mat.bcs for mat in filtered_matrices.matrices.itervalues()])

    # Write summary json
    combined_summary = matrix_summary.copy()
    combined_summary.update(summary)
    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(combined_summary), f, indent=4, sort_keys=True)

    # Write the filtered barcodes file
    write_filtered_barcodes(outs.filtered_barcodes, bcs_per_genome)

    return filtered_matrices

def write_filtered_barcodes(out_csv, bcs_per_genome):
    """ Args: bcs_per_genome - [genome1, genome2, ...]
              where genome1 = ['ACGT-1', ...] """
    with open(out_csv, 'w') as f:
        writer = csv.writer(f)
        for (genome, bcs) in bcs_per_genome.iteritems():
            for bc in bcs:
                writer.writerow([genome, bc])
