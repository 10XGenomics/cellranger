#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
from collections import OrderedDict, defaultdict
import csv
import enum
import json
import numpy as np
import pandas as pd
import random
import martian
import tenkit.safe_json as tk_safe_json
import cellranger.cell_calling as cr_cell
import cellranger.chemistry as cr_chem
import cellranger.matrix as cr_matrix
import cellranger.stats as cr_stats
import cellranger.constants as cr_constants
import cellranger.library_constants as lib_constants
import cellranger.rna.matrix as rna_matrix
import cellranger.rna.report_matrix as rna_report_mat
import cellranger.utils as cr_utils
import cellranger.feature.antibody.analysis as ab_utils

FILTER_BARCODES_MIN_MEM_GB = 2.0

__MRO__ = """
stage FILTER_BARCODES(
    in  string sample_id,
    in  h5     matrices_h5,
    in  json   raw_fastq_summary,
    in  json   attach_bcs_summary,
    in  int    recovered_cells,
    in  int    force_cells,
    in  h5     barcode_summary,
    in  csv    barcode_correction_csv,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    in  map    chemistry_def,
    in  json   cell_barcodes          "Cell barcode override",
    out json   summary,
    out csv    filtered_barcodes,
    out csv    aggregate_barcodes,
    out h5     filtered_matrices_h5,
    out path   filtered_matrices_mex,
    out csv    nonambient_calls,
    src py     "stages/counter/filter_barcodes",
) split using (
)
"""

class FilterMethod(enum.Enum):
    # Caller-provided list of barcodes
    MANUAL = 0
    # Take the top N barcodes by count
    TOP_N_BARCODES = 1
    # Take barcodes within an order of magnitude of the max by count
    ORDMAG = 2
    # The above, then find barcodes that differ from the ambient profile
    ORDMAG_NONAMBIENT = 3

def get_filter_method_name(fm):
    if fm == FilterMethod.MANUAL:
        return 'manual'
    elif fm == FilterMethod.TOP_N_BARCODES:
        return 'topn'
    elif fm == FilterMethod.ORDMAG:
        return 'ordmag'
    elif fm == FilterMethod.ORDMAG_NONAMBIENT:
        return 'ordmag_nonambient'
    else:
        raise ValueError('Unsupported filter method value %d' % fm)


def split(args):
    # We need to store one full copy of the matrix.
    mem_gb = 2 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrices_h5)
    mem_gb = max(mem_gb, FILTER_BARCODES_MIN_MEM_GB)

    return {
        'chunks': [],
        'join': {
            '__mem_gb': mem_gb,
        }
    }

def main(_args, _outs):
    martian.throw('main is not supposed to run.')


def join(args, outs, _chunk_defs, _chunk_outs):
    filtered_matrix = filter_barcodes(args, outs)

    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, args.gem_groups, cr_chem.get_description(args.chemistry_def))
    filtered_matrix.save_h5_file(outs.filtered_matrices_h5, extra_attrs=matrix_attrs)

    rna_matrix.save_mex(filtered_matrix,
                        outs.filtered_matrices_mex,
                        martian.get_pipelines_version())

def remove_bcs_with_high_umi_corrected_reads(correction_data, matrix):
    """ Given a CountMatrix and and csv file containing information about umi corrected reads,
        detect all barcodes with unusually high fraction of corrected reads (proobably aggregates),
        and remove them from the CoutMatrix """

    bcs_to_remove, reads_lost, removed_bcs_df = ab_utils.detect_aggregate_bcs(correction_data)
    bcs_to_remove = set(matrix.bc_to_int(bc) for bc in bcs_to_remove)
    # make sure filtered_bcs is in deterministic order or any later bootstrap sampling will not be deterministic
    filtered_bcs = [i for i in xrange(matrix.bcs_dim) if i not in bcs_to_remove]
    cleaned_matrix = matrix.select_barcodes(filtered_bcs)

    ### report how many aggregates were found, and the fraction of reads those accounted for
    metrics_to_report = {}
    metrics_to_report['ANTIBODY_number_highly_corrected_GEMs'] = len(bcs_to_remove)
    metrics_to_report['ANTIBODY_reads_lost_to_highly_corrected_GEMs'] = reads_lost

    return cleaned_matrix, metrics_to_report, removed_bcs_df

def filter_barcodes(args, outs):
    random.seed(0)
    np.random.seed(0)

    correction_data = pd.read_csv(args.barcode_correction_csv)
    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrices_h5)
    if np.isin('Antibody Capture', correction_data.library_type):
    	matrix, metrics_to_report, removed_bcs_df = remove_bcs_with_high_umi_corrected_reads(correction_data, raw_matrix)
    	### report all idenitified aggregate barcodes, together with their reads, umi corrected reads, fraction of corrected reads, and fraction of total reads
    	removed_bcs_df.to_csv(outs.aggregate_barcodes)
    	summary = metrics_to_report
    else:
        matrix = raw_matrix
        summary = {}

    if args.cell_barcodes is not None:
        method = FilterMethod.MANUAL
    elif args.force_cells is not None:
        method = FilterMethod.TOP_N_BARCODES
    else:
        method = FilterMethod.ORDMAG_NONAMBIENT

    summary['total_diversity'] = matrix.bcs_dim
    summary['filter_barcodes_method'] = get_filter_method_name(method)

    # Get unique gem groups
    unique_gem_groups = sorted(list(set(args.gem_groups)))

    # Get per-gem group cell load
    if args.recovered_cells is not None:
        gg_recovered_cells = int(float(args.recovered_cells) / float(len(unique_gem_groups)))
    else:
        gg_recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    if args.force_cells is not None:
        gg_force_cells = int(float(args.force_cells) / float(len(unique_gem_groups)))

    # Only use gene expression matrix for cell calling
    gex_matrix = matrix.view().select_features_by_type(lib_constants.GENE_EXPRESSION_LIBRARY_TYPE)

    # Make initial cell calls for each genome separately
    genomes = gex_matrix.get_genomes()

    # (gem_group, genome) => dict
    filtered_metrics_groups = OrderedDict()
    # (gem_group, genome) => list of barcode strings
    filtered_bcs_groups = OrderedDict()

    for genome in genomes:
        genome_matrix = gex_matrix.select_features_by_genome(genome)

        # Make initial cell calls for each gem group individually
        for gem_group in unique_gem_groups:

            gg_matrix = genome_matrix.select_barcodes_by_gem_group(gem_group)

            if method == FilterMethod.ORDMAG or \
               method == FilterMethod.ORDMAG_NONAMBIENT:
                gg_total_diversity = gg_matrix.bcs_dim
                gg_bc_counts = gg_matrix.get_counts_per_bc()
                gg_filtered_indices, gg_filtered_metrics, msg = cr_stats.filter_cellular_barcodes_ordmag(
                    gg_bc_counts, gg_recovered_cells, gg_total_diversity)
                gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

            elif method == FilterMethod.MANUAL:
                with(open(args.cell_barcodes)) as f:
                    cell_barcodes = json.load(f)
                gg_filtered_bcs, gg_filtered_metrics, msg = cr_stats.filter_cellular_barcodes_manual(
                    gg_matrix, cell_barcodes)

            elif method == FilterMethod.TOP_N_BARCODES:
                gg_bc_counts = gg_matrix.get_counts_per_bc()
                gg_filtered_indices, gg_filtered_metrics, msg = cr_stats.filter_cellular_barcodes_fixed_cutoff(
                    gg_bc_counts, gg_force_cells)
                gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

            else:
                martian.exit("Unsupported BC filtering method: %s" % method)

            if msg is not None:
                martian.log_info(msg)

            filtered_metrics_groups[(gem_group, genome)] = gg_filtered_metrics
            filtered_bcs_groups[(gem_group, genome)] = gg_filtered_bcs

    # Do additional cell calling
    outs.nonambient_calls = None

    if method == FilterMethod.ORDMAG_NONAMBIENT:
        # We need the full gene expression matrix instead of just a view
        full_gex_matrix = matrix.select_features_by_type(lib_constants.GENE_EXPRESSION_LIBRARY_TYPE)

        # Track these for recordkeeping
        eval_bcs_arrays = []
        umis_per_bc_arrays = []
        loglk_arrays = []
        pvalue_arrays = []
        pvalue_adj_arrays = []
        nonambient_arrays = []
        genome_call_arrays = []

        # Do it by gem group, but agnostic to genome
        for gg in unique_gem_groups:
            gg_matrix = full_gex_matrix.select_barcodes_by_gem_group(gg)

            # Take union of initial cell calls across genomes
            gg_bcs = sorted(list(reduce(set.union,
                                        [set(bcs) for group, bcs in filtered_bcs_groups.iteritems() if group[0] == gg])))

            result = cr_cell.find_nonambient_barcodes(gg_matrix, gg_bcs)
            if result is None:
                print 'Failed at attempt to call non-ambient barcodes in GEM group %s' % gg
                continue

            # Assign a genome to the cell calls by argmax genome counts
            genome_counts = []
            for genome in genomes:
                genome_counts.append(gg_matrix.view() \
                                     .select_features_by_genome(genome) \
                                     .select_barcodes(result.eval_bcs) \
                                     .get_counts_per_bc())
            genome_counts = np.column_stack(genome_counts)
            genome_calls = np.array(genomes)[np.argmax(genome_counts, axis=1)]

            umis_per_bc = gg_matrix.get_counts_per_bc()

            eval_bcs_arrays.append(np.array(gg_matrix.bcs)[result.eval_bcs])
            umis_per_bc_arrays.append(umis_per_bc[result.eval_bcs])
            loglk_arrays.append(result.log_likelihood)
            pvalue_arrays.append(result.pvalues)
            pvalue_adj_arrays.append(result.pvalues_adj)
            nonambient_arrays.append(result.is_nonambient)
            genome_call_arrays.append(genome_calls)

            # Update the lists of cell-associated barcodes
            for genome in genomes:
                eval_bc_strs = np.array(gg_matrix.bcs)[result.eval_bcs]
                filtered_bcs_groups[(gg, genome)].extend(eval_bc_strs[(genome_calls == genome) & (result.is_nonambient)])

        if len(eval_bcs_arrays) > 0:
            nonambient_summary = pd.DataFrame(OrderedDict([
                ('barcode', np.concatenate(eval_bcs_arrays)),
                ('umis', np.concatenate(umis_per_bc_arrays)),
                ('ambient_loglk', np.concatenate(loglk_arrays)),
                ('pvalue', np.concatenate(pvalue_arrays)),
                ('pvalue_adj', np.concatenate(pvalue_adj_arrays)),
                ('nonambient', np.concatenate(nonambient_arrays)),
                ('genome', np.concatenate(genome_call_arrays)),
            ]))
            nonambient_summary.to_csv(outs.nonambient_calls)

    # Record all filtered barcodes
    genome_filtered_bcs = defaultdict(set)
    filtered_bcs = set()
    for (gem_group, genome), bcs in filtered_bcs_groups.iteritems():
        genome_filtered_bcs[genome].update(bcs)
        filtered_bcs.update(bcs)

    # Combine initial-cell-calling metrics
    for genome in genomes:
        # Merge metrics over all gem groups for this genome
        txome_metrics = [v for k,v in filtered_metrics_groups.iteritems() if k[1] == genome]
        txome_summary = cr_stats.merge_filtered_metrics(txome_metrics)

        # Append method name to metrics
        summary.update({
            ('%s_%s_%s' % (genome,
                           key,
                           get_filter_method_name(method))): txome_summary[key] \
            for (key,_) in txome_summary.iteritems()})

        summary['%s_filtered_bcs' % genome] = len(genome_filtered_bcs[genome])

        # NOTE: This metric only applies to the initial cell calls
        summary['%s_filtered_bcs_cv' % genome] = txome_summary['filtered_bcs_cv']

    # Deduplicate and sort filtered barcode sequences
    # Sort by (gem_group, barcode_sequence)
    barcode_sort_key = lambda x: cr_utils.split_barcode_seq(x)[::-1]

    for genome, bcs in genome_filtered_bcs.iteritems():
        genome_filtered_bcs[genome] = sorted(list(set(bcs)), key=barcode_sort_key)
    filtered_bcs = sorted(list(set(filtered_bcs)), key=barcode_sort_key)

    # Re-compute various metrics on the filtered matrix
    reads_summary = cr_utils.merge_jsons_as_dict([args.raw_fastq_summary, args.attach_bcs_summary])
    matrix_summary = rna_report_mat.report_genomes(matrix,
                                                   reads_summary=reads_summary,
                                                   barcode_summary_h5_path=args.barcode_summary,
                                                   recovered_cells=args.recovered_cells,
                                                   cell_bc_seqs=genome_filtered_bcs)

    # Write metrics json
    combined_summary = matrix_summary.copy()
    combined_summary.update(summary)
    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(combined_summary), f, indent=4, sort_keys=True)

    # Write the filtered barcodes file
    write_filtered_barcodes(outs.filtered_barcodes, genome_filtered_bcs)

    # Select cell-associated barcodes
    filtered_matrix = matrix.select_barcodes_by_seq(filtered_bcs)

    return filtered_matrix

def write_filtered_barcodes(out_csv, bcs_per_genome):
    """ Args:
        bcs_per_genome (dict of str to list): Map each genome to its cell-associated barcodes
    """
    with open(out_csv, 'w') as f:
        writer = csv.writer(f)
        for (genome, bcs) in bcs_per_genome.iteritems():
            for bc in bcs:
                writer.writerow([genome, bc])
