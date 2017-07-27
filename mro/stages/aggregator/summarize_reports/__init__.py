#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import json
import martian
import numpy as np

import tenkit.stats as tk_stats

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import cellranger.molecule_counter as cr_mol_counter
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.data as cr_webshim_data
from cellranger.webshim.constants.shared import PIPELINE_AGGR

__MRO__ = """
stage SUMMARIZE_AGGREGATED_REPORTS(
    in  string aggregation_id,
    in  string aggregation_desc,
    in  map    gem_group_index,
    in  h5     filtered_matrices_h5,
    in  h5     barcode_summary_h5,
    in  path   analysis,
    in  json   normalize_depth_summary,
    in  json   count_genes_summary,
    in  json   analyze_matrices_summary,
    out json   summary,
    out html   web_summary,
    src py     "stages/aggregator/summarize_reports",
) split using (
)
"""

def split(args):
    matrix_mem_gb = cr_matrix.GeneBCMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrices_h5)
    chunks = [{
        '__mem_gb': matrix_mem_gb,
    }]
    return {'chunks': chunks}

def main(args, outs):
    summary = {}

    # add stats from matrices
    filtered_mats = cr_matrix.GeneBCMatrices.load_h5(args.filtered_matrices_h5)
    genomes = filtered_mats.get_genomes()
    cells_per_genome = {}
    for genome in genomes:
        matrix = filtered_mats.matrices[genome]
        cells_per_genome[genome] = matrix.bcs_dim
        median_gene_counts = np.median(matrix._sum(matrix.m >= cr_constants.MIN_READS_PER_GENE, axis=0))
        median_umi_counts = np.median(matrix._sum(matrix.m, axis=0))
        summary.update({
            '%s_filtered_bcs' % genome: cells_per_genome[genome],
            '%s_filtered_bcs_median_counts' % genome: median_umi_counts,
            '%s_filtered_bcs_median_unique_genes_detected' % genome: median_gene_counts,
        })
    del filtered_mats

    # get metrics from other summaries
    if args.analyze_matrices_summary:
        with open(args.analyze_matrices_summary) as reader:
            analysis_summary = json.load(reader)
        summary.update(analysis_summary)

    with open(args.normalize_depth_summary, 'r') as reader:
        data = json.load(reader)
        raw_conf_mapped_per_genome = data['raw_conf_mapped_per_genome']
        downsample_map = data['downsample_info']
        mol_counter_metrics = data['mol_counter_metrics']

    with open(args.count_genes_summary, 'r') as reader:
        data = json.load(reader)
        flt_conf_mapped_per_genome = data['flt_conf_mapped_per_genome']

    for genome in flt_conf_mapped_per_genome:
        frac_reads_in_cells = tk_stats.robust_divide(flt_conf_mapped_per_genome[genome], raw_conf_mapped_per_genome[genome])
        summary['%s_filtered_bcs_conf_mapped_barcoded_reads_cum_frac' % genome] = frac_reads_in_cells

    # Pass chemistry metrics through to output
    summary.update({k:v for k,v in mol_counter_metrics.iteritems() if k.startswith('chemistry_')})

    # Molecule counter metrics
    gem_groups = []
    total_reads_per_gem_group = []
    downsampled_reads_per_gem_group = []
    for (gg, submetrics) in mol_counter_metrics[cr_mol_counter.GEM_GROUPS_METRIC].iteritems():
        gem_groups.append(gg)
        total_reads = submetrics[cr_mol_counter.GG_TOTAL_READS_METRIC]
        total_reads_per_gem_group.append(total_reads)
        # If metric is missing, assume no downsampling was done
        downsampled = submetrics.get(cr_mol_counter.GG_DOWNSAMPLED_READS_METRIC,
                                     total_reads)
        downsampled_reads_per_gem_group.append(downsampled)
    total_reads = sum(total_reads_per_gem_group)
    downsampled_reads = sum(downsampled_reads_per_gem_group)
    total_cells = sum(cells_per_genome.values())
    mean_reads_per_cell = tk_stats.robust_divide(total_reads, total_cells)
    downsampled_mean_reads_per_cell = tk_stats.robust_divide(downsampled_reads, total_cells)
    summary.update({
        'pre_normalization_total_reads': total_reads,
        'post_normalization_total_reads': downsampled_reads,
        'filtered_bcs_transcriptome_union': total_cells,
        'pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc': mean_reads_per_cell,
        'post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc': downsampled_mean_reads_per_cell,
    })

    # Downsampling metrics
    gem_group_index = args.gem_group_index
    agg_batches = []
    lowest_frac_reads_kept = 1.0
    for (gg, rpg) in zip(gem_groups, total_reads_per_gem_group):
        dinfo = downsample_map[str(gg)]
        (library_id, old_gg) = gem_group_index[str(gg)]
        batch = library_id + ('-%d' % old_gg if old_gg > 1 else '')
        agg_batches.append(batch)
        # calc summary metrics
        frac_reads_kept = dinfo['frac_reads_kept']
        lowest_frac_reads_kept = min(lowest_frac_reads_kept, frac_reads_kept)
        summary['%s_frac_reads_kept' % batch] = frac_reads_kept
        summary['%s_pre_normalization_raw_reads_per_filtered_bc' % batch] = tk_stats.robust_divide(dinfo['total_reads'], dinfo['cells'])
        summary['%s_pre_normalization_cmb_reads_per_filtered_bc' % batch] = tk_stats.robust_divide(dinfo['cmb_reads'], dinfo['cells'])
        # this is an internal metric, so keep using gem group instead of batch
        summary['%s_total_reads_per_gem_group' % gg] = frac_reads_kept * rpg
    summary['lowest_frac_reads_kept'] = lowest_frac_reads_kept

    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)

    # build web summary
    sample_properties = cr_webshim.get_sample_properties(args.aggregation_id, args.aggregation_desc, genomes, version=martian.get_pipelines_version(), agg_batches=agg_batches)

    sample_data_paths = cr_webshim_data.SampleDataPaths(
        summary_path=outs.summary,
        barcode_summary_path=args.barcode_summary_h5,
        analysis_path=args.analysis,
    )

    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)
    cr_webshim.build_web_summary_html(outs.web_summary, sample_properties, sample_data, PIPELINE_AGGR)

def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]
    cr_utils.copy(chunk_out.summary, outs.summary)
    cr_utils.copy(chunk_out.web_summary, outs.web_summary)
