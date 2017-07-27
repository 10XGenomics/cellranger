#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import collections
import json
import numpy as np

import tenkit.safe_json as tk_json
import tenkit.stats as tk_stats

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils

__MRO__  = """
stage NORMALIZE_DEPTH(
    in  h5     molecules,
    in  map    detect_cells_gg_metrics,
    in  string normalization_mode,
    out h5     out_molecules,
    out json   summary,
    src py     "stages/aggregator/normalize_depth",
) split using (
    in  bool   downsample,
    in  int    chunk_start,
    in  int    chunk_len,
    in  map    downsample_map,
)
"""

def split(args):
    # default to downsampling by mapped reads
    downsample = True
    use_raw_reads = False

    if args.normalization_mode == cr_constants.NORM_MODE_RAW:
        use_raw_reads = True
    elif args.normalization_mode == cr_constants.NORM_MODE_NONE:
        downsample = False

    # compute downsample rates for each gem group
    downsample_map = args.detect_cells_gg_metrics
    with cr_mol_counter.MoleculeCounter.open(args.molecules, 'r') as mol_counter:
        for (gg, submetrics) in mol_counter.get_metric(cr_mol_counter.GEM_GROUPS_METRIC).iteritems():
            info = downsample_map[str(gg)]
            info['total_reads'] = submetrics[cr_mol_counter.GG_TOTAL_READS_METRIC]
            reads = info['total_reads'] if use_raw_reads else info['cmb_reads']
            cells = info['cells']
            info['rpc'] = tk_stats.robust_divide(reads, cells) if cells > 0 else 0.0

    lowest_rpc = min([gg['rpc'] for gg in downsample_map.values()])
    for gg, info in downsample_map.iteritems():
        if downsample and len(downsample_map) > 1:
            if lowest_rpc == 0:
                # one or more samples are empty. just do the naive thing for now.
                frac_reads_kept = 0.0
            else:
                frac_reads_kept = tk_stats.robust_divide(lowest_rpc, info['rpc'])
        else:
            frac_reads_kept = 1.0
        info['frac_reads_kept'] = frac_reads_kept

    # Split the molecule info h5 into equi-RAM chunks, preserving (barcode, gem_group) boundaries
    # Assumes the molecule_info is sorted by (barcode, gem_group)
    chunks = []
    with cr_mol_counter.MoleculeCounter.open(args.molecules, 'r') as mol_counter:
        for chunk_start, chunk_len in mol_counter.get_chunks(cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK, preserve_boundaries=False):
            chunks.append({
                'downsample': downsample,
                'downsample_map': downsample_map,
                'chunk_start': str(chunk_start),
                'chunk_len': str(chunk_len),
                '__mem_gb': cr_mol_counter.MoleculeCounter.estimate_mem_gb(chunk_len),
            })
    return {'chunks': chunks}

def main(args, outs):
    np.random.seed(0)

    with cr_mol_counter.MoleculeCounter.open(args.molecules, 'r', start=int(args.chunk_start), length=int(args.chunk_len)) as ctr_in:
        with cr_mol_counter.MoleculeCounter.open(outs.out_molecules, 'w') as ctr_out:
            metrics_in = ctr_in.get_all_metrics()
            metrics_out = metrics_in.copy()

    	    reads = ctr_in.get_column('reads')
      	    gem_groups = ctr_in.get_column('gem_group')

            if args.downsample and len(args.downsample_map) > 1:
                downsample_func = np.vectorize(lambda gem_group, read_count: np.random.binomial(read_count, args.downsample_map[str(gem_group)]['frac_reads_kept']))

                # downsample metrics
                for gg in metrics_out[cr_mol_counter.GEM_GROUPS_METRIC]:
                    frac_reads_kept = args.downsample_map[str(gg)]['frac_reads_kept']
                    total_reads_in = metrics_in[cr_mol_counter.GEM_GROUPS_METRIC][gg][cr_mol_counter.GG_TOTAL_READS_METRIC]
                    total_reads_out = round(frac_reads_kept * total_reads_in)
                    metrics_out[cr_mol_counter.GEM_GROUPS_METRIC][gg][cr_mol_counter.GG_DOWNSAMPLED_READS_METRIC] = total_reads_out

                ctr_out.set_all_metrics(metrics_out)

                # downsample molecule info
                subsampled_reads = downsample_func(gem_groups, reads)
                for col in cr_mol_counter.MOLECULE_INFO_COLUMNS:
                    if col == 'reads':
                        data = subsampled_reads
                    else:
                        data = ctr_in.get_column(col)
                    ctr_out.add_many(col, data)

                # pass reference info
                for col in cr_mol_counter.MOLECULE_REF_COLUMNS:
                    ctr_out.set_ref_column(col, ctr_in.get_ref_column(col))

            else:
                subsampled_reads = reads

            # collect summary stats
            genomes = ctr_in.get_ref_column('genome_ids')
            raw_conf_mapped_per_genome = {}
            if len(genomes) == 1:
                genome = genomes[0]
                raw_conf_mapped_per_genome[genome] = subsampled_reads.sum()
            else:
                genome_ids = ctr_in.get_column('genome')
                genome_index = cr_reference.get_genome_index(genomes)
                for genome in genomes:
                    genome_id = cr_reference.get_genome_id(genome, genome_index)
                    raw_conf_mapped_per_genome[genome] = subsampled_reads[genome_ids == genome_id].sum()

            summary = {'raw_conf_mapped_per_genome': raw_conf_mapped_per_genome, 'mol_counter_metrics': metrics_out}

            with open(outs.summary, 'w') as f:
                tk_json.dump_numpy(summary, f, pretty=True)

def join(args, outs, chunk_defs, chunk_outs):
    downsample = chunk_defs[0].downsample
    downsample_map = chunk_defs[0].downsample_map
    if downsample and len(downsample_map) > 1:
        input_h5_filenames = [chunk_out.out_molecules for chunk_out in chunk_outs]
        cr_mol_counter.MoleculeCounter.concatenate(outs.out_molecules, input_h5_filenames)
    else:
        # just copy input molecules
        cr_utils.copy(args.molecules, outs.out_molecules)

    # merge summaries
    summary = merge_summaries(chunk_outs)
    summary['downsample_info'] = downsample_map
    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)

def merge_summaries(chunk_outs):
    raw_conf_mapped_per_genome = collections.Counter()
    mol_counter_metrics = None
    for chunk_out in chunk_outs:
        with open(chunk_out.summary, 'r') as f:
            chunk_summary = json.load(f)
            raw_conf_mapped_per_genome += collections.Counter(chunk_summary['raw_conf_mapped_per_genome'])
            if mol_counter_metrics is None:
                mol_counter_metrics = chunk_summary['mol_counter_metrics']
    return {'raw_conf_mapped_per_genome': dict(raw_conf_mapped_per_genome), 'mol_counter_metrics': mol_counter_metrics}
