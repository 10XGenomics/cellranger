#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import numpy as np
import random
import csv

import cellranger.utils as cr_utils
import cellranger.stats as cr_stats
import cellranger.molecule_counter as cr_mol_counter
import cellranger.reference as cr_reference
import cellranger.constants as cr_constants

__MRO__ = """
stage DETECT_CELLS(
    in  h5[]  molecule_chunks,
    out csv   cell_barcodes,
    out map   gem_group_metrics,
    src py    "stages/aggregator/detect_cells",
) split using (
    in  h5    molecule_h5,
    in  int   chunk_start,
    in  int   chunk_len,
    in  int   gem_group,
    in  int   recovered_cells,
    in  int   force_cells,
)
"""

def split(args):
    chunks = []
    for molecule_h5 in args.molecule_chunks:
        with cr_mol_counter.MoleculeCounter.open(molecule_h5, 'r') as mol_counter:
            barcode_whitelist = mol_counter.get_barcode_whitelist()
            whitelist_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(barcode_whitelist, use_min=False)
            for gem_group, chunk_start, chunk_len in mol_counter.get_chunks_by_gem_group():
                mol_mem_gb = cr_mol_counter.MoleculeCounter.estimate_mem_gb(chunk_len)
                recovered_cells = mol_counter.get_metric(cr_mol_counter.GEM_GROUPS_METRIC)[gem_group].get(cr_mol_counter.GG_RECOVERED_CELLS_METRIC, None)
                force_cells = mol_counter.get_metric(cr_mol_counter.GEM_GROUPS_METRIC)[gem_group].get(cr_mol_counter.GG_FORCE_CELLS_METRIC, None)
                chunks.append({
                    'molecule_h5': molecule_h5,
                    'gem_group': str(gem_group),
                    'recovered_cells': recovered_cells,
                    'force_cells': force_cells,
                    'chunk_start': chunk_start,
                    'chunk_len': chunk_len,
                    '__mem_gb': whitelist_mem_gb + mol_mem_gb,
                })

    return {'chunks': chunks}

def main(args, outs):
    random.seed(0)
    np.random.seed(0)

    with cr_mol_counter.MoleculeCounter.open(args.molecule_h5, 'r', start=int(args.chunk_start), length=int(args.chunk_len)) as ctr_in:
        genome_ids = ctr_in.get_ref_column('genome_ids')
        gene_ids = ctr_in.get_ref_column('gene_ids')
        barcode_whitelist = cr_utils.load_barcode_whitelist(ctr_in.get_barcode_whitelist())

        # Estimate BC diversity and recovered cells per gem group
        gg_total_diversity = len(barcode_whitelist)

        bc_counts_per_genome = get_bc_counts(genome_ids, gene_ids, ctr_in)
        top_bcs_per_genome = {}
        total_conf_mapped_cell_reads = 0
        total_cells = 0
        recovered_cells = args.recovered_cells or cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP
        for genome, (barcodes, umi_counts, read_counts) in bc_counts_per_genome.iteritems():
            if args.force_cells is not None:
                top_bc_indices, filter_summary, _ = cr_stats.filter_cellular_barcodes_fixed_cutoff(umi_counts, args.force_cells)
            else:
                top_bc_indices, filter_summary, _ = cr_stats.filter_cellular_barcodes_ordmag(umi_counts, recovered_cells, gg_total_diversity)
            top_bcs_per_genome[genome] = barcodes[top_bc_indices]
            total_conf_mapped_cell_reads += read_counts[top_bc_indices].sum()
            total_cells += filter_summary['filtered_bcs']

        write_filtered_barcodes(outs.cell_barcodes, args.gem_group, ctr_in, top_bcs_per_genome)

        outs.gem_group_metrics = {'cells': int(total_cells), 'cmb_reads': int(total_conf_mapped_cell_reads)}

def join(args, outs, chunk_defs, chunk_outs):
    barcodes_csv = [chunk_out.cell_barcodes for chunk_out in chunk_outs]
    cr_utils.concatenate_files(outs.cell_barcodes, barcodes_csv)
    outs.gem_group_metrics = {cd.gem_group: co.gem_group_metrics for (cd, co) in zip(chunk_defs, chunk_outs)}

def get_bc_counts(genomes, genes, molecule_counter):
    genome_ids = molecule_counter.get_column('genome')
    genome_index = cr_reference.get_genome_index(genomes)
    conf_mapped_reads = molecule_counter.get_column('reads')
    barcodes = molecule_counter.get_column('barcode')

    bc_counts = {}
    for genome in genomes:
        genome_id = cr_reference.get_genome_id(genome, genome_index)
        genome_indices = genome_ids == genome_id
        if genome_indices.sum() == 0:
            # edge case - there's no data for this genome (e.g. empty sample, false barnyard sample, or nothing confidently mapped)
            continue
        bcs_for_genome = barcodes[genome_indices]
        # only count UMIs with at least one conf mapped read
        umi_conf_mapped_to_genome = conf_mapped_reads[genome_indices] > 0
        bc_breaks = bcs_for_genome[1:] - bcs_for_genome[:-1]
        bc_breaks = np.concatenate(([1], bc_breaks)) # first row is always a break
        bc_break_indices = np.nonzero(bc_breaks)[0]
        unique_bcs = bcs_for_genome[bc_break_indices]
        umis_per_bc = np.add.reduceat(umi_conf_mapped_to_genome, bc_break_indices)
        cmb_reads_per_bc = np.add.reduceat(conf_mapped_reads[genome_indices], bc_break_indices)
        bc_counts[genome] = (unique_bcs, umis_per_bc, cmb_reads_per_bc)

    return bc_counts

def write_filtered_barcodes(out_csv, gem_group, mol_counter, bcs_per_genome):
    with open(out_csv, 'wb') as f:
        writer = csv.writer(f)
        for (genome, bc_ids) in bcs_per_genome.iteritems():
            for bc_id in bc_ids:
                formatted_barcode = cr_utils.format_barcode_seq(
                    mol_counter.decompress_barcode_seq(bc_id),
                    gem_group)
                writer.writerow([genome, formatted_barcode])
