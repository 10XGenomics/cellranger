#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import collections
import json
import numpy as np
import tables

import tenkit.safe_json as tk_json

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import cellranger.molecule_counter as cr_mol_counter
import cellranger.reference as cr_reference

__MRO__ = """
stage COUNT_GENES_FROM_MOLECULES(
    in  map  gem_group_index,
    in  h5   raw_molecules,
    in  csv  filtered_barcodes,
    out h5   filtered_molecules,
    out h5   raw_matrices_h5,
    out h5   filtered_matrices_h5,
    out h5   barcode_summary_h5,
    out path raw_matrices_mex,
    out path filtered_matrices_mex,
    out json summary,
    src py   "stages/aggregator/count_genes",
) split using (
    in  int  chunk_start,
    in  int  chunk_len,
)
"""

def split(args):
    chunks = []
    with cr_mol_counter.MoleculeCounter.open(args.raw_molecules, 'r') as mol_counter:
        barcode_whitelist = mol_counter.get_barcode_whitelist()
        matrix_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(barcode_whitelist, use_min=False)
        for chunk_start, chunk_len in mol_counter.get_chunks(cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK):
            mol_mem_gb = cr_mol_counter.MoleculeCounter.estimate_mem_gb(chunk_len)
            tot_mem_gb = 4 * (matrix_mem_gb + mol_mem_gb)
            threads = max(1, round(tot_mem_gb / 8.0))
            chunks.append({
                'chunk_start': str(chunk_start),
                'chunk_len': str(chunk_len),
                '__mem_gb': tot_mem_gb,
                '__threads': threads,
            })
    join_mem_gb = max(matrix_mem_gb, cr_constants.MIN_MEM_GB)
    join = {'__mem_gb': join_mem_gb}
    return {'chunks': chunks, 'join': join}


def main(args, outs):
    molecule_counter = cr_mol_counter.MoleculeCounter.open(args.raw_molecules, 'r', start=int(args.chunk_start), length=int(args.chunk_len))

    filtered_bcs_per_genome = cr_utils.load_barcode_csv(args.filtered_barcodes)

    raw_matrices = cr_matrix.GeneBCMatrices.build_from_mol_counter(molecule_counter)
    raw_matrices.save_h5(outs.raw_matrices_h5)
    raw_matrices.save_mex(outs.raw_matrices_mex)
    raw_matrices.save_barcode_summary_h5(outs.barcode_summary_h5)

    filtered_matrices = raw_matrices.filter_barcodes(filtered_bcs_per_genome)
    filtered_matrices.save_h5(outs.filtered_matrices_h5)
    filtered_matrices.save_mex(outs.filtered_matrices_mex)

    genome_ids = molecule_counter.get_ref_column('genome_ids')

    with cr_mol_counter.MoleculeCounter.open(outs.filtered_molecules, 'w') as ctr_out:
        summary = write_filtered_molecules(molecule_counter, ctr_out, genome_ids, filtered_bcs_per_genome)

    with open(outs.summary, 'w') as f:
        tk_json.dump_numpy(summary, f, pretty=True)

def write_filtered_molecules(ctr_in, ctr_out, genomes, bcs_per_genome):
    ctr_out.set_all_metrics(ctr_in.get_all_metrics())

    filtered_bc_tuples = set()
    genome_ids = ctr_in.get_column('genome')
    genome_index = cr_reference.get_genome_index(genomes)
    for (genome, formatted_bcs) in bcs_per_genome.iteritems():
        genome_id = cr_reference.get_genome_id(genome, genome_index)
        for formatted_bc in formatted_bcs:
            (bc, gg) = cr_utils.split_barcode_seq(formatted_bc)
            cbc = cr_mol_counter.MoleculeCounter.compress_barcode_seq(bc)
            filtered_bc_tuples.add((genome_id, gg, cbc))

    def keep_molecule(genome_id, gem_group, barcode):
        tup = (genome_id, gem_group, barcode)
        return (tup in filtered_bc_tuples)

    filter_func = np.vectorize(keep_molecule)

    gem_groups = ctr_in.get_column('gem_group')
    barcodes = ctr_in.get_column('barcode')
    filter_index = filter_func(genome_ids, gem_groups, barcodes)

    for col in cr_mol_counter.MOLECULE_INFO_COLUMNS:
        data = ctr_in.get_column(col)
        filtered_data = data[filter_index]
        ctr_out.add_many(col, filtered_data)

    for col in cr_mol_counter.MOLECULE_REF_COLUMNS:
        ctr_out.set_ref_column(col, ctr_in.get_ref_column(col))

    # summarize filtered data
    genomes = ctr_out.get_ref_column('genome_ids')
    filtered_reads = ctr_out.get_column('reads')
    flt_conf_mapped_per_genome = {}
    if len(genomes) == 1:
        genome = genomes[0]
        flt_conf_mapped_per_genome[genome] = filtered_reads.sum()
    else:
        genome_ids = ctr_out.get_column('genome')
        genome_index = cr_reference.get_genome_index(genomes)
        for genome in genomes:
            genome_id = cr_reference.get_genome_id(genome, genome_index)
            flt_conf_mapped_per_genome[genome] = filtered_reads[genome_ids == genome_id].sum()
    summary = {'flt_conf_mapped_per_genome': flt_conf_mapped_per_genome}
    return summary

def join(args, outs, chunk_defs, chunk_outs):
    matrix_attrs = cr_matrix.make_matrix_attrs_aggr(args.gem_group_index, "Unknown")
    cr_matrix.concatenate_h5([chunk_out.raw_matrices_h5 for chunk_out in chunk_outs], outs.raw_matrices_h5, extra_attrs=matrix_attrs)
    cr_matrix.concatenate_h5([chunk_out.filtered_matrices_h5 for chunk_out in chunk_outs], outs.filtered_matrices_h5, extra_attrs=matrix_attrs)

    cr_matrix.concatenate_mex_dirs([chunk_out.raw_matrices_mex for chunk_out in chunk_outs], outs.raw_matrices_mex)
    cr_matrix.concatenate_mex_dirs([chunk_out.filtered_matrices_mex for chunk_out in chunk_outs], outs.filtered_matrices_mex)

    merged_molecules = [chunk_out.filtered_molecules for chunk_out in chunk_outs]
    cr_mol_counter.MoleculeCounter.concatenate(outs.filtered_molecules, merged_molecules)

    barcode_summaries = [chunk_out.barcode_summary_h5 for chunk_out in chunk_outs]
    merge_barcode_summaries(barcode_summaries, outs.barcode_summary_h5)

    # merge summaries
    summary = merge_summaries(chunk_outs)
    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)

def merge_summaries(chunk_outs):
    flt_conf_mapped_per_genome = collections.Counter()
    for chunk_out in chunk_outs:
        with open(chunk_out.summary, 'r') as f:
            chunk_summary = json.load(f)
            flt_conf_mapped_per_genome += collections.Counter(chunk_summary['flt_conf_mapped_per_genome'])
    return {'flt_conf_mapped_per_genome': dict(flt_conf_mapped_per_genome)}

def merge_barcode_summaries(input_files, output_file):
    # each chunk produces a barcode summary containing ALL barcodes from ALL gem groups, not just the ones being counted
    # in that chunk. so the datasets need to be squashed rather than concatenated
    with tables.open_file(output_file, mode = 'w') as fout:
        dsets = {}
        # init datasets using the first input
        if len(input_files) > 0:
            with tables.open_file(input_files[0], mode = 'r') as fin:
                for node in fin.walk_nodes('/', 'Array'):
                    dsets[node.name] = fout.create_carray('/', node.name, obj=node[:])
        # add data from the other inputs
        for input_file in input_files[1:]:
            with tables.open_file(input_file, mode = 'r') as fin:
                for (name, carray) in dsets.iteritems():
                    if name == cr_constants.H5_BC_SEQUENCE_COL:
                        continue # don't modify the barcode sequences
                    carray[:] += fin.get_node('/', name)[:]
