#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import itertools
import math
import numpy as np
import tenkit.bam as tk_bam
import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils

""" Report info on each detected molecule
"""

__MRO__ = '''
stage REPORT_MOLECULES(
    in  bam[]  inputs,
    in  path   reference_path,
    in  map    align,
    in  string barcode_whitelist,
    in  json   extract_reads_summary,
    in  json   attach_bcs_and_umis_summary,
    in  json   mark_duplicates_summary,
    in  csv    filtered_barcodes,
    in  int[]  gem_groups,
    in  int    recovered_cells,
    in  int    force_cells,
    out h5     output,
    src py     "stages/counter/report_molecules",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
'''

MAX_MEM_GB = 64

def split(args):
    chunk_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist)
    whitelist_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups, use_min=False)

    # Estimate the total number of rows in the final molecule info. Worst case.
    total_reads = cr_utils.get_metric_from_json(args.extract_reads_summary, 'total_reads')
    mol_info_rows = total_reads

    # Memory for sorting in MoleculeCounter.concatenate_sort:
    # N = total number of rows
    # 8*N bytes to store the sort indices
    # (8+8+8)*N bytes to load, concatenate, and index into a 64-bit data column
    mol_info_mem_gb = int(math.ceil((32 * mol_info_rows)/1e9))
    join_mem_gb = min(MAX_MEM_GB, max(cr_constants.MIN_MEM_GB, whitelist_mem_gb + mol_info_mem_gb))

    chunks = []
    for chunk_input in args.inputs:
        chunks.append({
            'chunk_input': chunk_input,
            '__mem_gb': chunk_mem_gb,
        })
    join = {
        '__mem_gb': join_mem_gb,
    }
    return {'chunks': chunks, 'join': join}

def main(args, outs):
    outs.coerce_strings()

    in_bam = tk_bam.create_bam_infile(args.chunk_input)

    counter = cr_mol_counter.MoleculeCounter.open(outs.output, mode='w')

    mol_data_keys = cr_mol_counter.MoleculeCounter.get_data_columns()
    mol_data_columns = {key:idx for idx, key in enumerate(mol_data_keys)}

    gene_index = cr_reference.GeneIndex.load_pickle(cr_utils.get_reference_genes_index(args.reference_path))
    genomes = cr_utils.get_reference_genomes(args.reference_path)
    genome_index = cr_reference.get_genome_index(genomes)
    none_gene_id = len(gene_index.get_genes())

    # store reference index columns
    # NOTE - these must be cast to str first, as unicode is not supported
    counter.set_ref_column('genome_ids', [str(genome) for genome in genomes])
    counter.set_ref_column('gene_ids', [str(gene.id) for gene in gene_index.genes])
    counter.set_ref_column('gene_names', [str(gene.name) for gene in gene_index.genes])

    filtered_bcs_per_genome = cr_utils.load_barcode_csv(args.filtered_barcodes)
    filtered_bcs = set()
    for _, bcs in filtered_bcs_per_genome.iteritems():
        filtered_bcs |= set(bcs)

    gg_metrics = collections.defaultdict(lambda: {cr_mol_counter.GG_CONF_MAPPED_FILTERED_BC_READS_METRIC: 0})

    for (gem_group, barcode, gene_ids), reads_iter in itertools.groupby(in_bam, key=cr_utils.barcode_sort_key):
        if barcode is None or gem_group is None:
            continue
        is_cell_barcode = cr_utils.format_barcode_seq(barcode, gem_group) in filtered_bcs
        molecules = collections.defaultdict(lambda: np.zeros(len(mol_data_columns), dtype=np.uint64))

        compressed_barcode = cr_mol_counter.MoleculeCounter.compress_barcode_seq(barcode)
        gem_group = cr_mol_counter.MoleculeCounter.compress_gem_group(gem_group)

        read_positions = collections.defaultdict(set)
        for read in reads_iter:
            umi = cr_utils.get_read_umi(read)
            # ignore read2 to avoid double-counting. the mapping + annotation should be equivalent.
            if read.is_secondary or umi is None or read.is_read2:
                continue

            raw_umi = cr_utils.get_read_raw_umi(read)
            raw_bc, raw_gg = cr_utils.split_barcode_seq(cr_utils.get_read_raw_barcode(read))
            proc_bc, proc_gg = cr_utils.split_barcode_seq(cr_utils.get_read_barcode(read))

            if cr_utils.is_read_conf_mapped_to_transcriptome(read, cr_utils.get_high_conf_mapq(args.align)):
                assert len(gene_ids) == 1

                mol_key, map_type = (umi, gene_index.gene_id_to_int(gene_ids[0])), 'reads'

                read_pos = (read.tid, read.pos)
                uniq_read_pos = read_pos not in read_positions[mol_key]
                read_positions[mol_key].add(read_pos)

                if is_cell_barcode: gg_metrics[int(gem_group)][cr_mol_counter.GG_CONF_MAPPED_FILTERED_BC_READS_METRIC] += 1

            elif read.is_unmapped:
                mol_key, map_type, uniq_read_pos = (umi, none_gene_id), 'unmapped_reads', False
            else:
                mol_key, map_type, uniq_read_pos = (umi, none_gene_id), 'nonconf_mapped_reads', False
            molecules[mol_key][mol_data_columns[map_type]] += 1
            molecules[mol_key][mol_data_columns['umi_corrected_reads']] += int(not raw_umi == umi)
            molecules[mol_key][mol_data_columns['barcode_corrected_reads']] += int(not raw_bc == proc_bc)
            molecules[mol_key][mol_data_columns['conf_mapped_uniq_read_pos']] += int(uniq_read_pos)

        for mol_key, molecule in sorted(molecules.items()):
            umi, gene_id = mol_key
            genome = cr_utils.get_genome_from_str(gene_index.int_to_gene_id(gene_id), genomes)
            genome_id = cr_reference.get_genome_id(genome, genome_index)
            counter.add(barcode = compressed_barcode,
                        gem_group = gem_group,
                        umi = cr_mol_counter.MoleculeCounter.compress_umi_seq(umi),
                        gene = gene_id,
                        genome = genome_id,
                        **{key:molecule[col_idx] for key, col_idx in mol_data_columns.iteritems()})

    in_bam.close()

    counter.set_metric(cr_mol_counter.GEM_GROUPS_METRIC, dict(gg_metrics))

    counter.save()

def join(args, outs, chunk_defs, chunk_outs):
    summary = cr_utils.merge_jsons_as_dict([
        args.extract_reads_summary,
        args.attach_bcs_and_umis_summary,
        args.mark_duplicates_summary,
    ])
    gem_groups = sorted(set(args.gem_groups))
    metrics = cr_mol_counter.MoleculeCounter.get_metrics_from_summary(summary, gem_groups, args.recovered_cells, args.force_cells)
    input_h5_filenames = [chunk_out.output for chunk_out in chunk_outs]
    # update with metrics that were computed in the chunks
    chunk_metric = cr_mol_counter.GG_CONF_MAPPED_FILTERED_BC_READS_METRIC
    for gg, count in cr_mol_counter.MoleculeCounter.sum_gem_group_metric(input_h5_filenames, chunk_metric).iteritems():
        metrics[cr_mol_counter.GEM_GROUPS_METRIC][gg][chunk_metric] = count
    # make sure to sort globally by gem group. since the input is a barcode-sorted BAM, we assume it's already sorted by barcode.
    sort_columns = ['gem_group']
    cr_mol_counter.MoleculeCounter.concatenate_sort(outs.output, input_h5_filenames, sort_columns, metrics=metrics)
