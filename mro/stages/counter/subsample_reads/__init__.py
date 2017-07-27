#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import numpy as np
import martian
import os
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
from cellranger.molecule_counter import MoleculeCounter
import cellranger.report as cr_report
import cellranger.utils as cr_utils

__MRO__ = """
stage SUBSAMPLE_READS(
    in  h5     molecule_info,
    in  csv    filtered_barcodes,
    out map[]  subsampled_matrices,
    out pickle chunked_reporter,
    out json   summary,
    src py     "stages/counter/subsample_reads",
) split using (
    in  int    chunk_start,
    in  int    chunk_len,
    in  map    subsample_info,
)
"""

# Limit memory usage in the split function
SPLIT_MEM_GB = 5

def split(args):
    # Get the cell count
    filtered_bcs_per_genome = cr_utils.load_barcode_csv(args.filtered_barcodes)
    filtered_bcs = set()
    for _, bcs in filtered_bcs_per_genome.iteritems():
        filtered_bcs |= set(bcs)
    n_cells = len(filtered_bcs)

    if n_cells == 0:
        return {'chunks': [{'chunk_start': 0, 'chunk_len': 0, 'subsample_info': {}}]}

    # Get required info from the mol info
    with MoleculeCounter.open(args.molecule_info, 'r') as mol_counter:
        n_molecule_info_entries = mol_counter.nrows()
        barcode_whitelist = mol_counter.get_barcode_whitelist()
        gem_groups = mol_counter.get_gem_groups()

        raw_reads = mol_counter.get_total_raw_reads()
        raw_rpc = tk_stats.robust_divide(raw_reads, n_cells)
        mapped_reads = mol_counter.get_total_conf_mapped_filtered_bc_reads()

    mapped_read_frac = tk_stats.robust_divide(mapped_reads, raw_reads)

    subsamplings = list() # track subsample info definitions

    # Calculate extra deciles to add in based on raw reads
    if raw_reads > 0:
        subsampling_deciles = [round(decile * raw_rpc) for decile in np.arange(0.1, 1.1, 0.1)]
    else:
        subsampling_deciles = []

    # All target depths
    target_rpcs = cr_constants.SUBSAMPLE_READS_PER_CELL + subsampling_deciles

    for subsample_type, rpc_multiplier in [(cr_constants.RAW_SUBSAMPLE_TYPE, mapped_read_frac),
                                           (cr_constants.MAPPED_SUBSAMPLE_TYPE, 1.0)]:
        # Generate subsampling definitions
        for target_rpc in target_rpcs:
            target_mapped_reads = int(float(target_rpc) * float(n_cells) * rpc_multiplier)

            subsample_rate = tk_stats.robust_divide(target_mapped_reads, mapped_reads)

            if subsample_rate > 1.0:
                continue

            subsamplings.append({'subsample_type': subsample_type,
                                 'target_rpc': target_rpc,
                                 'subsample_rate': subsample_rate,
                                 'all_target_rpc': target_rpcs,
            })

    # Each chunk needs to store the entire gene-bc matrix and a piece of the mol info h5
    matrix_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(barcode_whitelist,
                                                                       gem_groups)
    chunk_len = cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK
    chunk_mem_gb = matrix_mem_gb + MoleculeCounter.estimate_mem_gb(chunk_len)
    join_mem_gb = matrix_mem_gb

    # Split the molecule info h5 into equi-RAM chunks
    chunks = []
    for subsample_info in subsamplings:
        for chunk_start in xrange(0, n_molecule_info_entries, chunk_len):
            chunks.append({
                'chunk_start': str(chunk_start),
                'chunk_len': str(min(n_molecule_info_entries-chunk_start, chunk_len)),
                'subsample_info': subsample_info,
                '__mem_gb': chunk_mem_gb,
            })
    join = {
        '__mem_gb': join_mem_gb,
    }

    if len(chunks) == 0:
        chunks.append({
            'chunk_start': str(0),
            'chunk_len': str(0),
            'subsample_info': {},
        })

    return {'chunks': chunks, 'join': join}

def main(args, outs):
    np.random.seed(0)

    subsample_rate = args.subsample_info.get('subsample_rate')
    if subsample_rate is None:
        return

    mol_counter = MoleculeCounter.open(args.molecule_info, 'r',
                                       start=int(args.chunk_start),
                                       length=int(args.chunk_len))

    # Subsample the matrices
    subsample_result = {}
    subsampled_raw_mats = cr_matrix.GeneBCMatrices.build_from_mol_counter(mol_counter,
                                                                          subsample_rate=subsample_rate,
                                                                          subsample_result=subsample_result)

    # Filter the subsampled matrices
    filtered_bcs_per_genome = cr_utils.load_barcode_csv(args.filtered_barcodes)
    subsampled_filt_mats = subsampled_raw_mats.filter_barcodes(filtered_bcs_per_genome)

    # Calculations for subsampled duplication rate
    reporter = cr_report.Reporter(genomes = map(str, mol_counter.get_ref_column('genome_ids')),
                                  subsample_types = cr_constants.ALL_SUBSAMPLE_TYPES,
                                  subsample_depths = args.subsample_info['all_target_rpc'])

    reporter.subsampled_duplication_frac_cb(subsampled_raw_mats,
                                            mol_counter,
                                            args.subsample_info['subsample_rate'],
                                            args.subsample_info['subsample_type'],
                                            args.subsample_info['target_rpc'],
                                            subsample_result['mapped_reads'],
    )

    mol_counter.close()

    reporter.save(outs.chunked_reporter)

    outs.subsampled_matrices = {}
    outs.subsampled_matrices['raw_matrices'] = martian.make_path('raw_matrices.h5')
    outs.subsampled_matrices['filtered_matrices'] = martian.make_path('filtered_matrices.h5')

    subsampled_raw_mats.save_h5(outs.subsampled_matrices['raw_matrices'])
    subsampled_filt_mats.save_h5(outs.subsampled_matrices['filtered_matrices'])


def join(args, outs, chunk_defs, chunk_outs):
    # Summarize genes and UMI counts
    chunks = zip(chunk_defs, chunk_outs)

    # Check for an empty chunk
    if len(chunks) == 0 or chunk_defs[0].subsample_info.get('subsample_type') is None or chunk_defs[0].subsample_info.get('subsample_rate') is None:
        outs.summary = None
        return

    chunk_key = lambda chunk: (chunk[0].subsample_info['subsample_type'],
                               chunk[0].subsample_info['target_rpc'],
                               chunk[0].subsample_info['subsample_rate'])

    # Merge reporter objects from main
    reporter_file_names = [chunk_out.chunked_reporter for chunk_out in chunk_outs if os.path.isfile(chunk_out.chunked_reporter)]
    merged_reporter = cr_report.merge_reporters(reporter_file_names)

    outs.subsampled_matrices = []

    # Aggregate the molecule info chunks that belong together
    for chunk_group, (subsample_key, chunk_iter) in enumerate(itertools.groupby(sorted(chunks, key=chunk_key), chunk_key)):
        subsample_type, target_rpc, subsample_rate = subsample_key

        if subsample_type is None or subsample_rate is None:
            continue

        # Aggregate information over chunks with same key
        chunk_raw_h5s = []
        chunk_filtered_h5s = []
        all_subsample_types = cr_constants.ALL_SUBSAMPLE_TYPES
        all_target_rpc = None

        for chunk_def, chunk_out in chunk_iter:
            # List of target rpcs should be identical among all chunks
            assert all_target_rpc is None or all_target_rpc == chunk_def.subsample_info['all_target_rpc']
            all_target_rpc = chunk_def.subsample_info['all_target_rpc']

            chunk_raw_h5s.append(chunk_out.subsampled_matrices['raw_matrices'])
            chunk_filtered_h5s.append(chunk_out.subsampled_matrices['filtered_matrices'])

        raw_matrices = cr_matrix.merge_matrices(chunk_raw_h5s)
        filtered_matrices = cr_matrix.merge_matrices(chunk_filtered_h5s)

        # Compute metrics on subsampled matrices
        merged_reporter.summarize_subsampled_matrices_cb(filtered_matrices, subsample_type, target_rpc)

        # Write the merged matrices
        outs.subsampled_matrices.append({
            'subsample_type': subsample_type,
            'target_rpc': target_rpc,
            'subsample_rate': subsample_rate,
            'all_subsample_types': all_subsample_types,
            'all_target_rpc': all_target_rpc,
            'raw_matrices': martian.make_path('%s_%s_%s_raw_matrices.h5' % (subsample_type,
                                                                            target_rpc,
                                                                            chunk_group)),
            'filtered_matrices': martian.make_path('%s_%s_%s_filtered_matrices.h5' % (subsample_type,
                                                                                      target_rpc,
                                                                                      chunk_group)),
            })

        assert not os.path.exists(outs.subsampled_matrices[-1]['raw_matrices'])
        assert not os.path.exists(outs.subsampled_matrices[-1]['filtered_matrices'])

        raw_matrices.save_h5(outs.subsampled_matrices[-1]['raw_matrices'])
        filtered_matrices.save_h5(outs.subsampled_matrices[-1]['filtered_matrices'])

    merged_reporter.report_summary_json(filename=outs.summary)
