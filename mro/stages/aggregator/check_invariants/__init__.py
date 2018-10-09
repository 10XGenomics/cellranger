#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import json
import martian
import numpy as np

import cellranger.h5_constants as h5_constants
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
from cellranger.molecule_counter import MoleculeCounter
import tenkit.safe_json as tk_safe_json

__MRO__ = """
stage CHECK_INVARIANTS(
    in  map[] input_sample_defs,
    in  h5    merged_raw_gene_bc_matrices_h5,
    out json  summary,
    src py    "stages/aggregator/check_invariants",
) split (
)
"""

def split(args):
    matrix_mem_gb = int(1.1 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.merged_raw_gene_bc_matrices_h5))
    join_args = {
        '__mem_gb': max(matrix_mem_gb, h5_constants.MIN_MEM_GB),
    }
    return {'chunks': [], 'join': join_args}

def join(args, outs, chunk_defs, chunk_outs):
    # compute invariants on input data
    input_genomes = set()
    input_features = set()
    input_bc_counts = {}
    input_feature_counts = {}
    input_num_gem_groups = 0

    for sample_def in args.input_sample_defs:
        library_id = sample_def['library_id']
        with MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], 'r') as mc:
            input_genomes.update(mol_counter_genomes(mc))
            input_features.update(mol_counter_features_id_type(mc))
            gem_groups = mc.get_gem_groups()
            input_num_gem_groups += len(gem_groups)

            mol_gem_group = mc.get_column('gem_group')

            mol_barcode_idx = mc.get_column('barcode_idx')
            for gg in gem_groups:
                input_bc_counts[(library_id, gg)] = np.zeros(len(mc.get_ref_column('barcodes')))
                bc_idx, counts = np.unique(mol_barcode_idx[mol_gem_group == gg], return_counts=True)
                input_bc_counts[(library_id, gg)][bc_idx] = counts
            del mol_barcode_idx

            mol_feature_idx = mc.get_column('feature_idx')
            for gg in gem_groups:
                input_feature_counts[(library_id, gg)] = np.zeros(len(mc.feature_reference.feature_defs))
                feature_idx, counts = np.unique(mol_feature_idx[mol_gem_group == gg], return_counts=True)
                input_feature_counts[(library_id, gg)][feature_idx] = counts
            del mol_feature_idx

    # compute invariants on output
    output_matrix = cr_matrix.CountMatrix.load_h5_file(args.merged_raw_gene_bc_matrices_h5)
    output_genomes = set(output_matrix.get_genomes())
    output_features = set(count_matrix_features_id_type(output_matrix))
    output_bc_counts = {}
    output_feature_counts = {}
    output_gem_index = cr_matrix.get_gem_group_index(args.merged_raw_gene_bc_matrices_h5)
    output_num_gem_groups = len(output_gem_index)

    for gg in output_gem_index:
        library_id, old_gg = output_gem_index[gg]
        matrix_gg = output_matrix.select_barcodes_by_gem_group(gg)
        output_bc_counts[(library_id, old_gg)] = matrix_gg.get_counts_per_bc()
        output_feature_counts[(library_id, old_gg)] = matrix_gg.get_counts_per_feature()

    exit_message = ('An internal problem in the aggr pipeline has been detected '
                     'that might lead to incorrect results. Please report this '
                     'problem to support@10xgenomics.com.')

    if input_genomes != output_genomes:
        martian.log_info('Genomes differ between input molecule files and aggregated matrix')
        martian.exit(exit_message)
    if input_features != output_features:
        martian.log_info('Features differ between input molecule files and aggregated matrix')
        martian.exit(exit_message)
    if input_num_gem_groups != output_num_gem_groups:
        martian.log_info('Number of GEM groups differs between input molecule files and aggregated matrix')
        martian.exit(exit_message)
    for lib_gg in input_bc_counts.keys():
        if np.any(input_bc_counts[lib_gg] < output_bc_counts[lib_gg]):
            martian.log_info('Barcode(s) in library {}, GEM group {} have higher UMI counts '
                             'in aggregated output compared to inputs'
                             .format(lib_gg[0], lib_gg[1]))
            martian.exit(exit_message)
        if np.any(input_feature_counts[lib_gg] < output_feature_counts[lib_gg]):
            martian.log_info('Feature(s) in library {}, GEM group {} have higher UMI counts '
                             'in aggregated output compared to inputs'
                             .format(lib_gg[0], lib_gg[1]))
            martian.exit(exit_message)

    summary = {
        'genomes_present': list(input_genomes),
        'num_features_in_ref': len(input_features),
        'num_gem_groups': input_num_gem_groups,
    }

    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)

def mol_counter_genomes(mol_counter):
    return cr_matrix.CountMatrix._get_genomes_from_feature_ref(mol_counter.feature_reference)

def mol_counter_features_id_type(mol_counter):
    return ((f.id, f.feature_type) for f in mol_counter.feature_reference.feature_defs)

def count_matrix_features_id_type(count_matrix):
    return ((f.id, f.feature_type) for f in count_matrix.feature_ref.feature_defs)
