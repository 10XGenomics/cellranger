#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import numpy as np

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
from cellranger.molecule_counter import MoleculeCounter

__MRO__ = """
stage MERGE_MOLECULES(
    in  map[]  sample_defs,
    in  map[]  libraries,
    in  map    gem_group_index,
    out h5     merged_molecules,
    out h5[]   chunked_molecules,
    src py     "stages/aggregator/merge_molecules",
) split using (
    in  string aggr_id,
    in  h5     molecule_h5,
)
"""

def split(args):
    chunks = []
    for sample_def in args.sample_defs:
        with MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], 'r') as mol_counter:
            mem_gb = int(1.5 * MoleculeCounter.estimate_mem_gb(mol_counter.nrows()))
            chunks.append({
                'aggr_id': sample_def[cr_constants.AGG_ID_FIELD],
                'molecule_h5': sample_def[cr_constants.AGG_H5_FIELD],
                '__mem_gb': mem_gb,
            })
    return {'chunks': chunks, 'join': {'__mem_gb': 6}}

def main(args, outs):
    with MoleculeCounter.open(args.molecule_h5, 'r') as in_mc:
        # Get the gem group and library mappings
        gg_map, lib_idx_map = get_library_mapping(args.aggr_id, args.libraries)

        # FIXME: Handle hetergeneous barcode whitelists
        merged_barcodes = in_mc.get_barcodes()

        # FIXME: Handle heterogeneous feature references
        merged_feature_ref = in_mc.get_feature_ref()

        # Remap the barcode info
        old_barcode_info = in_mc.get_barcode_info()
        new_pass_filter = old_barcode_info.pass_filter
        new_pass_filter[:,1] = lib_idx_map[new_pass_filter[:,1]]

        new_barcode_info = cr_mol_counter.BarcodeInfo(
            pass_filter=new_pass_filter,
            genomes=old_barcode_info.genomes,
        )

        with MoleculeCounter.open(outs.molecule_h5, 'w',
                                  feature_ref=merged_feature_ref,
                                  barcodes=merged_barcodes,
                                  library_info=args.libraries,
                                  barcode_info=new_barcode_info,
        ) as out_mc:

            # Copy the datasets, rewriting the ones we remap
            for col, ds in in_mc.columns.iteritems():
                if col == 'gem_group':
                    old_gg = ds[:]
                    new_gg = gg_map[old_gg]
                    out_mc.append_column(col, new_gg)

                elif col == 'library_idx':
                    old_idx = ds[:]
                    new_idx = lib_idx_map[old_idx]
                    out_mc.append_column(col, new_idx)

                else:
                    out_mc.append_column(col, ds[:])

            # Copy over all standard metrics
            out_metrics = in_mc.get_all_metrics()

            # Remap the per-gem-group and per-library metrics
            old_gg_metrics = in_mc.get_metric(cr_mol_counter.GEM_GROUPS_METRIC)
            gg_metrics = {str(gg_map[int(og)]):m for og,m in old_gg_metrics.iteritems()}
            old_lib_metrics = in_mc.get_metric(cr_mol_counter.LIBRARIES_METRIC)
            lib_metrics = {str(lib_idx_map[int(ol)]):m for ol,m in old_lib_metrics.iteritems()}

            out_metrics[cr_mol_counter.GEM_GROUPS_METRIC] = gg_metrics
            out_metrics[cr_mol_counter.LIBRARIES_METRIC] = lib_metrics

            out_mc.set_all_metrics(out_metrics)

def get_library_mapping(aggr_id, libraries):
    """Get the mapping of gem groups and library indices to their new values.

    Args:
      aggr_id (str): The label given to a set of libraries in the aggr CSV file.
      libraries (list of dict): New library info.
    Returns:
      tuple of (gem_group_map, library_map) (np.array, np.array):
        gem_group_map maps the old gem group integer ro the new one
        library_map maps the old library index integer to the new one
    """
    for i, lib in enumerate(libraries):
        lib['index'] = i

    my_libs = [lib for lib in libraries if lib['aggr_id'] == aggr_id]
    max_old_gg = max(lib['gem_group'] for lib in my_libs)
    max_old_lib_idx = max(lib['old_library_index'] for lib in my_libs)

    gem_group_map = np.zeros(1+max_old_gg, dtype=MoleculeCounter.get_column_dtype('gem_group'))
    lib_idx_map = np.zeros(1+max_old_lib_idx, dtype=MoleculeCounter.get_column_dtype('library_idx'))
    for lib in my_libs:
        gem_group_map[lib['old_gem_group']] = lib['gem_group']
        lib_idx_map[lib['old_library_index']] = lib['index']

    return gem_group_map, lib_idx_map


def join(args, outs, chunk_defs, chunk_outs):
    molecules = [chunk_out.molecule_h5 for chunk_out in chunk_outs]
    metrics = MoleculeCounter.naive_concatenate_metrics(molecules)
    metrics[cr_mol_counter.IS_AGGREGATED_METRIC] = True
    MoleculeCounter.concatenate(outs.merged_molecules, molecules, metrics=metrics)
