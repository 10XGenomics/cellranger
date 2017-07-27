#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import numpy as np

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter

__MRO__ = """
stage MERGE_MOLECULES(
    in  map[]  sample_defs,
    in  map    gem_group_index,
    out h5     merged_molecules,
    out h5[]   chunked_molecules,
    src py     "stages/aggregator/merge_molecules",
) split using (
    in  string library_id,
    in  h5     molecule_h5,
)
"""

def split(args):
    chunks = []
    for sample_def in args.sample_defs:
        with cr_mol_counter.MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], 'r') as mol_counter:
            mem_gb = int(1.5 * cr_mol_counter.MoleculeCounter.estimate_mem_gb(mol_counter.nrows()))
            chunks.append({
                'library_id': sample_def[cr_constants.AGG_ID_FIELD],
                'molecule_h5': sample_def[cr_constants.AGG_H5_FIELD],
                '__mem_gb': mem_gb,
            })
    return {'chunks': chunks}

def main(args, outs):
    with cr_mol_counter.MoleculeCounter.open(args.molecule_h5, 'r') as in_mc:
        with cr_mol_counter.MoleculeCounter.open(outs.merged_molecules, 'w') as out_mc:
            remapped_gem_groups = remap_gems(in_mc.get_column('gem_group'), args.gem_group_index, args.library_id)
            sort_index = np.lexsort([remapped_gem_groups])

            for col in cr_mol_counter.MOLECULE_INFO_COLUMNS:
                if col == 'gem_group':
                    arr = remapped_gem_groups
                else:
                    arr = in_mc.get_column(col)
                out_mc.add_many(col, arr[sort_index])

            for col in cr_mol_counter.MOLECULE_REF_COLUMNS:
                array = in_mc.get_ref_column(col)
                out_mc.set_ref_column(col, array)

            out_metrics = in_mc.get_all_metrics()
            gg_metrics = {}
            for (gg, metrics) in in_mc.get_metric(cr_mol_counter.GEM_GROUPS_METRIC).iteritems():
                for ng, (sid, og) in args.gem_group_index.iteritems():
                    if sid == args.library_id and og == gg:
                        gg_metrics[int(ng)] = metrics

            out_metrics[cr_mol_counter.GEM_GROUPS_METRIC] = gg_metrics
            out_mc.set_all_metrics(out_metrics)

def remap_gems(old_gems, gem_group_index, library_id):
    new_gems = np.copy(old_gems)
    for ng, (sid, og) in gem_group_index.iteritems():
        if sid == library_id:
            new_gems[old_gems==og] = ng
    return new_gems

def join(args, outs, chunk_defs, chunk_outs):
    molecules = [chunk_out.merged_molecules for chunk_out in chunk_outs]
    outs.chunked_molecules = molecules
    metrics = cr_mol_counter.MoleculeCounter.naive_concatenate_metrics(molecules)
    metrics[cr_mol_counter.IS_AGGREGATED_METRIC] = True
    cr_mol_counter.MoleculeCounter.concatenate(outs.merged_molecules, molecules, metrics=metrics)
