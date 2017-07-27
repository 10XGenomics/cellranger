#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
import json

__MRO__ = '''
stage SETUP_SAMPLES(
    in  map[] sample_defs,
    out map   gem_group_index,
    out json  gem_group_index_json,
    src py    "stages/aggregator/setup_samples",
)
'''

def main(args, outs):
    # gem_group_index looks like: { new_gg : (library_id, old_gg) }
    current_gg = 1
    gg_index = {}
    for sample_def in args.sample_defs:
        # remap gem groups
        old_gem_groups = get_gem_groups(sample_def[cr_constants.AGG_H5_FIELD])
        for old_gg in sorted(old_gem_groups):
            gg_index[current_gg] = (sample_def[cr_constants.AGG_ID_FIELD], old_gg)
            current_gg += 1
    outs.gem_group_index = gg_index
    with open(outs.gem_group_index_json, 'w') as outfile:
        json.dump({"gem_group_index": gg_index}, outfile)

def get_gem_groups(molecule_h5):
    with cr_mol_counter.MoleculeCounter.open(molecule_h5, 'r') as ctr:
        gem_groups = ctr.get_gem_groups()
    return gem_groups
