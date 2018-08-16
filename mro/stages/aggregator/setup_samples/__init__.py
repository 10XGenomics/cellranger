#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

import copy
import cellranger.constants as cr_constants
from cellranger.molecule_counter import MoleculeCounter
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
    new_gg = 0
    gg_index = {}
    libraries = []

    ### Batch info
    # If a column 'batch' is given in sample_defs (read from input csv), that
    # column will be used as batch identifier, otherwise, aggr_id will be 
    # used as batch identifier.
    # Each batch will have a distinct batch_id, which is an increasing integer. 
    batch_name_to_id = {}

    for sample_def in args.sample_defs:
        seen_ggs = set()

        aggr_id = sample_def[cr_constants.AGG_ID_FIELD]
        batch_name = sample_def.get(cr_constants.AGG_BATCH_FIELD, aggr_id)
        if batch_name not in batch_name_to_id:
            batch_name_to_id[batch_name] = len(batch_name_to_id) 

        with MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], 'r') as mc:
            old_libraries = mc.get_library_info()

        for lib_idx, old_lib in enumerate(old_libraries):
            # Remap gem groups
            old_gg = old_lib['gem_group']

            # Increment gem group if this is a new one from the same input sample
            if old_gg not in seen_ggs:
                new_gg += 1

            gg_index[new_gg] = (aggr_id, old_gg)

            # Remap libraries
            new_lib = copy.deepcopy(old_lib)
            new_lib['gem_group'] = new_gg

            # Make the new library id unique
            new_lib['library_id'] += ".%d" % (new_gg)
            new_lib['old_library_index'] = lib_idx
            new_lib['old_gem_group'] = old_gg
            new_lib['aggr_id'] = sample_def[cr_constants.AGG_ID_FIELD]
            new_lib['batch_name'] = batch_name
            new_lib['batch_id'] = batch_name_to_id[batch_name]
            libraries.append(new_lib)

            # Track gem groups
            seen_ggs.add(old_gg)

    outs.libraries = libraries
    outs.gem_group_index = gg_index

    # Write the "gem group index" (a legacy structure) for Loupe
    with open(outs.gem_group_index_json, 'w') as outfile:
        json.dump({"gem_group_index": gg_index}, outfile)
