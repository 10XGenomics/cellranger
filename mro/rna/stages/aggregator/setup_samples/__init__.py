#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#


import copy
import json
from collections import OrderedDict

import martian

import cellranger.constants as cr_constants
import cellranger.rna.library as rna_library
from cellranger.molecule_counter import get_library_info

__MRO__ = """
stage SETUP_SAMPLES(
    in  map[] sample_defs,
    out map   gem_group_index,
    out json  gem_group_index_json,
    out map[] libraries,
    out bool  chemistry_batch_correction,
    src py    "stages/aggregator/setup_samples",
)
"""


def main(args, outs):
    new_gg = 0
    gg_index = {}
    libraries = []
    chemistry_batch_correction = False

    ### Batch info
    # If a column 'batch' is given in sample_defs (read from input csv), that
    # column will be used as batch identifier and chemistry_batch_correction will
    # be turned on. otherwise, aggr_id will be used as batch identifier.
    # Each batch will have a distinct batch_id, which is an increasing integer.
    batch_name_to_id = {}

    sample_defs = [] if args.sample_defs is None else args.sample_defs
    new_library_ids = OrderedDict()
    disable_crispr = True
    for sample_def in sample_defs:
        seen_ggs = set()

        aggr_id = sample_def[cr_constants.AGG_ID_FIELD]

        if cr_constants.AGG_BATCH_FIELD in sample_def:
            chemistry_batch_correction = True
            batch_name = sample_def[cr_constants.AGG_BATCH_FIELD]
        else:
            batch_name = aggr_id

        if batch_name not in batch_name_to_id:
            batch_name_to_id[batch_name] = len(batch_name_to_id)

        old_libraries = get_library_info(sample_def[cr_constants.AGG_H5_FIELD])

        for lib_idx, old_lib in enumerate(old_libraries):
            # Check for presence of crispr libraries
            if old_lib["library_type"] == rna_library.CRISPR_LIBRARY_TYPE:
                disable_crispr = False
            # Remap gem groups
            old_gg = int(old_lib["gem_group"])

            # Increment gem group if this is a new one from the same input sample
            if old_gg not in seen_ggs:
                new_gg += 1

            gg_index[new_gg] = (aggr_id, old_gg)

            # Remap libraries
            new_lib = copy.deepcopy(old_lib)
            new_lib["gem_group"] = new_gg

            # Make the new library id unique
            lib_key = (old_lib["library_id"], new_gg)
            new_lib_id = new_library_ids.setdefault(lib_key, len(new_library_ids))

            new_lib["library_id"] = new_lib_id
            new_lib["old_library_index"] = lib_idx
            new_lib["old_gem_group"] = old_gg
            new_lib["aggr_id"] = sample_def[cr_constants.AGG_ID_FIELD]
            new_lib["batch_name"] = batch_name
            new_lib["batch_id"] = batch_name_to_id[batch_name]
            libraries.append(new_lib)

            # Track gem groups
            seen_ggs.add(old_gg)

    if chemistry_batch_correction is True and len(batch_name_to_id) <= 1:
        chemistry_batch_correction = False
        martian.log_info(
            "Warning: only one batch sepecified in the input csv, chemistry_batch_correction is disabled."
        )

    outs.libraries = libraries
    outs.gem_group_index = gg_index
    outs.chemistry_batch_correction = chemistry_batch_correction

    outs.target_gene_ids_csv = None
    outs.disable_crispr_aggr = disable_crispr

    # Write the "gem group index" (a legacy structure) for Loupe
    with open(outs.gem_group_index_json, "w") as outfile:
        json.dump({"gem_group_index": gg_index}, outfile)
