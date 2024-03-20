#!/usr/bin/env python3
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved
#
"""Prepare input for _CRISPR_ANALYZER stage in aggr."""

import cellranger.molecule_counter as cr_mc
import cellranger.rna.library as rna_library

__MRO__ = """
stage CRISPR_AGGR_INPUT_PREP(
    in  h5   merged_molecules,
    out csv  feature_reference,
    src py   "stages/aggregator/crispr_aggr_input_prep",
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


def main(args, outs):
    with cr_mc.MoleculeCounter.open(args.merged_molecules, "r") as mc:
        # Create feature_reference csv file
        feature_ref = mc.get_feature_ref()
        with open(outs.feature_reference, "w") as file_handle:
            feature_ref.select_features_by_type(rna_library.CRISPR_LIBRARY_TYPE).to_csv(file_handle)
