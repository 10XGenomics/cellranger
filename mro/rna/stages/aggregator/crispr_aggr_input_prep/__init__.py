#!/usr/bin/env python3
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved
#
"""Prepare input for _CRISPR_ANALYZER stage in aggr."""

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

import cellranger.molecule_counter as cr_mc
import cellranger.rna.library as rna_library

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

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


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CRISPR_AGGR_INPUT_PREP stage inputs."""

    merged_molecules: mro_filetypes.h5


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CRISPR_AGGR_INPUT_PREP stage outputs."""

    feature_reference: mro_filetypes.csv


ChunkInputs = StageInputs
ChunkOutputs = StageOutputs


def main(args: StageInputs, outs: StageOutputs):
    with cr_mc.MoleculeCounter.open(args.merged_molecules, "r") as mc:
        # Create feature_reference csv file
        feature_ref = mc.get_feature_ref()
        with open(outs.feature_reference, "w") as file_handle:
            feature_ref.select_features_by_type(rna_library.CRISPR_LIBRARY_TYPE).to_csv(file_handle)
