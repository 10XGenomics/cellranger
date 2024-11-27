#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
"""Aggr preflight check + convert legacy molecule info h5 to current version."""

import os
import shutil

import martian

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
import cellranger.rna.library as rna_library
from cellranger.aggr.preprocessing import check_molecule_info_version_split
from cellranger.analysis.constants import CBC_MAX_NCELLS
from cellranger.feature.antigen.specificity import (
    BEAM_AB,
    BEAM_T,
    MHC_ALLELE,
    NO_ALLELE,
    TARGETING_ANTIGEN,
)
from cellranger.molecule_counter_converter import convert_v2_to_v4, upgrade_file

__MRO__ = """
stage CHECK_MOLECULE_INFO_VERSION(
    in  map[]       sample_defs,
    in  string      product_type,
    in  bool        is_pd,
    out map[]       updated_sample_defs,
    out bool        is_not_pd,
    out string      beam_mode,
    out bool        is_spatial,
    out map<string> antigen_specificity_controls,
    out csv         feature_reference,
    out bool        disable_antigen_aggr,
    src py          "stages/aggregator/check_molecule_info_version",
) split (
    in  int         mol_h5_version,
    in  map         sample_def,
    out map         updated_sample_def,
) using (
    volatile = strict,
)
"""
SP_PRODUCT = "sp"


def split(args):
    return check_molecule_info_version_split(args)


def main(args, outs):
    outs.updated_sample_def = args.sample_def.copy()

    if args.mol_h5_version == cr_mol_counter.CURR_FILE_VERSION:
        # avoid copying, pass it along
        return

    v2_mole_info_h5 = args.sample_def[cr_constants.AGG_H5_FIELD]
    v2_file_basename = os.path.basename(v2_mole_info_h5)
    out_mole_info_h5 = martian.make_path(v2_file_basename)

    if args.mol_h5_version == 2:
        convert_v2_to_v4(v2_mole_info_h5, out_mole_info_h5)
    else:
        shutil.copy(v2_mole_info_h5, out_mole_info_h5)

    try:
        with cr_mol_counter.MoleculeCounter.open(out_mole_info_h5, "r+") as mc:
            upgrade_file(mc)
    except ValueError as err:
        martian.exit(str(err))

    outs.updated_sample_def[cr_constants.AGG_H5_FIELD] = out_mole_info_h5


def join(args, outs, chunk_defs, chunk_outs):
    outs.is_not_pd = not args.is_pd
    outs.is_spatial = bool(args.product_type == SP_PRODUCT)
    outs.updated_sample_defs = [chunk_out.updated_sample_def for chunk_out in chunk_outs]

    outs.beam_mode = None
    outs.disable_antigen_aggr = True
    antigen_specificity_controls = {}  # {mhc_allele: control_feature_id}
    feature_ref = None
    for sample_def in outs.updated_sample_defs:
        with cr_mol_counter.MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], "r") as mc:
            if feature_ref is None:
                feature_ref = mc.feature_reference
            for lib in mc.library_info:
                if lib["library_type"] == rna_library.ANTIGEN_LIBRARY_TYPE:
                    outs.disable_antigen_aggr = False
                    antigen_feature_ref = mc.feature_reference.select_features(
                        mc.feature_reference.get_indices_for_type(rna_library.ANTIGEN_LIBRARY_TYPE)
                    )
                    # check if we have controls
                    if TARGETING_ANTIGEN in antigen_feature_ref.all_tag_keys:
                        # in case of beam_t
                        if MHC_ALLELE in antigen_feature_ref.all_tag_keys:
                            outs.beam_mode = BEAM_T
                            for feature in antigen_feature_ref.feature_defs:
                                if feature.tags[TARGETING_ANTIGEN] == "False":
                                    if feature.tags[MHC_ALLELE] not in antigen_specificity_controls:
                                        antigen_specificity_controls[feature.tags[MHC_ALLELE]] = (
                                            feature.id
                                        )
                                    else:
                                        # Each separate run has been checked to follow this assert
                                        # and we have checked that all feature references match.
                                        assert (
                                            antigen_specificity_controls[feature.tags[MHC_ALLELE]]
                                            == feature.id
                                        )
                        else:
                            outs.beam_mode = BEAM_AB
                            for feature in antigen_feature_ref.feature_defs:
                                if feature.tags[TARGETING_ANTIGEN] == "False":
                                    if NO_ALLELE not in antigen_specificity_controls:
                                        antigen_specificity_controls[NO_ALLELE] = feature.id
                                    else:
                                        # Each separate run has been checked to follow this assert
                                        # and we have checked that all feature references match.
                                        assert antigen_specificity_controls[NO_ALLELE] == feature.id
                    elif MHC_ALLELE in antigen_feature_ref.all_tag_keys:
                        outs.beam_mode = BEAM_T
                    else:
                        outs.beam_mode = BEAM_AB
    if feature_ref is not None:
        with open(outs.feature_reference, "w") as file_handle:
            feature_ref.to_csv(file_handle)
    outs.antigen_specificity_controls = antigen_specificity_controls
    if any("batch" in sample_def for sample_def in outs.updated_sample_defs):
        ncells = 0
        for sample_def in outs.updated_sample_defs:
            with cr_mol_counter.MoleculeCounter.open(
                sample_def[cr_constants.AGG_H5_FIELD], "r"
            ) as mc:
                ncells += mc.get_num_filtered_barcodes_for_library(0)

        if ncells > CBC_MAX_NCELLS:
            martian.exit(
                f"You provided {ncells:,} cells in total, but chemistry batch correction only supports up to {CBC_MAX_NCELLS:,} cells."
            )
