#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#

"""Martian stage AGGREGATOR_PREFLIGHT."""

from typing import NamedTuple

import martian

import cellranger.constants as cr_constants
import cellranger.hdf5 as cr_h5
import cellranger.molecule_counter as cr_mol_counter
from cellranger.chemistry import CHEMISTRY_SC3P_LT, SC3P_CHEMISTRIES, SC5P_CHEMISTRIES
from cellranger.molecule_counter import TARGETING_METHOD_METRIC
from cellranger.molecule_counter_converter import convert_v2_to_v4
from cellranger.rna.library import (
    ANTIBODY_LIBRARY_TYPE,
    FEATURE_LIBRARY_TYPES,
    GENE_EXPRESSION_LIBRARY_TYPE,
)
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC, TARGETING_METHOD_TL
from cellranger.utils import string_is_ascii

SC3P_CHEMISTRY_NAMES = {chemistry["name"] for chemistry in SC3P_CHEMISTRIES}
SC5P_CHEMISTRY_NAMES = {chemistry["name"] for chemistry in SC5P_CHEMISTRIES}

__MRO__ = """
stage AGGREGATOR_PREFLIGHT(
    in  map[]  sample_defs,
    in  string normalization_mode,
    in  bool   is_pd,
    src py     "stages/aggregator/aggregator_preflight",
)
"""


class AggregatorPreflightStageInputs:
    sample_defs: list[dict[str, str]]
    normalization_mode: str
    is_pd: bool


def incompat_msg(reason) -> str:
    return (
        f"The datasets you are trying to aggregate were created with different {reason}s, "
        f"but the 'aggr' command requires identical {reason}s in order to combine datasets. "
        f"Please re-run the original pipeline ('count' or 'multi', as the case may be) with "
        f"uniform {reason} in order to aggregate these data."
    )


NORM_MODES = ["mapped", cr_constants.NORM_MODE_NONE]


def _get_include_introns(analysis_parameters):
    """Get include_introns value from analysis_parameters."""
    if analysis_parameters is None:
        return cr_mol_counter.INTRON_MODE_HISTORIC_DEFAULT

    return analysis_parameters.get("include_introns", cr_mol_counter.INTRON_MODE_HISTORIC_DEFAULT)


def _get_filter_probes(analysis_parameters):
    """Return filter_probes value from analysis_parameters or !no_probe_filter if present for backward compatibility."""
    if analysis_parameters is None:
        return cr_mol_counter.FILTER_PROBES_PARAM_DEFAULT

    # Since CR 7.0 filter_probes is an Optional<bool>
    filter_probes = analysis_parameters.get(cr_mol_counter.FILTER_PROBES, "missing_key")
    if filter_probes != "missing_key":
        return filter_probes

    no_probe_filter = analysis_parameters.get(cr_mol_counter.NO_PROBE_FILTER)
    if no_probe_filter is not None:
        return not no_probe_filter

    return cr_mol_counter.FILTER_PROBES_PARAM_DEFAULT


def convert_v2_to_v4_if_needed(filename: str, is_pd: bool) -> str:
    """Convert v2 to v4 molecule_info.h5 if necessary, otherwise return the original filename."""
    _mc_h5, file_version = cr_mol_counter.get_h5py_file_and_version(filename)
    min_version = 2 if is_pd else 3
    if file_version < min_version:
        martian.exit(
            f"The molecule info HDF5 file (file: {filename}, format version {file_version}) "
            "was produced by an older software version. Reading these files is unsupported. "
            "Please use Cell Ranger 3.0 or newer to generate a supported molecule info file."
        )
    elif file_version == 2:
        assert is_pd
        # Convert v2 to v4, because MoleculeCounter.open does not support opening v2.
        v4_filename = "molecule_info.h5"
        convert_v2_to_v4(filename, v4_filename)
        return v4_filename
    else:
        # No conversion needed for versions v3+.
        return filename


# pylint: disable=too-many-branches,too-many-statements
# ruff: noqa: PLR0915, PLR0912
def main(args: AggregatorPreflightStageInputs, _outs: None):
    if args.normalization_mode is not None and args.normalization_mode not in NORM_MODES:
        martian.exit("Normalization mode must be one of: {}".format(", ".join(NORM_MODES)))

    global_fasta_hash = None
    global_gtf_hash = None
    global_feature_ref = None

    # TODO make assertions about the required metrics!
    # Gather library types for an aggr check
    library_types = set()

    # Gather target set info for targeted samples
    targeting_methods = set()
    target_sets = set()
    observed_ends = set()
    target_panel_hashes = set()

    class GexRtlFeatureCompatibility(NamedTuple):
        has_deprecated_probes: bool
        has_fb: bool

    gex_rtl_feature_compatibility: set[GexRtlFeatureCompatibility] = set()

    include_introns = None
    filter_probes = None

    # Gather antigen control info
    observed_ag_control_features = set()

    # check that all molecule files conform to spec
    libraries_seen = set()
    for sample in args.sample_defs:
        library_id = sample[cr_constants.AGG_ID_FIELD]
        if len(library_id) == 0:
            martian.exit(f"Library ID cannot be empty: {sample}")

        if not string_is_ascii(library_id):
            martian.exit(
                f"Library ID {library_id} contains unicode characters, only ASCII is allowed."
            )

        if cr_constants.AGG_BATCH_FIELD in sample:
            batch_name = sample[cr_constants.AGG_BATCH_FIELD]
            if not string_is_ascii(batch_name):
                martian.exit(
                    f"Batch ID {batch_name} contains unicode characters, only ASCII is allowed."
                )

        if library_id in libraries_seen:
            martian.exit(f"Same library ID is specified on multiple rows: {library_id}")
        else:
            libraries_seen.add(library_id)
        if not cr_h5.is_hdf5(sample[cr_constants.AGG_H5_FIELD]):
            martian.exit(
                f"Input molecule file is not a valid HDF5 file: {sample[cr_constants.AGG_H5_FIELD]}"
            )
        mol_h5 = sample[cr_constants.AGG_H5_FIELD]
        converted_mol_h5 = convert_v2_to_v4_if_needed(mol_h5, args.is_pd)
        try:
            counter = cr_mol_counter.MoleculeCounter.open(converted_mol_h5, "r")
        except ValueError as err:
            martian.exit(str(err))
        with counter:
            library_info = counter.get_library_info()
            sample_libs = frozenset(x["library_type"] for x in library_info)
            library_types.add(sample_libs)
            mol_cr_version = counter.get_metric(cr_constants.CELLRANGER_VERSION_KEY)
            # pre CR 7.0 there was no targeting_method metric and hybrid capture was the only targeting method
            targeting_method = counter.get_metric(
                TARGETING_METHOD_METRIC,
                (
                    TARGETING_METHOD_HC
                    if any(counter.is_targeted_library(x) for x in library_info)
                    else None
                ),
            )

            if targeting_method == TARGETING_METHOD_HC:
                martian.exit(
                    f"Hybrid capture targeted analysis is no longer supported by Cell Ranger; please use Cell Ranger 7.1 or earlier: {mol_h5}"
                )

            if not mol_cr_version:
                martian.exit(
                    f"Molecule file was produced with old cellranger version (missing version number): {mol_h5}"
                )

            molecule_info_type = counter.get_molecule_info_type()

            if molecule_info_type == cr_mol_counter.MOLECULE_INFO_TYPE_RAW and not args.is_pd:
                martian.exit(
                    "Raw molecule info files (raw_molecule_info.h5) from Cellranger Multi cannot be used with Cellranger Aggr. \n \
       Please provide only sample molecule info files from Cellranger Multi (sample_molecule_info.h5) or molecule infos produced by Cellranger Count (molecule_info.h5)."
                )

            mol_fasta_hash = counter.get_metric("reference_fasta_hash", "no-transcriptome")
            if global_fasta_hash is None:
                global_fasta_hash = mol_fasta_hash
            elif global_fasta_hash != mol_fasta_hash:
                martian.exit(
                    "{} (hashes: {} != {})".format(
                        incompat_msg("genome reference"), global_fasta_hash, mol_fasta_hash
                    )
                )

            mol_gtf_hash = counter.get_metric("reference_gtf_hash")
            if global_gtf_hash is None:
                global_gtf_hash = mol_gtf_hash
            elif global_gtf_hash != mol_gtf_hash:
                martian.exit(
                    "{} (hashes: {} != {})".format(
                        incompat_msg("annotation GTF"), global_gtf_hash, mol_gtf_hash
                    )
                )

            chemistry_name = counter.get_metric("chemistry_name")

            if chemistry_name == CHEMISTRY_SC3P_LT["name"]:
                martian.exit(
                    f"Library {library_id} was created with the {CHEMISTRY_SC3P_LT['name']} chemistry, which is no longer supported."
                )

            chemistry_endedness = counter.get_metric("chemistry_endedness")
            if chemistry_endedness is not None:
                observed_ends.add(chemistry_endedness)
            elif chemistry_name in SC3P_CHEMISTRY_NAMES:
                # Pre-cellranger 2.1 or so does not have chemistry_endedness, so infer based on chemistry_name
                observed_ends.add(cr_constants.THREE_PRIME)
            elif chemistry_name in SC5P_CHEMISTRY_NAMES:
                observed_ends.add(cr_constants.FIVE_PRIME)

            # check analysis parameters
            analysis_parameters = counter.get_metric("analysis_parameters")
            include_introns_sample = _get_include_introns(analysis_parameters)

            # Intron mode must match for any samples where it is relevant (not templated ligation)
            if not targeting_method == TARGETING_METHOD_TL and not args.is_pd:
                if include_introns is None:
                    include_introns = include_introns_sample
                elif include_introns != include_introns_sample:
                    martian.exit(
                        f"Molecule file was produced with a different include-introns setting than previous file: {mol_h5}.\n"
                        f" All samples must be run with include-introns set to the same value."
                    )

            filter_probes_sample = _get_filter_probes(analysis_parameters)
            if filter_probes is None:
                filter_probes = filter_probes_sample
            elif filter_probes_sample is None:
                filter_probes_sample = filter_probes
            elif filter_probes != filter_probes_sample:
                martian.exit(
                    "Molecule files provided were produced with different --filter-probes settings. All samples must be run with the same --filter-probes setting."
                )
            # Agg of Antibody Capture and Gene expression not allowed
            if {ANTIBODY_LIBRARY_TYPE} in library_types:
                if {ANTIBODY_LIBRARY_TYPE, GENE_EXPRESSION_LIBRARY_TYPE} in library_types:
                    martian.exit(
                        "Aggr with Antibody Capture and Gene Expression + Antibody Capture libraries is not supported."
                    )
            mol_feature_ref = counter.feature_reference
            assert mol_feature_ref is not None
            if global_feature_ref is None:
                global_feature_ref = mol_feature_ref
            elif not global_feature_ref.has_compatible_target_set(mol_feature_ref):
                martian.exit(
                    "Aggr of targeted Gene Expression datasets analyzed with different Target Panel or Probe set CSV files is not supported."
                )
            elif not global_feature_ref.equals_ignore_target_set(mol_feature_ref):
                martian.exit(incompat_msg("feature reference"))

            antigen_equals, antigen_equals_msg = global_feature_ref.equals_antigen_capture_content(
                mol_feature_ref, args.is_pd
            )
            if not antigen_equals:
                martian.exit(antigen_equals_msg)

            antigen_control_this_sample = mol_feature_ref.get_antigen_control()
            if antigen_control_this_sample is not None:
                observed_ag_control_features.add(antigen_control_this_sample)

            if counter.is_aggregated():
                martian.exit(
                    f"Molecule file was aggregated from multiple samples: {mol_h5}.\n"
                    f" Aggregated outputs cannot be re-aggregated, please pass each of the original samples instead."
                )

            if counter.nrows() == 0:
                martian.exit(
                    f"Cannot aggregate file because it contains no data: {mol_h5}.\n"
                    f" Please remove this file from the aggregation and try again."
                )
            if counter.get_num_filtered_barcodes() == 0:
                martian.exit(
                    f"Cannot aggregate file because it contains no cells: {mol_h5}.\n"
                    f" Please remove this file from the aggregation and try again."
                )

            for lib_key, metrics in (counter.get_metric(cr_mol_counter.LIBRARIES_METRIC)).items():
                lib_total_reads = metrics[cr_mol_counter.TOTAL_READS_METRIC]
                if lib_total_reads == 0:
                    lib_type = library_info[int(lib_key)]["library_type"]
                    martian.log_warn(
                        f"Library {lib_key} ({lib_type}) has zero reads in file: {mol_h5}\n"
                        "Barcodes from this library will have zero UMI counts."
                    )

            # Track targeting-specific fields
            targeting_methods.add(targeting_method)

            if (target_feature_indices := mol_feature_ref.get_target_feature_indices()) is not None:
                target_sets.add(frozenset(target_feature_indices))
                target_panel_hashes.add(counter.get_metric("target_panel_hash"))

            # Track relevant properties for RTL/GEX aggr feature compatibility
            has_fb = any(
                mol_feature_ref.has_feature_type(lib_type) for lib_type in FEATURE_LIBRARY_TYPES
            )
            has_deprecated_probes = mol_feature_ref.has_deprecated_probes()
            gex_rtl_feature_compatibility.add(
                GexRtlFeatureCompatibility(has_deprecated_probes, has_fb)
            )

    # Acceptable combinations of targeting_method
    # OK non-targeted + non-targeted
    # OK hybrid_capture + hybrid_capture
    # OK hybrid_capture + non-targeted
    # OK templated_ligation + templated_ligation
    # OK templated_ligation + non-targeted
    # NO templated_ligation + hybrid_capture
    fail = None
    if len(targeting_methods) == 1:
        fail = False
    elif len(targeting_methods) >= 3:
        fail = True
    elif targeting_methods == {None, TARGETING_METHOD_HC}:
        # aggr of hybrid_capture and its parent is permitted.
        fail = False
    elif targeting_methods == {None, TARGETING_METHOD_TL}:
        # aggr of templated_ligation and non-targeted is permitted.
        fail = False
    elif targeting_methods == {TARGETING_METHOD_HC, TARGETING_METHOD_TL}:
        # aggr of hybrid_capture and templated_ligation is not permitted.
        fail = True
    assert fail is not None
    if fail:
        martian.exit("Aggregating samples from different 10x products is not supported.")

    if len(target_sets) > 0:
        # Only allow a single target set (with or without any nontargeted)
        if len(target_sets) > 1:
            martian.exit(
                "Aggr of targeted Gene Expression datasets analyzed with different Target Panel or Probe Set CSV files is not supported."
            )

        # Disallow 3P + 5P together if there is a HC targeted sample
        if len(observed_ends) > 1 and TARGETING_METHOD_HC in targeting_methods:
            martian.exit(
                "Aggr of Hybrid Capture targeted Gene Expression datasets and other samples with different endedness (e.g. 5P Gene Expression vs 3P Gene Expression) is not supported."
            )

        if len(target_panel_hashes) > 1:
            martian.exit(
                "Aggr must only include samples run using the same target panel or probe set csv file."
            )

    # RTL + GEX when the RTL sample has deprecated probes + both have antibody will not work unless they
    # both have the same probe set specified (there is already a check they are the same above)
    if (
        GexRtlFeatureCompatibility(has_deprecated_probes=False, has_fb=True)
        in gex_rtl_feature_compatibility
        and GexRtlFeatureCompatibility(has_deprecated_probes=True, has_fb=True)
        in gex_rtl_feature_compatibility
    ):
        martian.exit(
            "Aggr of Flex and Gene Expression with Feature Barcode requires that all input analyses include the same probe-set. Please rerun the Gene Expression with Feature Barcode analyses with the same probe-set parameter as used for the Flex analyses."
        )

    if len(observed_ag_control_features) > 1:
        martian.exit(
            "The datasets you are trying to aggregate have incompatible control "
            "feature ids. Please re-run the original multi pipelines with uniform "
            "[antigen-specificity] sections."
        )
