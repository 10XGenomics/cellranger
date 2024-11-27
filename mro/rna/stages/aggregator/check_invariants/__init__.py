#!/usr/bin/env python3
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import math
from collections.abc import Iterable

import martian
import numpy as np
from six import ensure_binary

import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import tenkit.safe_json as tk_safe_json
from cellranger.molecule_counter import MoleculeCounter
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

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
    matrix_mem_gb = int(
        1.5 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.merged_raw_gene_bc_matrices_h5)
    )
    mol_info_mem_gb = (
        np.max(
            [
                int(
                    math.ceil(
                        MoleculeCounter.estimate_mem_gb(
                            MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], "r").nrows()
                        )
                    )
                )
                for sample_def in args.input_sample_defs
            ]
        )
        + 1
    )  # 1 for buffer for storing invariant data
    join_args = {
        "__mem_gb": 12 + int(0.8 * max(matrix_mem_gb, h5_constants.MIN_MEM_GB, mol_info_mem_gb)),
    }
    return {"chunks": [], "join": join_args}


def _update_barcode_counts(
    mc: MoleculeCounter,
    mol_gem_group: Iterable[int],
    gem_groups: list[int],
    library_id: bytes,
    input_bc_counts: dict[tuple[bytes, int], np.ndarray],
):
    # trim the barcodes by only retain barcodes with a count or pass filters
    mol_barcode_idx = mc.get_column("barcode_idx")
    mol_barcode_info = mc.get_barcode_info()
    pf_bc_idx = mol_barcode_info.pass_filter[:, 0]
    for gg in gem_groups:
        bc_idx, counts = np.unique(mol_barcode_idx[mol_gem_group == gg], return_counts=True)
        full_input_bc_counts = np.zeros(len(mc.get_ref_column("barcodes")))
        full_input_bc_counts[bc_idx] = counts
        trimmed_bc_idx = np.union1d(pf_bc_idx, bc_idx[counts > 0])
        input_bc_counts[(library_id, gg)] = full_input_bc_counts[trimmed_bc_idx]


def _update_feature_counts(
    mc: MoleculeCounter,
    mol_gem_group: Iterable[int],
    gem_groups: list[int],
    library_id: bytes,
    input_feature_counts: dict[tuple[bytes, int], np.ndarray],
):
    mol_feature_idx = mc.get_column("feature_idx")
    for gg in gem_groups:
        num_features = len(mc.feature_reference.feature_defs)
        input_feature_counts[(library_id, gg)] = np.zeros(num_features)
        feature_idx, counts = np.unique(mol_feature_idx[mol_gem_group == gg], return_counts=True)
        input_feature_counts[(library_id, gg)][feature_idx] = counts


def _update_counts(
    mc: MoleculeCounter,
    library_id: bytes,
    input_bc_counts: dict[tuple[bytes, int], np.ndarray],
    input_feature_counts: dict[tuple[bytes, int], np.ndarray],
) -> int:
    gem_groups = mc.get_gem_groups()
    mol_gem_group = mc.get_column("gem_group")

    _update_barcode_counts(mc, mol_gem_group, gem_groups, library_id, input_bc_counts)
    _update_feature_counts(mc, mol_gem_group, gem_groups, library_id, input_feature_counts)
    return len(gem_groups)


def join(args, outs, chunk_defs, chunk_outs):
    # compute invariants on input data
    input_genomes = set()
    input_features = set()
    input_bc_counts: dict[tuple[bytes, int], np.ndarray] = {}
    input_feature_counts: dict[tuple[bytes, int], np.ndarray] = {}
    input_num_gem_groups = 0

    # need to do this up front so we know to tally targeted features only
    target_gene_ids = load_target_feature_ids(args.input_sample_defs)
    if target_gene_ids is None:
        is_targeted_aggr = False
        input_target_features = None
    else:
        is_targeted_aggr = True
        input_target_features = {
            (gene_id, GENE_EXPRESSION_LIBRARY_TYPE) for gene_id in target_gene_ids
        }

    for sample_def in args.input_sample_defs:
        library_id = ensure_binary(sample_def["library_id"])
        with MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], "r") as mc:
            input_genomes.update(mc.get_genomes())
            input_features.update(mol_counter_features_id_type(mc))
            input_num_gem_groups += _update_counts(
                mc, library_id, input_bc_counts, input_feature_counts
            )

    # compute invariants on output
    output_matrix = cr_matrix.CountMatrix.load_h5_file(args.merged_raw_gene_bc_matrices_h5)
    output_genomes = set(output_matrix.get_genomes())
    if len(output_genomes) == 0:
        # FIXME this is the GEX-less case, many invariants below will fail
        return
    output_features = set(count_matrix_features_id_type(output_matrix))
    output_target_features = count_matrix_target_features(output_matrix)
    output_bc_counts = {}
    output_feature_counts = {}
    output_gem_index = cr_matrix.get_gem_group_index(args.merged_raw_gene_bc_matrices_h5)
    assert output_gem_index is not None
    output_num_gem_groups = len(output_gem_index)

    for gg, (library_id, old_gg) in output_gem_index.items():
        assert isinstance(library_id, bytes), type(library_id)
        matrix_gg = output_matrix.select_barcodes_by_gem_group(gg)
        output_bc_counts[(library_id, old_gg)] = matrix_gg.get_counts_per_bc()
        output_feature_counts[(library_id, old_gg)] = matrix_gg.get_counts_per_feature()

    exit_message = (
        "An internal problem in the aggr pipeline has been detected "
        "that might lead to incorrect results. Please report this "
        "problem to support@10xgenomics.com."
    )

    if input_genomes != output_genomes:
        martian.log_info("Genomes differ between input molecule files and aggregated matrix")
        martian.exit(exit_message)
    if input_features != output_features:
        martian.log_info("Features differ between input molecule files and aggregated matrix")
        martian.exit(exit_message)
    if is_targeted_aggr and input_target_features != output_target_features:
        martian.log_info(
            "Target features from input molecule files differ from features in aggregated matrix"
        )
        martian.exit(exit_message)
    if input_num_gem_groups != output_num_gem_groups:
        martian.log_info(
            "Number of GEM groups differs between input molecule files and aggregated matrix"
        )
        martian.exit(exit_message)
    for lib_gg, in_count in input_bc_counts.items():
        if len(in_count) != len(output_bc_counts[lib_gg]):
            martian.log_info(
                f"Barcode list for library {lib_gg[0]}, GEM group {lib_gg[1]} has different length "
                "in aggregated output compared to input."
            )
            martian.exit(exit_message)
        if np.any(in_count < output_bc_counts[lib_gg]):
            martian.log_info(
                f"Barcode(s) in library {lib_gg[0]}, GEM group {lib_gg[1]} have higher UMI counts "
                "in aggregated output compared to inputs"
            )
            martian.exit(exit_message)
        if len(input_feature_counts[lib_gg]) != len(output_feature_counts[lib_gg]):
            martian.log_info(
                f"Feature list for library {lib_gg[0]}, GEM group {lib_gg[1]} has different length "
                "in aggregated output compared to input."
            )
            if is_targeted_aggr:
                # Permit aggr of 5'/3' GEX with RTL. The RTL analysis has additional DEPRECATED
                # probe features after the gene features that GEX does not have.
                pass
            else:
                martian.exit(exit_message)
        elif np.any(input_feature_counts[lib_gg] < output_feature_counts[lib_gg]):
            martian.log_info(
                f"Feature(s) in library {lib_gg[0]}, GEM group {lib_gg[1]} have higher UMI counts "
                "in aggregated output compared to inputs"
            )
            martian.exit(exit_message)

    summary = {
        "genomes_present": list(input_genomes),
        "num_features_in_ref": len(input_features),
        "num_gem_groups": input_num_gem_groups,
    }

    with open(outs.summary, "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)


def mol_counter_features_id_type(mol_counter):
    return ((f.id, f.feature_type) for f in mol_counter.feature_reference.feature_defs)


def count_matrix_features_id_type(count_matrix):
    return ((f.id, f.feature_type) for f in count_matrix.feature_ref.feature_defs)


def count_matrix_target_features(count_matrix):
    target_feature_indices = count_matrix.feature_ref.get_target_feature_indices()
    if target_feature_indices is None:
        return None
    else:
        return {
            (f.id, f.feature_type)
            for f in count_matrix.feature_ref.feature_defs
            if f.index in target_feature_indices
        }


def load_target_feature_ids(sample_defs):
    """Load the ids of target features from the first targeted sample in sample_defs,.

    or else None if there is no targeted sample. For the result to be meaningful we
    assume all target sets present are the same (which holds for targeted aggr).
    """
    for sample_def in sample_defs:
        with MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], "r") as mc:
            feature_ref = mc.get_feature_ref()
            if feature_ref.has_target_features():
                return feature_ref.get_target_feature_ids()
    return None
