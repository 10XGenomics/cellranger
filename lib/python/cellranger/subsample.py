#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import pickle
from collections.abc import Collection, Iterable
from copy import copy
from typing import TYPE_CHECKING, TypedDict

import numpy as np
import pandas as pd
from six import ensure_binary

import cellranger.molecule_counter as cr_mc
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import tenkit.stats as tk_stats
from cellranger.constants import OFF_TARGET_SUBSAMPLE, ON_TARGET_SUBSAMPLE
from cellranger.fast_utils import (  # pylint: disable=no-name-in-module,unused-import
    FilteredBarcodes,
)

if TYPE_CHECKING:
    from cellranger.molecule_counter import MoleculeCounter

#####################################
# Subsampling types and target depths

# Fixed depths (reads per cell) for subsampling
#   - for targeted gene expression libraries
SUBSAMPLE_TARGETED_FIXED_DEPTHS = [
    100,
    250,
    500,
    1000,
    2500,
    3000,
    5000,
    10000,
    15000,
    20000,
    30000,
    40000,
    50000,
]

#  - for bulk comparisons
SUBSAMPLE_BULK_FIXED_DEPTHS = [
    int(_x)
    for _x in [
        1e4,
        5e4,
        1e5,
        2.5e5,
        5e5,
        1e6,
        2.5e6,
        5e6,
        7.5e6,
        1e7,
        5e7,
        1e8,
        1e9,
    ]
]

#  - for all other libraries
SUBSAMPLE_FIXED_DEPTHS = [1000, 3000, 5000, 10000, 20000, 30000, 50000]

# Number of additional quantile-based depth targets
SUBSAMPLE_NUM_ADDITIONAL_DEPTHS = 10


RAW_SUBSAMPLE_TYPE = "raw_rpc"
MAPPED_SUBSAMPLE_TYPE = "conf_mapped_barcoded_filtered_bc_rpc"
RAW_CELLS_SUBSAMPLE_TYPE = "raw_barcoded_filtered_bc_rpc"
BULK_SUBSAMPLE_TYPE = "raw_reads"

ALL_SUBSAMPLE_TYPES = [RAW_SUBSAMPLE_TYPE, MAPPED_SUBSAMPLE_TYPE]

SUBSAMPLE_TARGET_MODES = [ON_TARGET_SUBSAMPLE, OFF_TARGET_SUBSAMPLE, None]


def get_genomes(molecule_info: MoleculeCounter) -> list[str]:
    assert molecule_info.feature_reference is not None
    genomes = {genome: () for genome in molecule_info.feature_reference.get_genomes()}
    for genome in molecule_info.get_barcode_info().genomes:
        genomes[genome] = ()
    return sorted(genomes)


def get_num_cells_per_library(library_info, filtered_barcodes_csv):
    """Get the number of cell-associated (i.e., filtered) barcodes per library.

    Note:
        We assume cell barcodes are assigned per GEM group, not per library.

    Args:
        library_info (dict): library_info metadata, probably from a MoleculeCounter
        filtered_barcodes_csv (str): path to filtered_barcodes.csv file

    Returns:
        np.array of int: number of cell-associated barcodes per library
    """
    # get number of cells per GEM group
    num_cells_per_gg = FilteredBarcodes(filtered_barcodes_csv).cells_per_gem_group()
    # each library maps to a single gem group
    num_cells_per_lib = np.array(
        [num_cells_per_gg.get(lib["gem_group"], 0) for lib in library_info]
    )
    return num_cells_per_lib


def get_cell_associated_barcodes(genomes: Iterable[str], filtered_barcodes_csv):
    """Get cell-associated barcodes by genome.

    Args:
        genomes (list of str): Genome names.
        filtered_barcodes_csv (str): Path to CSV file.

    Returns:
        dict of (str, set): Map genome to list of cell-assoc barcodes.
            Empty-string key is for all genomes.
    """
    # Get all cell-assoc barcodes (ignoring genome) for the '' (blank) genome string
    cell_bcs = {
        genome: cr_utils.get_cell_associated_barcode_set(
            filtered_barcodes_csv, ensure_binary(genome)
        )
        for genome in genomes
    }
    # All cell-associated barcodes
    cell_bcs[""] = set.union(*cell_bcs.values())
    return cell_bcs


def compute_target_depths(
    max_target: float, num_targets: int
) -> np.ndarray[int, np.dtype[np.int_]]:
    """Construct a list of sorted, unique, integer-valued subsampling depths.

    Generally corresponding to target read pairs per cell.

    Args:
        max_target (float): the largest target depth
        num_targets (int): desired number of targets, including max_target. There
            will be fewer than this many targets in case num_targets > max_target.

    Returns:
        numpy array of int: target subsampling depths (sorted, distinct, nonzero)
    """
    distinct_targets = np.unique(
        np.linspace(start=0, stop=max_target, num=num_targets + 1, dtype=int)
    )
    return distinct_targets[distinct_targets > 0]


class SubsamplingDef(TypedDict):
    """The dictionary of metadata returned by make_subsamplings."""

    library_type: str
    subsample_type: str
    target_read_pairs_per_cell: int
    library_subsample_rates: list[float]


def _subsampling_for_depth(
    target_depth: int,
    subsample_type: str,
    library_type: str,
    num_cells_per_lib: np.ndarray[int, np.dtype[np.float64]],
    usable_frac_per_lib: np.ndarray[int, np.dtype[np.float64]],
    max_computed_depth: int | None,
    lib_indices: np.ndarray[int, np.dtype[np.int32]],
    raw_reads_per_lib: np.ndarray[int, np.dtype[np.float64]],
    usable_reads_per_lib: np.ndarray[int, np.dtype[np.float64]],
    library_count: int,
) -> SubsamplingDef:
    if subsample_type == BULK_SUBSAMPLE_TYPE:
        target_usable_reads_per_lib = np.full(
            fill_value=target_depth, shape=num_cells_per_lib.shape, dtype=float
        )
    elif subsample_type == MAPPED_SUBSAMPLE_TYPE:
        target_usable_reads_per_lib = target_depth * num_cells_per_lib
    else:
        # convert target raw read depth to usable read depth
        target_usable_reads_per_lib = target_depth * num_cells_per_lib * usable_frac_per_lib

    # compute subsampling rates (frac. of usable reads)
    subsample_rates: np.ndarray[int, np.dtype[np.float64]] = np.zeros(library_count, dtype=float)
    if subsample_type == BULK_SUBSAMPLE_TYPE:
        denominator = raw_reads_per_lib
    else:
        denominator = usable_reads_per_lib

    for index in lib_indices:
        if denominator[index] != 0.0:
            subsample_rates[index] = target_usable_reads_per_lib[index] / denominator[index]

    # for the largest computed (non-default) subsampling depth,
    # make sure we're subsampling the smallest library to rate=1.0
    if target_depth == max_computed_depth:
        max_rate = np.max(subsample_rates)
        if max_rate != 0.0:
            subsample_rates = subsample_rates / max_rate
    # zero out rates that are > 1
    # This can only apply to the the default subsampling targets,
    # for which we still want to run subsampling jobs with rate=0.
    subsample_rates[subsample_rates > 1.0] = 0.0

    return {
        "library_type": library_type,
        "subsample_type": subsample_type,
        "target_read_pairs_per_cell": int(target_depth),
        "library_subsample_rates": list(subsample_rates),
    }


def make_subsamplings(
    subsample_type: str,
    library_info,
    library_type: str,
    num_cells_per_lib: np.ndarray,
    raw_reads_per_lib: np.ndarray,
    usable_reads_per_lib: np.ndarray,
    fixed_depths: list[int],
    num_additional_depths: int,
) -> list[SubsamplingDef]:
    """Create metadata for subsampling jobs of a specified subsampling type.

    (raw or usable reads per cell) and for a specified library type.

    Args:
        subsample_type (str): subsample based on raw, usable, or bulk reads? (raw
            corresponds to raw rpc, usable to confidently mapped transcriptomic rpc,
            and bulk to bulk transcriptomic reads)
        library_info (dict): per-library metadata from MoleculeCounter
        library_type (str): library type to use for subsampling
        num_cells_per_lib (np.array of int): number of filtered barcodes per library
        raw_reads_per_lib (np.array of int): number of raw reads per library
        usable_reads_per_lib (np.array of int): number of usable reads per library
        fixed_depths (list of int): fixed subsampling depths (reads per cell)
            to include by default
        num_additional_depths (int): number of subsampling depths to use,
            in addition to the defaults
    Returns:
        list of dict: list of subsampling metadata, each of which is:
            {'library_type': <str>,
             'subsample_type': <str>,
             'target_read_pairs_per_cell': <int>,
             'library_subsample_rates': <np.array of floats>}
    """
    lib_indices = np.array(
        [i for i, lib in enumerate(library_info) if lib["library_type"] == library_type],
        dtype=np.int32,
    )

    # Do casts once up-front, as they're surprisingly expensive.
    num_cells_per_lib = num_cells_per_lib.astype(float)
    raw_reads_per_lib = raw_reads_per_lib.astype(float)
    raw_rppc_per_lib = raw_reads_per_lib / num_cells_per_lib
    usable_reads_per_lib = usable_reads_per_lib.astype(float)
    usable_rppc_per_lib = usable_reads_per_lib / num_cells_per_lib
    # fraction of usable reads per library
    usable_frac_per_lib = usable_reads_per_lib / raw_reads_per_lib

    # Pick a range of target depths that are feasible for all of the given libraries
    if subsample_type == BULK_SUBSAMPLE_TYPE:
        max_target_depth = np.min(raw_reads_per_lib[lib_indices])
    else:
        max_target_depth = np.min(
            (
                raw_rppc_per_lib
                if subsample_type in (RAW_SUBSAMPLE_TYPE, RAW_CELLS_SUBSAMPLE_TYPE)
                else usable_rppc_per_lib
            )[lib_indices]
        )
    computed_depths = compute_target_depths(max_target_depth, num_additional_depths)
    # make note of the maximum after rounding
    # N.B. computed_depths is empty if max_target_depth is less than 1 read pair per cell.
    # This might happen if there are very few usable reads for this library type.
    max_computed_depth = np.max(computed_depths) if len(computed_depths) > 0 else None
    target_depths: np.ndarray[int, np.dtype[np.int_]] = (
        np.concatenate(  # pylint: disable=unexpected-keyword-arg
            [computed_depths, fixed_depths], dtype=int
        )
    )
    target_depths.sort()
    target_depths = np.unique(target_depths)

    return [
        _subsampling_for_depth(
            target_depth,
            subsample_type,
            library_type,
            num_cells_per_lib,
            usable_frac_per_lib,
            max_computed_depth,
            lib_indices,
            raw_reads_per_lib,
            usable_reads_per_lib,
            len(library_info),
        )
        for target_depth in target_depths
    ]


def construct_all_subsamplings(
    molecule_info_h5,
    filtered_barcodes_csv,
    is_targeted: bool = False,
    include_bulk: bool = True,
) -> list[SubsamplingDef]:
    """Construct subsampling metadata for a range of target read depths,.

    both raw reads per cell and usable reads per cell.

    Args:
        molecule_info_h5 (str): path to molecule_info.h5 file
        filtered_barcodes_csv (str): path to filtered_barcodes.csv file
        is_targeted (bool, optional): when subsampling to usable reads per cell,
            also restrict to on-target reads. Defaults to False.
        include_bulk (bool, optional): include subsampling based on bulk read counts,
            ignoring the number of cells. Does not apply to targeted GEX libraries.
            Defaults to True.

    Returns:
        list of dict: metadata for subsampling job, produced by make_subsamplings
            and consumed by run_subsampling
    """
    subsamplings: list[SubsamplingDef] = []
    # Get required info from the mol info
    with cr_mc.MoleculeCounter.open(molecule_info_h5, "r") as mc:
        library_info = mc.library_info
        num_cells_per_lib = get_num_cells_per_library(library_info, filtered_barcodes_csv)

        if not is_targeted:
            usable_reads_per_lib = np.array(mc.get_usable_read_pairs_per_library())
        else:
            usable_reads_per_lib = np.array(mc.get_on_target_usable_read_pairs_per_library())

        transcriptomic_reads_per_lib = np.array(mc.get_transcriptomic_read_pairs_per_library())

        if num_cells_per_lib.sum() == 0:
            return subsamplings

        assert library_info is not None
        for library_type in rna_library.sorted_library_types(library_info):
            libraries = [l for l in library_info if l["library_type"] == library_type]
            is_targeted_gex = library_type == rna_library.GENE_EXPRESSION_LIBRARY_TYPE and any(
                rna_library.has_target_set(lib) for lib in libraries
            )
            if is_targeted_gex:
                subsample_types = copy(ALL_SUBSAMPLE_TYPES)
            else:
                subsample_types = copy(ALL_SUBSAMPLE_TYPES)
                if include_bulk:
                    subsample_types += [BULK_SUBSAMPLE_TYPE]

            # Only add raw_cells subsampling if the relevant metrics are available
            try:
                _ = np.array(mc.get_read_pairs_in_filtered_barcodes_per_library())
                subsample_types += [RAW_CELLS_SUBSAMPLE_TYPE]
            except ValueError:
                pass

            for subsample_type in subsample_types:
                if subsample_type == RAW_CELLS_SUBSAMPLE_TYPE:
                    raw_reads_per_lib = np.array(
                        mc.get_read_pairs_in_filtered_barcodes_per_library()
                    )
                else:
                    raw_reads_per_lib = np.array(mc.get_raw_read_pairs_per_library())

                if is_targeted_gex:
                    fixed_depths = SUBSAMPLE_TARGETED_FIXED_DEPTHS
                    subsampling_usable_reads_per_lib = usable_reads_per_lib
                elif subsample_type == BULK_SUBSAMPLE_TYPE:
                    fixed_depths = SUBSAMPLE_BULK_FIXED_DEPTHS
                    subsampling_usable_reads_per_lib = transcriptomic_reads_per_lib
                else:
                    fixed_depths = SUBSAMPLE_FIXED_DEPTHS
                    subsampling_usable_reads_per_lib = usable_reads_per_lib

                subsamplings.extend(
                    make_subsamplings(
                        subsample_type,
                        library_info,
                        library_type,
                        num_cells_per_lib,
                        raw_reads_per_lib,
                        subsampling_usable_reads_per_lib,
                        fixed_depths,
                        SUBSAMPLE_NUM_ADDITIONAL_DEPTHS,
                    )
                )

    return subsamplings


class SubsampleDataDict(TypedDict):
    """The dictionary returned by run_subsampling."""

    umis_per_bc: np.ndarray[tuple[int, int, int], np.dtype[np.int64]]
    features_det_per_bc: np.ndarray[tuple[int, int, int], np.dtype[np.int64]]
    read_pairs_per_bc: np.ndarray[tuple[int, int, int], np.dtype[np.int64]]
    read_pairs: np.ndarray[tuple[int, int], np.dtype[np.int64]]
    umis: np.ndarray[tuple[int, int], np.dtype[np.int64]]
    lib_type_genome_any_reads: np.ndarray[tuple[int, int], np.dtype[np.bool_]]
    total_features_det: np.ndarray[tuple[int, int, int], np.dtype[np.int64]]


def _make_group_keys(is_bulk_subsampling, mol_gem_group: np.ndarray, mol_barcode_idx: np.ndarray):
    """Compute tallies for each barcode."""
    if is_bulk_subsampling:
        # dummy gem groups and barcodes (all 0s) to use for bulk subsampling (one mega cell)
        return (
            np.zeros((mol_gem_group.shape), dtype=mol_gem_group.dtype),
            np.zeros((mol_barcode_idx.shape), dtype=mol_barcode_idx.dtype),
        )
    else:
        return (mol_gem_group, mol_barcode_idx)


def run_subsampling(
    molecule_info_h5,
    subsample_info: Collection[SubsamplingDef],
    filtered_barcodes_csv,
    feature_indices: Collection[int] | None,
    chunk_start,
    chunk_len,
) -> SubsampleDataDict:
    """Runs a subsampling chunk.

    Args:
        molecule_info_h5: Path to a MoleculeCounter file.
        subsample_info: A subsampling produced by make_subsamplings
        filtered_barcodes_csv: A CSV of filtered (cell) barcodes
        feature_indices: indices of filtered features
        chunk_start: integer chunk start
        chunk_len: integer chunk len

    Returns:
        dict: data dictionary
    """
    # NOTE: setting the random seed below, within each task
    with cr_mc.MoleculeCounter.open(molecule_info_h5, "r") as mc:
        # Get cell-associated barcodes
        genomes = get_genomes(mc)

        # Load chunk of relevant data from the mol_info
        chunk = slice(int(chunk_start), int(chunk_start) + int(chunk_len))
        mol_library_idx = mc.get_column_lazy("library_idx")[chunk]
        mol_read_pairs = mc.get_column_lazy("count")[chunk]
        mol_gem_group = mc.get_column_lazy("gem_group")[chunk]
        mol_barcode_idx = mc.get_column_lazy("barcode_idx")[chunk]
        mol_feature_idx: np.ndarray[int, np.dtype[np.int_]] = mc.get_column_lazy("feature_idx")[
            chunk
        ]

        barcodes = mc.get_ref_column("barcodes")

        if feature_indices is not None:
            # subset molecules to targeted panel for each library
            feature_indices = np.array(feature_indices, dtype=np.int64)
            feature_indices.sort()
            mask = np.isin(mol_feature_idx, feature_indices)
            mol_library_idx = mol_library_idx[:][mask]
            mol_read_pairs = mol_read_pairs[:][mask]
            mol_gem_group = mol_gem_group[:][mask]
            mol_barcode_idx = mol_barcode_idx[:][mask]
            mol_feature_idx = mol_feature_idx[:][mask]

        # Give each cell-associated barcode an integer index
        filtered_barcodes = FilteredBarcodes(filtered_barcodes_csv)

        # Give each genome an integer index
        genome_to_int = {g: i for i, g in enumerate(genomes)}
        assert mc.feature_reference is not None
        feature_int_to_genome_int = np.fromiter(
            (
                genome_to_int[f.tags.get(rna_feature_ref.GENOME_FEATURE_TAG, "")]
                for f in mc.feature_reference.feature_defs
            ),
            dtype=int,
        )
        mol_genome_idx: np.ndarray[int, np.dtype[np.int_]] = feature_int_to_genome_int[
            mol_feature_idx
        ]

        assert mc.library_info is not None
        # determine which (library type, genome) pairs have any associated reads
        lib_types = rna_library.sorted_library_types(mc.library_info)
        lib_type_to_int = {l: i for i, l in enumerate(lib_types)}
        lib_idx_to_lib_type_idx: np.ndarray[int, np.dtype[np.int64]] = np.fromiter(
            (lib_type_to_int[lib["library_type"]] for lib in mc.library_info), dtype=np.int64
        )

        lib_type_genome_any_reads: np.ndarray[tuple[int, int], np.dtype[np.bool_]] = np.zeros(
            (len(lib_types), len(genomes)), dtype=bool
        )
        for lib_idx, genome_idx in set(
            zip(mol_library_idx[mol_read_pairs > 0], mol_genome_idx[mol_read_pairs > 0])
        ):
            lib_type_idx = lib_idx_to_lib_type_idx[lib_idx]
            lib_type_genome_any_reads[lib_type_idx, genome_idx] = True

        # Run each subsampling task on this chunk of data
        n_tasks = len(subsample_info)
        n_genomes = len(genomes)
        n_cells = filtered_barcodes.num_cells()
        # total features summed across lib-types
        n_features: int = feature_int_to_genome_int.shape[0]

        umis_per_bc: np.ndarray[tuple[int, int, int], np.dtype[np.int64]] = np.zeros(
            (n_tasks, n_genomes, n_cells), np.int64
        )
        read_pairs_per_bc: np.ndarray[tuple[int, int, int], np.dtype[np.int64]] = np.zeros(
            (n_tasks, n_genomes, n_cells), np.int64
        )
        features_det_per_bc: np.ndarray[tuple[int, int, int], np.dtype[np.int64]] = np.zeros(
            (n_tasks, n_genomes, n_cells), dtype=np.int64
        )
        read_pairs_per_task: np.ndarray[tuple[int, int], np.dtype[np.int64]] = np.zeros(
            (n_tasks, n_genomes), dtype=np.int64
        )
        umis_per_task: np.ndarray[tuple[int, int], np.dtype[np.int64]] = np.zeros(
            (n_tasks, n_genomes), dtype=np.int64
        )
        features_det_per_task: np.ndarray[tuple[int, int, int], np.dtype[np.int64]] = np.zeros(
            (n_tasks, n_genomes, n_features), dtype=np.int64
        )

        for task_idx, task in enumerate(subsample_info):
            print(f"subsampling task: {task}")
            _run_subsample_task(
                task,
                n_features,
                barcodes,
                genomes,
                mol_library_idx,
                mol_feature_idx,
                mol_genome_idx,
                mol_read_pairs,
                mol_gem_group,
                mol_barcode_idx,
                filtered_barcodes,
                umis_per_bc[task_idx],
                read_pairs_per_bc[task_idx],
                features_det_per_bc[task_idx],
                features_det_per_task[task_idx],
                read_pairs_per_task[task_idx],
                umis_per_task[task_idx],
            )

    return {
        "umis_per_bc": umis_per_bc,
        "features_det_per_bc": features_det_per_bc,
        "read_pairs_per_bc": read_pairs_per_bc,
        "read_pairs": read_pairs_per_task,
        "umis": umis_per_task,
        "lib_type_genome_any_reads": lib_type_genome_any_reads,
        "total_features_det": features_det_per_task,
    }


def _run_subsample_task(
    task: SubsamplingDef,
    n_features: int,
    barcodes: np.ndarray[int, np.dtype[np.bytes_]],
    genomes: list[str],
    mol_library_idx: np.ndarray[int, np.dtype[np.int_]],
    mol_feature_idx: np.ndarray[int, np.dtype[np.int_]],
    mol_genome_idx: np.ndarray[int, np.dtype[np.int_]],
    mol_read_pairs: np.ndarray[int, np.dtype[np.int_]],
    mol_gem_group: np.ndarray[int, np.dtype[np.int_]],
    mol_barcode_idx: np.ndarray[int, np.dtype[np.int_]],
    filtered_barcodes: FilteredBarcodes,
    umis_per_bc: np.ndarray[tuple[int, int], np.dtype[np.int64]],
    read_pairs_per_bc: np.ndarray[tuple[int, int], np.dtype[np.int64]],
    features_det_per_bc: np.ndarray[tuple[int, int], np.dtype[np.int64]],
    features_det_per_task: np.ndarray[tuple[int, int], np.dtype[np.int64]],
    read_pairs_per_task: np.ndarray[int, np.dtype[np.int64]],
    umis_per_task: np.ndarray[int, np.dtype[np.int64]],
):
    # Set random seed
    np.random.seed(1)
    # Per-library subsampling rates
    rates_per_library = np.array(task["library_subsample_rates"], dtype=np.float64)
    is_bulk_subsampling = task["subsample_type"] == BULK_SUBSAMPLE_TYPE
    is_raw_cells_subsampling = task["subsample_type"] == RAW_CELLS_SUBSAMPLE_TYPE

    if np.count_nonzero(rates_per_library) == 0:
        return

    mol_rate = rates_per_library[mol_library_idx]
    if np.isnan(mol_rate).any():
        # subsampling rates contain NaNs i.e. library_subsample_rates contains Nones
        # subsampling cannot be done
        return

    if np.any(mol_rate < 0) or np.any(mol_rate > 1):
        raise ValueError("subsampling probabilities cannot be < 0 or > 1")

    for (gg, bc_idx), (feature_idx, genome_idx, read_pairs) in cr_utils.numpy_groupby(
        values=(
            mol_feature_idx,
            mol_genome_idx,
            np.random.binomial(mol_read_pairs, mol_rate),
        ),
        keys=_make_group_keys(is_bulk_subsampling, mol_gem_group, mol_barcode_idx),
    ):
        barcode = cr_utils.format_barcode_seq(barcodes[bc_idx], gg)

        for this_genome_idx, genome in enumerate(genomes):
            is_cell_barcode = filtered_barcodes.contains(barcode=barcode, genome=genome)

            # Skip entries from non-cell barcodes when doing raw cells subsampling
            if not is_cell_barcode and is_raw_cells_subsampling:
                continue

            umis = np.flatnonzero((read_pairs > 0) & (genome_idx == this_genome_idx))
            this_genome_read_pairs: int = np.sum(read_pairs[genome_idx == this_genome_idx])

            # Tally UMIs and number of features detected
            if is_bulk_subsampling:
                umis_per_bc[this_genome_idx, :] = len(umis)
                read_pairs_per_bc[this_genome_idx, :] = this_genome_read_pairs
                # this number doesn't make sense for bulk subsampling, set to 0
                features_det_per_bc[this_genome_idx, :] = 0
                features_det_per_task[this_genome_idx, :] = np.bincount(
                    feature_idx[umis], minlength=n_features
                )
            elif is_cell_barcode:
                # This is a cell-associated barcode for this genome
                cell_idx = filtered_barcodes.index_of_barcode(barcode)
                umis_per_bc[this_genome_idx, cell_idx] = len(umis)
                read_pairs_per_bc[this_genome_idx, cell_idx] = this_genome_read_pairs
                features_det_per_bc[this_genome_idx, cell_idx] = np.count_nonzero(
                    np.bincount(feature_idx[umis])
                )
                # add total features detected in cells
                features_det_per_task[this_genome_idx, :] += np.bincount(
                    feature_idx[umis], minlength=n_features
                )

            # Tally numbers for duplicate fraction
            read_pairs_per_task[this_genome_idx] += this_genome_read_pairs
            umis_per_task[this_genome_idx] += len(umis)


def make_metric_name(
    name,
    library_type,
    genome,
    ss_type,
    ss_depth,
    target_mode,
):
    lt_prefix = rna_library.get_library_type_metric_prefix(library_type)
    assert target_mode in SUBSAMPLE_TARGET_MODES
    target_suffix = "" if target_mode is None else ("_" + target_mode)
    return f"{lt_prefix}{genome}_{ss_type}_{ss_depth}_{name}{target_suffix}"


def compute_dup_frac(read_pairs, umis):
    return tk_stats.robust_divide(read_pairs - umis, read_pairs) if read_pairs > 0 else 0.0


def join_metrics(metrics):
    """Join together a list of metric dicts (each value is a numpy vector).

    :param metrics: list of metrics
    :return: joined dictionary
    """
    if len(metrics) == 0:
        return {}
    with open(metrics[0], "rb") as f:
        data = pickle.load(f)
    for m in metrics[1:]:
        with open(m, "rb") as f:
            chunk_data = pickle.load(f)
            for k, v in chunk_data.items():
                data[k] += v
    return data


def calculate_subsampling_metrics(
    data, molecule_info_h5, filtered_barcodes_csv, subsample_info, target_mode
):
    """Calculate subsampling metrics (summary) from a joined data structure from run_subsampling.

    :param data: A merged dictionary of data from run_subsampling
    :param molecule_info_h5: path to a MoleculeInfo file
    :param filtered_barcodes_csv: path to a list of cell barcodes
    :param subsample_info: subsampling info produced by construct_all_subsamplings
    :param target_mode: String of target mode for metrics suffix.
        Must be one of constants.SUBSAMPLE_TARGET_MODES
    :return: dict (JSON) metrics
    """
    # Compute metrics for each subsampling rate
    summary = {}

    with cr_mc.MoleculeCounter.open(molecule_info_h5, "r") as mc:
        genomes = get_genomes(mc)
        assert mc.library_info is not None
        lib_types = rna_library.sorted_library_types(mc.library_info)
        lib_type_map = {lt: idx for (idx, lt) in enumerate(lib_types)}
        assert mc.feature_reference is not None
        num_target_features = mc.feature_reference.count_target_feature_indices()

    filtered_barcodes = FilteredBarcodes(filtered_barcodes_csv)

    for i, task in enumerate(subsample_info):
        lib_type = task["library_type"]
        lib_type_idx = lib_type_map[lib_type]
        ss_type = task["subsample_type"]
        ss_depth = task["target_read_pairs_per_cell"]

        if rna_library.has_genomes(lib_type):
            genome_ints = list(range(data["umis_per_bc"].shape[1]))
        else:
            genome_ints = [0]

        # Per-genome metrics
        for g in genome_ints:
            if not data["lib_type_genome_any_reads"][lib_type_idx, g]:
                continue
            genome = genomes[g]

            # Only compute on cell-associated barcodes for this genome.
            # This only matters when there are multiple genomes present.
            cell_inds = filtered_barcodes.sorted_barcode_indices(genome)

            mean_reads_per_cell = np.mean(data["read_pairs_per_bc"][i, g, cell_inds])
            summary[
                make_metric_name(
                    "subsampled_filtered_bcs_mean_read_counts",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = mean_reads_per_cell

            median_reads_per_cell = np.median(data["read_pairs_per_bc"][i, g, cell_inds])
            summary[
                make_metric_name(
                    "subsampled_filtered_bcs_median_read_counts",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = median_reads_per_cell

            median_umis_per_cell = np.median(data["umis_per_bc"][i, g, cell_inds])
            summary[
                make_metric_name(
                    "subsampled_filtered_bcs_median_counts",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = median_umis_per_cell

            mean_umis_per_cell = np.mean(data["umis_per_bc"][i, g, cell_inds])
            summary[
                make_metric_name(
                    "subsampled_filtered_bcs_mean_counts",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = mean_umis_per_cell

            if ss_type == BULK_SUBSAMPLE_TYPE:
                # if bulk, count total number of genes detected at that ss_rate
                median_features_per_cell = np.count_nonzero(data["total_features_det"][i, g])
                mean_features_per_cell = median_features_per_cell
            else:
                median_features_per_cell = np.median(data["features_det_per_bc"][i, g, cell_inds])
                mean_features_per_cell = np.mean(data["features_det_per_bc"][i, g, cell_inds])
            summary[
                make_metric_name(
                    "subsampled_filtered_bcs_median_unique_genes_detected",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = median_features_per_cell
            summary[
                make_metric_name(
                    "subsampled_filtered_bcs_mean_unique_genes_detected",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = mean_features_per_cell

            dup_frac = compute_dup_frac(data["read_pairs"][i, g], data["umis"][i, g])
            summary[
                make_metric_name(
                    "subsampled_duplication_frac",
                    lib_type,
                    genome,
                    ss_type,
                    ss_depth,
                    target_mode,
                )
            ] = dup_frac

            # fraction of targeted genes detected
            if target_mode is not None and target_mode == ON_TARGET_SUBSAMPLE:
                num_target_features_detected = np.count_nonzero(data["total_features_det"][i, g])
                summary[
                    make_metric_name(
                        "subsampled_frac_targeted_genes_detected",
                        lib_type,
                        genome,
                        ss_type,
                        ss_depth,
                        target_mode,
                    )
                ] = np.divide(num_target_features_detected, float(num_target_features))

        # Whole-dataset duplication frac
        all_read_pairs = np.sum(data["read_pairs"][i, :])
        all_umis = np.sum(data["umis"][i, :])
        dup_frac = compute_dup_frac(all_read_pairs, all_umis)
        summary[
            make_metric_name(
                "subsampled_duplication_frac",
                lib_type,
                rna_library.MULTI_REFS_PREFIX,
                ss_type,
                ss_depth,
                target_mode,
            )
        ] = dup_frac
    return summary


def parse_subsample_summary(summary, suffix=None, prefix=None, include_mean=True) -> pd.DataFrame:
    """Convenience function.

    Args:
        summary (dict): loaded JSON, as produced by calculate_subsampling_metrics
        suffix (str, optional): subsampling metric suffix
        prefix (sre, optional): subsampling metric prefix
        include_mean (bool, optional): Include mean UMIs and genes per barcode. Defaults to True.

    Returns:
        pd.DataFrame: subsampling metrics
    """
    munged_data = []
    names = [
        "Mean reads per cell",
        "Median reads per cell",
        "Median UMIs per cell",
        "Median genes per cell",
    ]
    if include_mean:
        names.extend(["Mean UMIs per cell", "Mean genes per cell"])
    metrics = [
        "subsampled_filtered_bcs_mean_read_counts",
        "subsampled_filtered_bcs_median_read_counts",
        "subsampled_filtered_bcs_median_counts",
        "subsampled_filtered_bcs_median_unique_genes_detected",
    ]
    if include_mean:
        metrics.extend(
            [
                "subsampled_filtered_bcs_mean_counts",
                "subsampled_filtered_bcs_mean_unique_genes_detected",
            ]
        )
    if suffix is not None:
        metrics = [f"{x}_{suffix}" for x in metrics]
    for name, metric in zip(names, metrics):
        keys = [x for x in summary.keys() if metric in x and MAPPED_SUBSAMPLE_TYPE in x]
        if prefix is not None:
            keys = [x for x in keys if x.startswith(prefix)]
        ss_target = [int(x.split(metric)[0].rsplit("_", 2)[1]) for x in keys]
        data = [summary[x] for x in keys]
        for target, val in zip(ss_target, data):
            munged_data.append([target, val, name])
    df = pd.DataFrame(munged_data, columns=["target", "val", "name"])
    df = df.pivot(index="target", columns="name", values="val").astype(int)
    # remove rows that are all 0s
    df = df[(df != 0).any(axis=1)]
    # handle the case where all subsamplings were skipped due to super low coverage
    if df.shape == (0, 0):
        df = pd.DataFrame(columns=names)
    df.index.name = "Target mean reads per cell"
    if prefix is not None:
        df["genome"] = prefix
    return df


def get_selected_task_indices(
    ss_tasks, filter_library_type=None, filter_subsample_type=None, filter_target_read_pairs=None
) -> list[int]:
    """Parse through subsample info to extract task indices of interest that match filters.

    Args:
        ss_tasks (list): list of subsamplings, as generated by construct_all_subsamplings
        filter_library_type (list): list of library types to filter by (or None)
        filter_subsample_type (list): list of subsample types to filter by (or None)
        filter_target_read_pairs (list): list of subsample depths to filter by (or None)

    Returns:
        task_indices (list): list of indices that correspond to the tasks of interest
    """
    task_indices = []
    for task_idx, task_ss_info in enumerate(ss_tasks):
        if (
            filter_library_type is not None
            and task_ss_info["library_type"] not in filter_library_type
        ):
            continue
        if (
            filter_subsample_type is not None
            and task_ss_info["subsample_type"] not in filter_subsample_type
        ):
            continue
        if (
            filter_target_read_pairs is not None
            and task_ss_info["target_read_pairs_per_cell"] not in filter_target_read_pairs
        ):
            continue
        task_indices.append(task_idx)
    return task_indices
