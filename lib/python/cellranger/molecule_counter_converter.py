#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Function to convert V2 to V3 files.

Separated out from the main molecule_counter to avoid dependencies.
"""

from __future__ import annotations

import random
from collections import OrderedDict

import h5py
import numpy as np

import cellranger.barcodes.utils as bc_utils
import cellranger.cell_calling_helpers
import cellranger.rna.library as rna_library
import tenkit.seq as tk_seq
from cellranger import molecule_counter as mc
from cellranger.feature_ref import FeatureDef, FeatureReference
from cellranger.molecule_counter import (
    BARCODE_IDX_COL_NAME,
    COUNT_COL_NAME,
    FEATURE_IDX_COL_NAME,
    GEM_GROUP_COL_NAME,
    LIBRARY_IDX_COL_NAME,
    MOLECULE_INFO_COLUMNS,
    UMI_COL_NAME,
    UMI_TYPE_COL_NAME,
    V3_METRICS_GROUP_NAME,
    create_dataset,
    set_file_version,
)

BC_LENGTH_METRIC_V2 = "chemistry_barcode_read_length"
BC_WHITELIST_METRIC_V2 = "chemistry_barcode_whitelist"


def decompress_barcode_seq(x, barcode_length, bits=64) -> bytes:
    """Decompress a 2 bit encode barcode.

    Args:
        x:
        barcode_length:
        bits:
    """
    x = np.uint64(x)
    assert barcode_length <= (bits // 2 - 1)
    if x & (1 << (bits - 1)):
        return b"N" * barcode_length
    result = bytearray(barcode_length)
    mask = np.uint64(0b11)
    for i in range(barcode_length):
        result[(barcode_length - 1) - i] = ord(tk_seq.NUCS[x & mask])
        x = x >> np.uint64(2)
    return bytes(result)


def convert_v2_to_v4(v2_mole_info_h5, out_v4_mole_info_h5):
    """Given the input v2 molecule info h5 file, convert it into v4 file."""

    def build_feature_ref(gene_ids, gene_names, genome_index):
        """Build a feature ref.

        Args:
            gene_ids:
            gene_names:
            genome_index:
        """
        feature_defs = []
        if len(genome_index) == 1:
            genome = next(iter(genome_index.keys()))
            for idx, (gene_id, gene_name) in enumerate(zip(gene_ids, gene_names)):
                feature_defs.append(
                    FeatureDef(
                        index=idx,
                        id=gene_id,
                        name=gene_name,
                        feature_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
                        tags={"genome": genome},
                    )
                )
        else:
            for idx, (gene_id, gene_name) in enumerate(zip(gene_ids, gene_names)):
                genome = gene_id.split("_")[0]
                feature_defs.append(
                    FeatureDef(
                        index=idx,
                        id=gene_id,
                        name=gene_name,
                        feature_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
                        tags={"genome": genome},
                    )
                )

        return FeatureReference(feature_defs, ["genome"])

    random.seed(0)
    np.random.seed(0)

    v2_mc_in = h5py.File(v2_mole_info_h5, "r")
    v2_metrics = mc.get_v2_metrics(v2_mole_info_h5)

    v2_genome_ids = v2_mc_in["genome_ids"]
    v2_genome_name_to_index = {g: i for i, g in enumerate(v2_genome_ids)}

    # Feature Ref
    new_feature_ref = build_feature_ref(
        v2_mc_in["gene_ids"], v2_mc_in["gene_names"], v2_genome_name_to_index
    )

    # barcode whitelist
    barcode_length = v2_metrics[BC_LENGTH_METRIC_V2]
    barcode_whitelist = bc_utils.load_barcode_whitelist(v2_metrics[BC_WHITELIST_METRIC_V2])
    barcode_to_idx = OrderedDict((k, i) for i, k in enumerate(barcode_whitelist))

    # <-> genome information goes into feature_idx in v3
    v2_genomes = np.asarray(v2_mc_in["genome"], dtype=np.uint8)
    # <-> feature_idx in v3
    v2_gene = np.asarray(v2_mc_in["gene"], dtype=mc.MOLECULE_INFO_COLUMNS["feature_idx"])
    v2_conf_mapped_reads = np.asarray(
        v2_mc_in["reads"], dtype=mc.MOLECULE_INFO_COLUMNS["count"]
    )  # <-> count in v3
    # <-> transit into barcode_idx in v3
    v2_barcodes = np.asarray(v2_mc_in["barcode"], dtype=np.uint64)
    v2_umis = np.asarray(v2_mc_in["umi"], dtype=mc.MOLECULE_INFO_COLUMNS["umi"])  # <-> umi in v3

    # these are required to extrapolate transcriptomic reads
    v2_nonconf_mapped_reads = np.asarray(
        v2_mc_in["nonconf_mapped_reads"], dtype=mc.MOLECULE_INFO_COLUMNS["count"]
    )
    v2_unmapped_reads = np.asarray(
        v2_mc_in["unmapped_reads"], dtype=mc.MOLECULE_INFO_COLUMNS["count"]
    )
    v2_total_reads = v2_conf_mapped_reads + v2_nonconf_mapped_reads + v2_unmapped_reads
    del v2_nonconf_mapped_reads, v2_unmapped_reads

    library_info, chunks = mc.get_v2_library_info_and_chunks(v2_mc_in)
    barcode_info_genomes, barcode_info_pass_filter = [], []
    barcode_idx_list, feature_idx_list, library_idx_list = [], [], []
    gem_group_list, count_list, umi_list = [], [], []

    v2_metrics[mc.LIBRARIES_METRIC] = {}

    def convert_mc_chunk(gem_group, chunk_start, chunk_len, lib_idx):
        # extrapolated because total conf_mapped_reads are unavailable,
        #   only used in SUBSAMPLE_READS for internal metrics
        chunk_end = chunk_start + chunk_len
        total_reads_in_valid_barcodes = v2_total_reads[chunk_start:chunk_end].sum()
        conf_reads_in_valid_barcodes = v2_conf_mapped_reads[chunk_start:chunk_end].sum()
        total_reads_in_gem_group: int = v2_metrics[mc.GEM_GROUPS_METRIC][gem_group]["total_reads"]
        extrapolated_transcriptomic_reads = (
            float(conf_reads_in_valid_barcodes)
            / total_reads_in_valid_barcodes
            * total_reads_in_gem_group
        )
        extrapolated_transcriptomic_reads = int(round(extrapolated_transcriptomic_reads))

        # per library, raw_read_pairs and usable_read_pairs info
        v2_metrics[mc.LIBRARIES_METRIC][str(lib_idx)] = {
            mc.FEATURE_READS_METRIC: extrapolated_transcriptomic_reads,
            mc.USABLE_READS_METRIC: v2_metrics[mc.GEM_GROUPS_METRIC][gem_group][
                "conf_mapped_filtered_bc_reads"
            ],
            mc.TOTAL_READS_METRIC: total_reads_in_gem_group,
        }

        recovered_cells = v2_metrics[mc.GEM_GROUPS_METRIC][gem_group].get(
            mc.GG_RECOVERED_CELLS_METRIC, None
        )
        force_cells = v2_metrics[mc.GEM_GROUPS_METRIC][gem_group].get(
            mc.GG_FORCE_CELLS_METRIC, None
        )

        genomes_for_gem_group = v2_genomes[chunk_start:chunk_end]
        bcs_for_gem_group = v2_barcodes[chunk_start:chunk_end]
        reads_for_gem_group = v2_conf_mapped_reads[chunk_start:chunk_end]
        gene_for_gem_group = v2_gene[chunk_start:chunk_end]
        umis_for_gem_group = v2_umis[chunk_start:chunk_end]

        for genome_id in v2_genome_ids:
            g_idx = v2_genome_name_to_index[genome_id]
            genome_indices = genomes_for_gem_group == g_idx

            if genome_indices.sum() == 0:
                # edge case - there's no data for this genome (e.g. empty sample, false barnyard sample, or nothing confidently mapped)
                continue

            bcs_for_genome = bcs_for_gem_group[genome_indices]
            reads_for_genome = reads_for_gem_group[genome_indices]
            gene_for_genome = gene_for_gem_group[genome_indices]
            umis_for_genome = umis_for_gem_group[genome_indices]

            # only count UMIs with at least one conf mapped read
            umi_conf_mapped_to_genome = reads_for_genome > 0
            bc_breaks = bcs_for_genome[1:] - bcs_for_genome[:-1]
            # first row is always a break
            bc_breaks = np.concatenate(([1], bc_breaks))
            bc_break_indices = np.nonzero(bc_breaks)[0]
            unique_bcs = bcs_for_genome[bc_break_indices]
            umis_per_bc = np.add.reduceat(umi_conf_mapped_to_genome, bc_break_indices)

            if force_cells is not None:
                (
                    top_bc_indices,
                    _,
                    _,
                ) = cellranger.cell_calling_helpers.filter_cellular_barcodes_fixed_cutoff(
                    umis_per_bc, force_cells
                )
            else:
                (
                    top_bc_indices,
                    _,
                    _,
                ) = cellranger.cell_calling_helpers.filter_cellular_barcodes_ordmag(
                    umis_per_bc, recovered_cells
                )

            # barcode info
            barcode_seq_to_idx = {
                b: barcode_to_idx[decompress_barcode_seq(b, barcode_length)] for b in unique_bcs
            }

            barcode_info_genomes.append(genome_id)
            for b in unique_bcs[top_bc_indices]:
                barcode_info_pass_filter.append((barcode_seq_to_idx[b], lib_idx, g_idx))

            # data
            barcode_idx_list.append(np.vectorize(barcode_seq_to_idx.get)(bcs_for_genome))
            count_list.append(reads_for_genome)
            gem_group_list.append(
                np.full(
                    reads_for_genome.shape[0],
                    gem_group,
                    dtype=mc.MOLECULE_INFO_COLUMNS["gem_group"],
                )
            )
            library_idx_list.append(
                np.full(
                    reads_for_genome.shape[0],
                    lib_idx,
                    dtype=mc.MOLECULE_INFO_COLUMNS["library_idx"],
                )
            )
            feature_idx_list.append(gene_for_genome)
            umi_list.append(umis_for_genome)

    # Each gem group is a library, and we'll iterate over the chunks in it
    for gem_group, chunk_start, chunk_len, lib_idx in chunks:
        convert_mc_chunk(gem_group, chunk_start, chunk_len, lib_idx)

    new_barcode_info = mc.BarcodeInfo(
        pass_filter=np.array(barcode_info_pass_filter, dtype=mc.BARCODE_INFO_DTYPES["pass_filter"]),
        genomes=barcode_info_genomes,
    )

    with mc.MoleculeCounter.open(
        out_v4_mole_info_h5,
        "w",
        feature_ref=new_feature_ref,
        barcodes=barcode_whitelist,
        library_info=library_info,
        barcode_info=new_barcode_info,
    ) as out_mc:
        out_mc.append_column(BARCODE_IDX_COL_NAME, np.concatenate(barcode_idx_list))
        out_mc.append_column(COUNT_COL_NAME, np.concatenate(count_list))
        out_mc.append_column(FEATURE_IDX_COL_NAME, np.concatenate(feature_idx_list))
        out_mc.append_column(GEM_GROUP_COL_NAME, np.concatenate(gem_group_list))
        out_mc.append_column(UMI_COL_NAME, np.concatenate(umi_list))
        # library_idx is the same as gem_group_list
        out_mc.append_column(LIBRARY_IDX_COL_NAME, np.concatenate(library_idx_list))
        out_mc.set_all_metrics(v2_metrics)
        # Remove the umi_type column and set the version
        del out_mc.h5[UMI_TYPE_COL_NAME]
        set_file_version(out_mc, 4)


def convert_v3_to_v4(v3_molecule_info: mc.MoleculeCounter):
    """Converts a V3 or to V4.

    Args:
        v3_molecule_info: A V3 molecule info that we will convert by replacing the python
        pickle with a dictionary

    Returns: None
    """
    assert isinstance(v3_molecule_info, mc.MoleculeCounter)
    assert v3_molecule_info.file_version == 3
    metrics = v3_molecule_info.get_all_metrics()
    del v3_molecule_info.h5[V3_METRICS_GROUP_NAME]
    v3_molecule_info.set_all_metrics(metrics)
    set_file_version(v3_molecule_info, 4)


def convert_v4_to_v5(v4_molecule_info: mc.MoleculeCounter):
    """Converts a V4 file to v5.

    Adds a `umi_type` column,
    all previous files were txomic, e.g. 1
    """
    n = v4_molecule_info.nrows()
    create_dataset(v4_molecule_info, UMI_TYPE_COL_NAME)
    for _, chnk_len in v4_molecule_info.get_chunks(2000000, preserve_boundaries=False):
        v4_molecule_info.append_column(
            UMI_TYPE_COL_NAME, np.ones(chnk_len, dtype=MOLECULE_INFO_COLUMNS[UMI_TYPE_COL_NAME])
        )

    assert len(v4_molecule_info.columns[UMI_TYPE_COL_NAME]) == n, "Column `umi_type` does not match"
    set_file_version(v4_molecule_info, 5)


def upgrade_file(molecule_info: mc.MoleculeCounter):
    """Upgrades a molecule info to the latest format.

    Args:
        molecule_info: a molecule info
    """
    version = molecule_info.file_version
    # Versions <3 should be converted in check_molecule_info, and we don't trim in spatial
    assert version >= 3

    if version == 3:
        convert_v3_to_v4(molecule_info)
        version = molecule_info.file_version
    if version == 4:
        convert_v4_to_v5(molecule_info)
        version = molecule_info.file_version
    if version == 5:
        # Version 6 has an optional probeset and probe_idx data, but are not being placed in here.
        set_file_version(molecule_info, 6)
