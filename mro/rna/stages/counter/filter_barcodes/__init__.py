#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Identify partitions that contain cells."""

from __future__ import annotations

import dataclasses
import json
import os
from collections import defaultdict

import martian
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.cell_calling_helpers as helpers
import cellranger.chemistry as cr_chemistry
import cellranger.constants as cr_constants
import cellranger.feature.antibody.analysis as ab_utils
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
import cellranger.rna.matrix as rna_matrix
import tenkit.safe_json as tk_safe_json
from cellranger.cell_barcodes import InvalidCellBarcode
from cellranger.cell_calling_helpers import (
    FilterMethod,
    HighOccupancyGemSummary,
    MetricGroups,
    get_filter_method_from_string,
    get_filter_method_name,
)
from cellranger.csv_utils import (
    combine_csv,
    write_filtered_barcodes,
    write_isotype_normalization_csv,
)
from cellranger.fast_utils import (  # pylint: disable=no-name-in-module,unused-import
    MultiGraph,
    load_filtered_bcs_groups,
    save_filtered_bcs_groups,
)
from cellranger.feature.antibody.analysis import FRACTION_CORRECTED_READS, FRACTION_TOTAL_READS
from cellranger.metrics import BarcodeFilterResults
from cellranger.read_level_multiplexing import get_overhang_bc_defn, get_sample_tag_barcodes
from cellranger.reference_paths import get_reference_genomes
from cellranger.targeted.rtl_multiplexing import get_probe_bc_defn
from cellranger.targeted.simple_utils import determine_targeting_type_from_csv
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC

PROBE_BC_SAMPLE_ID = "id"
PROBE_BC_SEQS = "sequence"
PROBE_BC_OFFSET = "offset"
PROBE_BC_LEN = "length"

# FBc control tag name
FBC_CONTROL_TAG_NAME = "isotype_control"

# FBc Scale factors
FBC_SCALE_FACTOR = 10_000

# Number of FBc s needed in a count matrix to normalise
MINIMUM_FBCS_TO_NORMALIZE = 10

np.random.seed(0)


# TODO: We can get rid of all the gem-group logic here,
# because we run this only on data from a single gem-group

__MRO__ = """
struct ProbeBCDef(
    string   id,
    string[] sequence,
    int      offset,
    int      length,
)

struct CellCallingParam(
    int      per_gem_well,
    map<int> per_sample,
)


struct CellCalling(
    CellCallingParam recovered_cells,
    CellCallingParam force_cells,
    CellCallingParam emptydrops_minimum_umis,
    CellCallingParam global_minimum_umis,
    json             cell_barcodes,
    string           override_mode,
    string[]         override_library_types,
    bool             disable_ab_aggregate_detection,
    bool             disable_high_occupancy_gem_detection,
)

stage FILTER_BARCODES(
    in  map<ChemistryDef> chemistry_defs,
    in  string            sample_id,
    in  h5                matrices_h5,
    in  csv               barcode_correction_csv,
    in  bool              is_antibody_only,
    in  path              reference_path,
    in  int[]             gem_groups,
    in  CellCalling       config,
    in  csv               target_set,
    in  json              multi_graph,
    in  csv               per_barcode_metrics,
    in  bool              is_spatial,
    out json              summary,
    out csv               filtered_barcodes,
    out csv               aggregate_barcodes,
    out h5                filtered_matrices_h5,
    out path              filtered_matrices_mex,
    out csv               nonambient_calls,
    out csv               mitochondrial_summary,
    out csv               isotype_normalization_factors,
    src py                "stages/counter/filter_barcodes",
) split (
    in  ProbeBCDef        probe_bc_def,
    out json              filtered_metrics_groups,
    out bincode           filtered_bcs_groups,
    out csv               co_mitochondrial_summary,
) using (
    mem_gb   = 8,
    volatile = strict,
)
"""


def split(args):
    chemistry_def = cr_chemistry.get_primary_chemistry_def(args.chemistry_defs)
    chemistry_barcode = chemistry_def["barcode"]

    mem_gb = 2 + cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrices_h5, scale=6)

    if args.multi_graph is None:
        return {
            "chunks": [{"probe_bc_def": None, "__mem_gb": mem_gb}],
            "join": {
                "__mem_gb": mem_gb,
            },
        }
    else:
        # Read in the multi graph
        multi_graph = MultiGraph.from_path(args.multi_graph)
        if multi_graph.is_rtl_multiplexed():
            offset, length = get_probe_bc_defn(chemistry_barcode)
        elif multi_graph.is_oh_multiplexed():
            offset, length = get_overhang_bc_defn(chemistry_barcode)
        else:
            offset, length = None, None

        chunks = []
        # Perform cell-calling per sample when input is multiplexed
        if length is not None:
            probe_bcs = get_sample_tag_barcodes(multi_graph, chemistry_barcode)
            for sample_id, probe_bc_seqs in probe_bcs.items():
                chunks.append(
                    {
                        "probe_bc_def": {
                            "id": sample_id,
                            "sequence": probe_bc_seqs,
                            "offset": offset,
                            "length": length,
                        },
                        "__mem_gb": mem_gb,
                    }
                )
        else:
            chunks.append({"probe_bc_def": None, "__mem_gb": mem_gb})
        return {
            "chunks": chunks,
            "join": {
                "__mem_gb": mem_gb,
            },
        }


def main(args, outs):
    filter_barcodes(args, outs)


def join(args, outs, chunk_defs, chunk_outs):
    chemistry_def = cr_chemistry.get_primary_chemistry_def(args.chemistry_defs)

    filtered_metrics_groups = parse_chunked_metrics(chunk_outs)
    filtered_bcs, genome_filtered_bcs = parse_chunked_filtered_bcs(chunk_outs)
    summary = parse_chunked_summary(chunk_outs)

    # high occupancy GEM filtering to mitigate any GEMs that appear overloaded
    # in initial cell calls. Only applies if input has probe bcs and not disabled.
    config = helpers.CellCalling(args.config)
    probe_bc_offset, _ = get_probe_bc_defn(chemistry_def["barcode"])
    if probe_bc_offset and not config.disable_high_occupancy_gem_detection:
        (
            filtered_bcs,
            genome_filtered_bcs,
            high_occupancy_gem_metrics,
        ) = helpers.remove_bcs_from_high_occupancy_gems(
            filtered_bcs, genome_filtered_bcs, args.per_barcode_metrics, probe_bc_offset
        )
        high_occupancy_gem_metrics = dataclasses.asdict(high_occupancy_gem_metrics)
    else:
        high_occupancy_gem_metrics = dataclasses.asdict(HighOccupancyGemSummary())

    # Select cell-associated barcodes
    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrices_h5)
    filtered_matrix = raw_matrix.select_barcodes_by_seq(filtered_bcs)

    # If non-spatial targeted ensure all filtered barcodes have nonzero targeted UMI counts
    is_targeted = raw_matrix.feature_ref.has_target_features()
    if (not args.is_spatial) and is_targeted:
        filtered_bcs, genome_filtered_bcs = helpers.remove_cells_with_zero_targeted_counts(
            raw_matrix, filtered_bcs, genome_filtered_bcs, summary
        )
        filtered_matrix = raw_matrix.select_barcodes_by_seq(filtered_bcs)

    # subset the filtered matrix to only targeted genes
    if is_targeted:
        target_features = raw_matrix.feature_ref.get_target_feature_indices()
        filtered_matrix = filtered_matrix.remove_genes_not_on_list(target_features)

    # Write the filtered barcodes file
    write_filtered_barcodes(outs.filtered_barcodes, genome_filtered_bcs)

    # Normalise matrix on ABs in spatial
    libraries_in_h5 = filtered_matrix.get_library_types()
    # Check if there is a spatial FBc run in the matrix and if so normalise by the
    # control FBc counts
    normalisation_metrics = {"number_of_overflows_in_normalization": 0}
    outs.isotype_normalization_factors = None
    # This accounts for when the --no-libraries flag is passed in but a feature ref is defined.
    # An empty Antibody matrix will be passed in and we don't want to go through the waste of normalizing
    # when it's not needed.
    antibody_counts = sum(
        raw_matrix.view()
        .select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
        .get_counts_per_bc()
    )
    if (
        rna_library.ANTIBODY_LIBRARY_TYPE in libraries_in_h5
        and args.is_spatial
        and antibody_counts != 0
    ):
        # Get names of negative controls
        fbc_control_names = [
            x.id
            for x in filtered_matrix.view().feature_ref.feature_defs
            if x.tags.get(FBC_CONTROL_TAG_NAME, "FALSE") == "TRUE"
        ]
        number_of_fbcs_in_matrix = filtered_matrix.get_count_of_feature_type(
            rna_library.ANTIBODY_LIBRARY_TYPE
        )

        if fbc_control_names or number_of_fbcs_in_matrix > MINIMUM_FBCS_TO_NORMALIZE:
            # Normalise the filtered matrix
            normalisation_metrics["number_of_overflows_in_normalization"] = (
                filtered_matrix.normalise_library_type_by_control_features(
                    rna_library.ANTIBODY_LIBRARY_TYPE, FBC_SCALE_FACTOR, fbc_control_names
                )
            )

            # Compute normalisation factor for all barcodes in raw matrix
            if fbc_control_names:
                normalization_factor = FBC_SCALE_FACTOR / (
                    raw_matrix.view().select_features_by_ids(fbc_control_names).get_counts_per_bc()
                    + 1
                )
                normalisation_metrics["ANTIBODY_isotype_normalized"] = "On"
            else:
                normalization_factor = FBC_SCALE_FACTOR / (
                    raw_matrix.view()
                    .select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
                    .get_counts_per_bc()
                    + 1
                )
                normalisation_metrics["ANTIBODY_isotype_normalized"] = "Off"

            outs.isotype_normalization_factors = martian.make_path(
                "isotype_normalization_factors.csv"
            ).decode()
            filtered_bcs_set = set(filtered_bcs)
            in_tissue = [x in filtered_bcs_set for x in raw_matrix.bcs]
            write_isotype_normalization_csv(
                outs.isotype_normalization_factors, raw_matrix.bcs, in_tissue, normalization_factor
            )

    matrix_attrs = cr_matrix.make_matrix_attrs_count(
        args.sample_id, args.gem_groups, chemistry_def[cr_chemistry.CHEMISTRY_DESCRIPTION_FIELD]
    )
    filtered_matrix.save_h5_file(
        outs.filtered_matrices_h5,
        extra_attrs=matrix_attrs,
        sw_version=martian.get_pipelines_version(),
    )

    rna_matrix.save_mex(
        filtered_matrix, outs.filtered_matrices_mex, martian.get_pipelines_version()
    )

    parse_filtered_bcs_method(filtered_metrics_groups, summary)

    genomes = get_reference_genomes(args.reference_path, args.target_set)
    summary = helpers.combine_initial_metrics(
        genomes, filtered_metrics_groups, genome_filtered_bcs, summary
    )

    # Add keys in from high occupancy GEM filtering in if it was performed
    summary.update(high_occupancy_gem_metrics)

    # Add normalisation metrics
    summary.update(normalisation_metrics)

    # Write metrics json
    with open(ensure_binary(outs.summary), "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)

    chunk_aggregate_barcodes = [
        chunk.aggregate_barcodes
        for chunk in chunk_outs
        if chunk.aggregate_barcodes is not None and os.path.exists(chunk.aggregate_barcodes)
    ]
    if chunk_aggregate_barcodes:
        combine_csv(chunk_aggregate_barcodes, outs.aggregate_barcodes)
    else:
        outs.aggregate_barcodes = None

    chunk_nonambient_calls = [
        chunk.nonambient_calls
        for chunk in chunk_outs
        if chunk.nonambient_calls is not None and os.path.exists(chunk.nonambient_calls)
    ]
    combine_csv(chunk_nonambient_calls, outs.nonambient_calls)

    chunk_mitochondrial_summary = [
        x.co_mitochondrial_summary for x in chunk_outs if x.co_mitochondrial_summary is not None
    ]
    if chunk_mitochondrial_summary:
        combine_csv(chunk_mitochondrial_summary, outs.mitochondrial_summary)
    else:
        outs.mitochondrial_summary = None


def filter_barcodes(args, outs):  # pylint: disable=too-many-branches
    """Identify cell-associated partitions."""
    chemistry_def = cr_chemistry.get_primary_chemistry_def(args.chemistry_defs)
    chemistry_description = chemistry_def[cr_chemistry.CHEMISTRY_DESCRIPTION_FIELD]

    # Apply the cell calling config in args.config
    config = helpers.CellCalling(args.config)
    force_cells = helpers.get_force_cells(
        config.force_cells,
        args.probe_bc_def.get(PROBE_BC_SAMPLE_ID, None) if args.probe_bc_def else None,
    )
    recovered_cells = helpers.get_recovered_cells(
        config.recovered_cells,
        args.probe_bc_def.get(PROBE_BC_SAMPLE_ID, None) if args.probe_bc_def else None,
    )

    try:
        correction_data = pd.read_csv(
            ensure_str(args.barcode_correction_csv), converters={"barcode": ensure_binary}
        )
        correction_data["barcode"] = correction_data["barcode"].astype("S")
        ab_utils.augment_correction_table_with_corrected_reads_fraction(correction_data)
        ab_utils.augment_correction_table_with_read_fraction(correction_data)
    except pd.errors.EmptyDataError:
        correction_data = None  # will circumvent aggregate detection below

    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrices_h5)
    has_cmo_data = raw_matrix.feature_ref.has_feature_type(rna_library.MULTIPLEXING_LIBRARY_TYPE)

    # If there is a probe barcode def, restrict to data from this specific sample
    if args.probe_bc_def and args.probe_bc_def[PROBE_BC_SEQS]:
        probe_bc_seq = {seq.encode() for seq in args.probe_bc_def[PROBE_BC_SEQS]}
        offset = args.probe_bc_def[PROBE_BC_OFFSET]
        ending_idx = offset + args.probe_bc_def[PROBE_BC_LEN]
        mask = [bc[offset:ending_idx] in probe_bc_seq for bc in raw_matrix.bcs]
        bcs_idx = [i for i, x in enumerate(mask) if x]
        raw_matrix = raw_matrix.select_barcodes(bcs_idx)
        if correction_data is not None:
            correction_data = correction_data[
                correction_data["barcode"].str[offset:ending_idx].isin(probe_bc_seq)
            ]
    # Always do this whether multiplexed or not so the column is always present.
    # For singleplex data, this will ultimately be the same as FRACION_TOTAL_READS.
    if correction_data is not None:
        ab_utils.augment_correction_table_with_read_fraction(
            correction_data, ab_utils.FRACTION_SAMPLE_READS
        )

    is_targeted = raw_matrix.feature_ref.has_target_features()
    is_rtl = (
        determine_targeting_type_from_csv(args.target_set) != TARGETING_METHOD_HC
        if is_targeted
        else False
    )

    if config.cell_barcodes is not None:
        method = FilterMethod.MANUAL
    elif force_cells is not None:
        method = FilterMethod.TOP_N_BARCODES
    elif config.override_mode is not None:
        method = get_filter_method_from_string(config.override_mode)
    elif args.is_antibody_only:
        method = FilterMethod.ORDMAG
    elif is_targeted and not is_rtl:
        method = FilterMethod.TARGETED
    else:
        method = FilterMethod.ORDMAG_NONAMBIENT

    # Get unique gem groups
    unique_gem_groups = sorted(list(set(args.gem_groups)))

    if config.override_library_types is None:
        if args.is_antibody_only:
            feature_types = [rna_library.ANTIBODY_LIBRARY_TYPE]
        else:
            feature_types = [rna_library.GENE_EXPRESSION_LIBRARY_TYPE]
    else:
        feature_types = config.override_library_types

    # For targeted GEX, retrieve target gene indices for cell calling
    if is_targeted:
        target_features = raw_matrix.feature_ref.get_target_feature_indices()
        assert all(
            x in range(raw_matrix.features_dim) for x in target_features
        ), "Invalid index value found in target_features."
    else:
        target_features = None

    lib_types = raw_matrix.get_library_types()
    # In case of SC_MULTI determine library_types from input config.csv
    if args.multi_graph is not None:
        lib_types = MultiGraph.from_path(args.multi_graph).feature_types()

    num_probe_barcodes = len(args.probe_bc_def[PROBE_BC_SEQS]) if args.probe_bc_def else None
    if (
        method not in [FilterMethod.MANUAL]
        and correction_data is not None
        and any(
            x in correction_data.library_type.unique()
            for x in [rna_library.ANTIBODY_LIBRARY_TYPE, rna_library.ANTIGEN_LIBRARY_TYPE]
        )
    ):
        matrix, summary, removed_bcs_df = helpers.remove_antibody_antigen_aggregates(
            correction_data,
            raw_matrix,
            lib_types,
            config.disable_ab_aggregate_detection,
            num_probe_barcodes=num_probe_barcodes,
        )
        # report all identified aggregate barcodes, together with their umis,
        # umi corrected reads, fraction of corrected reads, and fraction of total reads
        removed_bcs_df = removed_bcs_df.round(
            {
                FRACTION_CORRECTED_READS: 3,
                FRACTION_TOTAL_READS: 3,
                ab_utils.FRACTION_SAMPLE_READS: 3,
            }
        )
        if len(removed_bcs_df) != 0:
            removed_bcs_df.to_csv(outs.aggregate_barcodes)
    else:
        matrix = raw_matrix
        summary = {}

    if config.cell_barcodes is not None:
        assert method == FilterMethod.MANUAL
    elif force_cells is not None:
        assert method == FilterMethod.TOP_N_BARCODES

    total_diversity_key = (
        args.probe_bc_def[PROBE_BC_SAMPLE_ID] + "_total_diversity"
        if args.probe_bc_def
        else "total_diversity"
    )
    summary[total_diversity_key] = matrix.bcs_dim

    genomes = get_reference_genomes(args.reference_path, args.target_set)

    # Get per-gem group cell load
    if recovered_cells is not None:
        gg_recovered_cells = int(float(recovered_cells) / float(len(unique_gem_groups)))
    elif method in (FilterMethod.ORDMAG, FilterMethod.ORDMAG_NONAMBIENT):
        gg_recovered_cells = None
    else:
        gg_recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    ### Get the initial cell calls
    probe_barcode_sample_id = args.probe_bc_def[PROBE_BC_SAMPLE_ID] if args.probe_bc_def else None
    try:
        filtered_metrics_groups, filtered_bcs_groups = helpers.call_initial_cells(
            matrix,
            genomes,
            probe_barcode_sample_id,
            unique_gem_groups,
            method,
            gg_recovered_cells,
            config.cell_barcodes,
            force_cells,
            feature_types,
            chemistry_description,
            target_features,
            has_cmo_data,
            num_probe_barcodes=num_probe_barcodes,
        )
    except InvalidCellBarcode as ex:
        # This is thrown deeper in the code, but caught here with a martian.exit call
        msg = (
            "Cell Barcodes did not match list of valid barcodes with observed reads.  "
            "If attempting manual cell calling, please be sure all input "
            f"barcodes are expected for that chemistry.\n\nError: {ex}"
        )
        martian.exit(msg)

    ### Do additional cell calling
    if method == FilterMethod.ORDMAG_NONAMBIENT:
        (
            filtered_bcs_groups,
            nonambient_summary,
            emptydrops_minimum_umis,
        ) = helpers.call_additional_cells(
            matrix,
            unique_gem_groups,
            genomes,
            filtered_bcs_groups,
            feature_types,
            chemistry_description,
            probe_barcode_sample_id,
            num_probe_barcodes=num_probe_barcodes,
            emptydrops_minimum_umis=helpers.get_emptydrops_minimum_umis(
                config.emptydrops_minimum_umis, probe_barcode_sample_id
            ),
        )
        if nonambient_summary.empty:
            outs.nonambient_calls = None
        else:
            nonambient_summary.to_csv(outs.nonambient_calls, index=False)
        for (genome, gg), value in emptydrops_minimum_umis.items():
            summary[
                f"{genome}_gem_group_{gg}_{get_filter_method_name(FilterMethod.ORDMAG_NONAMBIENT)}_threshold"
            ] = value
    else:
        outs.nonambient_calls = None

    ### Apply a minimum UMI threshold on cell calls
    if method not in [FilterMethod.MANUAL, FilterMethod.TOP_N_BARCODES]:
        filtered_bcs_groups = helpers.apply_global_minimum_umis_threshold(
            matrix,
            unique_gem_groups,
            genomes,
            filtered_bcs_groups,
            feature_types,
            global_minimum_umis=helpers.get_global_minimum_umis(
                config.global_minimum_umis, probe_barcode_sample_id
            ),
        )

        ### Apply a max mitochondrial read percentage threshold on cell calls
        max_mito_percent = helpers.get_max_mito_percent(
            config.max_mito_percent, probe_barcode_sample_id
        )
        filtered_bcs_groups, mitochondrial_summary = helpers.apply_mitochondrial_threshold(
            matrix,
            genomes,
            unique_gem_groups,
            filtered_bcs_groups,
            max_mito_percent,
        )
        mitochondrial_summary.to_csv(outs.co_mitochondrial_summary)
    else:
        outs.co_mitochondrial_summary = None

    def remap_keys_metrics(groups):
        return [{"key": k, "value": v.__dict__} for k, v in groups.items()]

    with open(outs.filtered_metrics_groups, "w") as f:
        tk_safe_json.dump_numpy(
            remap_keys_metrics(filtered_metrics_groups),
            f,
            indent=4,
            sort_keys=True,
        )

    save_filtered_bcs_groups(filtered_bcs_groups, outs.filtered_bcs_groups)

    with open(outs.summary, "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)


def parse_chunked_metrics(chunk_outs):
    filtered_metrics_groups_list = []
    for chunk_out in chunk_outs:
        if chunk_out.filtered_metrics_groups is not None:
            with open(chunk_out.filtered_metrics_groups) as infile:
                filtered_metrics_groups_list.append(json.load(infile))

    filtered_metrics_groups = defaultdict(set)
    for filtered_metrics in filtered_metrics_groups_list:
        for groups in filtered_metrics:
            key_tuple = MetricGroups(
                groups["key"][0], groups["key"][1], groups["key"][2], groups["key"][3]
            )
            result = BarcodeFilterResults(0)
            result.filtered_bcs = groups["value"]["filtered_bcs"]
            result.filtered_bcs_var = groups["value"]["filtered_bcs_var"]
            result.filtered_bcs_cv = groups["value"]["filtered_bcs_cv"]
            result.filtered_bcs_lb = groups["value"]["filtered_bcs_lb"]
            result.filtered_bcs_ub = groups["value"]["filtered_bcs_ub"]
            result.filtered_bcs_cutoff = groups["value"]["filtered_bcs_cutoff"]
            filtered_metrics_groups[key_tuple] = result
    return filtered_metrics_groups


def parse_chunked_filtered_bcs(chunk_outs):
    filtered_bcs_groups_list = [
        chunk_out.filtered_bcs_groups
        for chunk_out in chunk_outs
        if chunk_out.filtered_bcs_groups is not None
    ]

    return load_filtered_bcs_groups(filtered_bcs_groups_list)


def parse_chunked_summary(chunk_outs):
    summary = defaultdict(int)
    for chunk_out in chunk_outs:
        if chunk_out.summary is not None:
            with open(chunk_out.summary) as infile:
                per_chunk_summary = json.load(infile)
                for k, v in per_chunk_summary.items():
                    summary[k] += v
    return dict(summary)


def parse_filtered_bcs_method(filtered_metrics_groups, summary):
    for key_tuple in filtered_metrics_groups:
        if key_tuple.sample is not None:
            summary.update({key_tuple.sample + "_filter_barcodes_method": key_tuple.method})
        else:
            summary.update({"filter_barcodes_method": key_tuple.method})
