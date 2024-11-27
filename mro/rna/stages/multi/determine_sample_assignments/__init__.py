#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine sample assignments based on tag assignments for multiplexed data.

The input is the assignments of barcodes to tags,
and also the assignments of tags to samples.
The output is the sample barcodes json, which specifies assignments of barcodes to samples.
This stage also performs removal of multiplets and calculation of some multiplexing metrics
Important: this stage makes the assumption that all of the input is from a single gem well
"""

from __future__ import annotations

import itertools
import json
from copy import deepcopy

import martian

import cellranger.utils as cr_utils
from cellranger.analysis.jibes import (
    BLANK_FACTOR_NAME,
    MULTIPLETS_FACTOR_NAME,
    UNASSIGNED_FACTOR_NAME,
)
from cellranger.cr_io import hardlink_with_fallback
from cellranger.fast_utils import MultiGraph
from cellranger.feature.feature_assignments import (
    BarcodeAssignmentException,
    CellsPerFeature,
    SampleBarcodes,
    validate_force_sample_barcodes,
)
from cellranger.matrix import CountMatrix

__MRO__ = """
stage DETERMINE_SAMPLE_ASSIGNMENTS(
    in  json[]             barcodes_per_tag,
    in  BarcodeAssignments force_sample_barcodes,
    in  csv                filtered_barcodes,
    in  json               multi_graph,
    in  json               non_singlet_barcodes,
    in  int                gem_well,
    in  h5                 raw_feature_bc_matrix,
    out json               cells_per_tag,
    out json               sample_barcodes,
    out json               sample_cell_barcodes,
    out json               non_singlet_barcodes,
    out map<json>          sample_summaries,
    out json               summary,
    src py                 "stages/multi/determine_sample_assignments",
) using (
    mem_gb   = 8,
    volatile = strict,
)
"""


# Calculate per-cell multiplexing metrics
def calculate_sample_assignment_metrics(sample_barcodes, multiplets, filtered_barcodes_csv):
    """Undocumented.

    Args:
        sample_barcodes: SampleBarcodes: maps sample names as bytes to assigned barcodes as list of bytes
        multiplets: Set(bytes) holding the set of multiplet barcodes
        filtered_barcodes_csv: str csv file path holding filtered barcodes

    Returns:
        assignment_metrics: Dict[str] holding library-level assignment metrics
        sample_summaries: Dict[str,str] dict mapping sample names to JSONs holding sample-specific metrics
            most of the metrics are duplicates of the library-level metrics that are needed for the sample web summary
    """
    assignment_metrics = {}
    sample_summaries = {}

    # Count how many samples had singlets
    samples_with_any_singlets = set()
    for sample, barcodes in sample_barcodes.items():
        if len(barcodes) > 0:
            samples_with_any_singlets.add(sample)

    total_singlets = sum(len(bcs) for bcs in sample_barcodes.values())

    filtered_barcodes = cr_utils.load_barcode_csv(filtered_barcodes_csv)
    total_cell_associated_partitions = sum(len(barcodes) for barcodes in filtered_barcodes.values())

    assignment_metrics["cell_associated_partitions_not_assigned_any_samples"] = (
        total_cell_associated_partitions - total_singlets - len(multiplets)
    )
    assignment_metrics["cell_associated_partitions_identified_as_multiplets"] = len(multiplets)
    assignment_metrics["total_cell_associated_partitions"] = total_cell_associated_partitions
    assignment_metrics["samples_with_any_singlets"] = len(samples_with_any_singlets)

    assignment_metrics["total_singlets"] = total_singlets
    genomes = filtered_barcodes.keys()

    for sample_id in sample_barcodes:
        barcodes_for_sample = set(sample_barcodes[sample_id])

        sample_metrics = deepcopy(assignment_metrics)

        for genome in genomes:
            genome_barcodes = set(filtered_barcodes[genome])
            sample_metrics[genome.decode("utf8") + "_singlets_assigned_to_this_sample"] = len(
                genome_barcodes.intersection(barcodes_for_sample)
            )

        sample_metrics["singlets_assigned_to_this_sample"] = len(barcodes_for_sample)
        sample_metrics["singlets_assigned_to_other_samples"] = total_singlets - len(
            barcodes_for_sample
        )

        metrics_filename = "{}_summary.json".format(sample_id.decode("utf8"))
        metrics_path = martian.make_path(metrics_filename)

        with open(metrics_path, "w") as outf:
            json.dump(sample_metrics, outf, indent=4, sort_keys=True)

        sample_summaries[sample_id.decode("utf8")] = metrics_path

    return assignment_metrics, sample_summaries


def make_empty_non_singlet_bc(outs):
    """This file is used by COMPUTE_EXTRA_MULTIPLEXING_METRICS so we create an empty version of it if needed.

    Args:
        outs:
    """
    non_singlet_barcodes = SampleBarcodes()
    non_singlet_barcodes[MULTIPLETS_FACTOR_NAME.encode()] = []
    non_singlet_barcodes[UNASSIGNED_FACTOR_NAME.encode()] = []
    non_singlet_barcodes[BLANK_FACTOR_NAME.encode()] = []
    non_singlet_barcodes.save_to_file(outs.non_singlet_barcodes)


def _get_barcodes_per_tag(config: MultiGraph, args, outs):
    barcodes_per_tag_files = [x for x in args.barcodes_per_tag if x is not None]
    if config.is_multiplexed():
        assert len(barcodes_per_tag_files) <= 1
        # In case of CMO/HASHTAG multiplexing will be restricted to cell barcodes
        # however, in case of RTL & OH multiplexing will include cell and non-cell barcodes
        # CMO/HASHTAG multiplexing produces nothing if no cells are called
        barcodes_per_tag_file = (
            None if len(barcodes_per_tag_files) == 0 else barcodes_per_tag_files[0]
        )
        barcodes_per_tag = CellsPerFeature.load_from_file(barcodes_per_tag_file)
        if config.is_cmo_multiplexed() or config.is_hashtag_multiplexed():
            if barcodes_per_tag_file is not None:
                hardlink_with_fallback(barcodes_per_tag_file, outs.cells_per_tag)
            else:
                outs.cells_per_tag = None
                make_empty_non_singlet_bc(outs)
        return barcodes_per_tag
    else:
        return None


def _get_forced_mutiplets(
    force_sample_barcodes,
    args,
    outs,
    barcodes_per_tag: CellsPerFeature | None,
):
    assignment_metrics = {}
    cells_per_tag_fn = args.force_sample_barcodes["cells_per_tag"]
    non_singlet_bc_fn = args.force_sample_barcodes["non_singlet_barcodes"]
    try:
        validate_force_sample_barcodes(args.filtered_barcodes, cells_per_tag_fn, non_singlet_bc_fn)
    except BarcodeAssignmentException as ex:
        martian.exit(f"{ex}")
    hardlink_with_fallback(force_sample_barcodes, outs.sample_barcodes)
    hardlink_with_fallback(force_sample_barcodes, outs.sample_cell_barcodes)
    hardlink_with_fallback(non_singlet_bc_fn, outs.non_singlet_barcodes)
    hardlink_with_fallback(cells_per_tag_fn, outs.cells_per_tag)

    sample_barcodes = filtered_sample_barcodes = SampleBarcodes.load_from_file(
        force_sample_barcodes
    )
    multiplets = set(
        SampleBarcodes.load_from_file(outs.non_singlet_barcodes).get(
            MULTIPLETS_FACTOR_NAME.encode(), []
        )
    )
    cells_per_tag = CellsPerFeature.load_from_file(args.force_sample_barcodes["cells_per_tag"])
    for tag, barcodes in cells_per_tag.items():
        singlet_barcodes_for_tag = [bc for bc in barcodes if bc not in multiplets]
        metric_name = f"tag_{tag.decode()}_number_of_singlets"
        assignment_metrics[metric_name] = len(singlet_barcodes_for_tag)
    return (
        multiplets,
        assignment_metrics,
        force_sample_barcodes,
        sample_barcodes,
        filtered_sample_barcodes,
        barcodes_per_tag,
    )


def _get_singleplex_multiplets(config: MultiGraph, args, barcodes_per_tag: CellsPerFeature | None):
    assignment_metrics = {}
    filtered_sample_barcodes = SampleBarcodes()
    sample_barcodes = SampleBarcodes()
    sample_ids = list(config.sample_tag_ids().keys())
    assert len(sample_ids) == 1  # no muxing, multi config must have only one sample
    sample_id = sample_ids[0].encode()
    filtered_sample_barcodes[sample_id] = []
    for barcodes in cr_utils.load_barcode_csv(args.filtered_barcodes).values():
        filtered_sample_barcodes[sample_id] += barcodes
    filtered_sample_barcodes[sample_id].sort()
    sample_barcodes[sample_id] = CountMatrix.load_bcs_from_h5(args.raw_feature_bc_matrix)

    return (
        set(),
        assignment_metrics,
        None,
        sample_barcodes,
        filtered_sample_barcodes,
        barcodes_per_tag,
    )


def _get_multiplex_multiplets(config: MultiGraph, args, barcodes_per_tag: CellsPerFeature):
    assignment_metrics = {}
    # Determine which barcodes are multiplets
    multiplets = barcodes_per_tag.get_multiplets()

    filtered_barcodes = cr_utils.get_cell_associated_barcode_set(args.filtered_barcodes)

    # For each sample, for each tag add that tags barcodes to the sample
    # Filter out the multiplets
    filtered_sample_barcodes = SampleBarcodes()
    sample_barcodes = SampleBarcodes()

    for sample_id, sample_tags in config.sample_tag_ids().items():
        sample_id = sample_id.encode()
        filtered_sample_barcodes[sample_id] = []
        sample_barcodes[sample_id] = []

        for tag_name_str in sample_tags:
            tag_name = tag_name_str.encode()

            if tag_name not in barcodes_per_tag:
                # Tags with no cells assigned don't have an entry in cells_per_tag
                metric_name = f"tag_{tag_name_str}_number_of_singlets"
                assignment_metrics[metric_name] = 0
                continue

            # For RTL multiplexing cells per tag include noncell barcodes,
            # restrict "singlets" to bcs intersecting with filtered barcodes
            singlet_barcodes_for_tag = [
                bc
                for bc in barcodes_per_tag[tag_name]
                if bc not in multiplets and bc in filtered_barcodes
            ]

            # In RTL multiplexing the "tag" assignment is independent of
            # cell-calling and all valid barcodes can be unambiguously assigned to a sample.
            # For CMO-multiplexing this list is identical to the above, since singlets are
            # by necessity cell barcodes.
            all_barcodes_for_tag = [bc for bc in barcodes_per_tag[tag_name] if bc not in multiplets]

            metric_name = f"tag_{tag_name_str}_number_of_singlets"
            assignment_metrics[metric_name] = len(singlet_barcodes_for_tag)

            filtered_sample_barcodes[sample_id] += singlet_barcodes_for_tag
            sample_barcodes[sample_id] += all_barcodes_for_tag
    return (
        multiplets,
        assignment_metrics,
        None,
        sample_barcodes,
        filtered_sample_barcodes,
        barcodes_per_tag,
    )


def _get_multiplets(config: MultiGraph, args, outs, barcodes_per_tag: CellsPerFeature | None):
    # force_sample_barcodes overrides whatever the tags say
    # S.M. fixme currently this is valid only for CMO-multiplexing
    force_sample_barcodes = args.force_sample_barcodes["sample_barcodes"]
    if force_sample_barcodes is not None:
        return _get_forced_mutiplets(force_sample_barcodes, args, outs, barcodes_per_tag)

    # no multiplexing, so assign all barcodes to one sample
    elif force_sample_barcodes is None and not config.is_multiplexed():
        return _get_singleplex_multiplets(config, args, barcodes_per_tag)

    # multiplexed experiment
    else:
        assert barcodes_per_tag is not None
        return _get_multiplex_multiplets(config, args, barcodes_per_tag)


def main(args, outs):
    """Chunk phase.

    Args:
        args:
        outs:
    """
    if args.multi_graph is None:
        outs.non_singlet_barcodes = None
        outs.sample_cell_barcodes = None
        outs.sample_barcodes = None
        outs.sample_summaries = None
        outs.summary = None
        outs.cells_per_tag = None
        return

    config = MultiGraph.from_path(args.multi_graph)

    barcodes_per_tag = _get_barcodes_per_tag(config, args, outs)

    # Either CALL_TAGS_JIBES or CALL_TAGS_RTL is enabled.
    (
        multiplets,
        assignment_metrics,
        force_sample_barcodes,
        sample_barcodes,
        filtered_sample_barcodes,
        barcodes_per_tag,
    ) = _get_multiplets(config, args, outs, barcodes_per_tag)

    # Write sample barcodes JSON, unless we're using the presupplied inputs
    if force_sample_barcodes is None:
        sample_barcodes.save_to_file(outs.sample_barcodes)

        # In case of CMOs this file should be identical to the above sample_barcodes
        filtered_sample_barcodes.save_to_file(outs.sample_cell_barcodes)

        if args.non_singlet_barcodes is None:
            non_singlet_barcodes = SampleBarcodes()
            non_singlet_barcodes[MULTIPLETS_FACTOR_NAME.encode()] = list(multiplets)
            non_singlet_barcodes[UNASSIGNED_FACTOR_NAME.encode()] = []
            non_singlet_barcodes[BLANK_FACTOR_NAME.encode()] = []
            non_singlet_barcodes.save_to_file(outs.non_singlet_barcodes)
        else:
            hardlink_with_fallback(args.non_singlet_barcodes, outs.non_singlet_barcodes)

        if config.is_rtl_multiplexed() or config.is_oh_multiplexed():
            cells_per_tag = SampleBarcodes()
            cell_barcodes = set(itertools.chain(*filtered_sample_barcodes.values()))
            for tag in barcodes_per_tag:
                cells_per_tag[tag] = [bc for bc in barcodes_per_tag[tag] if bc in cell_barcodes]
            outs.cells_per_tag = martian.make_path("cells_per_tag.json")
            cells_per_tag.save_to_file(outs.cells_per_tag)
        elif not config.is_multiplexed():
            outs.cells_per_tag = None
    # can potentially generate relevant non-cell singlet metrics here too
    partial_summary, outs.sample_summaries = calculate_sample_assignment_metrics(
        filtered_sample_barcodes, multiplets, args.filtered_barcodes
    )

    assignment_metrics.update(partial_summary)
    with open(outs.summary, "w") as outf:
        json.dump(assignment_metrics, outf, indent=4, sort_keys=True)
