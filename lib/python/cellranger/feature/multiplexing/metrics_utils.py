#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Utils for multiplexing."""

from __future__ import annotations

# from collections import defaultdict
from six import ensure_binary

import cellranger.molecule_counter_extensions as cr_mce
import tenkit.stats as tk_stats

# pylint: disable=no-name-in-module,import-error
from cellranger.fast_utils import reads_per_feature
from cellranger.molecule_counter import (
    COUNT_COL_NAME,
    FEATURE_IDX_COL_NAME,
    LIBRARY_IDX_COL_NAME,
    MoleculeCounter,
)
from cellranger.rna.library import MULTIPLEXING_LIBRARY_TYPE

JIBES_METRIC_SUFFIX = "_jibes"


def add_jibes_suffix(metrics_dict):
    """Add a suffix to metrics json keys for JIBES metrics."""
    old_metric_keys = list(metrics_dict.keys())
    for k in old_metric_keys:
        new_k = k + JIBES_METRIC_SUFFIX
        metrics_dict[new_k] = metrics_dict.pop(k)

    return metrics_dict


def get_tag_read_counts(molecule_info):
    """Returns a dictionary of tag name to total read counts."""
    cnts, feat_id, feat_type, _ = reads_per_feature(molecule_info, False)
    return {
        ensure_binary(y): x
        for x, y, z in zip(cnts, feat_id, feat_type)
        if z == MULTIPLEXING_LIBRARY_TYPE
    }


def _get_tag_ids(molecule_info):
    with MoleculeCounter.open(molecule_info, "r") as mc:
        library_info = mc.get_library_info()
        feature_library_indices = {
            i
            for (i, l) in enumerate(library_info)
            if l["library_type"] == MULTIPLEXING_LIBRARY_TYPE
        }
        idxs = cr_mce.get_indices_for_values(
            mc, [LIBRARY_IDX_COL_NAME], list(feature_library_indices)
        )
        feature_indices = mc.get_column_with_indices(FEATURE_IDX_COL_NAME, idxs)
        counts = mc.get_column_with_indices(COUNT_COL_NAME, idxs)
        feature_ids = [feature_def.id for feature_def in mc.feature_reference.feature_defs]
    return [feature_ids[x] for x in feature_indices], counts


def get_frac_contaminant_tags(tag_assignment, molecule_info):
    """Calculate the fraction of reads from contaminant tags across all the barcodes.

    Args:
        tag_assignment (dict): from a TagAssigner object with the following structure:
            {feature_id : FeatureAssigner.FeatureAssignments Object}
        molecule_info (string): path to the molecule_info_h5 file

    Returns:
        float: a fractional value
    """
    contam_tag_reads, total_tag_reads = 0, 0
    tag_read_counts = get_tag_read_counts(molecule_info)
    # tag_ids, counts = _get_tag_ids(molecule_info)
    # tag_read_counts = defaultdict(int)
    # for i, tag_id in enumerate(tag_ids):
    #    tag_read_counts[tag_id] += counts[i]
    for tag_id, read_count in tag_read_counts.items():
        total_tag_reads += read_count
        if tag_id not in tag_assignment or tag_assignment[tag_id].is_contaminant:
            contam_tag_reads += read_count
    return tk_stats.robust_divide(contam_tag_reads, total_tag_reads)
