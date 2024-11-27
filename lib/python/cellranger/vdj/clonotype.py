# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.

"""Clonotype related functions."""

from collections import Counter


def extract_clonotype_id_from_name(clonotype_name: str):
    """Extract the clonotype id from a clonotype name string."""
    return int(clonotype_name.split("clonotype")[1])


def count_changes(initial_labels, final_labels):
    """Counts the number of changes needed to swap between two cluster labels.

    The changes are identified as splits in the initial cluster labels and
    joins among these split labels.
    """
    merged, start_edges, end_edges = Counter(), Counter(), Counter()
    if not len(initial_labels) == len(final_labels):
        raise ValueError("Both input label arrays must be the same size.")
    for pair in zip(initial_labels, final_labels):
        (
            init,
            last,
        ) = pair
        if pair not in merged:
            start_edges[init] += 1
            end_edges[last] += 1
        merged[pair] += 1
    splits = sum(start_edges.values()) - len(start_edges)
    joins = sum(end_edges.values()) - len(end_edges)
    return splits, joins
