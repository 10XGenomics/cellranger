#!/usr/bin/env python3
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Utils for feature-barcoding technology."""

# Do not add new things to this module.
# Instead, either find or create a module with a name that better describes
# the functionality implemented by the methods or classes you want to add.

import os

import numpy as np
import pandas as pd

import tenkit.safe_json as tk_safe_json


def get_gex_cell_list(filtered_barcodes):
    if not all_files_present([filtered_barcodes]):
        raise ValueError("Filtered barcodes file not present")

    gex_cell_list = set()
    with open(filtered_barcodes, "rb") as f:
        for line in f:
            gex_cell_list.add(line.split(b",")[-1].strip())

    # In case of multiple-genomes, multiplets can lead to duplicates in gex_cell_list, because the genome is stripped away
    gex_cell_list = list(gex_cell_list)
    gex_cell_list.sort()
    return gex_cell_list


def get_feature_counts_as_df(feature_counts_matrix, get_transpose=False):
    """Given a CountMatrix object with feature counts, return a pandas dataframe."""
    feature_defs = feature_counts_matrix.feature_ref.feature_defs
    if get_transpose:
        feature_counts_df = pd.DataFrame(
            feature_counts_matrix.m.T.toarray(),
            columns=[f.id for f in feature_defs],
            index=feature_counts_matrix.bcs,
        )
    else:
        feature_counts_df = pd.DataFrame(
            feature_counts_matrix.m.toarray(),
            index=[f.id for f in feature_defs],
            columns=feature_counts_matrix.bcs,
        )
    return feature_counts_df


def check_if_none_or_empty(matrix):
    return matrix is None or matrix.get_shape()[0] == 0 or matrix.get_shape()[1] == 0


def write_json_from_dict(input_dict, out_file_name):
    with open(out_file_name, "w") as f:
        tk_safe_json.dump_numpy(
            input_dict,
            f,
            indent=4,
            sort_keys=True,
            separators=(",", ": "),
        )


def get_depth_string(num_reads_per_cell):
    return str(np.round(float(num_reads_per_cell) / 1000, 1)) + "k"


def all_files_present(list_file_paths):
    if list_file_paths is None:
        return False

    if any(fpath is None for fpath in list_file_paths):
        return False

    if not all(os.path.isfile(fpath) for fpath in list_file_paths):
        return False

    return True
