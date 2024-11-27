#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#


import csv
import os
import shutil

import martian

__MRO__ = """
stage PARSE_PARAM_CSV(
    in  csv   params_csv,
    out csv   params_csv,
    out int   num_analysis_bcs,
    out int   random_seed,
    out int   num_pca_bcs,
    out int   num_pca_genes,
    out int   num_principal_comps,
    out int   cbc_knn,
    out float cbc_alpha,
    out float cbc_sigma,
    out bool  cbc_realign_panorama,
    out int   max_clusters,
    out int   graphclust_neighbors,
    out float neighbor_a,
    out float neighbor_b,
    out int   tsne_perplexity,
    out int   tsne_input_pcs,
    out int   tsne_max_dims,
    out int   tsne_max_iter,
    out int   tsne_stop_lying_iter,
    out int   tsne_mom_switch_iter,
    out float tsne_theta,
    out int   umap_n_neighbors,
    out int   umap_input_pcs,
    out int   umap_max_dims,
    out float umap_min_dist,
    out string umap_metric,
    src py    "stages/analyzer/parse_csv",
)
"""

ANALYSIS_PARAMS = {
    "num_analysis_bcs": int,
    "random_seed": int,
    "num_pca_bcs": int,
    "num_pca_genes": int,
    "num_principal_comps": int,
    "cbc_knn": int,
    "cbc_alpha": float,
    "cbc_sigma": float,
    "cbc_realign_panorama": bool,
    "max_clusters": int,
    "graphclust_neighbors": int,
    "neighbor_a": float,
    "neighbor_b": float,
    "tsne_perplexity": int,
    "tsne_input_pcs": int,
    "tsne_max_dims": int,
    "tsne_max_iter": int,
    "tsne_stop_lying_iter": int,
    "tsne_mom_switch_iter": int,
    "tsne_theta": float,
    "umap_n_neighbors": int,
    "umap_input_pcs": int,
    "umap_max_dims": int,
    "umap_min_dist": float,
    "umap_metric": str,
}


def main(args, outs):
    parsed = parse_parameters(args.params_csv)
    for param in ANALYSIS_PARAMS:
        if param in parsed:
            setattr(outs, param, parsed[param])
        else:
            setattr(outs, param, None)

    if args.params_csv is not None:
        shutil.copy(args.params_csv, outs.params_csv)


def parse_parameters(filename):
    if filename is None:
        return {}

    if not os.path.exists(filename):
        martian.exit(f"Parameters file does not exist: {filename}")

    if not os.access(filename, os.R_OK):
        martian.exit(f"Parameters file is not readable, please check file permissions: {filename}")

    params = {}
    with open(filename) as f:
        # skip comment lines
        reader = csv.reader(row for row in f if not row.startswith("#"))
        for i, row in enumerate(reader, start=1):
            if len(row) != 2:
                martian.exit("Row %d is incorrectly formatted (must have exactly 2 columns)" % i)
            name = row[0].strip().lower()
            value = row[1].strip()
            if name not in ANALYSIS_PARAMS:
                martian.exit(f"Unrecognized parameter: {name}")
            if name in params:
                martian.exit(f"Cannot specify the same parameter twice: {name}")
            required_type = ANALYSIS_PARAMS[name]
            try:
                if required_type == bool:
                    org_value = value
                    value = value.lower()
                    if value == "true":
                        value = 1
                    elif value == "false":
                        value = 0
                    else:
                        msg = f"Parameter {name} must be set to 'true' or 'false', not {org_value}."
                        martian.exit(msg)
                cast_value = required_type(value)
                params[name] = cast_value
            except ValueError:
                martian.exit(
                    f"Parameter {name} could not be cast to the required type: {required_type}"
                )

    return params
