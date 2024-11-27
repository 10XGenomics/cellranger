#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Postprocess the cell type diff expression generated.

Rewrite the header of the cell types diff expression csv to have the celltypes 
and generate the figure for websummary.
"""

import json
import os
from dataclasses import asdict, dataclass

import martian

from cellranger.analysis.analysis_types import DifferentialExpressionWithFeatures
from cellranger.websummary.analysis_tab_core import diffexp_table
from cellranger.websummary.react_components import ReactComponentEncoder

__MRO__ = """
stage TIDY_CELLTYPE_DIFFEXP(
    in  path diffexp_csv,
    in  json cell_types_map,
    out csv  cas_diffexp_csv,
    out json diffexp,
    src py   "stages/cas_cell_typing/tidy_celltype_diffexp",
)
"""


@dataclass
class DifferentialExpressionTable:
    """Differential expression table for a single bin size and clustering algorithm."""

    table: dict


def replace_header(header, cluster_mapping):
    """Replaces substrings in the header string based on a provided mapping dictionary.

    Parameters:
        header (str): The header line of a CSV file, where each column name is
                      separated by commas.
        cluster_mapping (dict): A dictionary where each key is a new string to
                                replace the existing 'Cluster X' substring and each
                                value is the integer X part of the 'Cluster X' pattern
                                to be replaced in the header.

    Returns:
        str: A new header string with the specified replacements made.
    """
    parts = header.split(",")
    new_parts = []
    for part in parts:
        part_split = part.split()
        if len(part_split) >= 2:
            cluster_name = " ".join(part_split[:2])
            for key, value in cluster_mapping.items():
                if f"Cluster {value}" == cluster_name:
                    part = part.replace(f"Cluster {value}", key)
        new_parts.append(part)
    return ",".join(new_parts)


def main(args, outs):
    if not args.cell_types_map or not args.diffexp_csv or not args.filtered_matrix:
        martian.clear(outs)
        return

    diffexp_csv_path = os.path.join(
        args.diffexp_csv, "gene_expression_graphclust", "differential_expression.csv"
    )
    with open(args.cell_types_map) as file:
        cell_types_map = json.load(file)

    with open(diffexp_csv_path) as infile:
        with open(martian.make_path(outs.cas_diffexp_csv), "w") as outfile:
            # Replace header with celltypes instead of clusters
            header = infile.readline().strip()
            new_header = replace_header(header, cell_types_map)
            outfile.write(new_header + "\n")
            # Write the rest of the file as is
            for line in infile:
                outfile.write(line)

    diffexp_with_features = DifferentialExpressionWithFeatures.from_diffexp_csv(diffexp_csv_path)
    cell_types_list = sorted(cell_types_map, key=lambda x: cell_types_map[x])

    diffexp_table_dict = DifferentialExpressionTable(
        table=json.loads(
            json.dumps(
                diffexp_table(
                    diffexp_with_features=diffexp_with_features,
                    cluster_names=cell_types_list,
                ),
                cls=ReactComponentEncoder,
            )
        )
    )

    with open(outs.diffexp, "w") as f:
        json.dump(asdict(diffexp_table_dict), f)
