#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Compute extra metrics for multiplexing from the uber molecule info.

and the sets of multiplet and singlet barcodes

- median CMO UMIs over all singlets
- fraction of CMO reads lost to multiplets
- fraction reads in cells on a per-tag basis
"""

from __future__ import annotations

import json

import numpy as np

import cellranger.feature_ref as feature_ref
import cellranger.rna.library as rna_library
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
from cellranger import matrix as cr_matrix
from cellranger.fast_utils import MultiGraph, tag_read_counts
from cellranger.feature.feature_assignments import SampleBarcodes

__MRO__ = """
stage COMPUTE_EXTRA_MULTIPLEXING_METRICS(
    in  h5     molecule_info,
    in  h5     filtered_feature_counts_matrix,
    in  json   multi_graph,
    in  json   sample_cell_barcodes,
    in  json   non_singlet_barcodes,
    in  string multiplexing_method,
    out json   summary,
    src py     "stages/multi/compute_extra_multiplexing_metrics",
) split (
) using (
    volatile = strict,
)
"""


# calculate median CMO UMIs per singlet
def calculate_median_cmo_umis_per_singlet(fbm_fname, singlet_barcodes_json, library_type):
    """Calculate the median CMO or Hashtag UMIs per singlet.

    Args:
        fbm_fname: File name of filtered feature_barcode_matrix
        singlet_barcodes_json: the SampleBarcodes in JSON format mapping sample names to list of barcodes assigned to it
        library_type: Multiplexing Capture (for CMO) or Antibody Capture (for Hashatag)

    Returns:
        dict containing the MULTIPLEXING_median_cmo_umis_per_singlet metric.
    """
    fbm = cr_matrix.CountMatrix.load_h5_file(fbm_fname)

    # Subset features to multiplexing types
    if library_type == rna_library.MULTIPLEXING_LIBRARY_TYPE:
        fbm = fbm.select_features_by_type(library_type)
    else:
        fbm = fbm.select_features_by_type_and_tag(library_type, feature_ref.HASHTAG_TAG)

    # If not multiplexing features, return None
    if fbm.features_dim == 0:
        return None
    # Subset down to singlet barcodes
    sample_barcodes = SampleBarcodes.load_from_file(singlet_barcodes_json)
    flattened_bcs = [bc for tag in sample_barcodes.values() for bc in tag]
    del sample_barcodes
    fbm = fbm.select_barcodes_by_seq(flattened_bcs)
    del flattened_bcs
    umis_per_singlet = np.ravel(fbm.m.sum(0))
    del fbm
    median_umis_per_singlet = np.median(umis_per_singlet)
    return median_umis_per_singlet


def split(args):
    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_feature_counts_matrix)
    num_cells = matrix_dims[1]
    nnz = matrix_dims[2]

    mem_gb = 2.0 + 2.0e-08 * nnz + 6.3e-06 * num_cells
    vmem_gb = 3.0 + 1.4e-08 * nnz + 2e-05 * num_cells

    if vmem_gb < mem_gb + 3.0:
        vmem_gb = mem_gb + 3.0

    return {
        "chunks": [],
        "join": {"__mem_gb": mem_gb, "__threads": 1, "__vmem_gb": vmem_gb},
    }


def main(args, _outs):
    pass


# read in molecule info for CMO reads only,
# figure out the indices for singlet and multiplet barcodes inside the molecule info,
# calculate these 3 metrics:
# - median CMO UMIs over all singlets
# - fraction of CMO reads lost to multiplets
# - fraction reads in cells on a per-tag basis
def join(args, outs, chunk_defs, chunk_outs):
    # pylint: disable=too-many-locals
    if args.multi_graph is None:
        outs.summary = None
        return

    config = MultiGraph.from_path(args.multi_graph)

    multiplexing_method = rna_library.BarcodeMultiplexingType(args.multiplexing_method)
    assert multiplexing_method.is_cell_multiplexed(), "Unsupported multiplexing method!"
    library_type = multiplexing_method.multiplexing_library_type()

    metrics = {}
    lib_prefix = rna_library.metric_prefix_map[library_type]
    singlet_metric_name = f"{lib_prefix}_median_cmo_umis_per_singlet"
    frac_reads_metric_name = f"{lib_prefix}_frac_reads_from_multiplets"

    def tag_frac_metric_name_maker(tag):
        return f"tag_{tag}_frac_reads_in_cells"

    # Caclulate median CMO UMI per singlet
    median_umi_per_singlet = calculate_median_cmo_umis_per_singlet(
        args.filtered_feature_counts_matrix, args.sample_cell_barcodes, library_type
    )

    if median_umi_per_singlet is None:  # Occurs if no multiplexing features
        metrics[singlet_metric_name] = 0.0
        metrics[frac_reads_metric_name] = float("nan")
        metrics[tag_frac_metric_name_maker("None")] = float("nan")
    else:
        # Calculate metrics when possible
        metrics[singlet_metric_name] = median_umi_per_singlet
        mol_info = args.molecule_info
        # Get read counts from the molecule info
        (
            tag_names,
            _singlet_counts,
            multiplet_counts,
            cell_associated_counts,
            total_counts,
        ) = tag_read_counts(
            mol_info,
            args.sample_cell_barcodes,
            args.non_singlet_barcodes,
            json.dumps(args.multiplexing_method),
        )

        # get the list of tags from multi graph
        tags = set(tag for names in config.sample_tag_ids().values() for tag in names)

        # Calculate tag reads in cells for each tag
        for i, tag in enumerate(tag_names):
            if tag in tags:
                metric_name = tag_frac_metric_name_maker(tag)
                frac_tag_reads_in_cells = tk_stats.robust_divide(
                    cell_associated_counts[i], total_counts[i]
                )
                metrics[metric_name] = frac_tag_reads_in_cells

        # How many reads went to barcodes assigned as multiplets?
        metrics[frac_reads_metric_name] = tk_stats.robust_divide(
            np.sum(multiplet_counts), np.sum(total_counts)
        )

    with open(outs.summary, "w") as outf:
        tk_safe_json.dump_numpy(metrics, outf, indent=4, sort_keys=True)
