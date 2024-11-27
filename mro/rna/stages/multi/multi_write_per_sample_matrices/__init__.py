#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Slice the filtered matrix into per sample matrices."""

from __future__ import annotations

import csv
import gc
import json

import martian

import cellranger.cr_io as cr_io
import cellranger.matrix as cr_matrix
import cellranger.rna.matrix as rna_matrix
import cellranger.utils as cr_utils
from cellranger.fast_utils import MultiGraph

__MRO__ = """
struct SampleMatrices(
    string sample,
    h5     filtered_matrix_h5,
    path   filtered_matrix_mex,
    h5     raw_matrix_h5,
    h5     raw_probe_bc_matrix,
    path   raw_matrix_mex,
    csv    filtered_barcodes,
    csv    aggregate_barcodes,
    csv    per_probe_metrics,
)

stage MULTI_WRITE_PER_SAMPLE_MATRICES(
    in  h5               matrix_h5,
    in  h5               raw_matrix_h5,
    in  csv              filtered_barcodes,
    in  csv              aggregate_barcodes,
    in  json             sample_barcodes,
    in  json             sample_cell_barcodes,
    in  json             multi_graph,
    in  map<h5>          sample_raw_probe_bc_matrices,
    in  map<csv>         samples_per_probe_metrics,
    out SampleMatrices[] sample_matrices,
    src py               "stages/multi/multi_write_per_sample_matrices",
) split (
    in  string           sample,
    out SampleMatrices   matrices,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


def split(args):
    mem_gib = 4 + cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.raw_matrix_h5, scale=1.5)
    return {
        "chunks": [
            {
                "sample": sample_id,
                "__mem_gb": mem_gib,
            }
            for sample_id in MultiGraph.from_path(args.multi_graph).sample_ids()
        ],
        "join": {
            "__mem_gb": 1,
        },
    }


def main(args, outs):  # pylint: disable=too-many-locals
    def load_sample_barcodes(barcodes_file: str) -> list[str]:
        with open(barcodes_file) as f:
            return json.load(f)[args.sample]

    barcodes = {x.encode() for x in load_sample_barcodes(args.sample_barcodes)}
    cell_barcodes = {x.encode() for x in load_sample_barcodes(args.sample_cell_barcodes)}

    chemistry = cr_matrix.CountMatrix.load_chemistry_from_h5(args.matrix_h5)
    gem_groups = list({cr_utils.split_barcode_seq(bc)[1] for bc in barcodes})
    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample, gem_groups, chemistry)

    # make paths for new per-sample (all genes) matrices with all barcodes,
    # for RTL-type multiplexing this includes non-cell barcodes
    sample_raw_matrix_h5 = martian.make_path(f"{args.sample}_raw_feature_barcode_matrix.h5").decode(
        "utf8"
    )
    sample_raw_matrix_mex = martian.make_path(f"{args.sample}_raw_feature_barcode_matrix").decode(
        "utf8"
    )

    chemistry = cr_matrix.CountMatrix.load_chemistry_from_h5(args.matrix_h5)
    gem_groups = list({cr_utils.split_barcode_seq(bc)[1] for bc in barcodes})
    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample, gem_groups, chemistry)

    filter_matrix(
        args.raw_matrix_h5,
        sample_raw_matrix_h5,
        sample_raw_matrix_mex,
        barcodes,
        matrix_attrs,
    )
    # force garbage collection to ensure we free the memory used to load the raw matrix
    gc.collect()

    # make paths for new per-sample matrices (these are already filtered for target panel in case of targeted)
    sample_filtered_matrix_h5 = martian.make_path(
        f"{args.sample}_filtered_feature_barcode_matrix.h5"
    ).decode("utf8")
    sample_filtered_matrix_mex = martian.make_path(
        f"{args.sample}_filtered_feature_barcode_matrix"
    ).decode("utf8")

    filter_matrix(
        args.matrix_h5,
        sample_filtered_matrix_h5,
        sample_filtered_matrix_mex,
        cell_barcodes,
        matrix_attrs,
    )

    sample_raw_probe_bc_matrix = hardlink_sample_file(
        args.sample, args.sample_raw_probe_bc_matrices, "_raw_probe_bc_matrix.h5"
    )
    sample_per_probe_metrics = hardlink_sample_file(
        args.sample, args.samples_per_probe_metrics, "_raw_probe_bc_matrix.csv"
    )

    # write the filtered barcodes file for sample
    sample_filtered_barcodes_csv = martian.make_path(f"{args.sample}_filtered_barcodes.csv").decode(
        "utf8"
    )
    filter_csv(
        args.filtered_barcodes,
        sample_filtered_barcodes_csv,
        cell_barcodes,
        barcode_col=1,
        has_header=False,
    )

    multi_graph = MultiGraph.from_path(args.multi_graph)

    if (
        args.aggregate_barcodes is not None
        and not multi_graph.is_cmo_multiplexed()
        and not multi_graph.is_hashtag_multiplexed()
    ):
        sample_aggregate_barcodes_csv = martian.make_path(
            f"{args.sample}_aggregate_barcodes.csv"
        ).decode("utf8")
        filter_csv(
            args.aggregate_barcodes,
            sample_aggregate_barcodes_csv,
            barcodes,
            barcode_col=0,
            has_header=True,
        )
    else:
        sample_aggregate_barcodes_csv = None

    outs.matrices = {
        "sample": args.sample,
        "filtered_matrix_h5": sample_filtered_matrix_h5,
        "filtered_matrix_mex": sample_filtered_matrix_mex,
        "raw_matrix_h5": sample_raw_matrix_h5,
        "raw_matrix_mex": sample_raw_matrix_mex,
        "raw_probe_bc_matrix": sample_raw_probe_bc_matrix,
        "per_probe_metrics": sample_per_probe_metrics,
        "filtered_barcodes": sample_filtered_barcodes_csv,
        "aggregate_barcodes": sample_aggregate_barcodes_csv,
    }


def join(args, outs, chunk_defs, chunk_outs):
    outs.sample_matrices = [c.matrices for c in chunk_outs]


def filter_matrix(h5_in, h5_out, mex_out, barcodes: set[bytes], matrix_attrs):
    """Filter the provided count matrix by the provided barcodes.

    Write out filtered h5 and mex files.
    """
    matrix = cr_matrix.CountMatrix.load_h5_file(h5_in)

    # filter matrix by sample barcodes, maintain original order
    sample_bcs = [b for b in matrix.bcs if b in barcodes]
    matrix = matrix.select_barcodes_by_seq(sample_bcs)

    matrix.save_h5_file(
        h5_out,
        extra_attrs=matrix_attrs,
        sw_version=martian.get_pipelines_version(),
    )

    rna_matrix.save_mex(matrix, mex_out, martian.get_pipelines_version())


def filter_csv(
    in_path: str, out_path: str, barcodes: set[bytes], barcode_col: int, has_header: bool
):
    """Filter the provided CSV by the provided barcode set.

    Write out the new CSV to out_path.
    """
    with open(in_path, newline="") as csvfile, open(out_path, "w", newline="") as outfile:
        reader = csv.reader(csvfile)
        writer = csv.writer(outfile)
        if has_header:
            writer.writerow(next(reader))
        for row in reader:
            if row[barcode_col].encode() in barcodes:
                writer.writerow(row)


def hardlink_sample_file(
    sample: str, per_sample_files: dict[str, str] | None, filename_suffix: str
) -> str | None:
    """Hardlink a sample-wise file into local outs.

    The sample name is prepended to filename_suffix to form the filename.
    """
    if per_sample_files is None:
        return None
    sample_file = per_sample_files.get(sample, None)
    if sample_file is None:
        return None
    path = martian.make_path(sample + filename_suffix).decode("utf8")
    # copy the raw probe bc matrix if it exists. Not doing this screws things
    # up when vdrmode=strict
    cr_io.hardlink_with_fallback(sample_file, path)
    return path
