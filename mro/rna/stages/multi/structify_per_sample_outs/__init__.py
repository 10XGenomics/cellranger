#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to put all of the per sample outputs together into a struct."""

import os
from shutil import copyfile
from typing import Any

import martian

import cellranger.cr_io as cr_io
from cellranger.fast_utils import MultiGraph

__MRO__ = """
struct SampleBamFile(
    string  sample,
    bam     bam_file,
    bam.bai bai_index_file,
    bam.csi csi_index_file,
)

struct SampleMetrics(
    string sample,
    json   summary,
    csv    per_barcode_metrics,
)

struct SampleMoleculeInfo(
    string sample,
    h5     h5_file,
    json   summary,
)

struct SampleMatrices(
    string sample,
    h5     filtered_matrix_h5,
    path   filtered_matrix_mex,
    h5     raw_matrix_h5, # per-sample matrix but not filtered by target panel
    h5     raw_probe_bc_matrix,
    path   raw_matrix_mex,
    csv    filtered_barcodes,
    csv    aggregate_barcodes,
    csv    per_probe_metrics,
)

struct SampleSlfeOuts(
    string  sample,
    bam     bam_file,
    bam.bai bai_index_file,
    bam.csi csi_index_file,
    json    metrics_summary,
    csv     per_barcode_metrics,
    h5      molecule_info,
    h5      filtered_matrix_h5,
    path    filtered_matrix_mex,
    h5      raw_matrix_h5,
    h5      raw_probe_bc_matrix,
    path    raw_matrix_mex,
    csv     filtered_barcodes,
    csv     aggregate_barcodes,
    csv     feature_reference,
    csv     target_panel,
    csv     probe_set,
    csv     per_probe_metrics,
)

stage STRUCTIFY_PER_SAMPLE_OUTS(
    in  SampleBamFile[]      sample_bams,
    in  SampleMetrics[]      sample_metrics,
    in  SampleMoleculeInfo[] sample_molecule_infos,
    in  SampleMatrices[]     sample_matrices,
    in  json                 multi_graph,
    in  csv                  feature_reference,
    in  csv                  target_panel,
    in  csv                  probe_set,
    out map<SampleSlfeOuts>  sample_outs,
    out bam                  unassigned_alignments,
    out bam.bai              unassigned_alignments_bai_index,
    out bam.csi              unassigned_alignments_csi_index,
    src py                   "stages/multi/structify_per_sample_outs",
) using (
    volatile = false,
)
"""


def mappify(sample_structs):
    """Convert list of structs to map of structs keyed by sample name."""
    structs = {}
    for item in sample_structs:
        structs[item["sample"]] = item

    return structs


def hard_link(f, relative_path=None):
    """Make a new hard link in this stage directory to the file f, defaulting to the basename of f."""
    if relative_path is None:
        relative_path = os.path.basename(f)

    new_path = martian.make_path(relative_path).decode("utf8")
    cr_io.hardlink_with_fallback(f, new_path)

    return new_path


def hard_link_index_files(sample_bam):
    """Either a CSI or BAI index file will be present, we check for which one exists, and link it.

    Args:
        sample_bam: a SampleBamFile
    """
    bai_index = sample_bam["bai_index_file"]
    csi_index = sample_bam["csi_index_file"]
    csi_file = None
    bai_file = None
    if bai_index is not None and os.path.exists(bai_index):
        bai_file = hard_link(bai_index)
    if csi_index is not None and os.path.exists(csi_index):
        csi_file = hard_link(csi_index)
    return (bai_file, csi_file)


def main(args, outs):
    if args.multi_graph is None:
        return

    multi_graph = MultiGraph.from_path(args.multi_graph)

    # running into problems with Martian when disabling this stage so
    # instead we soft-disable it based on the input.
    if (
        args.sample_metrics is None
        or args.sample_molecule_infos is None
        or args.sample_matrices is None
    ):
        outs.sample_outs = {}
        for sample_id in multi_graph.sample_ids():
            outs.sample_outs[sample_id] = None
        outs.unassigned_alignments = None
        outs.unassigned_alignments_bai_index = None
        outs.unassigned_alignments_csi_index = None
        return

    sample_metrics = mappify(args.sample_metrics)
    sample_molecule_infos = mappify(args.sample_molecule_infos)
    sample_matrices = mappify(args.sample_matrices)

    # we can't just put the old file paths into the struct,
    # we have to hard-link all the files into this directory so that Martian doesn't destroy them
    if args.sample_bams is None:
        sample_bams = None
        outs.unassigned_alignments = None
        outs.unassigned_alignments_index = None
        outs.unassigned_alignments_bai_index = None
        outs.unassigned_alignments_csi_index = None
    else:
        sample_bams = mappify(args.sample_bams)
        unassigned = sample_bams["unassigned"]
        outs.unassigned_alignments = hard_link(unassigned["bam_file"])
        bai_file, csi_file = hard_link_index_files(unassigned)
        outs.unassigned_alignments_bai_index = bai_file
        outs.unassigned_alignments_csi_index = csi_file

    sample_outs = {}
    for sample, molecule_info in sample_molecule_infos.items():
        assert sample not in ["Unassigned", "unassigned"]

        sample_outs[sample] = link_sample_files(
            args,
            sample,
            molecule_info,
            sample_bams,
            sample_matrices[sample],
            sample_metrics[sample],
        )

    outs.sample_outs = sample_outs


def link_sample_files(
    args, sample, molecule_info, sample_bams, sample_data, sample_metrics
) -> dict[str, Any]:
    """Copy or hardlink all sample-level files into stage outs."""
    sample_files = {}

    for key, file in [
        ("feature_reference", args.feature_reference),
        ("target_panel", args.target_panel),
        ("probe_set", args.probe_set),
    ]:
        if file is None:
            sample_files[key] = None
        else:
            sample_files[key] = martian.make_path(sample + "_" + key)
            copyfile(file, sample_files[key])

    sample_files["sample"] = sample
    if sample_bams is None:
        sample_files["bam_file"] = None
        sample_files["bai_index_file"] = None
        sample_files["csi_index_file"] = None
    else:
        sample_bam = sample_bams[sample]
        sample_files["bam_file"] = hard_link(sample_bam["bam_file"])
        bai_file, csi_file = hard_link_index_files(sample_bam)
        sample_files["bai_index_file"] = bai_file
        sample_files["csi_index_file"] = csi_file
    sample_files["molecule_info"] = hard_link(molecule_info["h5_file"])

    def hard_link_sample_data(key: str):
        if sample_data[key] is None:
            sample_files[key] = None
        else:
            sample_files[key] = hard_link(sample_data[key])

    hard_link_sample_data("filtered_matrix_h5")
    hard_link_sample_data("raw_matrix_h5")
    hard_link_sample_data("filtered_barcodes")
    hard_link_sample_data("aggregate_barcodes")
    hard_link_sample_data("per_probe_metrics")
    hard_link_sample_data("raw_probe_bc_matrix")

    # can't hard link a directory, have to hardlink the contents
    # do this for both the matrix MEX and the all_genes version of the matrix mex
    for prefix in ["filtered_", "raw_"]:
        key = prefix + "matrix_mex"
        matrix_mex_dir = os.path.basename(sample_data[key])
        sample_files[key] = martian.make_path(matrix_mex_dir).decode("utf8")
        os.makedirs(sample_files[key])
        hard_link(
            os.path.join(sample_data[key], "barcodes.tsv.gz"),
            relative_path=os.path.join(matrix_mex_dir, "barcodes.tsv.gz"),
        )
        hard_link(
            os.path.join(sample_data[key], "features.tsv.gz"),
            relative_path=os.path.join(matrix_mex_dir, "features.tsv.gz"),
        )
        hard_link(
            os.path.join(sample_data[key], "matrix.mtx.gz"),
            relative_path=os.path.join(matrix_mex_dir, "matrix.mtx.gz"),
        )

    # these two files need to have a <samplename>_ prefix added to not clash
    sample_files["metrics_summary"] = hard_link(
        sample_metrics["summary"],
        relative_path=sample + "_" + os.path.basename(sample_metrics["summary"]),
    )
    if sample_metrics["per_barcode_metrics"] is not None:
        sample_files["per_barcode_metrics"] = hard_link(
            sample_metrics["per_barcode_metrics"],
            relative_path=sample + "_" + os.path.basename(sample_metrics["per_barcode_metrics"]),
        )
    else:
        sample_files["per_barcode_metrics"] = None

    return sample_files
