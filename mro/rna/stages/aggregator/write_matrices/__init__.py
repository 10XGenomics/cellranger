#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#


from __future__ import annotations

import json

import martian
import numpy as np

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
import cellranger.rna.matrix as rna_matrix
import cellranger.spatial.spatial_aggr_files as sa_files
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
from cellranger.matrix import CountMatrix
from cellranger.molecule_counter import MoleculeCounter

__MRO__ = """
stage WRITE_MATRICES(
    in  map[] sample_defs,
    in  map   gem_group_index,
    in  h5    molecules,
    in  h5[]  raw_matrices_h5,
    in  int   raw_nnz,
    in  h5[]  filtered_matrices_h5,
    in  int   filtered_nnz,
    in  json  summary,
    in  bool  is_pd,
    out h5    raw_matrix_h5,
    out h5    filtered_matrix_h5,
    out path  filtered_matrix_mex,
    out h5    barcode_summary_h5,
    out json  summary,
    src py    "stages/aggregator/write_matrices",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    with MoleculeCounter.open(args.molecules, "r") as mc:
        num_barcodes = mc.get_barcode_list_size()

    # Worst case number of nonzero elements in final matrix
    if args.is_pd:
        num_nonzero = args.raw_nnz
    else:
        num_nonzero = args.filtered_nnz
    mem_gib = 3 + CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, num_nonzero, scale=1.41)
    print(f"{args.raw_nnz=},{args.filtered_nnz=},{num_nonzero=},{num_barcodes=},{mem_gib=}")
    return {
        "chunks": [],
        "join": {"__mem_gb": mem_gib},
    }


def join(args, outs, chunk_defs, chunk_outs):
    version = martian.get_pipelines_version()
    with open(args.summary) as f:
        summary = json.load(f)

    with MoleculeCounter.open(args.molecules, "r") as mc:
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()
        barcode_seqs = mc.get_barcodes()

    lib_types = rna_library.sorted_library_types(library_info)
    is_antibody_only = False
    if (
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE not in lib_types
        and rna_library.ANTIBODY_LIBRARY_TYPE in lib_types
    ):  # No GEX features found
        is_antibody_only = True

    # make attrs for user-added columns in aggr csv
    extra_attrs = get_custom_aggr_columns(args.sample_defs)
    # track original library/gem info
    library_map = cr_matrix.make_library_map_aggr(args.gem_group_index)
    extra_attrs.update(library_map)
    sw_version = martian.get_pipelines_version()

    # Merge raw matrix for PD purposes
    if args.is_pd:
        martian.log_info("Starting to merge raw matrices.")
        cr_matrix.create_merged_matrix_from_col_concat(
            args.raw_matrices_h5, outs.raw_matrix_h5, extra_attrs, sw_version
        )
        if is_antibody_only:
            raw_matrix = cr_matrix.CountMatrix.load_h5_file(outs.raw_matrix_h5)
            raw_matrix = raw_matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
            raw_matrix.save_h5_file(
                outs.raw_matrix_h5,
                extra_attrs=extra_attrs,
                sw_version=sw_version,
            )
            del raw_matrix
        martian.log_info("Finished merging raw matrices.")

    # Merge filtered matrix
    cr_matrix.create_merged_matrix_from_col_concat(
        args.filtered_matrices_h5, outs.filtered_matrix_h5, extra_attrs, sw_version
    )
    filt_mat = cr_matrix.CountMatrix.load_h5_file(outs.filtered_matrix_h5)
    if is_antibody_only:
        filt_mat = filt_mat.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
        filt_mat.save_h5_file(
            outs.filtered_matrix_h5,
            extra_attrs=extra_attrs,
            sw_version=martian.get_pipelines_version(),
        )
    martian.log_info("Finished merging filtered matrices.")

    # Summarize the matrix across library types and genomes
    for lib_type in lib_types:
        libtype_prefix = rna_library.get_library_type_metric_prefix(lib_type)

        if rna_library.has_genomes(lib_type):
            genomes = filt_mat.get_genomes()
        else:
            genomes = [None]

        mat_lib = filt_mat.view().select_features_by_type(lib_type)

        for genome in genomes:
            if genome is None:
                mat = mat_lib
                genome_idx = None
            else:
                mat = mat_lib.select_features_by_genome(genome)
                genome_idx = barcode_info.genomes.index(genome)

            # Select barcodes passing filter for this (lib_type, genome)
            filtered_bcs = MoleculeCounter.get_filtered_barcodes(
                barcode_info,
                library_info,
                barcode_seqs,
                genome_idx=genome_idx,
                library_type=lib_type,
            )
            mat = mat.select_barcodes_by_seq(filtered_bcs)

            median_features = np.median(
                mat.count_ge(axis=0, threshold=cr_constants.MIN_COUNTS_PER_GENE)
            )
            median_counts = np.median(mat.sum(axis=0))
            genome_prefix = genome if genome is not None else rna_library.MULTI_REFS_PREFIX

            prefixes = (libtype_prefix, genome_prefix)
            if genome is not None:
                flt_reads = summary["{}{}_flt_mapped_reads".format(*prefixes)]
                raw_reads = summary["{}{}_raw_mapped_reads".format(*prefixes)]
                frac_reads_in_cells = tk_stats.robust_divide(flt_reads, raw_reads)

                summary[
                    "{}{}_filtered_bcs_conf_mapped_barcoded_reads_cum_frac".format(*prefixes)
                ] = frac_reads_in_cells
            martian.log_info("Finished collecting metrics")

            summary.update(
                {
                    "{}{}_filtered_bcs_median_counts".format(*prefixes): median_counts,
                    "{}{}_filtered_bcs_median_unique_genes_detected".format(
                        *prefixes
                    ): median_features,
                }
            )

        # Compute frac reads in cells across all genomes
        prefixes = [(libtype_prefix, g) for g in genomes if g is not None]
        if len(prefixes) == 0:
            prefixes = [(libtype_prefix, rna_library.MULTI_REFS_PREFIX)]
        flt_reads = sum(summary["{}{}_flt_mapped_reads".format(*p)] for p in prefixes)
        raw_reads = sum(summary["{}{}_raw_mapped_reads".format(*p)] for p in prefixes)

        frac_reads_in_cells = tk_stats.robust_divide(flt_reads, raw_reads)
        summary[
            f"{libtype_prefix}{rna_library.MULTI_REFS_PREFIX}_filtered_bcs_conf_mapped_barcoded_reads_cum_frac"
        ] = frac_reads_in_cells
    if args.is_pd:
        outs.filtered_matrix_mex = None
    else:
        martian.log_info("Writing MEX File")
        # Write MEX format (do it last because it converts the matrices to COO)
        rna_matrix.save_mex(filt_mat, outs.filtered_matrix_mex, version)
        martian.log_info("Finished writing MEX file")

    with open(outs.summary, "w") as f:
        tk_safe_json.dump_numpy(summary, f, indent=4, sort_keys=True)


def get_custom_aggr_columns(sample_defs):
    custom_attrs = {}
    n_samples = len(sample_defs)

    for i, sample_def in enumerate(sample_defs):
        for key, val in sample_def.items():
            if key not in [
                cr_constants.AGG_ID_FIELD,
                cr_constants.AGG_H5_FIELD,
                cr_constants.AGG_BATCH_FIELD,
                cr_constants.AGG_CLOUPE_FIELD,
                sa_files.AGG_SPATIAL_FIELD,
            ]:
                if key not in custom_attrs:
                    custom_attrs[key] = [None] * n_samples
                custom_attrs[key][i] = val

    return custom_attrs
