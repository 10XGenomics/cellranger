#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import h5py
import json
import martian
import numpy as np
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats

import cellranger.constants as cr_constants
import cellranger.library_constants as lib_constants
import cellranger.matrix as cr_matrix
from cellranger.matrix import CountMatrix
from cellranger.molecule_counter import MoleculeCounter
import cellranger.io as cr_io
import cellranger.rna.matrix as rna_matrix
import cellranger.rna.library as rna_library

def split(args):

    with MoleculeCounter.open(args.molecules, 'r') as mc:
        library_info = mc.get_library_info()
        # For memory request calculation
        num_gem_groups = len(set(lib['gem_group'] for lib in library_info))

        # Number of barcodes in the full matrix
        num_barcodes = mc.get_ref_column_lazy('barcodes').shape[0]

    # Worst case number of nonzero elements in final matrix
    num_nonzero = args.raw_nnz
    join_mem_gb = CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes*num_gem_groups,
                                                                    num_nonzero)

    return {
        'chunks': [],
        'join': {
            '__mem_gb': join_mem_gb,
            '__threads': 2
        }
    }

def main(args, outs):
    assert "No main"

def join(args, outs, chunk_defs, chunk_outs):

    version = martian.get_pipelines_version()

    with open(args.summary) as f:
        summary = json.load(f)

    with MoleculeCounter.open(args.molecules, 'r') as mc:
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()
        barcode_seqs = mc.get_barcodes()

    lib_types = sorted(set(lib['library_type'] for lib in library_info))

    # make attrs for user-added columns in aggr csv
    extra_attrs = get_custom_aggr_columns(args.sample_defs)
    # track original library/gem info
    library_map = cr_matrix.make_library_map_aggr(args.gem_group_index)
    extra_attrs.update(library_map)

    # Merge raw matrix
    raw_matrix = cr_matrix.merge_matrices(args.raw_matrices_h5)
    raw_matrix.save_h5_file(outs.raw_matrix_h5, extra_attrs=extra_attrs)

    genomes = raw_matrix.get_genomes()

    # Create barcode summary HDF5 file w/ GEX data for the barcode rank plot
    with h5py.File(outs.barcode_summary_h5, 'w') as f:
        cr_io.create_hdf5_string_dataset(f, cr_constants.H5_BC_SEQUENCE_COL, raw_matrix.bcs)

        gex_bc_counts = raw_matrix.view().select_features_by_type(lib_constants.GENE_EXPRESSION_LIBRARY_TYPE).sum(axis=0).astype('uint64')
        genome_key = genomes[0] if len(genomes) == 1 else lib_constants.MULTI_REFS_PREFIX
        f.create_dataset('_%s_transcriptome_conf_mapped_deduped_barcoded_reads' % genome_key,
                         data=gex_bc_counts)

    rna_matrix.save_mex(raw_matrix,outs.raw_matrix_mex, version)
    del raw_matrix

    # Merge filtered matrix
    filt_mat = cr_matrix.merge_matrices(args.filtered_matrices_h5)
    filt_mat.save_h5_file(outs.filtered_matrix_h5, extra_attrs=extra_attrs)

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
            filtered_bcs = MoleculeCounter.get_filtered_barcodes(barcode_info,
                                                                 library_info,
                                                                 barcode_seqs,
                                                                 genome_idx=genome_idx,
                                                                 library_type=lib_type)
            mat = mat.select_barcodes_by_seq(filtered_bcs)

            median_features = np.median(mat.count_ge(axis=0,
                                                     threshold=cr_constants.MIN_COUNTS_PER_GENE))
            median_counts = np.median(mat.sum(axis=0))
            genome_prefix = genome if genome is not None else lib_constants.MULTI_REFS_PREFIX

            prefixes = (libtype_prefix, genome_prefix)
            if genome is not None:
                flt_reads = summary['%s%s_flt_mapped_reads' % prefixes]
                raw_reads = summary['%s%s_raw_mapped_reads' % prefixes]
                frac_reads_in_cells = tk_stats.robust_divide(flt_reads, raw_reads)

                summary['%s%s_filtered_bcs_conf_mapped_barcoded_reads_cum_frac' % prefixes] =  frac_reads_in_cells

            summary.update({
                '%s%s_filtered_bcs_median_counts' % prefixes: median_counts,
                '%s%s_filtered_bcs_median_unique_genes_detected' % prefixes: median_features,
            })

        # Compute frac reads in cells across all genomes
        prefixes = [(libtype_prefix, g) for g in genomes if g is not None]
        if len(prefixes) == 0:
            prefixes = [(libtype_prefix, lib_constants.MULTI_REFS_PREFIX)]
        flt_reads = sum(summary['%s%s_flt_mapped_reads' % p] for p in prefixes)
        raw_reads = sum(summary['%s%s_raw_mapped_reads' % p] for p in prefixes)

        frac_reads_in_cells = tk_stats.robust_divide(flt_reads, raw_reads)
        summary['%s%s_filtered_bcs_conf_mapped_barcoded_reads_cum_frac' % (
            libtype_prefix, lib_constants.MULTI_REFS_PREFIX)] = frac_reads_in_cells


    # Write MEX format (do it last because it converts the matrices to COO)
    rna_matrix.save_mex(filt_mat, outs.filtered_matrix_mex, version)

    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)

def get_custom_aggr_columns(sample_defs):
    custom_attrs = {}
    n_samples = len(sample_defs)

    for i, sample_def in enumerate(sample_defs):
        for key, val in sample_def.items():
            if key not in cr_constants.AGG_METADATA_FIELDS:
                if key not in custom_attrs:
                    custom_attrs[key] = [None] * n_samples
                custom_attrs[key][i] = val

    return custom_attrs
