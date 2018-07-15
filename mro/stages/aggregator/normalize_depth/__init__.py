#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
from collections import defaultdict, OrderedDict
import copy
import h5py
import itertools
import json
import numpy as np
import scipy.sparse as sp_sparse

import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats

import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.library_constants as lib_constants
import cellranger.matrix as cr_matrix
import cellranger.molecule_counter as cr_mol_counter
from cellranger.molecule_counter import MoleculeCounter
import cellranger.io as cr_io
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils

__MRO__  = """
stage NORMALIZE_DEPTH(
    in  h5     molecules,
    in  map    detect_cells_gg_metrics,
    in  string normalization_mode,
    out h5     out_molecules,
    out json   summary,
    src py     "stages/aggregator/normalize_depth",
) split using (
    in  float[] frac_reads_kept,
    in  int[]   num_cells,
    in  int     chunk_start,
    in  int     chunk_len,
    out json    chunk_summary,
    out h5      raw_matrix_h5,
    out h5      filtered_matrix_h5,
)
"""

def split(args):
    # default to downsampling by mapped reads
    downsample = True

    if args.normalization_mode == cr_constants.NORM_MODE_NONE:
        downsample = False

    # compute downsample rates for each library
    with MoleculeCounter.open(args.molecules, 'r') as mc:
        library_info = mc.get_library_info()
        usable_reads = mc.get_usable_read_pairs_per_library()
        cells = np.array([mc.get_num_filtered_barcodes_for_library(lib_idx) for lib_idx in xrange(len(library_info))])

        print "Libraries: %s" % library_info
        print "Usable reads: %s" % usable_reads
        print "Cells: %s" % cells

        usable_rpc = np.zeros(len(library_info), dtype=float)
        for i in xrange(len(library_info)):
            usable_rpc[i] = tk_stats.robust_divide(usable_reads[i], cells[i]) if cells[i] > 0 else 0.0

    # Determine lowest depth for each library type
    lt_rpcs = defaultdict(list)
    for lib, rpc in itertools.izip(library_info, usable_rpc):
        lt_rpcs[lib['library_type']].append(rpc)
    min_rpc_by_lt = {lt: min(rpcs) for lt, rpcs in lt_rpcs.iteritems()}

    for lib_idx in xrange(len(library_info)):
        lib_type = library_info[lib_idx]['library_type']
        print "%s Usable read pairs per cell: %s" % (lib_type, usable_rpc[lib_idx])
        print "%s Minimum read pairs usable per cell: %d" % (lib_type, min_rpc_by_lt[lib_type])

    if not downsample:
        frac_reads_kept = np.ones(len(library_info), dtype=float)
    else:
        frac_reads_kept = np.zeros(len(library_info), dtype=float)
        for i in xrange(len(library_info)):
            lib_type = library_info[i]['library_type']
            min_rpc = min_rpc_by_lt[lib_type]
            if min_rpc == 0:
                frac_reads_kept[i] = 0
            else:
                frac_reads_kept[i] = tk_stats.robust_divide(min_rpc, usable_rpc[i])

    # Split the molecule info h5 into equi-RAM chunks, preserving (barcode, gem_group) boundaries
    # Assumes the molecule_info is sorted by (gem_group, barcode)
    tgt_chunk_len = cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK

    chunks = []
    with MoleculeCounter.open(args.molecules, 'r') as mc:
        for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
            chunks.append({
                'frac_reads_kept': list(frac_reads_kept),
                'num_cells': list(cells),
                'chunk_start': chunk_start,
                'chunk_len': chunk_len,
                # Request enough for two copies
                '__mem_gb': MoleculeCounter.estimate_mem_gb(chunk_len, scale=2.0),
            })
    return {
        'chunks': chunks,
        'join': {
            '__mem_gb': h5_constants.MIN_MEM_GB
        }
    }

def summarize_read_matrix(matrix, library_info, barcode_info, barcode_seqs):
    """Summarize matrix of read-pair counts"""
    lib_types = sorted(set(lib['library_type'] for lib in library_info))

    view = matrix.view()
    summary = {}

    for lib_type in lib_types:
        if rna_library.has_genomes(lib_type):
            sum_genomes = map(str, barcode_info.genomes)
        else:
            sum_genomes = [lib_constants.MULTI_REFS_PREFIX]

        for genome in sum_genomes:
            m = view.select_features_by_type(lib_type)
            if rna_library.has_genomes(lib_type):
                m = m.select_features_by_genome(genome)
                genome_idx = barcode_info.genomes.index(genome)
            else:
                genome_idx = None

            prefix = '%s%s' % (rna_library.get_library_type_metric_prefix(lib_type),
                               genome)
            summary['%s_raw_mapped_reads' % prefix] = m.sum()

            filtered_bcs = MoleculeCounter.get_filtered_barcodes(barcode_info,
                                                                 library_info,
                                                                 barcode_seqs,
                                                                 genome_idx=genome_idx,
                                                                 library_type=lib_type)
            filtered_m = m.select_barcodes_by_seq(filtered_bcs)
            summary['%s_flt_mapped_reads' % prefix] = filtered_m.sum()

            summary['%s_filtered_bcs' % prefix] = len(filtered_bcs)
    return summary

def main(args, outs):
    np.random.seed(0)

    with MoleculeCounter.open(args.molecules, 'r') as mc:
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()

        metrics_in = mc.get_all_metrics()
        metrics_out = copy.deepcopy(metrics_in)

        # Compute subsampling rate and approximate new total readpair count
        frac_reads_kept = np.array(args.frac_reads_kept, dtype=float)
        total_reads_in = mc.get_raw_read_pairs_per_library()
        total_reads_out = total_reads_in * frac_reads_kept

        for lib_idx, _ in enumerate(library_info):
            metrics_out[cr_mol_counter.LIBRARIES_METRIC][str(lib_idx)][cr_mol_counter.DOWNSAMPLED_READS_METRIC] = total_reads_out[lib_idx]

        # downsample molecule info
        chunk = slice(args.chunk_start, args.chunk_start + args.chunk_len)
        mol_library_idx = mc.get_column_lazy('library_idx')[chunk]
        mol_read_pairs = mc.get_column_lazy('count')[chunk]

        mol_rate = frac_reads_kept[mol_library_idx]
        del mol_library_idx

        new_read_pairs = np.random.binomial(mol_read_pairs, mol_rate)
        del mol_read_pairs
        del mol_rate

        keep_mol = np.flatnonzero(new_read_pairs)
        new_read_pairs = new_read_pairs[keep_mol]

        mol_gem_group = mc.get_column_lazy('gem_group')[chunk][keep_mol]
        mol_barcode_idx = mc.get_column_lazy('barcode_idx')[chunk][keep_mol]
        mol_feature_idx = mc.get_column_lazy('feature_idx')[chunk][keep_mol]

        # Assert that gem groups start at 1 and are contiguous
        gem_groups = sorted(set(lib['gem_group'] for lib in library_info))
        assert(min(gem_groups) == 1 and \
               np.all(np.diff(np.array(gem_groups,dtype=int)) == 1))

        feature_ref = mc.get_feature_ref()

        # FIXME: assumes one barcode whitelist
        whitelist_size = mc.h5['barcodes'].shape[0]
        num_bcs = whitelist_size * len(gem_groups)
        num_features = feature_ref.get_num_features()

        # Convert molecule barcode indices into matrix barcode indices
        # FIXME: assumes one barcode whitelist
        mol_barcode_idx += (mol_gem_group-1)*whitelist_size
        ones = np.ones(len(mol_barcode_idx), dtype=cr_matrix.DEFAULT_DATA_DTYPE)
        umi_matrix = sp_sparse.coo_matrix((ones,
                                           (mol_feature_idx, mol_barcode_idx)),
                                          shape=(num_features, num_bcs))

        # Create a read-count matrix so we can summarize reads per barcode
        read_matrix = sp_sparse.coo_matrix((new_read_pairs,
                                           (mol_feature_idx, mol_barcode_idx)),
                                          shape=(num_features, num_bcs))
        del ones
        del mol_feature_idx
        del mol_barcode_idx
        del new_read_pairs

        # Get all barcodes strings for the raw matrix
        # FIXME: assumes one barcode whitelist
        barcode_seqs = mc.get_barcodes()
        barcodes = [cr_utils.format_barcode_seq(bc, gg) for gg, bc in \
                    itertools.product(gem_groups, barcode_seqs)]

        # Get mapped reads per barcode per library,genome
        read_summary = {}
        read_matrix = cr_matrix.CountMatrix(feature_ref, barcodes, read_matrix)
        read_matrix.m = read_matrix.m.tocsc(copy=True)
        read_summary = summarize_read_matrix(read_matrix, library_info, barcode_info,
                                             barcode_seqs)
        del read_matrix

        # Construct the raw UMI matrix
        raw_umi_matrix = cr_matrix.CountMatrix(feature_ref, barcodes, umi_matrix)
        raw_umi_matrix.save_h5_file(outs.raw_matrix_h5)

        # Construct the filtered UMI matrix
        filtered_bcs = MoleculeCounter.get_filtered_barcodes(barcode_info,
                                                             library_info,
                                                             barcode_seqs)
        filtered_umi_matrix = raw_umi_matrix.select_barcodes_by_seq(filtered_bcs)
        filtered_umi_matrix.save_h5_file(outs.filtered_matrix_h5)

        summary = {
            'read_summary': read_summary,
            'mol_metrics': metrics_out,
        }

        with open(outs.chunk_summary, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)

def join(args, outs, chunk_defs, chunk_outs):
    with MoleculeCounter.open(args.molecules, 'r') as mc:
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()
        barcode_seqs = mc.get_barcodes()

    lib_types = sorted(set(lib['library_type'] for lib in library_info))

    summary = {
        'frac_reads_kept': chunk_defs[0].frac_reads_kept,
        'num_cells_by_library': chunk_defs[0].num_cells,
    }

    # Merge read summary metrics
    read_summary = defaultdict(int)
    for filename in [co.chunk_summary for co in chunk_outs]:
        with open(filename) as f:
            d = json.load(f)
            for k in d['read_summary'].iterkeys():
                read_summary[k] += d['read_summary'][k]
    summary.update(read_summary)

    # Get summary metrics
    with open(chunk_outs[0].chunk_summary) as f:
        mol_metrics = json.load(f)['mol_metrics']
    chem_keys = [k for k in mol_metrics.iterkeys() if k.startswith('chemistry')]
    for k in chem_keys:
        summary[k] = mol_metrics[k]
    print json.dumps(mol_metrics, indent=4, sort_keys=True)

    # Merge raw matrix
    raw_chunks = [co.raw_matrix_h5 for co in chunk_outs]
    raw_matrix = cr_matrix.merge_matrices(raw_chunks)
    raw_matrix.save_h5_file(outs.raw_matrices_h5)

    genomes = raw_matrix.get_genomes()

    # Create barcode summary HDF5 file w/ GEX data for the barcode rank plot
    with h5py.File(outs.barcode_summary_h5, 'w') as f:
        cr_io.create_hdf5_string_dataset(f, cr_constants.H5_BC_SEQUENCE_COL, raw_matrix.bcs)

        gex_bc_counts = raw_matrix.view().select_features_by_type(lib_constants.GENE_EXPRESSION_LIBRARY_TYPE).sum(axis=0).astype('uint64')
        genome_key = genomes[0] if len(genomes) == 1 else lib_constants.MULTI_REFS_PREFIX
        f.create_dataset('_%s_transcriptome_conf_mapped_deduped_barcoded_reads' % genome_key,
                         data=gex_bc_counts)


    # Merge filtered matrix
    filt_chunks = [co.filtered_matrix_h5 for co in chunk_outs]
    filt_mat = cr_matrix.merge_matrices(filt_chunks)
    filt_mat.save_h5_file(outs.filtered_matrices_h5)


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

    # Report normalization metrics
    all_batches = OrderedDict()

    # These are all per-library-type
    min_frac_reads_kept = np.ones(len(lib_types), dtype='float')
    total_raw_read_pairs = np.zeros(len(lib_types), dtype='uint64')
    total_ds_raw_read_pairs = np.zeros(len(lib_types), dtype='uint64')
    total_cells = np.zeros(len(lib_types), dtype='uint64')

    for lib_type_idx, lib_type in enumerate(lib_types):
        lib_inds = [i for i, lib in enumerate(library_info) if lib['library_type'] == lib_type]
        for lib_idx in lib_inds:
            aggr_id = library_info[lib_idx]['aggr_id']
            old_gg = library_info[lib_idx]['old_gem_group']
            batch = aggr_id + ('-%d' % old_gg if old_gg > 1 else '')
            all_batches[batch] = None

            n_cells = summary['num_cells_by_library'][lib_idx]
            total_cells[lib_type_idx] += n_cells

            lib_metrics = mol_metrics[cr_mol_counter.LIBRARIES_METRIC][str(lib_idx)]
            raw_read_pairs = lib_metrics[cr_mol_counter.TOTAL_READS_METRIC]
            mapped_read_pairs = lib_metrics[cr_mol_counter.USABLE_READS_METRIC]
            ds_read_pairs = lib_metrics[cr_mol_counter.DOWNSAMPLED_READS_METRIC]

            total_raw_read_pairs[lib_type_idx] += raw_read_pairs
            total_ds_raw_read_pairs[lib_type_idx] += ds_read_pairs

            frac_reads_kept = summary['frac_reads_kept'][lib_idx]
            min_frac_reads_kept[lib_type_idx] = min(min_frac_reads_kept[lib_type_idx],
                                                    frac_reads_kept)

            pre_norm_raw_rppc = tk_stats.robust_divide(raw_read_pairs, n_cells)
            pre_norm_mapped_rppc = tk_stats.robust_divide(mapped_read_pairs, n_cells)

            # Prefix with library type and batch
            p = (rna_library.get_library_type_metric_prefix(lib_type),
                 batch)
            summary.update({
                '%s%s_frac_reads_kept' % p: frac_reads_kept,
                '%s%s_pre_normalization_raw_reads_per_filtered_bc' % p: pre_norm_raw_rppc,
                '%s%s_pre_normalization_cmb_reads_per_filtered_bc' % p: pre_norm_mapped_rppc,
            })
    summary['batches'] = all_batches.keys()

    for lib_type_idx, lib_type in enumerate(lib_types):
        mean_rppc = tk_stats.robust_divide(total_raw_read_pairs[lib_type_idx],
                                           total_cells[lib_type_idx])
        ds_mean_rppc = tk_stats.robust_divide(total_ds_raw_read_pairs[lib_type_idx],
                                              total_cells[lib_type_idx])

        p = rna_library.get_library_type_metric_prefix(lib_type)
        summary.update({
            '%spre_normalization_total_reads' % p: total_raw_read_pairs[lib_type_idx],
            '%spost_normalization_total_reads' % p: total_ds_raw_read_pairs[lib_type_idx],
            '%sfiltered_bcs_transcriptome_union' % p: total_cells[lib_type_idx],
            '%spre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc' % p: mean_rppc,
            '%spost_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc' % p: ds_mean_rppc,
            '%slowest_frac_reads_kept' % p: min_frac_reads_kept[lib_type_idx],
        })

    # Write MEX format (do it last because it converts the matrices to COO)
    raw_matrix.save_mex(outs.raw_matrices_mex,
                        rna_feature_ref.save_features_tsv)
    filt_mat.save_mex(outs.filtered_matrices_mex,
                      rna_feature_ref.save_features_tsv)

    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)
