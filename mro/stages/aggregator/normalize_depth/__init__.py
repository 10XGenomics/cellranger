#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
from collections import defaultdict, OrderedDict
import copy
import itertools
import json
import numpy as np
import scipy.sparse as sp_sparse
from cellranger.logperf import LogPerf
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats

import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.library_constants as lib_constants
import cellranger.matrix as cr_matrix
from cellranger.matrix import CountMatrix
import cellranger.molecule_counter as cr_mol_counter
from cellranger.molecule_counter import MoleculeCounter
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils

__MRO__  = """
stage NORMALIZE_DEPTH(
    in  map     gem_group_index,
    in  h5      molecules,
    in  string  normalization_mode,
    in  map     gem_group_barcode_ranges,
    out h5[]    raw_matrices_h5,
    out int     raw_nnz,
    out h5[]    filtered_matrices_h5,
    out int     filtered_nnz,
    src py      "stages/aggregator/normalize_depth",
) split (
    in  float[] frac_reads_kept,
    in  int[]   num_cells,
    in  int     chunk_start,
    in  int     chunk_len,
    out json    chunk_summary,
    out h5      raw_matrix_h5,
    out h5      filtered_matrix_h5,
    out int     raw_nnz,
    out int     filtered_nnz,
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

    # For memory request calculation
    num_gem_groups = len(set(lib['gem_group'] for lib in library_info))

    with MoleculeCounter.open(args.molecules, 'r') as mc:
        # Number of barcodes in the full matrix
        num_barcodes = mc.get_ref_column_lazy('barcodes').shape[0]

        for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
            mol_mem_gb = MoleculeCounter.estimate_mem_gb(chunk_len, scale=2.0, cap=False)
            print 'molecule_info mem_gb = %d' % mol_mem_gb

            # Worst case number of nonzero elements in chunk matrix
            num_nonzero = chunk_len
            matrix_mem_gb = CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes*num_gem_groups,
                                                                     num_nonzero)
            print 'matrix mem_gb = %d' % matrix_mem_gb

            mem_gb = max(h5_constants.MIN_MEM_GB, matrix_mem_gb + mol_mem_gb)

            chunks.append({
                'frac_reads_kept': list(frac_reads_kept),
                'num_cells': list(cells),
                'chunk_start': chunk_start,
                'chunk_len': chunk_len,
                # Request enough for two copies
                '__mem_gb': mem_gb,
            })

    # Join is not loading the merged matrix, so it doesn't need much memory.
    # WRITE_MATRICES will use the precise nnz counts to make an appropriate mem request.
    return {
        'chunks': chunks,
        'join': {
            '__mem_gb': 3,
            '__threads': 2
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

    LogPerf.mem()

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

        # Compute matrix dimensions
        # Get the range of possible barcode indices for each gem group.
        gg_barcode_idx_start = np.zeros(1+len(gem_groups), dtype=int)
        gg_barcode_idx_len = np.zeros(1+len(gem_groups), dtype=int)
        for gg_str, idx_range in sorted(args.gem_group_barcode_ranges.iteritems(),
                                        key=lambda kv: int(kv[0])):
            gg = int(gg_str)
            gg_barcode_idx_start[gg] = idx_range[0]
            gg_barcode_idx_len[gg] = idx_range[1] - idx_range[0]

        num_bcs = gg_barcode_idx_len.sum()
        num_features = feature_ref.get_num_features()

        print 'downsampled'
        LogPerf.mem()

        # Convert molecule barcode indices into matrix barcode indices
        # The molecule info barcode_idx is in this space:
        #  [W_0, W_1, ...] where W_i is distinct original whitelist i.
        # The matrix is in, e.g., this space:
        #  [w_0-1, w_1-2, w_0-3, ...] where w_i-j is a copy of whitelist i for gem group j.

        # Return to the original whitelist index
        mol_barcode_idx -= gg_barcode_idx_start.astype(np.uint64)[mol_gem_group]

        # Offset by the cumulative whitelist length up to a barcode's gem group
        gg_barcode_matrix_start = np.cumsum(gg_barcode_idx_len).astype(np.uint64)
        mol_barcode_idx += gg_barcode_matrix_start[mol_gem_group - 1]

        ones = np.ones(len(mol_barcode_idx), dtype=cr_matrix.DEFAULT_DATA_DTYPE)
        umi_matrix = sp_sparse.coo_matrix((ones,
                                           (mol_feature_idx, mol_barcode_idx)),
                                          shape=(num_features, num_bcs))
        print 'created umi matrix'
        LogPerf.mem()

        # Create a read-count matrix so we can summarize reads per barcode
        read_matrix = sp_sparse.coo_matrix((new_read_pairs,
                                           (mol_feature_idx, mol_barcode_idx)),
                                          shape=(num_features, num_bcs))
        del ones
        del mol_feature_idx
        del mol_barcode_idx
        del new_read_pairs

        # Get all barcodes strings for the raw matrix
        barcode_seqs = mc.get_barcodes()

        print len(barcode_seqs), len(gem_groups)
        print 'creating barcode strings'
        LogPerf.mem()

        barcodes = []
        for gg in gem_groups:
            idx_start = gg_barcode_idx_start[gg]
            idx_end = idx_start + gg_barcode_idx_len[gg]
            gg_bcs = np.array([cr_utils.format_barcode_seq(bc, gg) for bc in barcode_seqs[idx_start:idx_end]])
            barcodes.append(gg_bcs)
        barcodes = np.concatenate(barcodes)
        barcodes.flags.writeable = False

        print 'created barcode strings'
        LogPerf.mem()

        # Get mapped reads per barcode per library,genome
        read_summary = {}
        read_matrix = CountMatrix(feature_ref, barcodes, read_matrix)
        read_matrix.m = read_matrix.m.tocsc(copy=True)
        read_summary = summarize_read_matrix(read_matrix, library_info, barcode_info,
                                             barcode_seqs)
        del read_matrix

        print 'created read matrix'
        LogPerf.mem()
        # Construct the raw UMI matrix
        raw_umi_matrix = CountMatrix(feature_ref, barcodes, umi_matrix)
        raw_umi_matrix.save_h5_file(outs.raw_matrix_h5)
        outs.raw_nnz = raw_umi_matrix.m.nnz

        # Construct the filtered UMI matrix
        filtered_bcs = MoleculeCounter.get_filtered_barcodes(barcode_info,
                                                             library_info,
                                                             barcode_seqs)
        filtered_umi_matrix = raw_umi_matrix.select_barcodes_by_seq(filtered_bcs)
        filtered_umi_matrix.save_h5_file(outs.filtered_matrix_h5)
        outs.filtered_nnz = filtered_umi_matrix.m.nnz

        print 'created filtered umi matrix'
        LogPerf.mem()

        summary = {
            'read_summary': read_summary,
            'mol_metrics': metrics_out,
        }

        with open(outs.chunk_summary, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)

    # Don't write MEX from chunks.
    outs.raw_matrices_mex = None
    outs.filtered_matrices_mex = None


def join(args, outs, chunk_defs, chunk_outs):

    # Pass through the matrix chunks and nnz counts
    outs.raw_matrices_h5 = [o.raw_matrix_h5 for o in chunk_outs]
    outs.raw_nnz = sum(o.raw_nnz for o in chunk_outs)
    outs.filtered_matrices_h5 = [o.filtered_matrix_h5 for o in chunk_outs]
    outs.filted_nnz = sum(o.filtered_nnz for o in chunk_outs)

    with MoleculeCounter.open(args.molecules, 'r') as mc:
        library_info = mc.get_library_info()

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

            # Prefix with batch and library type
            if lib_type.lower().startswith(rna_library.CUSTOM_LIBRARY_TYPE_PREFIX.lower()):
                lib_prefix = rna_library.CUSTOM_LIBRARY_TYPE_PREFIX + '_'
            else:
                lib_prefix = rna_library.get_library_type_metric_prefix(lib_type)

            p = (batch, lib_prefix)
            summary.update({
                '%s_%sfrac_reads_kept' % p: frac_reads_kept,
                '%s_%spre_normalization_raw_reads_per_filtered_bc' % p: pre_norm_raw_rppc,
                '%s_%spre_normalization_cmb_reads_per_filtered_bc' % p: pre_norm_mapped_rppc,
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

    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)
