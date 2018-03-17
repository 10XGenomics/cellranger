#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import h5py as h5
import itertools
import json
import numpy as np
import operator
import os
import pandas as pd
import shutil
import tables
import scipy.io as sp_io
import scipy.sparse as sp_sparse
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.stats as cr_stats
import cellranger.utils as cr_utils
from numpy.compat import asbytes

MATRIX_H5_FILETYPE = 'matrix'

class GeneBCMatrix:
    def __init__(self, genes, bcs, dtype='int32'):
        self.genes = list(genes)
        self.genes_dim = len(self.genes)
        self.gene_ids_map = {gene.id:i for i, gene in enumerate(self.genes)}

        self.bcs = list(bcs)
        self.bcs_dim = len(self.bcs)
        self.bcs_map = {bc:i for i, bc in enumerate(self.bcs)}

        self.dtype = dtype
        self.m = sp_sparse.lil_matrix((self.genes_dim, self.bcs_dim), dtype=dtype)

    def gene_id_to_int(self, gene_id):
        if gene_id not in self.gene_ids_map:
            raise KeyError("Specified gene not found in matrix: %s" % gene_id)
        return self.gene_ids_map[gene_id]

    def gene_ids_to_ints(self, gene_ids):
        return sorted([self.gene_id_to_int(gene_id) for gene_id in gene_ids])

    def gene_id_to_name(self, gene_id):
        i = self.gene_id_to_int(gene_id)
        return self.int_to_gene_name(i)

    def int_to_gene_id(self, i):
        return self.genes[i].id

    def int_to_gene_name(self, i):
        return self.genes[i].name

    def bc_to_int(self, bc):
        if bc not in self.bcs_map:
            raise KeyError("Specified barcode not found in matrix: %s" % bc)
        return self.bcs_map[bc]

    def bcs_to_ints(self, bcs):
        return sorted([self.bc_to_int(bc) for bc in bcs])

    def int_to_bc(self, j):
        return self.bcs[j]

    def ints_to_bcs(self, jj):
        return [self.int_to_bc(j) for j in jj]

    def add(self, gene_id, bc, value=1):
        i, j = self.gene_id_to_int(gene_id), self.bc_to_int(bc)
        self.m[i,j] += value

    def get(self, gene_id, bc):
        i, j = self.gene_id_to_int(gene_id), self.bc_to_int(bc)
        return self.m[i,j]

    def get_nonzero(self):
        i_array, j_array = self.m.nonzero()
        return [(self.genes[i], self.bcs[j], self.m[i, j])
                 for i, j in itertools.izip(i_array, j_array)]

    def merge(self, group):
        data = group.data.read()
        indices = group.indices.read()
        indptr = group.indptr.read()
        shape = group.shape.read()

        self.m += sp_sparse.csc_matrix((data, indices, indptr), shape=shape).tocoo()

    def save_dense_csv(self, filename):
        dense_gbm = pd.DataFrame(self.m.toarray(), index=[gene.id for gene in self.genes], columns=self.bcs)
        dense_gbm.to_csv(filename, index=True, header=True)

    def _save_h5(self, f, group, attr, arr):
        atom = tables.Atom.from_dtype(arr.dtype)
        if arr.size > 0:
            ds = f.create_carray(group, attr, atom, arr.shape)
        else:
            ds = f.create_earray(group, attr, atom, arr.shape)
        ds[:] = arr

    def save_h5(self, f, group):
        self._save_h5(f, group, cr_constants.H5_GENE_IDS_ATTR, np.array([gene.id for gene in self.genes]))
        self._save_h5(f, group, cr_constants.H5_GENE_NAMES_ATTR, np.array([gene.name for gene in self.genes]))
        self._save_h5(f, group, cr_constants.H5_BCS_ATTR, np.array(self.bcs))

        for attr, dtype in cr_constants.H5_MATRIX_ATTRS.iteritems():
            arr = np.array(getattr(self.m, attr), dtype=dtype)
            self._save_h5(f, group, attr, arr)

    def save_mex(self, target):
        """
        Writes the sparse matrix `self.m` to Matrix Market file `target`.
        Parameters
        ----------
        target : str
            Matrix Market filename (extension .mtx).
        """

        # Add the extension if necessary
        if target[-4:] != '.mtx':
            target = target + '.mtx'

        # Supports only integers currently
        assert( self.m.dtype in ['uint32', 'int32', 'uint64', 'int64'] )
        # Ensure that the matrix is in the COO format
        assert( type(self.m) == sp_sparse.coo.coo_matrix )

        rows, cols = self.m.shape
        # Header fields in the file
        rep = 'coordinate'
        field = 'integer'
        symmetry = 'general'
        comment=''

        with open(target, 'wb') as stream:
            # write initial header line
            stream.write(asbytes('%%MatrixMarket matrix {0} {1} {2}\n'.format(rep, field, symmetry)))

            # write comments
            for line in comment.split('\n'):
                stream.write(asbytes('%%%s\n' % (line)))

            # write shape spec
            stream.write(asbytes('%i %i %i\n' % (rows, cols, self.m.nnz)))
            # write row, col, val in 1-based indexing
            for r, c, d in itertools.izip(self.m.row+1, self.m.col+1, self.m.data):
                stream.write(asbytes(("%i %i %i\n" % (r, c, d))))

    @staticmethod
    def preprocess_matrix(matrix, num_bcs=None, use_bcs=None, use_genes=None, exclude_genes=None, force_cells=None):
        if force_cells is not None:
            bc_counts = matrix.get_reads_per_bc()
            bc_indices, _, _ = cr_stats.filter_cellular_barcodes_fixed_cutoff(bc_counts, force_cells)
            matrix = matrix.select_barcodes(bc_indices)
        elif use_bcs is not None:
            bc_seqs = cr_utils.load_csv_rownames(use_bcs)
            bc_indices = matrix.bcs_to_ints(bc_seqs)
            matrix = matrix.select_barcodes(bc_indices)
        elif num_bcs is not None and num_bcs < matrix.bcs_dim:
            bc_indices = np.sort(np.random.choice(np.arange(matrix.bcs_dim), size=num_bcs, replace=False))
            matrix = matrix.select_barcodes(bc_indices)

        include_indices = list(range(matrix.genes_dim))
        if use_genes is not None:
            include_ids = cr_utils.load_csv_rownames(use_genes)
            include_indices = matrix.gene_ids_to_ints(include_ids)

        exclude_indices = []
        if exclude_genes is not None:
            exclude_ids = cr_utils.load_csv_rownames(exclude_genes)
            exclude_indices = matrix.gene_ids_to_ints(exclude_ids)

        gene_indices = np.array(sorted(list(set(include_indices) - set(exclude_indices))), dtype=int)

        matrix = matrix.select_genes(gene_indices)

        matrix, _, _ = matrix.select_nonzero_axes()
        return matrix

    @staticmethod
    def load_dims(group):
        '''Load the shape of the matrix, as well as the number of nonzero entries.'''
        (genes, bcs) = getattr(group, cr_constants.H5_MATRIX_SHAPE_ATTR).read()
        entries = getattr(group, cr_constants.H5_MATRIX_DATA_ATTR).nrows
        return (genes, bcs, entries)

    @staticmethod
    def get_mem_gb_from_matrix_dim(nonzero_entries):
        ''' Estimate memory usage of loading a matrix. '''
        matrix_mem_gb = round(np.ceil(1.0 * nonzero_entries / cr_constants.NUM_MATRIX_ENTRIES_PER_MEM_GB))
        return cr_constants.MATRIX_MEM_GB_MULTIPLIER * matrix_mem_gb

    @staticmethod
    def get_mem_gb_from_matrix_h5(matrix_h5):
        matrix_dims = GeneBCMatrices.load_dims_from_h5(matrix_h5)
        (genes_dim, bcs_dim, nonzero_entries) = matrix_dims.values()[0]
        return GeneBCMatrix.get_mem_gb_from_matrix_dim(nonzero_entries)

    @staticmethod
    def get_mem_gb_from_group(group):
        _, _, nonzero_entries = GeneBCMatrix.load_dims(group)
        return GeneBCMatrix.get_mem_gb_from_matrix_dim(nonzero_entries)

    @staticmethod
    def load_mtx(genome_dir):
        barcodes_tsv = os.path.join(genome_dir, "barcodes.tsv")
        genes_tsv = os.path.join(genome_dir, "genes.tsv")
        matrix_mtx = os.path.join(genome_dir, "matrix.mtx")
        for filepath in [barcodes_tsv, genes_tsv, matrix_mtx]:
            if not os.path.exists(filepath):
                raise IOError("Required file not found: %s" % filepath)
        barcodes = pd.read_csv(barcodes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        genes = pd.read_csv(genes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        genes = [cr_constants.Gene(gene_id, None, None, None, None) for gene_id in genes]
        matrix = sp_io.mmread(matrix_mtx)
        gbm = GeneBCMatrix(genes, barcodes)
        gbm.m = matrix
        return gbm

    @staticmethod
    def load(group):
        gene_ids = list(getattr(group, cr_constants.H5_GENE_IDS_ATTR).read())

        if hasattr(group, cr_constants.H5_GENE_NAMES_ATTR):
            gene_names = list(getattr(group, cr_constants.H5_GENE_NAMES_ATTR).read())
        else:
            gene_names = gene_ids

        assert len(gene_ids) == len(gene_names)
        genes = [cr_constants.Gene(id, name, None, None, None) for id, name in itertools.izip(gene_ids, gene_names)]
        bcs = list(getattr(group, cr_constants.H5_BCS_ATTR).read())
        matrix = GeneBCMatrix(genes, bcs)

        shape = getattr(group, cr_constants.H5_MATRIX_SHAPE_ATTR).read()
        data = getattr(group, cr_constants.H5_MATRIX_DATA_ATTR).read()
        indices = getattr(group, cr_constants.H5_MATRIX_INDICES_ATTR).read()
        indptr = getattr(group, cr_constants.H5_MATRIX_INDPTR_ATTR).read()

        # quick check to make sure indptr increases monotonically (to catch overflow bugs)
        assert np.all(np.diff(indptr)>=0)

        matrix.m = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)

        return matrix

    @staticmethod
    def load_chunk(group, col_start, col_end):
        ''' Load a submatrix specified by the given column (barcode) range from an h5 group
        Args: col_start, col_end - half-open interval of column indices to load'''
        # Check bounds
        shape = getattr(group, cr_constants.H5_MATRIX_SHAPE_ATTR).read()
        assert col_start >= 0 and col_start < shape[1]
        assert col_end >= 0 and col_end <= shape[1]

        # Load genes and barcodes
        genes = GeneBCMatrix.load_genes_from_h5_group(group)
        bcs = GeneBCMatrix.load_bcs_from_h5_group(group)[col_start:col_end]
        matrix = GeneBCMatrix(genes, bcs)

        # Get views into full matrix
        data = getattr(group, cr_constants.H5_MATRIX_DATA_ATTR)
        indices = getattr(group, cr_constants.H5_MATRIX_INDICES_ATTR)
        indptr = getattr(group, cr_constants.H5_MATRIX_INDPTR_ATTR)

        # Determine extents of selected columns
        ind_start = indptr[col_start]
        if col_end < len(indptr)-1:
            # Last index (end-exclusive) is the start of the next column
            ind_end = indptr[col_end]
        else:
            # Last index is the last index in the matrix
            ind_end = len(data)

        chunk_data = data[ind_start:ind_end]
        chunk_indices = indices[ind_start:ind_end]
        chunk_indptr = np.append(indptr[col_start:col_end], ind_end) - ind_start
        chunk_shape = (shape[0], col_end - col_start)

        matrix.m = sp_sparse.csc_matrix((chunk_data, chunk_indices, chunk_indptr), shape=chunk_shape)

        return matrix

    @staticmethod
    def load_h5(filename):
        """ Load just the first GeneBCMatrix from a GeneBCMatrices file. """
        matrices = GeneBCMatrices.load_h5(filename).matrices.values()
        if len(matrices) > 0:
            return matrices[0]
        else:
            return None

    @staticmethod
    def load_bcs_from_h5_group(group):
        """ Load just the barcode sequences from an h5 """
        return list(getattr(group, cr_constants.H5_BCS_ATTR).read())

    @staticmethod
    def load_genes_from_h5_group(group):
        """ Load just the genes from an h5 """
        gene_ids = list(getattr(group, cr_constants.H5_GENE_IDS_ATTR).read())

        if hasattr(group, cr_constants.H5_GENE_NAMES_ATTR):
            gene_names = list(getattr(group, cr_constants.H5_GENE_NAMES_ATTR).read())
        else:
            gene_names = gene_ids

        assert len(gene_ids) == len(gene_names)
        genes = [cr_constants.Gene(id, name, None, None, None) for id, name in itertools.izip(gene_ids, gene_names)]

        return genes

    @staticmethod
    def hstack(matrix_list, strict=True):
        if not matrix_list:
            return None

        # Compile genes for new stacked matrix
        if strict:
            # All input matrices must have equivalent gene lists, this also becomes output list
            genes = list(matrix_list[0].genes)
            for matrix in matrix_list[1:]:
                if (len(genes) != len(matrix.genes) or any(genes[i] != matrix.genes[i]
                                                           for i in xrange(len(genes)))):
                    raise ValueError('Stacking GeneBCMatrix objects with incompatible genes')
        else:
            # Output gene list is sorted union of input gene lists
            genes = sorted(set().update(m.genes for m in matrix_list))

        # Compile barcodes for new stacked matrix
        bcs = list()
        bcs_group_offset = 0
        bcs_map = dict()
        for i, matrix in enumerate(matrix_list):
            last_bc = matrix.bcs[-1]
            if cr_constants.GEM_GROUP_SEP in last_bc:
                # matrix has GEM group tagged barcodes, so strip these and renumber
                n_bcs_groups = int(last_bc.split(cr_constants.GEM_GROUP_SEP)[1])
                for bc in matrix.bcs:
                    bc_seq, old_group = bc.split(cr_constants.GEM_GROUP_SEP)
                    new_bc = cr_constants.GEM_GROUP_SEP.join([bc_seq,
                                                              str(int(old_group)+bcs_group_offset)])
                    bcs.append(new_bc)
                    bcs_map[(i, bc)] = new_bc
                bcs_group_offset += n_bcs_groups
            else:
                # matrix does not have GEM group tagged barcodes, so give all same new number
                for bc in matrix.bcs:
                    new_bc = cr_constants.GEM_GROUP_SEP.join([bc, str(1+bcs_group_offset)])
                    bcs.append(new_bc)
                    bcs_map[(i, bc)] = new_bc
                bcs_group_offset += 1

        # Set up new matrix and add counts
        new_mat = GeneBCMatrix(genes, bcs)
        if strict:
            new_mat.m = sp_sparse.hstack([matrix.m for matrix in matrix_list])
        else:
            for i, matrix in enumerate(matrix_list):
                for gene, bc, value in matrix.get_nonzero():
                    new_mat.add(gene.id, bcs_map[(i, bc)], value)

        return new_mat

    def tolil(self):
        if type(self.m) is not sp_sparse.lil_matrix:
            self.m = self.m.tolil()

    def tocoo(self):
        if type(self.m) is not sp_sparse.coo_matrix:
            self.m = self.m.tocoo()

    def tocsc(self):
        # Convert from lil to csc matrix for efficiency when analyzing data
        if type(self.m) is not sp_sparse.csc_matrix:
            self.m = self.m.tocsc()

    def _sum(self, matrix, axis=0):
        return np.reshape(np.asarray(matrix.sum(axis=axis)), matrix.shape[1-axis])

    def _topN(self, array, topN=cr_constants.TOP_N):
        if array.shape:
            items = [(i, count) for i, count in enumerate(array)]
        else:
            items = [(0, array.item())]
        sorted_items = sorted(items, key=operator.itemgetter(1), reverse=True)
        return sorted_items[:topN]

    def select_nonzero_axes(self):
        new_mat = GeneBCMatrix(list(self.genes), list(self.bcs))
        new_mat.m = self.m

        nonzero_bcs = np.flatnonzero(new_mat.get_reads_per_bc())
        if new_mat.bcs_dim > len(nonzero_bcs):
            new_mat = new_mat.select_barcodes(nonzero_bcs)

        nonzero_genes = np.flatnonzero(new_mat.get_reads_per_gene())
        if new_mat.genes_dim > len(nonzero_genes):
            new_mat = new_mat.select_genes(nonzero_genes)

        return new_mat, nonzero_bcs, nonzero_genes

    """ Select a subset of barcodes and return the resulting GeneBCMatrix """
    def select_barcodes(self, barcode_indices):
        new_mat = GeneBCMatrix(list(self.genes), np.array(self.bcs)[barcode_indices])
        new_mat.m = self.m[:, barcode_indices]
        return new_mat

    def select_barcodes_by_seq(self, barcode_seqs):
        return self.select_barcodes([self.bc_to_int(bc) for bc in barcode_seqs])

    def select_barcodes_by_gem_group(self, gem_group):
        return self.select_barcodes_by_seq(
            [bc for bc in self.bcs if gem_group == cr_utils.split_barcode_seq(bc)[1]])

    """ Select a subset of genes and return the resulting GeneBCMatrix """
    def select_genes(self, gene_indices):
        new_genes = [cr_constants.Gene(gene[0], gene[1], None, None, None) for \
                     gene in np.array(self.genes)[gene_indices]]
        new_mat = GeneBCMatrix(new_genes, list(self.bcs))
        new_mat.m = self.m[gene_indices,:]
        return new_mat

    def get_unique_genes_per_bc(self):
        return self._sum(self.m > 0, axis=0)

    def get_reads_per_bc(self):
        return self._sum(self.m, axis=0)

    def get_reads_per_gene(self):
        return self._sum(self.m, axis=1)

    def get_top_bcs(self, cutoff):
        reads_per_bc = self.get_reads_per_bc()
        index = max(0, min(reads_per_bc.size, cutoff) - 1)
        value = sorted(reads_per_bc, reverse=True)[index]
        return np.nonzero(reads_per_bc >= value)[0]

    def report(self, genome, barcode_summary_h5, recovered_cells, cell_bc_seqs):
        d = {}

        filtered_mat = self.select_barcodes_by_seq(cell_bc_seqs)
        cell_bc_indices = self.bcs_to_ints(cell_bc_seqs)
        n_cell_bcs = len(cell_bc_seqs)

        # Don't compute metrics if no cells detected
        if n_cell_bcs == 0:
            return d

        # Compute matrix density
        d['filtered_gene_bc_matrix_density'] = tk_stats.robust_divide(filtered_mat.m.getnnz(), filtered_mat.m.shape[0]*filtered_mat.m.shape[1])

        reads_per_gene = filtered_mat.get_reads_per_gene()
        top_genes_with_reads = {filtered_mat.int_to_gene_id(i): int(count) for i, count in filtered_mat._topN(reads_per_gene)}
        d['filtered_bcs_top_genes_with_reads'] = top_genes_with_reads

        unique_bcs_per_gene = filtered_mat._sum(filtered_mat.m >= cr_constants.MIN_READS_PER_BARCODE, axis=1)
        top_genes_with_unique_bcs = {filtered_mat.int_to_gene_id(i): int(count) for i, count in filtered_mat._topN(unique_bcs_per_gene)}
        d['filtered_bcs_top_genes_with_unique_bcs'] = top_genes_with_unique_bcs

        # Total genes and counts
        total_genes_detected = np.count_nonzero(reads_per_gene)
        total_counts = int(reads_per_gene.sum())
        d['filtered_bcs_total_unique_genes_detected'] = total_genes_detected
        d['filtered_bcs_total_counts'] = total_counts

        def _summarize_per_barcode(a):
            mean = np.mean(a)
            stddev = np.std(a)
            return {
                'mean': mean,
                'median': np.median(a),
                'cv': tk_stats.robust_divide(float(stddev), float(mean)),
                'iqr': np.percentile(a, 75) - np.percentile(a, 25),
            }

        # Unique genes per bc
        unique_genes_per_bc = filtered_mat._sum(filtered_mat.m >= cr_constants.MIN_READS_PER_GENE, axis=0)
        unique_genes_stats = _summarize_per_barcode(unique_genes_per_bc)
        for stat, value in unique_genes_stats.iteritems():
            d['filtered_bcs_%s_unique_genes_detected' % stat] = value

        # Counts per bc
        counts_per_bc_stats = _summarize_per_barcode(filtered_mat._sum(filtered_mat.m, axis=0))
        for stat, value in counts_per_bc_stats.iteritems():
            d['filtered_bcs_%s_counts' % stat] = value

        # Cumulative fraction of counts going to top bcs
        d['filtered_bcs_cum_frac'] = tk_stats.robust_divide(filtered_mat.m.sum(), self.m.sum())

        # cDNA PCR duplication in top bcs
        dupe_candidate_h5_key = cr_utils.format_barcode_summary_h5_key(genome, cr_constants.TRANSCRIPTOME_REGION, cr_constants.CONF_MAPPED_BC_READ_TYPE)
        if dupe_candidate_h5_key in barcode_summary_h5:
            n_reads = barcode_summary_h5[dupe_candidate_h5_key][list(cell_bc_indices)].sum()
            n_deduped_reads = filtered_mat.m.sum()
        else:
            n_reads = 0
            n_deduped_reads = 0
        d['filtered_bcs_%s_dupe_reads_frac' % cr_constants.CDNA_PCR_DUPE_TYPE] = 1 - tk_stats.robust_divide(n_deduped_reads, n_reads)

        # Reads per top bc for the various read types (computed over top bcs)
        for read_type in cr_constants.MATRIX_REPORT_READ_TYPES:
            # Compute (n_reads)/(n_bcs) over all bcs and over top bcs
            per_bc_metric = 'filtered_bcs_%s_reads_per_filtered_bc' % read_type

            # Cumulative fraction of reads going to top bcs
            frac_metric = 'filtered_bcs_%s_reads_cum_frac' % read_type

            if read_type in cr_constants.MATRIX_USE_MATRIX_FOR_READ_TYPE:
                n_reads = filtered_mat.m.sum()
                n_all_reads = self.m.sum()
            else:
                h5_key = cr_utils.format_barcode_summary_h5_key(genome, cr_constants.TRANSCRIPTOME_REGION, read_type)
                if h5_key in barcode_summary_h5:
                    n_reads = barcode_summary_h5[h5_key][list(cell_bc_indices)].sum()
                    n_all_reads = barcode_summary_h5[h5_key][()].sum()
                else:
                    n_reads = 0
                    n_all_reads = 0
            d[per_bc_metric] = tk_stats.robust_divide(n_reads, n_cell_bcs)
            d[frac_metric] = tk_stats.robust_divide(n_reads, n_all_reads)
        return d

    def correct_for_saturation(self, diversity, genome):
        m_coo = self.m.tocoo()
        for i, j, old_count in itertools.izip(m_coo.row, m_coo.col, m_coo.data):
            new_count = cr_stats.convert_umi_count_to_transcript_count(old_count, diversity)
            self.m[i,j] = new_count

class GeneBCMatrices:
    def __init__(self, genomes=[], genes=[], barcode_whitelist=None, dtype='int32'):
        self.matrices = {}
        for genome, g in zip(genomes, genes):
            self.matrices[genome] = GeneBCMatrix(g, barcode_whitelist, dtype=dtype)

    def get_genomes(self):
        return self.matrices.keys()

    def get_matrix(self, genome):
        return self.matrices.get(genome)

    def add(self, genome, gene_id, bc):
        self.matrices[genome].add(gene_id, bc)

    def tolil(self):
        for matrix in self.matrices.itervalues():
            matrix.tolil()

    def tocoo(self):
        for matrix in self.matrices.itervalues():
            matrix.tocoo()

    def tocsc(self):
        for matrix in self.matrices.itervalues():
            matrix.tocsc()

    def save_mex(self, base_dir):
        """ Save matrices in Matrix Market Exchange format """
        self.tocoo()
        for name, matrix in self.matrices.iteritems():
            mex_dir = os.path.join(base_dir, name)
            cr_utils.makedirs(mex_dir, allow_existing=True)

            out_matrix_prefix = os.path.join(mex_dir, 'matrix')
            out_genes_fn = os.path.join(mex_dir, 'genes.tsv')
            out_barcodes_fn = os.path.join(mex_dir, 'barcodes.tsv')

            matrix.save_mex(out_matrix_prefix)
            with open(out_genes_fn, 'w') as f:
                for gene in matrix.genes:
                    f.write(gene.id + '\t' + gene.name + '\n')
            with open(out_barcodes_fn, 'w') as f:
                for bc in matrix.bcs:
                    f.write(bc + '\n')

    def save_dense_csv(self, basename):
        for genome, gbm in self.matrices.iteritems():
            gbm.save_dense_csv("%s_%s.csv" % (basename, genome))

    def save_h5(self, filename, extra_attrs={}):
        self.tocsc()
        filters = tables.Filters(complevel = cr_constants.H5_COMPRESSION_LEVEL)
        with tables.open_file(filename, 'w', filters = filters) as f:
            f.set_node_attr('/', cr_constants.H5_FILETYPE_KEY, MATRIX_H5_FILETYPE)
            # set optional top-level attributes
            for (k,v) in extra_attrs.iteritems():
                f.set_node_attr('/', k, v)
            for genome, matrix in self.matrices.iteritems():
                group = f.create_group(f.root, genome)
                matrix.save_h5(f, group)

    def save_barcode_summary_h5(self, filename):
        """ Generate a minimal barcode summary h5 without going through the reporter.
        NOTE: only use this if all genomes have the same set of barcodes, i.e. a raw matrix.
        """
        bc_sequences = None
        bc_table_cols = {}
        total_conf_mapped_deduped_reads = None

        for (genome, matrix) in self.matrices.iteritems():
            if bc_sequences is None:
                bc_sequences = np.array(matrix.bcs)
                bc_table_cols[cr_constants.H5_BC_SEQUENCE_COL] = bc_sequences
            conf_mapped_deduped_reads_key = cr_utils.format_barcode_summary_h5_key(genome,
                cr_constants.TRANSCRIPTOME_REGION, cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)
            conf_mapped_deduped_reads = matrix.get_reads_per_bc()

            if len(bc_sequences) != len(conf_mapped_deduped_reads):
                raise ValueError('Cannot write barcode summary since different genomes have different number of barcodes!')
            bc_table_cols[conf_mapped_deduped_reads_key] = conf_mapped_deduped_reads

            # Track total counts (across genomes)
            if total_conf_mapped_deduped_reads is None:
                total_conf_mapped_deduped_reads = conf_mapped_deduped_reads.copy()
            else:
                total_conf_mapped_deduped_reads += conf_mapped_deduped_reads

        # Record the 'multi'-prefixed (aka total) counts
        # for the web summary to display the barcode rank plot.
        key = cr_utils.format_barcode_summary_h5_key(cr_constants.MULTI_REFS_PREFIX,
                                                     cr_constants.TRANSCRIPTOME_REGION,
                                                     cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)
        bc_table_cols[key] = total_conf_mapped_deduped_reads

        cr_utils.write_h5(filename, bc_table_cols)

    def merge(self, filename):
        self.tocoo()
        with tables.open_file(filename, 'r') as f:
            for group in f.list_nodes(f.root):
                genome = group._v_name
                if genome in self.matrices:
                    self.matrices[genome].merge(group)
                else:
                    self.matrices[genome] = GeneBCMatrix.load(group)

    @staticmethod
    def load_mtx(basedir):
        matrices = GeneBCMatrices()
        for genome in os.listdir(basedir):
            matrices.matrices[genome] = GeneBCMatrix.load_mtx(os.path.join(basedir,genome))
        return matrices

    @staticmethod
    def load_h5(filename):
        matrices = GeneBCMatrices()
        with tables.open_file(filename, 'r') as f:
            for group in f.list_nodes(f.root):
                genome = group._v_name
                matrices.matrices[genome] = GeneBCMatrix.load(group)
        return matrices

    @staticmethod
    def load_dims_from_h5(filename):
        dims = {}
        with tables.open_file(filename, 'r') as f:
            for group in f.list_nodes(f.root):
                genome = group._v_name
                dims[genome] = GeneBCMatrix.load_dims(group)
        return dims

    @staticmethod
    def count_cells_from_h5(filename):
        # NOTE - this double-counts doublets.
        total_cells = 0
        for (genes, bcs, entries) in GeneBCMatrices.load_dims_from_h5(filename).values():
            total_cells += bcs
        return total_cells

    @staticmethod
    def load_genomes_from_h5(filename):
        genomes = []
        with tables.open_file(filename, 'r') as f:
            for group in f.list_nodes(f.root):
                genome = group._v_name
                genomes.append(genome)
        return genomes

    @staticmethod
    def load_chemistry_from_h5(filename):
        with tables.open_file(filename, 'r') as f:
            try:
                chemistry = f.get_node_attr('/', cr_constants.H5_CHEMISTRY_DESC_KEY)
            except AttributeError:
                chemistry = "Unknown"
        return chemistry

    def _get_reads(self, func):
        reads = None
        for matrix in self.matrices.itervalues():
            genome_reads = getattr(matrix, func)()
            if reads is not None:
                reads += genome_reads
            else:
                reads = genome_reads
        return reads

    def get_unique_genes_per_bc(self):
        return self._get_reads('get_unique_genes_per_bc')

    def get_reads_per_bc(self):
        return self._get_reads('get_reads_per_bc')

    def get_reads_per_gene(self):
        return self._get_reads('get_reads_per_gene')

    def _get_stacked_reads_per_bc(self, genomes, barcode_seqs):
        """ Return gene-bc-matrices vertically stacked across genomes """
        # Ensure barcode indices are equivalent across matrices
        bcs = self.matrices.values()[0].bcs
        for matrix in self.matrices.values():
            assert matrix.bcs == bcs

        reads_per_bc = None
        for matrix in [self.matrices[k] for k in genomes]:
            genome_reads_per_bc = matrix.select_barcodes_by_seq(barcode_seqs).get_reads_per_bc()
            if reads_per_bc is not None:
                reads_per_bc = np.vstack((reads_per_bc, genome_reads_per_bc))
            else:
                reads_per_bc = genome_reads_per_bc
        return reads_per_bc

    def union_barcodes(self, cell_bc_seqs, genomes=None):
        """ Take the union of all barcodes present in genome-specific matrices;
        this is used to get a total cell count """
        if genomes is None:
            genomes = self.matrices.keys()
        assert len(cell_bc_seqs) == len(genomes)
        barcodes = set()
        for txome_cell_bc_seqs in cell_bc_seqs:
            barcodes.update(txome_cell_bc_seqs)
        return sorted(barcodes)

    def _report_genome_agnostic_metrics(self, summary_json_paths, barcode_summary_h5, recovered_cells,
                                        cell_bc_seqs):
        """ Report metrics that are computed across all barcodes and all genomes """
        d = {}

        # Get total_reads and *_conf_mapped_reads_frac
        merged_jsons = cr_utils.merge_jsons_as_dict(summary_json_paths)
        total_reads = int(merged_jsons['total_reads'])
        conf_mapped_metrics = ['_'.join([ref,
                                         cr_constants.TRANSCRIPTOME_REGION,
                                         cr_constants.CONF_MAPPED_READ_TYPE,
                                         'reads_frac']) for ref in self.matrices.keys()]
        total_conf_mapped_reads = sum(float(merged_jsons.get(metric, 0)) * float(total_reads) for metric in conf_mapped_metrics)

        # Get number of cell bcs across all genomes
        cell_bcs_union = self.union_barcodes(cell_bc_seqs)
        n_cell_bcs_union = len(cell_bcs_union)
        d['filtered_bcs_transcriptome_union'] = n_cell_bcs_union
        d['%s_filtered_bcs' % cr_constants.MULTI_REFS_PREFIX] = n_cell_bcs_union

        # Report reads/cell across all genomes
        d['%s_%s_total_raw_reads_per_filtered_bc' % (cr_constants.MULTI_REFS_PREFIX, cr_constants.TRANSCRIPTOME_REGION)] = tk_stats.robust_divide(total_reads, n_cell_bcs_union)
        d['%s_%s_total_conf_mapped_reads_per_filtered_bc' % (cr_constants.MULTI_REFS_PREFIX, cr_constants.TRANSCRIPTOME_REGION)] = tk_stats.robust_divide(total_conf_mapped_reads, n_cell_bcs_union)

        # Total UMI counts across all matrices and all filtered barcodes
        total_umi_counts = 0
        for mat in self.matrices.values():
            total_umi_counts += mat.select_barcodes_by_seq(cell_bcs_union).m.sum()


        # Deviation from cell load
        if recovered_cells is None:
            d['%s_filtered_bcs_difference_from_recovered_cells' % cr_constants.MULTI_REFS_PREFIX] = 0
            d['%s_filtered_bcs_relative_difference_from_recovered_cells' % cr_constants.MULTI_REFS_PREFIX] = 0
        else:
            d['%s_filtered_bcs_difference_from_recovered_cells' % cr_constants.MULTI_REFS_PREFIX] = int(n_cell_bcs_union) - int(recovered_cells)
            d['%s_filtered_bcs_relative_difference_from_recovered_cells' % cr_constants.MULTI_REFS_PREFIX] = tk_stats.robust_divide(n_cell_bcs_union - recovered_cells, recovered_cells)

        # Duplicate these metrics across genomes for backwards-compat
        for genome in self.matrices.keys():
            d['%s_total_raw_reads_per_filtered_bc' % genome] = tk_stats.robust_divide(total_reads, n_cell_bcs_union)
            d['%s_total_conf_mapped_reads_per_filtered_bc' % genome] = tk_stats.robust_divide(total_conf_mapped_reads, n_cell_bcs_union)

            for read_type in cr_constants.MATRIX_REPORT_READ_TYPES:
                metric = '%s_total_%s_reads_per_filtered_bc' % (genome, read_type)
                if read_type in cr_constants.MATRIX_USE_MATRIX_FOR_READ_TYPE:
                    n_reads = total_umi_counts
                else:
                    h5_keys = ['%s_%s_%s_reads' % (txome, cr_constants.TRANSCRIPTOME_REGION, read_type) for txome in self.matrices.keys()]
                    h5_keys = [x for x in h5_keys if x in barcode_summary_h5]
                    n_reads = sum(np.array(barcode_summary_h5[h5_key]).sum() for h5_key in h5_keys)
                d[metric] = tk_stats.robust_divide(n_reads, n_cell_bcs_union)

        # Report frac reads in cells across all genomes
        total_conf_mapped_reads_in_cells = 0
        total_conf_mapped_barcoded_reads = 0

        for txome, matrix in self.matrices.iteritems():
            h5_key = '%s_%s_%s_reads' % (txome, cr_constants.TRANSCRIPTOME_REGION,
                                      cr_constants.CONF_MAPPED_BC_READ_TYPE)
            cmb_reads = barcode_summary_h5[h5_key]
            cell_bc_indices = matrix.bcs_to_ints(cell_bcs_union)
            total_conf_mapped_reads_in_cells += cmb_reads[list(cell_bc_indices)].sum() if cell_bc_indices else 0
            total_conf_mapped_barcoded_reads += cmb_reads[()].sum()
        d['multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac'] = tk_stats.robust_divide(total_conf_mapped_reads_in_cells, total_conf_mapped_barcoded_reads)


        # Compute fraction of reads usable (conf mapped, barcoded, filtered barcode)
        unique_barcodes = set(cell_bcs_union)
        in_unique_barcodes_vectorized = np.vectorize(lambda x: x in unique_barcodes)
        filtered_bc_h5_row = in_unique_barcodes_vectorized(np.array(barcode_summary_h5['bc_sequence']))

        usable_reads = 0

        for txome in self.matrices.keys():
            h5_key = '%s_%s_%s_reads' % (txome,
                                                 cr_constants.TRANSCRIPTOME_REGION,
                                                 cr_constants.CONF_MAPPED_BC_READ_TYPE)

            if h5_key not in barcode_summary_h5:
                continue

            usable_reads += (filtered_bc_h5_row * np.array(barcode_summary_h5[h5_key])).sum()

        d['%s_transcriptome_usable_reads_frac' % cr_constants.MULTI_REFS_PREFIX] = tk_stats.robust_divide(usable_reads, total_reads)


        # Compute matrix density across all genomes
        total_nonzero_entries, total_entries = 0, 0
        for matrix in self.matrices.values():
            filtered_mat = matrix.select_barcodes_by_seq(cell_bcs_union)
            total_nonzero_entries += filtered_mat.m.getnnz()
            total_entries += filtered_mat.m.shape[0] * filtered_mat.m.shape[1]
        d['%s_filtered_gene_bc_matrix_density' % cr_constants.MULTI_REFS_PREFIX] = tk_stats.robust_divide(total_nonzero_entries, total_entries)

        return d

    def report(self, summary_json_paths, barcode_summary_h5_path, recovered_cells, cell_bc_seqs):
        assert len(cell_bc_seqs) == len(self.matrices)

        barcode_summary_h5 = h5.File(barcode_summary_h5_path, 'r')

        d = {}

        d.update(self._report_genome_agnostic_metrics(
            summary_json_paths, barcode_summary_h5, recovered_cells, cell_bc_seqs))

        # Compute genome-specific metrics
        for i, (genome, matrix) in enumerate(self.matrices.iteritems()):
            for key, value in matrix.report(genome,
                                            barcode_summary_h5,
                                            recovered_cells,
                                            cell_bc_seqs=cell_bc_seqs[i],
                                        ).iteritems():
                key = '_'.join([genome, key])
                d[key] = value
        return d

    def correct_for_saturation(self, summary_json_paths):
        merged_jsons = cr_utils.merge_jsons_as_dict(summary_json_paths)

        for genome, matrix in self.matrices.iteritems():
            effective_umi_diversity = merged_jsons.get('%s_conf_mapped_effective_umi_diversity' % genome, 0)
            matrix.correct_for_saturation(float(effective_umi_diversity), genome)

    @staticmethod
    def build_from_mol_counter(molecule_counter, subsample_rate=1.0,
                               subsample_result=None):
        """ Construct a GeneBCMatrices object from a MoleculeCounter.
            Args: subsample_result (dict) - Return some metrics results into this dict. """

        # Reconstruct all barcode sequences in the original matrices
        barcode_whitelist = cr_utils.load_barcode_whitelist(molecule_counter.get_barcode_whitelist())
        barcode_length = molecule_counter.get_barcode_length() or len(barcode_whitelist[0])

        gem_groups = molecule_counter.get_gem_groups()
        barcode_seqs = cr_utils.format_barcode_seqs(barcode_whitelist, gem_groups)

        # Reconstruct Gene tuples from the molecule info ref columns
        gene_ids = molecule_counter.get_ref_column('gene_ids')
        genome_ids = molecule_counter.get_ref_column('genome_ids')
        gene_names = molecule_counter.get_ref_column('gene_names')
        gene_tuples = [cr_constants.Gene(gid, gname, None, None, None) for (gid, gname) in itertools.izip(gene_ids, gene_names)]
        genes = cr_utils.split_genes_by_genomes(gene_tuples, genome_ids)

        matrices = GeneBCMatrices(genome_ids, genes, barcode_seqs)

        # Track results of subsampling
        reads = 0

        for mol in molecule_counter.get_molecule_iter(barcode_length, subsample_rate=subsample_rate):
            matrices.add(mol.genome, mol.gene_id, mol.barcode)
            reads += mol.reads

        if subsample_result is not None:
            subsample_result['mapped_reads'] = reads

        return matrices

    def filter_barcodes(self, bcs_per_genome):
        """ Return GeneBCMatrices containing only the specified barcodes.
        Args: bcs_per_genome - dict of {genome: barcodes,...}
        Returns: A new GeneBCMatrices object w/ the specified barcodes """

        filtered_matrices = GeneBCMatrices()
        for (genome, matrix) in self.matrices.iteritems():
            filtered_bcs = bcs_per_genome[genome]
            filtered_matrix = matrix.select_barcodes_by_seq(filtered_bcs)
            filtered_matrices.matrices[genome] = filtered_matrix
        return filtered_matrices

    def report_summary_json(self, filename, summary_json_paths, barcode_summary_h5_path,
                            recovered_cells, cell_bc_seqs):
        """ summary_json_paths: paths to summary jsons containing total_reads and *_conf_mapped_reads_frac
            barcode_summary_h5_path: path to barcode summary h5 file
        """
        d = self.report(summary_json_paths,
                        barcode_summary_h5_path=barcode_summary_h5_path,
                        recovered_cells=recovered_cells,
                        cell_bc_seqs=cell_bc_seqs)
        with open(filename, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(d), f, indent=4, sort_keys=True)

    @staticmethod
    def h5_path(base_path):
        return os.path.join(base_path, "hdf5", "matrices.hdf5")

def merge_matrices(h5_filenames):
    matrices = None
    for h5_filename in h5_filenames:
        if matrices is None:
            matrices = GeneBCMatrices.load_h5(h5_filename)
        else:
            matrices.merge(h5_filename)
    if matrices is not None:
        matrices.tocsc()
    return matrices

def concatenate_mex_dirs(mex_dir_list, out_mex_dir):
    if len(mex_dir_list) == 0:
        return

    # copy tree structure of first dir (assuming genomes, bcs and genes are the same across all chunks)
    cr_utils.copytree(mex_dir_list[0], out_mex_dir)

    # concatenate mtx for each genome
    for genome_dir in os.listdir(out_mex_dir):
        mtx_list = [os.path.join(base_dir, genome_dir, "matrix.mtx") for base_dir in mex_dir_list]
        out_mtx = os.path.join(out_mex_dir, genome_dir, "matrix.mtx")
        concatenate_mtx(mtx_list, out_mtx)

def concatenate_mtx(mtx_list, out_mtx):
    if len(mtx_list) == 0:
        return

    with open(out_mtx, 'w') as out_file:
        # write header
        with open(mtx_list[0], 'r') as in_file:
            out_file.write(in_file.readline())
            out_file.write(in_file.readline())
            (genes, bcs, data) = map(int, in_file.readline().rstrip().split())
        for in_mtx in mtx_list[1:]:
            with open(in_mtx, 'r') as in_file:
                in_file.readline()
                in_file.readline()
                (_, _, mo_data) = map(int, in_file.readline().rstrip().split())
                data += mo_data
        out_file.write(' '.join(map(str, [genes, bcs, data])) + '\n')

        # write data
        for in_mtx in mtx_list:
            with open(in_mtx, 'r') as in_file:
                for i in range(3):
                    in_file.readline()
                shutil.copyfileobj(in_file, out_file)

def concatenate_h5(h5_list, out_h5, extra_attrs={}):
    if len(h5_list) == 0:
        return

    filters = tables.Filters(complevel = cr_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(out_h5, mode = 'w', filters = filters) as f:
        f.set_node_attr('/', cr_constants.H5_FILETYPE_KEY, MATRIX_H5_FILETYPE)

        # set optional top-level attributes
        for (k,v) in extra_attrs.iteritems():
            f.set_node_attr('/', k, v)

        # setup using first file as template
        group_dsets = {}
        group_shapes = {}

        dtypes = cr_constants.H5_MATRIX_ATTRS
        ext_attrs = ['data', 'indices']

        with tables.open_file(h5_list[0], mode = 'r') as fin:
            for in_group in fin.list_nodes(fin.root):
                genome = in_group._v_name
                out_group = f.create_group(f.root, genome)

                dsets = {}

                # copy static-size datasets
                barcodes = getattr(in_group, 'barcodes').read()
                dsets['barcodes'] = f.create_carray(out_group, 'barcodes', obj=barcodes)
                genes = getattr(in_group, 'genes').read()
                dsets['genes'] = f.create_carray(out_group, 'genes', obj=genes)
                gene_names = getattr(in_group, 'gene_names').read()
                dsets['gene_names'] = f.create_carray(out_group, 'gene_names', obj=gene_names)
                indptr = getattr(in_group, 'indptr').read().astype(dtypes['indptr'])
                dsets['indptr'] = f.create_carray(out_group, 'indptr', obj=indptr)
                shape = [len(genes), len(barcodes), 0]

                # initialize extendable datasets
                for name in ext_attrs:
                    atom = tables.Atom.from_dtype(np.dtype(dtypes[name]))
                    dsets[name] = f.create_earray(out_group, name, atom, (0,))

                group_dsets[genome] = dsets
                group_shapes[genome] = shape

        first_mat = True
        for h5_file in h5_list:
            with tables.open_file(h5_file, mode = 'r') as fin:
                for in_group in fin.list_nodes(fin.root):
                    genome = in_group._v_name
                    shape = group_shapes[genome]
                    dsets = group_dsets[genome]

                    data = getattr(in_group, 'data').read()
                    indices = getattr(in_group, 'indices').read()
                    indptr = getattr(in_group, 'indptr').read()

                    dsets['data'].append(data)
                    dsets['indices'].append(indices)

                    if not first_mat:
                        # combine column pointers
                        # since the set of nonzero columns in each chunk is disjoint, this is simple
                        old_indptr = dsets['indptr'][:]
                        dsets['indptr'][:] = old_indptr + indptr

                    shape[2] += len(data)

                first_mat = False

        # set shape
        for (genome, shape) in group_shapes.iteritems():
            genes_bcs = np.array(shape[0:2], dtype=dtypes['shape'])
            dsets = group_dsets[genome]
            dsets['shape'] = f.create_carray('/' + genome, 'shape', obj=genes_bcs)

        # sanity check dimensions
        expected_cols = shape[1]
        expected_nnz = shape[2]
        assert dsets['indptr'].nrows == expected_cols + 1
        assert dsets['indptr'][-1] == expected_nnz
        assert dsets['data'].nrows == expected_nnz
        assert dsets['indices'].nrows == expected_nnz

def make_matrix_attrs_count(sample_id, gem_groups, chemistry):
    matrix_attrs = make_library_map_count(sample_id, gem_groups)
    matrix_attrs[cr_constants.H5_CHEMISTRY_DESC_KEY] = chemistry
    return matrix_attrs

def make_matrix_attrs_aggr(gem_group_index, chemistry):
    matrix_attrs = make_library_map_aggr(gem_group_index)
    matrix_attrs[cr_constants.H5_CHEMISTRY_DESC_KEY] = chemistry
    return matrix_attrs

def get_matrix_attrs(matrix_h5):
    attrs = {}
    with tables.open_file(matrix_h5, 'r') as f:
        for key in cr_constants.H5_METADATA_ATTRS:
            try:
                val = f.get_node_attr('/', key)
                attrs[key] = val
            except AttributeError:
                pass
    return attrs

def make_library_map_count(sample_id, gem_groups):
    # store gem group mapping for use by Cell Loupe
    unique_gem_groups = sorted(set(gem_groups))
    library_map = {
        cr_constants.H5_LIBRARY_ID_MAPPING_KEY: np.array([sample_id]*len(unique_gem_groups), dtype=str),
        cr_constants.H5_ORIG_GEM_GROUP_MAPPING_KEY: np.array(unique_gem_groups, dtype=int)
    }
    return library_map

def make_library_map_aggr(gem_group_index):
    # store gem group mapping for use by Cell Loupe
    library_ids = []
    original_gem_groups = []
    # Sort numerically by new gem group
    for ng, (lid, og) in sorted(gem_group_index.iteritems(), key=lambda pair: int(pair[0])):
        library_ids.append(lid)
        original_gem_groups.append(og)
    library_map = {
        cr_constants.H5_LIBRARY_ID_MAPPING_KEY: np.array(library_ids, dtype=str),
        cr_constants.H5_ORIG_GEM_GROUP_MAPPING_KEY: np.array(original_gem_groups, dtype=int)
    }
    return library_map

def get_gem_group_index(matrix_h5):
    with tables.open_file(matrix_h5, mode = 'r') as f:
        try:
            library_ids = f.get_node_attr('/', cr_constants.H5_LIBRARY_ID_MAPPING_KEY)
            original_gem_groups = f.get_node_attr('/', cr_constants.H5_ORIG_GEM_GROUP_MAPPING_KEY)
        except AttributeError:
            return None
    library_map = {}
    for ng, (lid, og) in enumerate(zip(library_ids, original_gem_groups), start=1):
        library_map[ng] = (lid, og)
    return library_map
