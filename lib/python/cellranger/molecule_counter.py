#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
from collections import namedtuple, defaultdict
import itertools
import multiprocessing
import numpy as np
import tables
import cellranger.utils as cr_utils
import cellranger.constants as cr_constants

MOLECULE_H5_FILETYPE = 'molecule'

FILE_VERSION_KEY = 'file_version'
CURR_FILE_VERSION = 2

METRICS_GROUP_NAME = 'metrics'
METRICS_GROUP_PATH = '/metrics'

MOLECULE_INFO_COLUMNS = {
    'barcode'                    : np.uint64,   # Up to 31-mers
    'gem_group'                  : np.uint16,   # Up to 65k (256 is not enough for mega-scale datasets)
    'umi'                        : np.uint32,   # Up to 15-mers
    'gene'                       : np.uint32,   # Up to 4e9 genes; set to 1+max_genes if not conf
    'genome'                     : np.uint8,    # Up to 255 genomes
    'reads'                      : np.uint32,   # Up to 4e9 reads/mol
    'nonconf_mapped_reads'       : np.uint32,   # Up to 4e9 reads/mol
    'unmapped_reads'             : np.uint32,   # Up to 4e9 reads/mol
    'umi_corrected_reads'        : np.uint32,   # Reads where UMI sequence was corrected
    'barcode_corrected_reads'    : np.uint32,   # Reads where BC sequence was corrected
    'conf_mapped_uniq_read_pos'  : np.uint32,   # Confidently mapped reads where aligned position is unique
}

MOLECULE_REF_COLUMNS = ['genome_ids', 'gene_ids', 'gene_names']

# NOTE - key columns should be sorted in this order (i.e. gem group column is most significant, etc.)
MOLECULE_KEY_COLUMNS = ['gem_group', 'barcode', 'genome', 'gene', 'umi']

# Preserve contiguity of these when chunking a MoleculeCounter
CHUNK_COLUMNS = ['barcode', 'gem_group']

Molecule = namedtuple('Molecule',
                      ['barcode', # Barcode string, e.g. ACGT-1
                      'genome',  # String genome name, e.g. GRCh38
                      'gene_id', # String genome id, e.g., ENSG..
                      'reads',   # Num confidently mapped reads
                       ]
)

# Common top-level metrics
CHEMISTRY_DESC_METRIC = 'chemistry_description'
BC_LENGTH_METRIC = 'chemistry_barcode_read_length'
BC_WHITELIST_METRIC = 'chemistry_barcode_whitelist'
GEM_GROUPS_METRIC = 'gem_groups'
IS_AGGREGATED_METRIC = 'is_aggregated'

# Per-gem-group metrics
GG_TOTAL_READS_METRIC = 'total_reads'
GG_DOWNSAMPLED_READS_METRIC = 'downsampled_reads'
GG_CONF_MAPPED_FILTERED_BC_READS_METRIC = 'conf_mapped_filtered_bc_reads'
GG_RECOVERED_CELLS_METRIC = 'recovered_cells'
GG_FORCE_CELLS_METRIC = 'force_cells'

class MoleculeCounter:
    """ Streams a list of tuples w/named elements to or from an h5 file """
    def __init__(self):
        self.file_version = None
        self.h5 = None
        self.columns = {}
        self.ref_columns = {}

    @staticmethod
    def get_umi_bits():
        return np.dtype(MOLECULE_INFO_COLUMNS['umi']).itemsize * 8

    @staticmethod
    def get_barcode_bits():
        return np.dtype(MOLECULE_INFO_COLUMNS['barcode']).itemsize * 8

    def get_barcode_length(self):
        return self.get_metric(BC_LENGTH_METRIC)

    def get_barcode_whitelist(self):
        return self.get_metric(BC_WHITELIST_METRIC)

    def get_chemistry_description(self):
        return self.get_metric(CHEMISTRY_DESC_METRIC)

    def get_gem_groups(self):
        return self.get_metric(GEM_GROUPS_METRIC).keys()

    def is_aggregated(self):
        ret = self.get_metric(IS_AGGREGATED_METRIC)
        return ret if ret is not None else False

    @staticmethod
    def get_key_columns():
        return MOLECULE_KEY_COLUMNS

    @staticmethod
    def get_data_columns():
        return set(MOLECULE_INFO_COLUMNS.keys()) - set(MOLECULE_KEY_COLUMNS)

    @staticmethod
    def get_record_bytes():
        return sum([np.dtype(x).itemsize for x in MOLECULE_INFO_COLUMNS.values()])

    @staticmethod
    def estimate_mem_gb(chunk_len):
        mol_entries_per_gb = int(1e9 / MoleculeCounter.get_record_bytes())
        return max(cr_constants.MIN_MEM_GB, chunk_len / mol_entries_per_gb)

    @staticmethod
    def open(filename, mode, start=None, length=None):
        assert mode == 'r' or mode == 'w'

        mc = MoleculeCounter()

        if mode == 'w':
            assert start is None
            assert length is None
            filters = tables.Filters(complevel = cr_constants.H5_COMPRESSION_LEVEL)
            mc.h5 = tables.open_file(filename, mode = 'w', title = '10X', filters = filters)
            mc.h5.set_node_attr('/', FILE_VERSION_KEY, CURR_FILE_VERSION)
            mc.h5.set_node_attr('/', cr_constants.H5_FILETYPE_KEY, MOLECULE_H5_FILETYPE)
            mc.h5.create_group('/', METRICS_GROUP_NAME)

            for name, col_type in MOLECULE_INFO_COLUMNS.iteritems():
                atom = tables.Atom.from_dtype(np.dtype(col_type))
                # Create an (array, element_buffer) tuple
                # where element_buffer is a len=1 numpy array
                # designed to avoid excess allocations
                mc.columns[name] = (mc.h5.create_earray(mc.h5.root, name, atom, (0,)),
                                      np.array([0], dtype=np.dtype(col_type)))

        elif mode == 'r':
            mc.h5 = tables.open_file(filename, mode = 'r')
            try:
                mc.file_version = mc.h5.get_node_attr('/', FILE_VERSION_KEY)
            except AttributeError:
                mc.file_version = 1 # V1 doesn't have version field

            for node in mc.h5.walk_nodes('/', 'Array'):
                if node.name in MOLECULE_INFO_COLUMNS:
                    if start is None:
                        assert length is None
                        mc.columns[node.name] = (node, None)
                    else:
                        assert length is not None
                        mc.columns[node.name] = (node[start:(start+length)], None)
                elif node.name in MOLECULE_REF_COLUMNS:
                    mc.ref_columns[node.name] = node
                else:
                    raise AttributeError("Illegal column: %s" % node.name)

        return mc

    def nrows(self):
        col = self.get_column_lazy(MoleculeCounter.get_key_columns()[0])
        if isinstance(col, np.ndarray):
            # handle case where MoleculeCounter is loaded in chunks, so data is not lazy-loaded
            return col.shape[0]
        else:
            return col.nrows

    def get_chunk_key(self, idx):
        return tuple(self.get_column_lazy(col)[idx] for col in CHUNK_COLUMNS)

    def set_metric(self, key, value):
        self.h5.set_node_attr(METRICS_GROUP_PATH, key, value)

    def get_metric(self, key):
        try:
            value = self.h5.get_node_attr(METRICS_GROUP_PATH, key)
        except AttributeError:
            value = None
        return value

    def set_all_metrics(self, metrics):
        for (k,v) in metrics.iteritems():
            self.set_metric(k, v)

    def get_all_metrics(self):
        # this is the only way to iterate over attribues in pytables :(
        # NOTE: this may be slow, as there are potentially thousands of metrics
        group = self.h5.get_node(METRICS_GROUP_PATH)
        attrset = group._v_attrs
        metrics = {k: attrset[k] for k in attrset._f_list()}
        return metrics

    def add(self, **kwargs):
        """ Add a tuple spread out across the column arrays """
        for key, value in kwargs.iteritems():
            array, element_buf = self.columns[key]
            element_buf[0] = value
            array.append(element_buf)

    def add_many(self, col_name, values):
        """ Append an array of values to a column """
        array, _ = self.columns[col_name]
        array.append(values)

    def get_column_lazy(self, col_name):
        """ Retrieve column. Depending on how the file was opened,
        this may only be a file view instead of a full array. """
        array, _ = self.columns[col_name]
        return array

    def get_column_chunked(self, col_name, chunk_size):
        data = self.get_column_lazy(col_name)
        for i in xrange(0, len(data), chunk_size):
            yield data[i:(i+chunk_size)]

    def get_column(self, col_name):
        """ Retrieve an entire column of data. """
        array = self.get_column_lazy(col_name)
        return array[:]

    def set_ref_column(self, col_name, values):
        assert col_name in MOLECULE_REF_COLUMNS
        self.ref_columns[col_name] = self.h5.create_carray(self.h5.root, col_name, obj=np.array(values))

    def get_ref_column(self, col_name):
        array = self.ref_columns[col_name]
        return array[:]

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        self.h5.close()

    def save(self):
        self.h5.close()

    @staticmethod
    def concatenate(out_filename, in_filenames, metrics=None):
        # Append each column from each input h5 to the output h5
        out_mc = MoleculeCounter.open(out_filename, mode='w')
        ref_set = False
        for in_filename in in_filenames:
            in_mc = MoleculeCounter.open(in_filename, mode='r')
            # if no metrics specified, copy them from the first file
            if metrics is None:
                metrics = in_mc.get_all_metrics()
            for name, array_tuple in in_mc.columns.iteritems():
                h5_array, _ = array_tuple
                out_mc.add_many(name, h5_array[:])
            if not ref_set: # only set once
                for name, h5_array in in_mc.ref_columns.iteritems():
                    out_mc.set_ref_column(name, h5_array[:])
                ref_set = True
            in_mc.close()
        out_mc.set_all_metrics(metrics)
        out_mc.save()

    @staticmethod
    def concatenate_sort(out_filename, in_filenames, sort_cols, metrics=None):
        in_mcs = [MoleculeCounter.open(f, 'r') for f in in_filenames]
        out_mc = MoleculeCounter.open(out_filename, mode='w')
        if metrics is None:
            metrics = in_mcs[0].get_all_metrics()
        out_mc.set_all_metrics(metrics)
        for col, array in in_mcs[0].ref_columns.iteritems():
            out_mc.set_ref_column(col, array[:])
        sort_array = []
        # reverse sort columns so they get sorted in the right order
        for col in reversed(sort_cols):
            sort_array.append(np.concatenate([mc.get_column(col) for mc in in_mcs]))
        sort_index = np.lexsort(sort_array)
        for col in MOLECULE_INFO_COLUMNS:
            col_sorted = np.concatenate([mc.get_column(col) for mc in in_mcs])[sort_index]
            out_mc.add_many(col, col_sorted)
        for mc in in_mcs:
            mc.close()
        out_mc.save()

    def find_last_occurrence_of_chunk_key(self, from_row):
        num_rows = self.nrows()
        initial_chunk_key = self.get_chunk_key(from_row)
        for i in xrange(from_row, num_rows):
             chunk_key = self.get_chunk_key(i)
             if not chunk_key == initial_chunk_key:
                 return i - 1
        return num_rows - 1

    def bisect(self, query, key_func):
        """ Performs a binary search to find the leftmost insertion point of query.
        Takes a key function, where key_func(i) = the value to compare to at index i."""
        num_rows = self.nrows()
        lo = 0
        hi = num_rows
        exists = True
        while True:
            i = (hi + lo) / 2
            curr = key_func(i)
            if curr == query:
                break
            elif hi - lo <= 1:
                # non-matching case
                exists = False
                break
            elif curr < query:
                lo = i
            else:
                hi = i

        if exists:
            # backtrack to first occurrence
            for j in xrange(i, -1, -1):
                 curr = key_func(j)
                 if curr != query:
                     return j + 1
        return 0

    def get_chunks_from_partition(self, values, key_func):
        """ Get chunks by partitioning on the specified values."""
        starts = [0] + [self.bisect(val, key_func) for val in values[1:]]
        n = len(starts)
        for i in xrange(n):
            chunk_start = starts[i]
            chunk_end = starts[i+1] if i+1 < n else self.nrows()
            yield (chunk_start, chunk_end - chunk_start)

    def get_chunks_by_gem_group(self):
        """ Return exactly one chunk per gem group."""
        gem_group_arr = self.get_column('gem_group')
        # verify gem groups are sorted
        assert np.all(np.diff(gem_group_arr)>=0)
        unique_ggs = np.unique(gem_group_arr)
        gg_key = lambda i: gem_group_arr[i]
        chunk_iter = self.get_chunks_from_partition(unique_ggs, gg_key)
        for (gg, chunk) in zip(unique_ggs, chunk_iter):
            yield (gg, chunk[0], chunk[1])

    def get_chunks(self, target_chunk_len, preserve_boundaries=True):
        """ Get chunks, optionally preserving boundaries defined by get_chunk_key().
            Yields (chunk_start, chunk_len) which are closed intervals """
        num_rows = self.nrows()
        chunk_start, chunk_end = 0, 0
        while chunk_end < (num_rows - 1):
            target_chunk_end = min(num_rows - 1, chunk_start + target_chunk_len - 1)
            chunk_end = self.find_last_occurrence_of_chunk_key(target_chunk_end) if preserve_boundaries else target_chunk_end
            chunk_len = 1 + chunk_end - chunk_start
            yield (chunk_start, chunk_len)
            chunk_start = 1 + chunk_end

    @staticmethod
    def compress_barcode_seq(x):
        return cr_utils.compress_seq(x, MoleculeCounter.get_barcode_bits())

    @staticmethod
    def compress_gem_group(x):
        return MOLECULE_INFO_COLUMNS['gem_group'](x)

    @staticmethod
    def compress_umi_seq(x):
        return cr_utils.compress_seq(x, MoleculeCounter.get_umi_bits())

    def get_cdna_mol_counts_per_gene(self, gene_index, remove_none_gene=True):
        mol_genes = self.get_column('gene')

        num_genes = len(gene_index.get_genes())
        gene_counts = np.bincount(mol_genes, minlength=num_genes + 1)
        if remove_none_gene:
            gene_counts = gene_counts[:num_genes]

        return gene_counts

    @staticmethod
    def get_metrics_from_summary(summary, gem_groups, total_recovered_cells=None, total_force_cells=None):
        """ Extract relevant metrics from a summary dict."""
        mol_metrics = {}

        version_metrics = ['cellranger_version', 'reference_mkref_version', 'reference_fasta_hash', 'reference_gtf_hash']
        for m in version_metrics:
            mol_metrics[m] = summary[m]

        chemistry_metrics = [m for m in summary if m.startswith('chemistry')]
        for m in chemistry_metrics:
            mol_metrics[m] = summary[m]

        gem_group_metrics = {}
        for gg in gem_groups:
            total_read_metric = 'total_reads' if len(gem_groups) == 1 else '%s_total_reads_per_gem_group' % gg
            recovered_cells = total_recovered_cells / len(gem_groups) if total_recovered_cells is not None else None
            force_cells = total_force_cells / len(gem_groups) if total_force_cells is not None else None
            gem_group_metrics[gg] = {
                GG_TOTAL_READS_METRIC: summary[total_read_metric],
                GG_RECOVERED_CELLS_METRIC: recovered_cells,
                GG_FORCE_CELLS_METRIC: force_cells,
            }

        mol_metrics[GEM_GROUPS_METRIC] = gem_group_metrics
        return mol_metrics

    @staticmethod
    def naive_concatenate_metrics(mol_h5_list):
        combined_metrics = None
        for mol_h5 in mol_h5_list:
            with MoleculeCounter.open(mol_h5, mode='r') as counter:
                single_metrics = counter.get_all_metrics()
                if combined_metrics is None:
                    combined_metrics = single_metrics
                else:
                    # concatenate new gem groups to the metrics. if it collides with an existing
                    # gem group, the old one will be overwritten.
                    new_gg_metrics = counter.get_metric(GEM_GROUPS_METRIC)
                    combined_metrics[GEM_GROUPS_METRIC].update(new_gg_metrics)

        return combined_metrics

    def get_molecule_iter(self, barcode_length, subsample_rate=1.0):
        """ Return an iterator on Molecule tuples """
        assert subsample_rate >= 0 and subsample_rate <= 1.0

        # Store the previous compressed barcode so we don't have to decompress every single row
        prev_compressed_bc = None
        prev_gem_group = None
        prev_bc = None

        # Load the molecule data
        mol_barcodes = self.get_column('barcode')
        mol_gem_groups = self.get_column('gem_group')
        mol_genome_ints = self.get_column('genome')
        mol_gene_ints = self.get_column('gene')
        mol_reads = self.get_column('reads')

        gene_ids = self.get_ref_column('gene_ids')
        genome_ids = self.get_ref_column('genome_ids')

        if subsample_rate < 1.0:
            mol_reads = np.random.binomial(mol_reads, subsample_rate)

        for compressed_bc, gem_group, genome_int, gene_int, reads in itertools.izip(mol_barcodes,
                                                                                    mol_gem_groups,
                                                                                    mol_genome_ints,
                                                                                    mol_gene_ints,
                                                                                    mol_reads):
                if reads == 0:
                    continue

                # Decompress the cell barcode if necessary
                if compressed_bc == prev_compressed_bc and gem_group == prev_gem_group:
                    bc = prev_bc
                else:
                    bc = cr_utils.format_barcode_seq(self.decompress_barcode_seq(compressed_bc, barcode_length=barcode_length),
                                                     gem_group)
                yield Molecule(barcode=bc,
                               genome=genome_ids[genome_int],
                               gene_id=gene_ids[gene_int],
                               reads=reads)

    @staticmethod
    def get_compressed_bc_iter(barcodes):
        """ Yields compressed barcode tuples that can be compared against
            a MoleculeCounter's data. Useful for filtering a MoleculeCounter by barcode.
        Args: barcodes (iterable) - list of barcode strings (e.g., ACGT-1)
        Yields: (compressed_bc, compressed_gem_group) tuples """

        for barcode in barcodes:
            barcode_seq, gem_group = cr_utils.split_barcode_seq(barcode)
            compressed_bc = MoleculeCounter.compress_barcode_seq(barcode_seq)
            compressed_gg = MoleculeCounter.compress_gem_group(gem_group)
            yield compressed_bc, compressed_gg

    def get_total_raw_reads(self):
        """ Sum 'total_reads' across all gem groups """
        raw_reads = 0
        for _, gg_metrics in self.get_metric(GEM_GROUPS_METRIC).iteritems():
            raw_reads += gg_metrics[GG_TOTAL_READS_METRIC]
        return raw_reads

    def get_total_conf_mapped_filtered_bc_reads(self):
        conf_reads = 0
        for _, gg_metrics in self.get_metric(GEM_GROUPS_METRIC).iteritems():
            conf_reads += gg_metrics[GG_CONF_MAPPED_FILTERED_BC_READS_METRIC]
        return conf_reads

    @staticmethod
    def sum_gem_group_metric(mol_h5_list, metric_name):
        """ Combine a gemgroup-level integer metric across multiple h5 files """
        combined = defaultdict(int)
        for mol_h5 in mol_h5_list:
            with MoleculeCounter.open(mol_h5, mode='r') as counter:
                for gg, metrics in counter.get_metric(GEM_GROUPS_METRIC).iteritems():
                    combined[gg] += metrics[metric_name]
        return combined

    @staticmethod
    def get_total_conf_mapped_reads_in_cells_chunk(filename, filtered_bcs_set, start, length, queue):
        total_mapped_reads = 0
        with MoleculeCounter.open(filename, 'r', start, length) as mc:
            for barcode, gem_group, reads in itertools.izip(mc.get_column('barcode'),
                                                            mc.get_column('gem_group'),
                                                            mc.get_column('reads')):
                if reads < 1:
                    continue
                if (barcode, gem_group) not in filtered_bcs_set:
                    continue
                total_mapped_reads += reads
        queue.put(total_mapped_reads)

    @staticmethod
    def get_total_conf_mapped_reads_in_cells(filename, filtered_barcodes, mem_gb):
        """ Number of confidently mapped reads w/ valid, filtered barcodes.
            Because this is called from a 'split' function, we must stay within the given mem limit.
            NOTE: We re-open the file for each chunk IN ISOLATED PROCESSES
                  due to a possible memory leak in h5py. Tests show the mem usage is nondeterministic, too.
                  https://github.com/h5py/h5py/issues/763 (among many others)
        Args: filtered_barcodes (set) - set of barcode strings (e.g., ACGT-1)
              filename (str) - path to molecule info HDF5 file
              mem_gb (int) - limit memory usage to this value """

        filtered_bcs_set = set(MoleculeCounter.get_compressed_bc_iter(filtered_barcodes))

        entries_per_chunk = int(np.floor(float(mem_gb*1e9)) / MoleculeCounter.get_record_bytes())
        print 'Entries per chunk: %d' % entries_per_chunk

        with MoleculeCounter.open(filename, 'r') as mc:
            num_entries = mc.nrows()

        total_mapped_reads = 0
        for start in xrange(0, num_entries, entries_per_chunk):
            queue = multiprocessing.Queue()
            p = multiprocessing.Process(target=MoleculeCounter.get_total_conf_mapped_reads_in_cells_chunk,
                                        args=(filename, filtered_bcs_set, start, entries_per_chunk, queue))
            p.start()
            p.join()
            total_mapped_reads += queue.get()

        return total_mapped_reads

    def decompress_barcode_seq(self, bc_int, barcode_length=None):
        bc_len = barcode_length or self.get_barcode_length()
        return cr_utils.decompress_seq(bc_int,
                                       bc_len,
                                       MoleculeCounter.get_barcode_bits())
