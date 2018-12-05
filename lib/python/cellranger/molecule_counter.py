#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict, OrderedDict, namedtuple
import cPickle
import h5py
import tables
import itertools
import json
import math
import random
import numpy as np
from cellranger.feature_ref import FeatureReference, FeatureDef
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import cellranger.stats as cr_stats
import cellranger.h5_constants as h5_constants
import cellranger.library_constants as lib_constants
import cellranger.io as cr_io
import tenkit.seq as tk_seq

MOLECULE_H5_FILETYPE = 'molecule'

FILE_VERSION_KEY = 'file_version'
CURR_FILE_VERSION = 3

# Group whose attributes store various metadata
METRICS_GROUP_NAME = 'metrics'

# Group that tracks which barcodes passed filters (usually means they are
# cell-associated)
BARCODE_INFO_GROUP_NAME = 'barcode_info'

MOLECULE_INFO_COLUMNS = OrderedDict([
    ('gem_group'                  , np.uint16),   # Up to 65k
    ('barcode_idx'                , np.uint64),
    ('feature_idx'                , np.uint32),   # Up to 4e9 features
    ('library_idx'                , np.uint16),   # Up to 65k
    ('umi'                        , np.uint32),   # Up to 16-mers
    ('count'                      , np.uint32),   # Up to 4e9 readpairs/mol
])

MOLECULE_REF_COLUMNS = ['barcodes', 'library_info']

# Preserve contiguity of these when chunking a MoleculeCounter
CHUNK_COLUMNS = ['gem_group', 'barcode_idx']

# Common top-level metrics
CHEMISTRY_DESC_METRIC = 'chemistry_description'
BC_LENGTH_METRIC = 'chemistry_barcode_read_length'
BC_WHITELIST_METRIC = 'chemistry_barcode_whitelist'
GEM_GROUPS_METRIC = 'gem_groups'
LIBRARIES_METRIC = 'libraries'
IS_AGGREGATED_METRIC = 'is_aggregated'

# Per-library metrics
TOTAL_READS_METRIC = 'raw_read_pairs'
DOWNSAMPLED_READS_METRIC = 'downsampled_reads'
USABLE_READS_METRIC = 'usable_read_pairs'
GG_RECOVERED_CELLS_METRIC = 'recovered_cells'
GG_FORCE_CELLS_METRIC = 'force_cells'

HDF5_COMPRESSION = 'gzip'
# Number of elements per HDF5 chunk;
#   here, 1 MiB/(32 bytes per element)
HDF5_CHUNK_SIZE = 32768

# Per-barcode metadata. Sparse (not every barcode is listed)
BarcodeInfo = namedtuple('BarcodeInfo', [
    'pass_filter',         # Array-ized list of (barcode_idx, library_idx, genome_idx)
                           # tuples where presence indicates passing filter.
                           # This is a binary 3-d sparse matrix in COO format.
    'genomes',             # Genome strings for genome 0..j
])

BARCODE_INFO_DTYPES =  {
    'pass_filter': 'uint64',
    'genomes': 'str',
}


class MoleculeCounter:
    """ Streams a list of tuples w/named elements to or from an h5 file """
    def __init__(self):
        self.file_version = None
        self.h5 = None
        self.columns = OrderedDict()
        self.ref_columns = OrderedDict()
        self.library_info = None

    def get_barcode_whitelist(self):
        return self.get_metric(BC_WHITELIST_METRIC)

    def get_gem_groups(self):
        return map(int, self.get_metric(GEM_GROUPS_METRIC).keys())

    def is_aggregated(self):
        ret = self.get_metric(IS_AGGREGATED_METRIC)
        return ret if ret is not None else False

    @staticmethod
    def get_column_dtype(k):
        return np.dtype(MOLECULE_INFO_COLUMNS[k])

    @staticmethod
    def get_record_bytes():
        return sum([np.dtype(x).itemsize for x in MOLECULE_INFO_COLUMNS.values()])

    @staticmethod
    def estimate_mem_gb(chunk_len, scale=1.0, cap=True):
        """ Estimate memory usage of this object given a number of records. """
        mol_entries_per_gb = int(1e9 / MoleculeCounter.get_record_bytes())
        mem_gb = round(math.ceil(scale * chunk_len / mol_entries_per_gb))
        if cap:
            return max(h5_constants.MIN_MEM_GB, mem_gb)
        else:
            return mem_gb

    @staticmethod
    def build_barcode_info(filtered_barcodes_by_genome, library_info, barcodes):
        """Generate numpy arrays for per-barcode info
        Args:
          filtered_barcodes_by_genome (dict of str:list(str)): Keys are genomes, values are lists of filtered barcode strings.
          library_info (list of dict): Per-library metadata.
          barcodes (list of str): All barcode sequences (e.g. ['ACGT', ...]
        Returns:
          BarcodeInfo object
        """
        # Replace a genome string with its lexicographical rank
        genome_to_idx = {g:i for i, g in \
                         enumerate(sorted(filtered_barcodes_by_genome.keys()))}

        libraries_for_gem_group = defaultdict(list)
        for lib_idx, lib in enumerate(library_info):
            libraries_for_gem_group[lib['gem_group']].append(lib_idx)

        # Map a barcode sequence to its index into the MoleculeCounter
        #  'barcodes' array
        bc_seq_to_idx = {bc:i for i, bc in enumerate(barcodes)}

        # Populate the "pass filter" array of tuples
        pf_tuples = []
        for genome, bcs in filtered_barcodes_by_genome.iteritems():
            genome_idx = genome_to_idx[genome]
            for bc_str in bcs:
                seq, gg = cr_utils.split_barcode_seq(bc_str)
                barcode_idx = bc_seq_to_idx[seq]

                # FIXME: Assumes no per-library filtering, just per-gem-group
                library_inds = libraries_for_gem_group[gg]
                for library_idx in library_inds:
                    pf_tuples.append((barcode_idx, library_idx, genome_idx))

        if len(pf_tuples) > 0:
            pass_filter = np.array(pf_tuples, dtype=BARCODE_INFO_DTYPES['pass_filter'])
        else:
            pass_filter = np.zeros((0,3), dtype=BARCODE_INFO_DTYPES['pass_filter'])

        assert pass_filter.shape[0] == len(pf_tuples)
        assert pass_filter.shape[1] == 3

        # Sort by barcode index
        pass_filter = pass_filter[np.argsort(pass_filter[:,0]), :]

        return BarcodeInfo(
            pass_filter,
            genomes=sorted(filtered_barcodes_by_genome.keys()),
        )

    @staticmethod
    def get_filtered_barcodes(barcode_info, library_info, barcodes,
                              genome_idx=None, library_type=None):
        """Get a list of filtered barcode strings e.g. ['ACGT-1',...]
        Args:
          barcode_info (BarcodeInfo): Barcode info object.
          library_info (list of dict): Library info.
          barcodes (np.array): Barcode sequences.
          genome_idx (int): Restrict passing definition to this genome. None for no restriction.
          library_type (str): Restrict passing definition to this library type. None for no restriction.
        Returns:
          list of str
        """

        # Without restrictions, assumes passing filter in a single library or genome is sufficient
        # for a barcode to be passing filter overall.

        pass_filter = barcode_info.pass_filter

        pf_barcode_idx = pass_filter[:,0]
        pf_library_idx = pass_filter[:,1]
        pf_genome_idx = pass_filter[:,2]

        mask = np.ones(pass_filter.shape[0], dtype=bool)
        if genome_idx is not None:
            mask &= pf_genome_idx == genome_idx

        if library_type is not None:
            library_inds = np.array([i for i,lib in enumerate(library_info) if lib['library_type'] == library_type],
                                    dtype=MOLECULE_INFO_COLUMNS['library_idx'])
            mask &= np.isin(pf_library_idx, library_inds)
        inds = np.flatnonzero(mask)

        lib_to_gg = np.array([lib['gem_group'] for lib in library_info], dtype='uint64')

        pf_gem_group = lib_to_gg[pf_library_idx[inds]]

        # Take unique, sorted barcodes (sorted by (gem_group, barcode_idx))
        gg_bcs = np.unique(np.column_stack((pf_gem_group, pf_barcode_idx[inds])), axis=0)

        # Create barcode strings
        return [cr_utils.format_barcode_seq(barcodes[gg_bcs[i, 1]],
                                            gg_bcs[i, 0]) for i in xrange(gg_bcs.shape[0])]

    @staticmethod
    def save_barcode_info(bc_info, group):
        """Save barcode info to HDF5.
        Args:
          barcode_info (BarcodeInfo): Data.
          group (h5py.Group): Output group.
        """
        group.create_dataset('pass_filter', data=bc_info.pass_filter,
                             maxshape=(None, bc_info.pass_filter.shape[1]),
                             compression=HDF5_COMPRESSION,
                             shuffle=True)
        cr_io.create_hdf5_string_dataset(group, 'genomes', bc_info.genomes,
                                         compression=HDF5_COMPRESSION,
                                         shuffle=True)

    @staticmethod
    def load_barcode_info(group):
        """Load barcode info from an HDF5 group.
        Args:
          group (h5py.Group): Input group.
        Returns:
          BarcodeInfo object
        """
        return BarcodeInfo(
            pass_filter=group['pass_filter'][:],
            genomes=cr_io.read_hdf5_string_dataset(group['genomes']),
        )

    def get_barcode_info(self):
        return MoleculeCounter.load_barcode_info(self.h5[BARCODE_INFO_GROUP_NAME])

    @staticmethod
    def open(filename, mode, feature_ref=None, barcodes=None, library_info=None,
             barcode_info=None):
        """Open a molecule info object.

        Args:
          filename (str): Filename to open or create
          mode (str): 'r' for reading, 'w' for writing.
          feature_ref (FeatureReference): Required when mode is 'w'.
          barcodes (list of str): All possible barcode sequences. Required when mode is 'w'.
          library_info (list of dict): Library metadata. Required when mode is 'w'.
          barcode_info (BarcodeInfo): Per-barcode metadata.
        Returns:
          MoleculeInfo: A new object
        """
        assert mode == 'r' or mode == 'w'

        mc = MoleculeCounter()

        if mode == 'w':
            if feature_ref is None:
                raise ValueError('Feature reference must be specified when opening a molecule info object for writing')
            if barcodes is None:
                raise ValueError('Barcodes must be specified when opening a molecule info object for writing')
            if library_info is None:
                raise ValueError('Library info must be specified when opening a molecule info object for writing')
            if barcode_info is None:
                raise ValueError('Barcode info must be specified when opening a molecule info object for writing')

            mc.h5 = h5py.File(filename, 'w')
            cr_io.set_hdf5_attr(mc.h5, FILE_VERSION_KEY, CURR_FILE_VERSION)
            cr_io.set_hdf5_attr(mc.h5, h5_constants.H5_FILETYPE_KEY, MOLECULE_H5_FILETYPE)
            cr_io.set_hdf5_attr(mc.h5, FILE_VERSION_KEY, CURR_FILE_VERSION)

            mc.h5.create_group(METRICS_GROUP_NAME)

            # Write feature reference
            fref_group = mc.h5.create_group(h5_constants.H5_FEATURE_REF_ATTR)
            feature_ref.to_hdf5(fref_group)

            # Write barcodes
            # If there are multiple barcode lengths, use the largest for the numpy dtype.
            max_barcode_len = np.max(map(len, barcodes))
            barcode_dtype = np.dtype('S%d' % max_barcode_len)
            mc.h5.create_dataset('barcodes', data=np.fromiter(barcodes, barcode_dtype, count=len(barcodes)), compression=HDF5_COMPRESSION)

            # Write library info
            lib_info_json = json.dumps(library_info, indent=4, sort_keys=True)
            cr_io.create_hdf5_string_dataset(mc.h5, 'library_info', [lib_info_json])

            # Write barcode info
            g = mc.h5.create_group(BARCODE_INFO_GROUP_NAME)
            MoleculeCounter.save_barcode_info(barcode_info, g)

            # Create empty per-molecule datasets
            for name, col_type in MOLECULE_INFO_COLUMNS.iteritems():
                mc.columns[name] = mc.h5.create_dataset(name, (0,),
                                                        maxshape=(None,),
                                                        dtype=col_type,
                                                        compression=HDF5_COMPRESSION,
                                                        chunks=(HDF5_CHUNK_SIZE,))

        elif mode == 'r':
            mc.h5 = h5py.File(filename, 'r')

            try:
                mc.file_version = mc.h5.attrs[FILE_VERSION_KEY]
            except AttributeError:
                mc.file_version = 1 # V1 doesn't have version field

            if mc.file_version < CURR_FILE_VERSION:
                raise ValueError('The molecule info HDF5 file (format version %d) was produced by an older version of Cell Ranger. Reading these files is unsupported.' % mc.file_version)
            if mc.file_version > CURR_FILE_VERSION:
                raise ValueError('The molecule info HDF5 file (format version %d) was produced by an newer version of Cell Ranger. Reading these files is unsupported.' % mc.file_version)

            for key in mc.h5.keys():
                if key in MOLECULE_INFO_COLUMNS:
                    mc.columns[key] = mc.h5[key]
                elif key in MOLECULE_REF_COLUMNS:
                    mc.ref_columns[key] = mc.h5[key]
                elif key == h5_constants.H5_FEATURE_REF_ATTR:
                    mc.feature_reference = FeatureReference.from_hdf5(mc.h5[key])
                elif key == METRICS_GROUP_NAME \
                     or key == BARCODE_INFO_GROUP_NAME:
                    pass
                else:
                    raise AttributeError("Unrecognized dataset key: %s" % key)

            # Load library info
            mc.library_info = json.loads(cr_io.read_hdf5_string_dataset(mc.h5['library_info'])[0])

        return mc

    def nrows(self):
        return self.get_column_lazy(MOLECULE_INFO_COLUMNS.keys()[0]).shape[0]

    def get_chunk_key(self, idx):
        return tuple(self.get_column_lazy(col)[idx] for col in CHUNK_COLUMNS)

    def set_metric(self, key, value):
        """Set a metric. Serialize to Pickle."""
        self.h5[METRICS_GROUP_NAME].attrs[key] = cPickle.dumps(value)

    def get_metric(self, key):
        """Get a metric."""
        try:
            value = cPickle.loads(self.h5[METRICS_GROUP_NAME].attrs[key])
        except KeyError:
            value = None
        return value

    def set_all_metrics(self, metrics):
        for (k,v) in metrics.iteritems():
            self.set_metric(k, v)

    def get_all_metrics(self):
        return {k:cPickle.loads(v) for k,v in self.h5[METRICS_GROUP_NAME].attrs.iteritems()}

    def append_column(self, name, values):
        """Append an array of values to a column."""
        ds = self.columns[name]
        start = len(ds)
        end = start + len(values)
        ds.resize((end,))
        ds[start:end] = values

    def get_column_lazy(self, col_name):
        """ Retrieve column. Depending on how the file was opened,
        this may only be a file view instead of a full array. """
        return self.columns[col_name]

    def get_column(self, col_name):
        """Load an entire column of data into memory"""
        return self.get_column_lazy(col_name)[:]

    def set_ref_column(self, col_name, values):
        assert col_name in MOLECULE_REF_COLUMNS
        self.ref_columns[col_name] = self.h5.create_carray(self.h5.root, col_name, obj=np.array(values))

    def get_ref_column(self, col_name):
        """Load a reference array into memory as a numpy array"""
        return self.get_ref_column_lazy(col_name)[:]

    def get_ref_column_lazy(self, col_name):
        """Get a reference array as a lazy h5py Dataset"""
        return self.ref_columns[col_name]

    def get_feature_ref(self):
        return FeatureReference.from_hdf5(self.h5[h5_constants.H5_FEATURE_REF_ATTR])

    def get_barcodes(self):
        return self.h5['barcodes'][:]

    def get_num_filtered_barcodes_for_library(self, library_idx):
        """Count the number of barcodes passing filter for a library.
        Args:
          library_idx (int): Index of library to count.
        Returns:
          int: Number of filtered barcodes for this library.
        """
        pass_filter = self.h5[BARCODE_INFO_GROUP_NAME]['pass_filter'][:]
        this_lib = np.flatnonzero(pass_filter[:,1] == library_idx)
        barcode_inds = pass_filter[this_lib, 0]
        return len(np.unique(barcode_inds))

    def get_library_info(self):
        return json.loads(self.h5['library_info'][0])

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        self.h5.close()

    def save(self):
        self.h5.close()

    @staticmethod
    def merge_barcode_infos(bc_infos):
        """Merge a BarcodeInfo into another BarcodeInfo.
        Args:
          src_bc_infos (list of BarcodeInfo): Input BarcodeInfos.
        Returns:
          BarcodeInfo"""
        assert len(bc_infos) > 0
        genomes = bc_infos[0].genomes

        # Total number of barcodes with any information
        pfs = []
        for bc_info in bc_infos:
            assert bc_info.pass_filter.shape[1] == 3
            assert bc_info.genomes == genomes
            pfs.append(bc_info.pass_filter)

        new_pf = np.concatenate(pfs, axis=0)

        # Deduplicate the tuples. Unique throws an error on a zero-row array.
        if new_pf.shape[0] > 0:
            new_pf = np.unique(new_pf, axis=0)

        return BarcodeInfo(
            pass_filter=new_pf,
            genomes=genomes,
        )

    @staticmethod
    def concatenate(out_filename, in_filenames, metrics=None):
        """Concatenate MoleculeCounter HDF5 files
        Args:
          out_filename (str): Output HDF5 filename
          in_filenames (list of str): Input HDF5 filenames
          metrics (dict): Metrics to write
        """
        # Load reference info from first file
        first_mc = MoleculeCounter.open(in_filenames[0], 'r')
        feature_ref = first_mc.get_feature_ref()
        barcodes = first_mc.get_barcodes()
        library_info = first_mc.get_library_info()

        feature_ids = [f.id for f in feature_ref.feature_defs]

        # print 'Merging barcode info'
        bc_infos = []
        for filename in in_filenames:
            with MoleculeCounter.open(filename, 'r') as mc:
                bc_infos.append(mc.get_barcode_info())
        merged_bc_info = MoleculeCounter.merge_barcode_infos(bc_infos)

        # print 'Concatenating molecule info files'
        out_mc = MoleculeCounter.open(out_filename, mode='w',
                                      feature_ref=feature_ref,
                                      barcodes=barcodes,
                                      library_info=library_info,
                                      barcode_info=merged_bc_info)

        for filename in in_filenames:
            with MoleculeCounter.open(filename, mode='r') as in_mc:
                # Assert that these data are compatible
                assert in_mc.get_library_info() == library_info
                assert np.array_equal(in_mc.get_barcodes(), barcodes)
                fref = in_mc.get_feature_ref()
                assert [f.id for f in fref.feature_defs] == feature_ids

                # if no metrics specified, copy them from the first file
                if metrics is None:
                    metrics = in_mc.get_all_metrics()

                # Concatenate per-molecule datasets
                for name, ds in in_mc.columns.iteritems():
                    out_mc.append_column(name, ds[:])

        out_mc.set_all_metrics(metrics)
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
        return MoleculeCounter.bisect_static(self.nrows(), query, key_func)

    @staticmethod
    def bisect_static(num_rows, query, key_func):
        """ Performs a binary search to find the leftmost insertion point of query.
        Takes a key function, where key_func(i) = the value to compare to at index i."""
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
        return MoleculeCounter.get_chunks_from_partition_static(self.nrows(), values, key_func)

    @staticmethod
    def get_chunks_from_partition_static(num_rows, values, key_func):
        """ Get chunks by partitioning on the specified values."""
        starts = [0] + [MoleculeCounter.bisect_static(num_rows, val, key_func) for val in values[1:]]
        n = len(starts)
        for i in xrange(n):
            chunk_start = starts[i]
            chunk_end = starts[i+1] if i+1 < n else num_rows
            yield (chunk_start, chunk_end - chunk_start)

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
    def compress_gem_group(x):
        return MOLECULE_INFO_COLUMNS['gem_group'](x)

    @staticmethod
    def compress_umi_seq(x, umi_bits):
        return cr_utils.compress_seq(x, umi_bits)

    @staticmethod
    def get_metrics_from_summary(summary, libraries, total_recovered_cells=None, total_force_cells=None):
        """ Extract relevant metrics from a summary dict."""
        mol_metrics = {}

        version_metrics = ['cellranger_version', 'reference_mkref_version', 'reference_fasta_hash', 'reference_gtf_hash']
        for m in version_metrics:
            mol_metrics[m] = summary[m]

        chemistry_metrics = [m for m in summary if m.startswith('chemistry')]
        for m in chemistry_metrics:
            mol_metrics[m] = summary[m]

        # Per-library values
        lib_metrics = {}
        for lib_idx, lib in enumerate(libraries):
            lib_type_prefix = rna_library.get_library_type_metric_prefix(lib['library_type'])
            summary_name = '%s%s_total_read_pairs_per_library' % (lib_type_prefix, lib_idx)
            lib_metrics[str(lib_idx)] = {
                TOTAL_READS_METRIC: summary[summary_name],
            }

        # Per-gem-group values
        gg_metrics = {}
        gem_groups = sorted([lib['gem_group'] for lib in libraries])
        for gg in gem_groups:
            # Distribute the toplevel expected and forced cells parameters
            #   evenly among the gem groups.
            recovered_cells = total_recovered_cells / len(gem_groups) if total_recovered_cells is not None else None
            force_cells = total_force_cells / len(gem_groups) if total_force_cells is not None else None
            gg_metrics[str(gg)] = {
                GG_RECOVERED_CELLS_METRIC: recovered_cells,
                GG_FORCE_CELLS_METRIC: force_cells,
            }

        mol_metrics[LIBRARIES_METRIC] = lib_metrics
        mol_metrics[GEM_GROUPS_METRIC] = gg_metrics
        return mol_metrics

    @staticmethod
    def naive_concatenate_metrics(mol_h5_list):
        combined_metrics = None
        gg_metrics = {}
        lib_metrics = {}
        for mol_h5 in mol_h5_list:
            with MoleculeCounter.open(mol_h5, mode='r') as counter:
                single_metrics = counter.get_all_metrics()
                if combined_metrics is None:
                    combined_metrics = single_metrics
                    gg_metrics = counter.get_metric(GEM_GROUPS_METRIC)
                    lib_metrics = counter.get_metric(LIBRARIES_METRIC)
                else:
                    # concatenate new gem groups to the metrics. if it collides with an existing
                    # gem group, the old one will be overwritten.
                    new_gg_metrics = counter.get_metric(GEM_GROUPS_METRIC)
                    new_lib_metrics = counter.get_metric(LIBRARIES_METRIC)
                    gg_metrics.update(new_gg_metrics)
                    lib_metrics.update(new_lib_metrics)

        combined_metrics[GEM_GROUPS_METRIC] = gg_metrics
        combined_metrics[LIBRARIES_METRIC] = lib_metrics
        return combined_metrics

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

    def get_raw_read_pairs_per_library(self):
        """ Get raw read pairs per library.
        Returns:
          list of int: Order is by library index
        """
        return [self.get_metric(LIBRARIES_METRIC)[str(li)][TOTAL_READS_METRIC] for li,_ in enumerate(self.library_info)]

    def get_usable_read_pairs_per_library(self):
        """ Get usable read pairs per library.
        Returns:
          list of int: Order is by library index
        """
        return [self.get_metric(LIBRARIES_METRIC)[str(li)][USABLE_READS_METRIC] for li,_ in enumerate(self.library_info)]

    @staticmethod
    def _sum_metric(mol_h5_list, metric_name, metric_type):
        """ Combine a library- or gemgroup- level integer metric across multiple h5 files """
        assert metric_type is LIBRARIES_METRIC or \
            metric_type is GEM_GROUPS_METRIC
        combined = defaultdict(int)
        for mol_h5 in mol_h5_list:
            with MoleculeCounter.open(mol_h5, mode='r') as counter:
                for key, metrics in counter.get_metric(metric_type).iteritems():
                    combined[key] += metrics[metric_name]
        return combined

    @staticmethod
    def sum_library_metric(mol_h5_list, metric_name):
        return MoleculeCounter._sum_metric(mol_h5_list, metric_name, LIBRARIES_METRIC)

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
    def convert_v2_to_v3(v2_mole_info_h5, out_v3_mole_info_h5):
        """
        Given the input v2 molecule info h5 file, convert it into v3 file.
        """
        def get_v2_metrics(h5_file):
            group = tables.open_file(h5_file, 'r').get_node('/metrics')
            attrset = group._v_attrs
            return {k: attrset[k] for k in attrset._f_list()}

        def decompress_barcode_seq(x, barcode_length, bits=64):
            x = np.uint64(x)
            assert barcode_length <= (bits/2 - 1)
            if x & (1L << (bits-1)):
                return 'N' * barcode_length
            result = bytearray(barcode_length)
            for i in xrange(barcode_length):
                result[(barcode_length-1)-i] = tk_seq.NUCS[x & np.uint64(0b11)]
                x = x >> np.uint64(2)
            return str(result)
        
        def build_feature_ref(gene_ids, gene_names, genome_index):
            feature_defs = []
            if len(genome_index) == 1:
                genome = genome_index.keys()[0]
                for idx, (gene_id, gene_name) in enumerate(zip(gene_ids, gene_names)):
                    feature_defs.append(FeatureDef(index=idx,
                                                   id=gene_id,
                                                   name=gene_name,
                                                   feature_type=lib_constants.GENE_EXPRESSION_LIBRARY_TYPE,
                                                   tags={'genome': genome}))
            else:
                for idx, (gene_id, gene_name) in enumerate(zip(gene_ids, gene_names)):
                    genome = gene_id.split('_')[0]
                    feature_defs.append(FeatureDef(index=idx,
                                                   id=gene_id,
                                                   name=gene_name,
                                                   feature_type=lib_constants.GENE_EXPRESSION_LIBRARY_TYPE,
                                                   tags={'genome': genome}))

            return FeatureReference(feature_defs, ['genome'])
        
        def get_chunks_by_gem_group(gem_group_arr):
            """ Return exactly one chunk per gem group."""
            # verify gem groups are sorted
            assert np.all(np.diff(gem_group_arr)>=0)
            num_rows = gem_group_arr.shape[0]
            unique_ggs = np.unique(gem_group_arr)
            gg_key = lambda i: gem_group_arr[i]
            chunk_iter = MoleculeCounter.get_chunks_from_partition_static(num_rows, unique_ggs, gg_key)
            for (gg, chunk) in zip(unique_ggs, chunk_iter):
                yield (gg, chunk[0], chunk[1])
        
        random.seed(0)
        np.random.seed(0)

        v2_mc_in = h5py.File(v2_mole_info_h5, 'r')
        v2_metrics = get_v2_metrics(v2_mole_info_h5)

        v2_genome_ids = v2_mc_in['genome_ids']
        v2_genome_name_to_index = {g:i for i, g in enumerate(v2_genome_ids)}
        
        # Feature Ref
        new_feature_ref = build_feature_ref(v2_mc_in['gene_ids'], v2_mc_in['gene_names'], v2_genome_name_to_index)
        
        # barcode whitelist
        barcode_length = v2_metrics[BC_LENGTH_METRIC]
        barcode_whitelist = cr_utils.load_barcode_whitelist(v2_metrics[BC_WHITELIST_METRIC])
        barcode_to_idx = OrderedDict((k, i) for i,k in enumerate(barcode_whitelist))
        gg_total_diversity = len(barcode_whitelist)
        
        v2_genomes = np.asarray(v2_mc_in['genome'], dtype=np.uint8)  # <-> genome information goes into feature_idx in v3
        v2_gene = np.asarray(v2_mc_in['gene'], dtype=MOLECULE_INFO_COLUMNS['feature_idx'])# <-> feature_idx in v3
        v2_conf_mapped_reads = np.asarray(v2_mc_in['reads'], dtype=MOLECULE_INFO_COLUMNS['count']) # <-> count in v3
        v2_barcodes = np.asarray(v2_mc_in['barcode'], dtype=np.uint64) # <-> transit into barcode_idx in v3
        v2_umis = np.asarray(v2_mc_in['umi'], dtype=MOLECULE_INFO_COLUMNS['umi']) # <-> umi in v3
        v2_gem_groups = np.asarray(v2_mc_in['gem_group'], dtype=MOLECULE_INFO_COLUMNS['gem_group']) # <-> gem_group in v3 

        library_info = []
        barcode_info_genomes, barcode_info_pass_filter = [], []
        barcode_idx_list, feature_idx_list, library_idx_list = [], [], [] 
        gem_group_list, count_list, umi_list = [], [], []

        v2_metrics[LIBRARIES_METRIC] = {}
        # each gem_group is a library
        for lib_idx, (gem_group, chunk_start, chunk_len) in enumerate(get_chunks_by_gem_group(v2_gem_groups)):
            library_info.append({
                'gem_group': int(gem_group), 
                'library_id': str(lib_idx), 
                'library_type': lib_constants.GENE_EXPRESSION_LIBRARY_TYPE
            })

            # per library, raw_read_pairs and usable_read_pairs info
            v2_metrics[LIBRARIES_METRIC][str(lib_idx)] = {
                USABLE_READS_METRIC : v2_metrics[GEM_GROUPS_METRIC][gem_group]['conf_mapped_filtered_bc_reads'], 
                TOTAL_READS_METRIC : v2_metrics[GEM_GROUPS_METRIC][gem_group]['total_reads']
            }

            recovered_cells = v2_metrics[GEM_GROUPS_METRIC][gem_group].get(GG_RECOVERED_CELLS_METRIC, None)
            force_cells = v2_metrics[GEM_GROUPS_METRIC][gem_group].get(GG_FORCE_CELLS_METRIC, None)

            chunk_end = chunk_start + chunk_len
            genomes_for_gem_group = v2_genomes[chunk_start:chunk_end]
            bcs_for_gem_group = v2_barcodes[chunk_start:chunk_end]
            reads_for_gem_group = v2_conf_mapped_reads[chunk_start:chunk_end]
            gene_for_gem_group = v2_gene[chunk_start:chunk_end]
            umis_for_gem_group = v2_umis[chunk_start:chunk_end]
            
            for genome_id in v2_genome_ids:
                g_idx = v2_genome_name_to_index[genome_id]
                genome_indices = genomes_for_gem_group == g_idx

                if genome_indices.sum() == 0:
                    # edge case - there's no data for this genome (e.g. empty sample, false barnyard sample, or nothing confidently mapped)
                    continue

                bcs_for_genome = bcs_for_gem_group[genome_indices]
                reads_for_genome = reads_for_gem_group[genome_indices]
                gene_for_genome = gene_for_gem_group[genome_indices]
                umis_for_genome = umis_for_gem_group[genome_indices]

                # only count UMIs with at least one conf mapped read
                umi_conf_mapped_to_genome = reads_for_genome > 0
                bc_breaks = bcs_for_genome[1:] - bcs_for_genome[:-1]
                bc_breaks = np.concatenate(([1], bc_breaks)) # first row is always a break
                bc_break_indices = np.nonzero(bc_breaks)[0]
                unique_bcs = bcs_for_genome[bc_break_indices]
                umis_per_bc = np.add.reduceat(umi_conf_mapped_to_genome, bc_break_indices)
                
                if force_cells is not None:
                    top_bc_indices, _, _ = cr_stats.filter_cellular_barcodes_fixed_cutoff(umis_per_bc, force_cells)
                else:
                    top_bc_indices, _, _ = cr_stats.filter_cellular_barcodes_ordmag(umis_per_bc, recovered_cells, gg_total_diversity)
                
                # barcode info
                barcode_seq_to_idx = {b:barcode_to_idx[decompress_barcode_seq(b, barcode_length)] for b in unique_bcs}

                barcode_info_genomes.append(genome_id)
                for b in unique_bcs[top_bc_indices]:
                    barcode_info_pass_filter.append((barcode_seq_to_idx[b], lib_idx, g_idx))  

                # data
                barcode_idx_list.append(np.vectorize(barcode_seq_to_idx.get)(bcs_for_genome))
                count_list.append(reads_for_genome)
                gem_group_list.append(np.full(reads_for_genome.shape[0], gem_group, dtype=MOLECULE_INFO_COLUMNS['gem_group']))
                library_idx_list.append(np.full(reads_for_genome.shape[0], lib_idx, dtype=MOLECULE_INFO_COLUMNS['library_idx']))
                feature_idx_list.append(gene_for_genome)
                umi_list.append(umis_for_genome)
        
        new_barcode_info = BarcodeInfo(
            pass_filter=np.array(barcode_info_pass_filter, dtype=BARCODE_INFO_DTYPES['pass_filter']),
            genomes=barcode_info_genomes,
        )

        with MoleculeCounter.open(out_v3_mole_info_h5, 'w',
                                  feature_ref=new_feature_ref,
                                  barcodes=barcode_whitelist,
                                  library_info=library_info,
                                  barcode_info=new_barcode_info,
        ) as out_mc:
            out_mc.append_column('barcode_idx', np.concatenate(barcode_idx_list))
            out_mc.append_column('count', np.concatenate(count_list))
            out_mc.append_column('feature_idx', np.concatenate(feature_idx_list))
            out_mc.append_column('gem_group', np.concatenate(gem_group_list))
            out_mc.append_column('umi', np.concatenate(umi_list))
            # library_idx is the same as gem_group_list 
            out_mc.append_column('library_idx', np.concatenate(library_idx_list))
            out_mc.set_all_metrics(v2_metrics)
        
        return
