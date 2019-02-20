#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
cimport cython
import collections
import cPickle
import h5py as h5
import itertools
import json
import math
import numpy as np
import operator
import random
import tenkit.fasta as tk_fasta
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.library_constants as lib_constants
import cellranger.reference as cr_reference
import cellranger.rna.library as rna_library
import cellranger.stats as cr_stats
import cellranger.utils as cr_utils
import cellranger.io as cr_io

# Metrics
class Metric:
    def __init__(self, report_type=cr_constants.DEFAULT_REPORT_TYPE, **kwargs):
        self.report_type = report_type
        self.active = kwargs.get('always_active', False)

    def add(self, elem):
        self.active = True

    def merge(self, metric):
        raise NotImplementedError

    def report(self):
        raise NotImplementedError

# Dictionary metrics i.e. histograms, element counts
class DictionaryMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = collections.Counter()

    def add(self, elem, value=1):
        Metric.add(self, elem)
        self.d[elem] += value

    def merge(self, metric):
        for elem in metric.d:
            self.d[elem] += metric.d[elem]

    def report(self):
        return self.d

class PercentDictionaryMetric(DictionaryMetric):
    def report(self):
        total = sum(self.d.values())

        d = {}
        for key, value in self.d.iteritems():
            d[key] = tk_stats.robust_divide(float(value), float(total))
        return d

class CountDistinctIntegerMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.counts = None

    def add(self, value):
        raise NotImplementedError

    def add_many(self, counts):
        """ Takes an array of integer counts """
        if self.counts is None:
            self.counts = np.copy(counts)
        else:
            self.counts += counts
        self.active = True

    def merge(self, metric):
        if metric.counts is not None:
            self.add_many(metric.counts)

    def report(self):
        if self.counts is None:
            return 0
        else:
            return np.count_nonzero(self.counts)


class HistogramMetric(DictionaryMetric):
    def __init__(self, cutoffs, **kwargs):
        DictionaryMetric.__init__(self, **kwargs)
        self.cutoffs = cutoffs
        for cutoff in self.cutoffs:
            self.d[cutoff] = 0

    def add(self, elem, value=1):
        Metric.add(self, elem)
        for cutoff in reversed(self.cutoffs):
            if elem >= cutoff:
                self.d[cutoff] += value
                break

class NumpyHistogramMetric(Metric):
    """ A histogram with bins of size 1 and a final bin containing elements > max_value """
    def __init__(self, max_value=10, **kwargs):
        Metric.__init__(self, **kwargs)
        self.counts = np.zeros(2+max_value, dtype=np.int_)
        self.max_value = max_value

    def add(self, elem):
        assert elem >= 0
        Metric.add(self, elem)
        if elem > self.max_value:
            self.counts[-1] += 1
        else:
            self.counts[elem] += 1

    def add_many(self, elems):
        self.active = True
        elems = np.copy(elems).astype(np.int_)
        elems[elems > self.max_value] = 1 + self.max_value
        self.counts += np.bincount(elems, minlength=len(self.counts))

    def merge(self, metric):
        assert self.max_value == metric.max_value
        self.counts += metric.counts

    def report(self):
        d = {str(k):int(v) for k, v in itertools.izip(xrange(0, 1 + self.max_value), self.counts)}
        d[">%d" % self.max_value] = int(self.counts[-1])
        return d

class SequenceDistributionMetric(Metric):
    """ Count base frequencies at k positions """
    def __init__(self, k, **kwargs):
        Metric.__init__(self, **kwargs)
        assert k > 0
        self.k = k
        self.counts = np.zeros((len(tk_seq.NUCS), self.k), dtype=np.int_)

    def add(self, kmer, value=1):
        Metric.add(self, kmer)
        self.counts[[tk_seq.NUCS_INVERSE[nuc] for nuc in kmer], np.arange(0, self.k)] += int(value)

    def merge(self, metric):
        self.counts += metric.counts

    def report(self):
        d = collections.defaultdict(lambda: collections.defaultdict(int))
        for col, pos in enumerate(xrange(-self.k/2, self.k/2)):
            for nuc, row in tk_seq.NUCS_INVERSE.iteritems():
                d[str(pos)][nuc] = self.counts[row, col]
        return d

class MedianMetric(DictionaryMetric):
    def report(self):
        return cr_stats.compute_median_from_distribution(self.d)

class IQRMetric(DictionaryMetric):
    def report(self):
        return cr_stats.compute_iqr_from_distribution(self.d)

class TopNMetric(DictionaryMetric):
    def __init__(self, topN, **kwargs):
        DictionaryMetric.__init__(self, **kwargs)
        self.topN = topN

    def report(self):
        items = sorted(self.d.items(), key=operator.itemgetter(1), reverse=True)
        topN = {key: value for key, value in items[:self.topN]}
        return topN

class SubsampledTopNMetric(TopNMetric):
    def __init__(self, topN, sample_rate, **kwargs):
        TopNMetric.__init__(self, topN, **kwargs)
        self.sample_rate = sample_rate

    def add(self, elem, value=1):
        if cr_utils.downsample(self.sample_rate):
            TopNMetric.add(self, elem, value)

class EffectiveDiversityMetric(DictionaryMetric):
    def report(self):
        counts = np.array(self.d.values())
        return cr_stats.effective_diversity(counts)

class SetMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = set()

    def add(self, elem):
        Metric.add(self, elem)
        self.d.add(elem)

    def merge(self, metric):
        self.d.update(metric.d)

    def report(self):
        return len(self.d)

# Stats metrics i.e. mean, total count, percent, sum, stddev
class StatsMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.m1 = 0
        self.m0 = 0

    def merge(self, metric):
        self.m1 += metric.m1
        self.m0 += metric.m0

class PercentMetric(StatsMetric):
    def add(self, elem, filter=False):
        Metric.add(self, elem)
        self.m0 += elem
        if filter:
            self.m1 += elem

    def set_value(self, numerator, denominator):
        self.active = True
        self.m0 = denominator
        self.m1 = numerator

    def report(self):
        return float(self.m1) / float(self.m0) if self.m0 > 0 else 0

class RateMetric(StatsMetric):
    def add(self, elem, numerator=False):
        Metric.add(self, elem)
        if numerator:
            self.m1 += elem
        else:
            self.m0 += elem

    def report(self):
        return float(self.m1) / float(self.m0) if self.m0 > 0 else 0

class CountMetric(StatsMetric):
    def add(self, elem):
        Metric.add(self, elem)
        self.m0 += 1

    def set_value(self, value):
        self.active = True
        self.m0 = value

    def report(self):
        return self.m0

class MeanMetric(StatsMetric):
    def add(self, elem, weight=1):
        Metric.add(self, elem)
        self.m1 += elem * weight
        self.m0 += weight

    def add_many(self, elems):
        self.active = True
        self.m1 += elems.sum()
        self.m0 += len(elems)

    def report(self):
        return float(self.m1) / float(self.m0) if self.m0 > 0 else 0.0

class VarianceMetric(StatsMetric):
    """ Chan et. al (1979) parallel variance """
    def __init__(self, **kwargs):
        StatsMetric.__init__(self, **kwargs)
        self.m2 = 0

    def add(self, elem):
        Metric.add(self, elem)
        self.m0 += 1
        delta = float(elem - self.m1)
        self.m1 += delta / float(self.m0)
        self.m2 += delta * (float(elem - self.m1))

    def add_many(self, elems):
        self.active = True
        elems_mean = elems.mean()
        elems_m2 = np.square(elems - elems_mean).sum()
        delta = elems_mean - self.m1
        self.m1 = (self.m0 * self.m1 + len(elems) * elems_mean) / float(self.m0 + len(elems))
        self.m2 += elems_m2 + delta * delta * float(self.m0) * float(len(elems)) / (self.m0 + len(elems))
        self.m0 += len(elems)

    def merge(self, metric):
        delta = metric.m1 - self.m1
        self.mean = (self.m0 * self.m1 + metric.m0 * metric.m1) / float(self.m0 + metric.m0)
        self.m2 += metric.m2 + delta * delta * float(self.m0) * float(metric.m0) / (self.m0 + metric.m0)
        self.m0 += metric.m0

    def report(self):
        if self.m0 < 2:
            return float('nan')
        else:
            return tk_stats.robust_divide(self.m2, self.m0 - 1)

class CVMetric(VarianceMetric):
    """ Coeffcient of variation: sd/mean """
    def report(self):
        if self.m0 < 2:
            return float('nan')
        else:
            sd = math.sqrt(VarianceMetric.report(self))
            return tk_stats.robust_divide(sd, self.m1)

class DispersionMetric(VarianceMetric):
    """ Method of moments estimator for negative binomial 1/r """
    def report(self):
        if self.m0 < 2:
            return float('nan')
        else:
            var = VarianceMetric.report(self)
            return (var - self.m1) / (self.m1 * self.m1)

# For embedding metrics within one another
class EmbeddedMetric(Metric):
    def __init__(self, metric_cls, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = {}
        self.metric_cls = metric_cls

    def add(self, elem, *args):
        Metric.add(self, elem)
        if elem not in self.d:
            self.d[elem] = self.metric_cls(report_type=self.report_type)
        self.d[elem].add(*args)

    def merge(self, metric):
        for elem in metric.d:
            if elem not in self.d:
                self.d[elem] = metric.d[elem]
            else:
                self.d[elem].merge(metric.d[elem])

    def report(self):
        return {key: metric.report() for key, metric in self.d.iteritems()}

# For metrics binned by a numeric value
class BinnedMetric(EmbeddedMetric):
    def __init__(self, metric_cls, cutoffs, **kwargs):
        EmbeddedMetric.__init__(self, metric_cls, **kwargs)
        self.cutoffs = cutoffs
        for cutoff in cutoffs:
            self.d[cutoff] = self.metric_cls(report_type=self.report_type)

    def add(self, elem, *args):
        Metric.add(self, elem)
        for cutoff in reversed(self.cutoffs):
            if elem >= cutoff:
                self.d[cutoff].add(*args)
                break

METRICS = {
    # Raw fastq metrics
    'total_reads':               (CountMetric, {}),
    'total_read_pairs_per_library': (CountMetric, {'prefixes': ['library_indices']}),
    'total_read_pairs':          (CountMetric, {}),
    'read_bases_with_q30_frac':  (PercentMetric, {}),
    'read_N_bases_frac':         (PercentMetric, {}),
    'top_read_prefixes':         (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),
    'read2_bases_with_q30_frac': (PercentMetric, {}),
    'read2_N_bases_frac':        (PercentMetric, {}),
    'top_read2_prefixes':        (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),

    # Primer fastq metrics
    'read_with_perfect_primer_or_homopolymer_frac':  (PercentMetric, {}),
    'perfect_homopolymer_frac':                      (PercentMetric, {'prefixes': ['nucs']}),
    'perfect_primer_frac':                           (PercentMetric, {'prefixes': ['primer_names']}),

    # Aligned bam metrics
    'good_bc_frac':            (PercentMetric, {}),
    'barcodes_detected':       (SetMetric, {}),
    'corrected_bc_frac':       (PercentMetric, {}),
    'bc_bases_with_q30_frac':  (PercentMetric, {}),
    'bc_N_bases_frac':         (PercentMetric, {}),

    'sample_index_bases_with_q30_frac':   (PercentMetric, {}),
    'sample_index_N_bases_frac':          (PercentMetric, {}),

    'good_umi_frac':            (PercentMetric, {}),
    'corrected_umi_frac':       (PercentMetric, {}),
    'umi_bases_with_q30_frac':  (PercentMetric, {}),
    'umi_N_bases_frac':         (PercentMetric, {}),
    'failed_rescue_frac':       (PercentMetric, {}),

    # Transcript-based metrics
    'read_distance_to_transcript_three_prime': (HistogramMetric, {'args': [cr_constants.READ_POSITION_CUTOFFS]}),
    'read_distance_to_transcript_five_prime': (HistogramMetric, {'args': [cr_constants.READ_POSITION_CUTOFFS]}),

    # Hierarchical (exclusive) sequence filters
    'umi_filter_frac':            (PercentMetric, {'prefixes': ['umi_properties']}),
    'barcode_filter_frac':        (PercentMetric, {'prefixes': ['barcode_properties']}),

    # Non-hierarchical (non-exclusive) sequence properties
    'umi_property_frac':          (PercentMetric, {'prefixes': ['umi_properties']}),
    'barcode_property_frac':      (PercentMetric, {'prefixes': ['barcode_properties']}),

    'median_insert_size':    (MedianMetric, {'prefixes': ['references']}),
    'insert_size_histogram': (HistogramMetric, {'args': [cr_constants.INSERT_SIZE_CUTOFFS], 'prefixes': ['references']}),
    'iqr_insert_size':       (IQRMetric, {'prefixes': ['references']}),

    'reads_frac': (PercentMetric, {'prefixes': ['references', 'regions', 'read_types']}),
    'unmapped_reads_frac': (PercentMetric, {}),
    'antisense_reads_frac': (PercentMetric, {'prefixes': ['references']}),
    'discordant_pairs_frac': (PercentMetric, {'prefixes': ['references']}),

    'dupe_reads_frac': (PercentMetric, {'prefixes': ['references', 'dupe_types']}),
    'dupe_reads_rate': (RateMetric,{'prefixes': ['references', 'dupe_types']}),

    'umi_hamming_distance_per_dupe_group_histogram': (DictionaryMetric, {'prefixes': ['references', 'dupe_types']}),
    'reads_per_dupe_group_histogram': (NumpyHistogramMetric, {'args': [cr_constants.PER_DUPE_GROUP_MAX], 'prefixes': ['references', 'dupe_types']}),
    'umis_per_dupe_group_histogram': (NumpyHistogramMetric, {'args': [cr_constants.PER_DUPE_GROUP_MAX], 'prefixes': ['references', 'dupe_types']}),
    'reads_per_molecule_histogram': (NumpyHistogramMetric, {'args': [cr_constants.PER_DUPE_GROUP_MAX], 'prefixes': ['references', 'dupe_types']}),

    'top_raw_umis':       (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),
    'top_processed_umis': (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_PROCESSED_SEQUENCE_SAMPLE_RATE], 'prefixes': ['references', 'read_types']}),
    'effective_umi_diversity': (EffectiveDiversityMetric, {'prefixes': ['references', 'read_types']}),
    'low_support_umi_reads_frac': (PercentMetric, {}),

    'top_raw_barcodes':            (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),
    'effective_barcode_diversity': (EffectiveDiversityMetric, {'prefixes': ['references', 'read_types']}),

    'transcriptome_conf_mapped_barcoded_reads': (DictionaryMetric, {'kwargs': {'report_type': 'barcodes', 'always_active': True}, 'prefixes': ['library_types', 'references']}),
    'transcriptome_conf_mapped_deduped_barcoded_reads': (DictionaryMetric, {'kwargs': {'report_type': 'barcodes', 'always_active': True}, 'prefixes': ['library_types', 'references']}),

    ## Subsampled metrics
    'subsampled_duplication_frac': (PercentMetric, {'kwargs': {'always_active': True}, 'prefixes': ['references', 'subsample_types', 'subsample_depths']}),

    # Subsampled matrix summary metrics
    'subsampled_filtered_bcs_median_unique_genes_detected': (CountMetric, {'kwargs': {'always_active': True}, 'prefixes': ['references', 'subsample_types', 'subsample_depths']}),
    'subsampled_filtered_bcs_median_counts': (CountMetric, {'kwargs': {'always_active': True}, 'prefixes': ['references', 'subsample_types', 'subsample_depths']}),


}

cdef class _RawFastqSeqMetrics:
    cdef int total_bases
    cdef int called_bases
    cdef int q30_bases
    cdef int n_bases

    def __init__(self):
        self.total_bases = 0
        self.called_bases = 0
        self.q30_bases = 0
        self.n_bases = 0


cdef class _RawFastqMetricsCache:
    cdef int total_reads
    cdef int total_read_pairs
    cdef int primer_or_homopolymer_reads
    cdef dict homopolymer_reads            # nucleotide -> count
    cdef dict primer_reads                 # primer_name -> count
    cdef dict total_read_pairs_per_library # library_index -> count
    cdef dict _seq_types                   # seq_type -> counts


    def __init__(self):
        self.total_reads = 0 # total number of paired-end inserts
        self.total_read_pairs = 0 # total number of inserts with non-empty template sequence on both ends (after stripping BC + UMI)
        self.primer_or_homopolymer_reads = 0
        self.homopolymer_reads = {}
        self.primer_reads = {}
        self.total_read_pairs_per_library = {}
        self._seq_types = {}

    def finalize(self, reporter):
        get_metric = reporter._get_metric_attr
        total_reads = float(self.total_reads)
        total_read_pairs = float(self.total_read_pairs)
        get_metric('total_reads').set_value(total_reads)
        get_metric('total_read_pairs').set_value(total_read_pairs)

        for library_idx, count in self.total_read_pairs_per_library.iteritems():
            get_metric('total_read_pairs_per_library', library_idx).set_value(count)

        cdef _RawFastqSeqMetrics metrics
        for seq_type, metrics in self._seq_types.iteritems():
            total_bases = float(metrics.total_bases)
            called_bases = float(metrics.called_bases)
            q30_bases = float(metrics.q30_bases)
            n_bases = float(metrics.n_bases)
            get_metric('%s_bases_with_q30_frac' % seq_type).set_value(q30_bases, called_bases)
            get_metric('%s_N_bases_frac' % seq_type).set_value(n_bases, total_bases)

        for nuc, count in self.homopolymer_reads.iteritems():
            get_metric('perfect_homopolymer_frac', nuc).set_value(count, total_reads)

        for primer_name, count in self.primer_reads.iteritems():
            get_metric('perfect_primer_frac', primer_name).set_value(count, total_reads)

        primer_or_homopolymer_reads = float(self.primer_or_homopolymer_reads)
        get_metric('read_with_perfect_primer_or_homopolymer_frac').set_value(primer_or_homopolymer_reads, total_reads)

    cdef set_seq_qual_metrics(cache,
                              bytes seq,
                              const unsigned char[:] qvs,
                              bytes seq_type):
        cdef _RawFastqSeqMetrics metrics
        if seq_type in cache._seq_types:
            metrics = cache._seq_types[seq_type]
        else:
            metrics = _RawFastqSeqMetrics()
            cache._seq_types[seq_type] = metrics
        cdef unsigned char bases_q = 2+tk_fasta.ILLUMINA_QUAL_OFFSET
        cdef unsigned char q30_q
        with nogil:
            q30_q = bases_q+28
            for i in range(qvs.shape[0]):
                with cython.boundscheck(False):
                    with cython.wraparound(False):
                        if qvs[i] > bases_q:
                            metrics.called_bases += 1
                            if qvs[i] > q30_q:
                                metrics.q30_bases += 1

        metrics.total_bases += len(seq)
        metrics.n_bases += seq.count('N')


class Reporter:
    def __init__(self, umi_length=None, primers=None,
                 reference_path=None, high_conf_mapq=None, chroms=None,
                 barcode_whitelist=None, barcode_summary=None, barcode_dist=None,
                 gem_groups=None,
                 metrics_dict=METRICS,
                 gene_index=None,
                 umi_min_qual_threshold=0,
                 subsample_types=None,
                 subsample_depths=None,
                 genomes=None,
                 library_types=None,
                 num_libraries=None):
        self.chroms = chroms
        self.umi_length = umi_length
        self.high_conf_mapq = high_conf_mapq
        self.barcode_whitelist = barcode_whitelist
        self.whitelist_set = set(barcode_whitelist) if barcode_whitelist else None
        self.barcode_summary = barcode_summary
        self.barcode_dist = barcode_dist
        self.gem_groups = gem_groups
        if num_libraries is  None:
            self.library_indices = None
        else:
            self.library_indices = list(range(num_libraries))

        self.gene_index = gene_index
        self.umi_min_qual_threshold = umi_min_qual_threshold
        random.seed(0)

        self.barcode_types = cr_constants.BARCODE_TYPES
        self.molecule_types = cr_constants.MOLECULE_TYPES
        self.read_types = cr_constants.READ_TYPES
        self.dupe_types = cr_constants.DUPE_TYPES
        self.umi_properties = cr_constants.UMI_PROPERTIES
        self.barcode_properties = cr_constants.BARCODE_PROPERTIES

        self.regions = cr_constants.REGIONS

        assert not (reference_path is not None and genomes is not None)

        if reference_path is not None:
            self.genomes = cr_utils.get_reference_genomes(reference_path)
            self.references = self.genomes + [lib_constants.MULTI_REFS_PREFIX]
            self.has_multiple_genomes = len(self.genomes) > 1
        elif genomes is not None:
            self.genomes = genomes
            self.references = self.genomes + [lib_constants.MULTI_REFS_PREFIX]
            self.has_multiple_genomes = len(self.genomes) > 1
        else:
            self.genomes = []
            self.references = []
            self.has_multiple_genomes = False

        if library_types:
            self.library_types = map(rna_library.get_library_type_metric_prefix,
                                     library_types)
        else:
            self.library_types = []

        self.subsample_types = subsample_types
        self.subsample_depths = subsample_depths

        if primers is not None:
            self.primers = {}
            self.primer_names = [primer.name for primer in primers]
            for primer_name, primer in zip(self.primer_names, primers):
                seq = primer.seq
                seq_rc = tk_seq.get_rev_comp(primer.seq)
                self.primers[primer_name] = {
                    'seq': seq,
                    'seq_rc': seq_rc,
                }
        else:
            self.primers = {}
            self.primer_names = []

        if umi_length is not None:
            self.umi_primers = {}
            self.poly_t = "T" * min(self.umi_length, cr_constants.UMI_POLYT_SUFFIX_LENGTH)
            for name, seq in self.primers.iteritems():
                self.umi_primers[name] = {
                    'seq': seq['seq'][:self.umi_length],
                    'seq_rc': seq['seq_rc'][:self.umi_length],
                }
        else:
            self.umi_primers = {}
            self.poly_t = None

        self.homopolymers = {}
        self.nucs = tk_seq.NUCS
        for nuc in self.nucs:
            self.homopolymers[nuc] = nuc * cr_constants.HOMOPOLYMER_LENGTH

        self.metrics_dict = metrics_dict
        self.metrics_data = {}
        prefixes = []
        for base, (metric_cls, metric_dict) in self.metrics_dict.iteritems():
            # Track prefixes for passing to the webshim
            prefixes += metric_dict.get('prefixes', [])
            args = metric_dict.get('args', [])
            kwargs = metric_dict.get('kwargs', {})
            for key in self._get_metric_keys(base):
                self.metrics_data[key] = metric_cls(*args, **kwargs)

        self.prefixes = {}
        for prefix in prefixes:
            self.prefixes[prefix] = self.__dict__.get(prefix, [])

        self.reference_path = reference_path
        self.metadata = {}


    def _get_metric_keys(self, name):
        metric_cls, metric_dict = self.metrics_dict[name]
        prefixes = metric_dict.get('prefixes', [])
        kwargs = metric_dict.get('kwargs', {})

        always_active = kwargs.get('always_active', False)

        parts = [[name]]
        for prefix in prefixes:
            prefix = getattr(self, prefix)

            if prefix:
                parts.append(prefix)

        # Check to make sure all specified metrics are present for metrics that are always active
        if always_active and len(parts) != len(prefixes) + 1:
            return []

        # Return the set of keys
        keys = set(itertools.product(*parts))

        # Add bare keys
        keys.add((name,))

        return keys

    def _get_metric_attr(self, *args):
        return self.metrics_data[tuple(args)]

    def get_all_prefixes(self):
        return self.prefixes

    def raw_fastq_cb(self, read1, read2, bc_read, si_read, umi_read, library_idx,
                     skip_metrics=False):
        cdef _RawFastqMetricsCache cache = self.raw_fastq_cache
        get_metric = self._get_metric_attr

        name, seq, qual = read1
        _, seq2, qual2 = read2
        _, bc_seq, bc_qual = bc_read
        _, si_seq, si_qual = si_read
        _, umi_seq, umi_qual = umi_read

        # These both mean the same thing
        cache.total_reads += 1
        cache.total_read_pairs += 1

        if library_idx in cache.total_read_pairs_per_library:
            cache.total_read_pairs_per_library[library_idx] += 1
        else:
            cache.total_read_pairs_per_library[library_idx] = 1
        cache.set_seq_qual_metrics(seq, qual, b'read')

        if skip_metrics:
            return
        # Metrics below this line are for human consumption only

        if seq2:
            cache.set_seq_qual_metrics(seq2, qual2, b'read2')

        if bc_seq:
            cache.set_seq_qual_metrics(bc_seq, bc_qual, b'bc')

        if si_seq:
            cache.set_seq_qual_metrics(si_seq, si_qual, b'sample_index')

        if umi_seq:
            cache.set_seq_qual_metrics(umi_seq, umi_qual, b'umi')

        seq_prefix = seq[:cr_constants.READ_PREFIX_LENGTH]
        get_metric('top_read_prefixes').add(seq_prefix)
        if seq2:
            seq2_prefix = seq2[:cr_constants.READ_PREFIX_LENGTH]
            get_metric('top_read2_prefixes').add(seq2_prefix)

        read_with_perfect_primer_or_homopolymer = False
        for nuc, homopolymer_seq in self.homopolymers.iteritems():
            perfect_homopolymer = homopolymer_seq in seq
            if seq2:
                perfect_homopolymer = perfect_homopolymer or homopolymer_seq in seq2
            read_with_perfect_primer_or_homopolymer = (
                read_with_perfect_primer_or_homopolymer or perfect_homopolymer)

            if perfect_homopolymer:
                if nuc in cache.homopolymer_reads:
                    cache.homopolymer_reads[nuc] += 1
                else:
                    cache.homopolymer_reads[nuc] = 1

        for primer_name, primer in self.primers.iteritems():
            primer_seq = primer['seq']
            primer_seq_rc = primer['seq_rc']
            present = primer_seq in seq or primer_seq_rc in seq
            if seq2:
                present = present or primer_seq in seq2 or primer_seq_rc in seq2

            perfect_primer = present > 0
            read_with_perfect_primer_or_homopolymer = (
                read_with_perfect_primer_or_homopolymer or perfect_primer)

            if present:
                if primer_name in cache.primer_reads:
                    cache.primer_reads[primer_name] += 1
                else:
                    cache.primer_reads[primer_name] = 1

        if read_with_perfect_primer_or_homopolymer:
            cache.primer_or_homopolymer_reads += 1

    def raw_umi_cb(self, seq, qual):
        """ Returns the processed sequence """
        cache = self.raw_umi_cache

        low_min_qual = cr_utils.min_qual_below(qual, self.umi_min_qual_threshold)
        has_n = 'N' in seq
        is_homopolymer = cr_utils.is_homopolymer_seq(seq)
        has_primer = cr_utils.find_any_primers(seq, self.umi_primers)
        has_polyt = seq[:-cr_constants.UMI_POLYT_SUFFIX_LENGTH] == self.poly_t

        cache.total_umis += 1
        if low_min_qual: cache.low_min_qual_umis += 1
        if has_n: cache.has_n_umis += 1
        if is_homopolymer: cache.homopolymer_umis += 1
        if has_primer: cache.primer_umis += 1
        if has_polyt: cache.polyt_umis += 1

        valid = True

        if low_min_qual: cache.filtered_low_min_qual_umis += 1
        valid = valid and not low_min_qual

        if valid and has_n: cache.filtered_has_n_umis += 1
        valid = valid and not has_n

        if valid and is_homopolymer: cache.filtered_homopolymer_umis += 1
        valid = valid and not is_homopolymer

        return seq if valid else None

    def raw_barcode_cb(self, seq, qual):
        """ Returns the processed sequence """
        cache = self.raw_barcode_cache

        missed_whitelist = not cr_utils.is_barcode_on_whitelist(seq, self.whitelist_set)
        has_n = 'N' in seq
        is_homopolymer = cr_utils.is_homopolymer_seq(seq)
        low_min_qual = cr_utils.min_qual_below(qual, cr_constants.BARCODE_MIN_QUAL_THRESHOLD)

        cache.total_bcs += 1
        if missed_whitelist: cache.miss_whitelist_bcs += 1
        if has_n: cache.has_n_bcs += 1
        if is_homopolymer: cache.homopolymer_bcs += 1
        if low_min_qual: cache.low_min_qual_bcs += 1

        valid = True

        if missed_whitelist: cache.filtered_miss_whitelist_bcs += 1
        valid = valid and not missed_whitelist

        return seq if valid else None

    def extract_reads_init(self):
        self.raw_fastq_cache = _RawFastqMetricsCache()

    def extract_reads_finalize(self):
        self.raw_fastq_cache.finalize(self)

    def _set_dupe_metrics(self, read, dupe_type, genome):
        dupe_reads_frac = self._get_metric_attr('dupe_reads_frac', genome, dupe_type)
        dupe_reads_rate = self._get_metric_attr('dupe_reads_rate', genome, dupe_type)
        dupe_reads_frac.add(1, filter=read.is_duplicate)
        dupe_reads_rate.add(1, numerator=read.is_duplicate)

    def mark_dupes_bam_cb(self, read, dupe_type):
        assert self.high_conf_mapq
        if not cr_utils.is_read_dupe_candidate(read, self.high_conf_mapq):
            return

        genome = cr_utils.get_genome_from_read(read, self.chroms, self.genomes)
        self._set_dupe_metrics(read, dupe_type, genome)
        self._set_dupe_metrics(read, dupe_type, lib_constants.MULTI_REFS_PREFIX)

    def mark_dupes_group_cb(self, gene_id, umis, dupe_type):
        total_counts = sum(umis.values())
        total_umis = len(umis)
        if any([count > 1 for count in umis.itervalues()]):
            umi_hamming_distance = 0
        else:
            umi_hamming_distance = cr_utils.get_kmers_hamming_distance(umis.keys())

        for reference in [cr_utils.get_genome_from_str(gene_id, self.genomes), lib_constants.MULTI_REFS_PREFIX]:
            if total_counts > 0:
                reads_per_dupe_group_histogram = self._get_metric_attr(
                    'reads_per_dupe_group_histogram', reference, dupe_type)
                reads_per_dupe_group_histogram.add(total_counts)

            if total_umis > 0:
                umis_per_dupe_group_histogram = self._get_metric_attr(
                    'umis_per_dupe_group_histogram', reference, dupe_type)
                umis_per_dupe_group_histogram.add(total_umis)

            reads_per_molecule_histogram = self._get_metric_attr(
                'reads_per_molecule_histogram', reference, dupe_type)
            for count in umis.itervalues():
                reads_per_molecule_histogram.add(count)

            if umi_hamming_distance is not None:
                umi_hamming_distance_per_dupe_group_histogram = self._get_metric_attr(
                    'umi_hamming_distance_per_dupe_group_histogram', reference, dupe_type)
                umi_hamming_distance_per_dupe_group_histogram.add(umi_hamming_distance)

    def mark_dupes_corrected_cb(self, read):
        if read.is_secondary:
            return

        raw_umi_seq = cr_utils.get_read_raw_umi(read)
        processed_umi_seq = cr_utils.get_read_umi(read)
        self._get_metric_attr('corrected_umi_frac').add(1, filter=cr_utils.is_umi_corrected(raw_umi_seq, processed_umi_seq))

    def count_genes_bam_cb(self, records, library_info, library_prefixes, use_umis=True):
        """Determine whether to generate a UMI count; compute read-level metrics
        Args:
          records (iterable of pyam.AlignedSegment): Records for a single qname
          library_info (list of dict): Library metadata.
          library_prefixes (list of str): List of metric library prefixes, one per library in the BAM
          use_umis (bool): Use UMIs to count.
        Returns:
          A tuple as per the code below
        """
        assert self.high_conf_mapq is not None

        read1, read2, feature1, feature2, bc, umi = None, None, None, None, None, None
        for rec in records:
            if not cr_utils.is_read_dupe_candidate(rec, self.high_conf_mapq,
                                                   use_umis=use_umis):
                continue

            rec_bc = cr_utils.get_read_barcode(rec)
            rec_umi = cr_utils.get_read_umi(rec)
            feature_ids = cr_utils.get_read_gene_ids(rec)

            # We shouldn't see differing barcodes or UMIs
            # among records with the same qname.
            assert bc is None or bc == rec_bc
            assert umi is None or umi == rec_umi

            if rec.is_read2:
                # Normally we wouldn't see >1 primary record for a QNAME
                # but we can if the sequencing data were intentionally duplicated.
                # Prioritize the record that is not a duplicate
                if read2 is not None and cr_utils.is_read_umi_count(read2):
                    continue
                read2, feature2, bc, umi = rec, feature_ids[0], rec_bc, rec_umi
            else:
                if read1 is not None and cr_utils.is_read_umi_count(read1):
                    continue
                read1, feature1, bc, umi = rec, feature_ids[0], rec_bc, rec_umi

        # We didn't find any suitable records
        if read1 is None and read2 is None:
            return False, None, None, None

        # We shouldn't see discordant features between read1 and read2
        assert feature1 is None or feature2 is None or feature1 == feature2

        # Prefer the primary alignment for Read1
        read = read1 if read1 is not None else read2
        feature_id = feature1 if feature1 is not None else feature2

        library_idx = cr_utils.get_read_library_index(read)
        library_prefix = library_prefixes[library_idx]

        if rna_library.has_genomes(library_info[library_idx]['library_type']):
            genome = cr_utils.get_genome_from_read(read, self.chroms, self.genomes)
        else:
            genome = None

        # Update read counts that will end up in the barcode_summary.h5
        if genome is not None:
            conf_mapped_barcode_reads = self._get_metric_attr('transcriptome_conf_mapped_barcoded_reads',
                                                              library_prefix,
                                                              genome)
            conf_mapped_barcode_reads.add(bc)

        # Only report barcode_reads for multi_* if there are multiple genomes
        # or the library type is genomeless
        if self.has_multiple_genomes or genome is None:
            conf_mapped_barcode_reads = self._get_metric_attr('transcriptome_conf_mapped_barcoded_reads',
                                                              library_prefix,
                                                              lib_constants.MULTI_REFS_PREFIX)
            conf_mapped_barcode_reads.add(bc)

        # Count this read
        if use_umis:
            conf_mapped_deduped = cr_utils.is_read_umi_count(read)
        else:
            conf_mapped_deduped = not read.is_duplicate

        # Update UMI counts that will end up in the barcode_summary.h5
        if conf_mapped_deduped:
            if genome is not None:
                conf_mapped_deduped_barcode_reads = self._get_metric_attr(
                    'transcriptome_conf_mapped_deduped_barcoded_reads',
                    library_prefix, genome)
                conf_mapped_deduped_barcode_reads.add(bc)

            # Only report barcode_reads for multi_* if there are multiple genomes
            # or the library type is genomeless
            if self.has_multiple_genomes or genome is None:
                conf_mapped_deduped_barcode_reads = self._get_metric_attr(
                    'transcriptome_conf_mapped_deduped_barcoded_reads',
                    library_prefix, lib_constants.MULTI_REFS_PREFIX)
                conf_mapped_deduped_barcode_reads.add(bc)

            return True, genome, feature_id, bc

        return False, None, None, None

    def _set_insert_size_metrics(self, read, insert_size):
        genome = cr_utils.get_genome_from_read(read, self.chroms, self.genomes)

        if insert_size is not None:
            for reference in [genome, lib_constants.MULTI_REFS_PREFIX]:
                median_insert_size = self._get_metric_attr('median_insert_size', reference)
                insert_size_histogram = self._get_metric_attr('insert_size_histogram', reference)
                median_insert_size.add(insert_size)
                insert_size_histogram.add(insert_size)

                iqr_insert_size = self._get_metric_attr('iqr_insert_size', reference)
                iqr_insert_size.add(insert_size)

    def store_pipeline_metadata(self, version):
        if self.metadata is None:
            self.metadata = {}
        self.metadata[cr_constants.CELLRANGER_VERSION_KEY] = version

    def store_reference_metadata(self, reference_path, ref_type, metric_prefix):
        """ ref_type - string e.g., 'Transcriptome'
            metric_prefix - string e.g., 'vdj' """

        if self.metadata is None:
            self.metadata = {}

        ref_metadata = cr_utils._load_reference_metadata_file(reference_path)

        for key in cr_constants.REFERENCE_METADATA_KEYS:
            value = ref_metadata.get(key, '')
            if value is None:
                value = ''

            # Backward compatibility with old reference metadata jsons that don't contain the type field
            if key == cr_constants.REFERENCE_TYPE_KEY and value == '':
                self.metadata['%s%s' % (metric_prefix, cr_constants.REFERENCE_TYPE_KEY)] = ref_type
                continue

            if np.isscalar(value):
                self.metadata['%s%s' % (metric_prefix, key)] = value
            elif key == cr_constants.REFERENCE_GENOMES_KEY:
                # Special case for genome key
                self.metadata['%s%s' % (metric_prefix, key)] = cr_reference.get_ref_name_from_genomes(value)
            else:
                self.metadata['%s%s' % (metric_prefix, key)] = ', '.join(str(x) for x in value)


    def store_chemistry_metadata(self, chemistry_def):
        """ Store the chemistry definition as metrics in the summary json """
        if self.metadata is None:
            self.metadata = {}
        for key, value in chemistry_def.iteritems():
            self.metadata['chemistry_%s' % key] = value

    def merge(self, reporter):
        for key, metric2 in reporter.metrics_data.iteritems():
            metric1 = self.metrics_data.get(key)
            if isinstance(metric2, Metric):
                if metric1 and metric1.active and metric2.active:
                    metric1.merge(metric2)
                elif metric2.active:
                    self.metrics_data[key] = metric2

    def report(self, report_type):
        summary = {}
        for key, metric in self.metrics_data.iteritems():
            if metric.active and metric.report_type == report_type:

                name = '_'.join([str(x) for x in key[1:] + (key[0],)]) # ensure keys are strings before output
                summary[name] = metric.report()

        if self.metadata is not None and report_type == cr_constants.DEFAULT_REPORT_TYPE:
            summary.update(self.metadata)

        return summary

    def to_json(self):
        '''Return data as safe-for-json dictionary'''
        data = self.report(cr_constants.DEFAULT_REPORT_TYPE)
        return tk_safe_json.json_sanitize(data)

    def report_summary_json(self, filename):
        '''Write to JSON file'''
        with open(filename, 'w') as f:
            tk_safe_json.dump_numpy(self.to_json(), f, pretty=True)

    def report_barcodes_csv(self, filename):
        barcodes = sorted(self._get_metric_attr('barcodes_detected').d)
        with open(filename, 'w') as f:
            for bc in barcodes:
                f.write(bc + '\n')

    def report_barcodes_h5(self, filename):
        data = self.report('barcodes')

        if self.barcode_whitelist:
            bc_sequences = cr_utils.format_barcode_seqs(self.barcode_whitelist, self.gem_groups)
        elif self.barcode_summary is not None:
            bc_sequences = self.barcode_summary
        else:
            # Get all observed bc sequences
            bc_sequences = sorted(list(reduce(lambda x, y: x.union(y.keys()), data.values(), set())))

        # Build the columns for the table
        bc_table_cols = {cr_constants.H5_BC_SEQUENCE_COL: bc_sequences}

        for metric_name, metric_data in data.iteritems():
            counts = np.array([metric_data.get(bc, 0) for bc in bc_sequences], dtype=np.uint32)
            bc_table_cols[metric_name] = counts

        cr_io.write_h5(filename, bc_table_cols)

    def save(self, filename):
        with open(filename, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(filename):
        with open(filename, 'rb') as f:
            return cPickle.load(f)

def merge_reporters(filenames):
    reporter = None
    for filename in filenames:
        tmp_reporter = Reporter.load(filename)
        if reporter is None:
            reporter = tmp_reporter
        else:
            reporter.merge(tmp_reporter)

    return reporter

def merge_jsons(in_filenames, out_filename, dicts=[]):
    """ Merge a list of json files and optional dicts """
    d = cr_utils.merge_jsons_as_dict(in_filenames)
    for data in dicts:
        cr_utils.update_require_unique_key(d, data)

    with open(out_filename, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(d), f, indent=4, sort_keys=True)

def merge_h5(in_filenames, out_filename):
    """ Merge a list of h5 files """
    out_h5 = h5.File(out_filename, 'a')
    for filename in in_filenames:
        if filename is None:
            continue
        in_h5 = h5.File(filename, 'r')
        for name in in_h5.keys():
            # If the dataset already exists,
            # They must be equal or one must be all-zero.
            if name in out_h5.keys():
                src_data, dest_data = in_h5[name][()], out_h5[name][()]
                if src_data.dtype.kind != 'S' and dest_data.dtype.kind != 'S':
                    # Both numeric
                    if not np.any(src_data):
                        # Source is all zero. Do nothing.
                        continue
                    elif not np.any(dest_data):
                        # Dest is all zero. Overwrite.
                        del out_h5[name]
                        h5.h5o.copy(in_h5.id, name, out_h5.id, name)
                    else:
                        # Both non-zero. Assert equality and do nothing.
                        assert np.array_equal(src_data, dest_data)
                else:
                    # Either are non-numeric. Assert equality and do nothing.
                    assert np.array_equal(src_data, dest_data)
            else:
                # Only exists in src. Copy to dest.
                h5.h5o.copy(in_h5.id, name, out_h5.id, name)

    out_h5.flush()
    out_h5.close()
