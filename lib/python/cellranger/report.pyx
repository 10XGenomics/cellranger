#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import collections
import itertools
import json
import math
import operator
import os.path
import pickle
import random
from collections.abc import Iterable
from sys import intern

import numpy as np

import cellranger.constants as cr_constants
import cellranger.rna.library as rna_library
import cellranger.stats as cr_stats
import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
from cellranger.reference_paths import get_ref_name_from_genomes, get_reference_genomes

DEFAULT_REPORT_TYPE = "summary"


# Metrics
class Metric:
    def __init__(self, report_type: str = DEFAULT_REPORT_TYPE, **kwargs):
        self.report_type: str = report_type
        self.active: bool = kwargs.get(_ALWAYS_ACTIVE, False)

    def add(self, elem):
        self.active = True

    def merge(self, metric: "Metric") -> None:
        raise NotImplementedError

    def report(self) -> int | float | dict:
        raise NotImplementedError


# Dictionary metrics i.e. histograms, element counts
class DictionaryMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = {}  # type: dict

    def add(self, elem: int | bytes | str, value: int | float = 1):
        Metric.add(self, elem)
        self.d[elem] = self.d.get(elem, 0) + value

    def merge(self, metric: "DictionaryMetric"):
        for elem, value in metric.d.items():
            self.add(elem, value)

    def report(self) -> dict[int | bytes | str, int | float]:
        return self.d


class PercentDictionaryMetric(DictionaryMetric):
    def report(self):
        total = sum(self.d.values())

        d = {}
        for key, value in self.d.items():
            d[key] = tk_stats.robust_divide(float(value), float(total))
        return d


class CountDistinctIntegerMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.counts = None

    def add(self, value):
        raise NotImplementedError

    def add_many(self, counts):
        """Takes an array of integer counts."""
        if self.counts is None:
            self.counts = np.copy(counts)
        else:
            self.counts += counts
        self.active = True

    def merge(self, metric) -> None:
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
    """A histogram with bins of size 1 and a final bin containing elements > max_value."""

    def __init__(self, max_value=10, **kwargs):
        Metric.__init__(self, **kwargs)
        self.counts = np.zeros(2 + max_value, dtype=np.int_)
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

    def merge(self, metric) -> None:
        assert self.max_value == metric.max_value
        self.counts += metric.counts

    def report(self) -> dict[str, int]:
        d = {str(k): int(v) for k, v in zip(range(1 + self.max_value), self.counts)}
        d[">%d" % self.max_value] = int(self.counts[-1])
        return d


class SequenceDistributionMetric(Metric):
    """Count base frequencies at k positions."""

    def __init__(self, k, **kwargs):
        Metric.__init__(self, **kwargs)
        assert k > 0
        self.k = k
        self.counts = np.zeros((len(tk_seq.NUCS), self.k), dtype=np.int_)

    def add(self, kmer: bytes, value=1):
        Metric.add(self, kmer)
        self.counts[
            [tk_seq.NUCS_INVERSE[kmer[i : i + 1]] for i in range(len(kmer))], np.arange(0, self.k)
        ] += int(value)

    def merge(self, metric: "SequenceDistributionMetric") -> None:
        self.counts += metric.counts

    def report(self) -> dict[str, dict[bytes, int]]:
        d = collections.defaultdict(
            lambda: collections.defaultdict(int)
        )  # type: Dict[str, Dict[bytes, int]]
        for col, pos in enumerate(range(-self.k // 2, self.k // 2)):
            for nuc, row in tk_seq.NUCS_INVERSE.items():
                d[str(pos)][nuc] = self.counts[row, col]
        return d


def compute_percentile_from_distribution(
    counter: dict[int | bytes | str, int | float], percentile: float
) -> float:
    """Takes a Counter object (or value:frequency dict) and computes a single percentile.

    Uses Type 7 interpolation from:
      Hyndman, R.J.; Fan, Y. (1996). "Sample Quantiles in Statistical Packages".
    """
    assert 0 <= percentile <= 100

    n: int = np.sum(list(counter.values())).item()
    h: float = (n - 1) * (percentile / 100.0)
    h_floor: float = np.floor(h).item()
    h_ceil: float = np.ceil(h).item()
    lower_value = None  # Optional[int]

    cum_sum = 0
    for value, freq in sorted(counter.items()):
        cum_sum += freq
        if cum_sum > h_floor and lower_value is None:
            lower_value = value
        if cum_sum > h_ceil:
            return lower_value + (h - h_floor) * (value - lower_value)
    raise AssertionError("Search failed")


# Test for compute_percentile_from_distribution()
# def test_percentile(x, p):
#    c = Counter()
#    for xi in x:
#        c[xi] += 1
#    my_res = np.array([compute_percentile_from_distribution(c, p_i) for p_i in p], dtype=float)
#    numpy_res = np.percentile(x, p)
#    print np.sum(np.abs(numpy_res - my_res))


class MedianMetric(DictionaryMetric):
    def report(self):
        return compute_percentile_from_distribution(self.d, 50)


def compute_iqr_from_distribution(counter: dict[int | bytes | str, int | float]) -> float:
    p25 = compute_percentile_from_distribution(counter, 25)
    p75 = compute_percentile_from_distribution(counter, 75)
    return p75 - p25


class IQRMetric(DictionaryMetric):
    def report(self):
        return compute_iqr_from_distribution(self.d)


class TopNMetric(DictionaryMetric):
    def __init__(self, topN, **kwargs):
        DictionaryMetric.__init__(self, **kwargs)
        self.topN = topN

    def report(self):
        items = sorted(self.d.items(), key=operator.itemgetter(1), reverse=True)
        topN = {key: value for key, value in items[: self.topN]}
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

    def merge(self, metric) -> None:
        self.d.update(metric.d)

    def report(self):
        return len(self.d)


# Stats metrics i.e. mean, total count, percent, sum, stddev
class StatsMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.m1 = 0
        self.m0 = 0

    def merge(self, metric) -> None:
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

    def add_value(self, numerator, denominator):
        self.active = True
        self.m0 += denominator
        self.m1 += numerator

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
    """Chan et. al (1979) parallel variance."""

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
        self.m2 += elems_m2 + delta * delta * float(self.m0) * float(len(elems)) / (
            self.m0 + len(elems)
        )
        self.m0 += len(elems)

    def merge(self, metric) -> None:
        delta = metric.m1 - self.m1
        self.mean = (self.m0 * self.m1 + metric.m0 * metric.m1) / float(self.m0 + metric.m0)
        self.m2 += metric.m2 + delta * delta * float(self.m0) * float(metric.m0) / (
            self.m0 + metric.m0
        )
        self.m0 += metric.m0

    def report(self):
        if self.m0 < 2:
            return float("nan")
        else:
            return tk_stats.robust_divide(self.m2, self.m0 - 1)


class CVMetric(VarianceMetric):
    """Coeffcient of variation: sd/mean."""

    def report(self):
        if self.m0 < 2:
            return float("nan")
        else:
            sd = math.sqrt(VarianceMetric.report(self))
            return tk_stats.robust_divide(sd, self.m1)


class DispersionMetric(VarianceMetric):
    """Method of moments estimator for negative binomial 1/r."""

    def report(self):
        if self.m0 < 2:
            return float("nan")
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

    def merge(self, metric) -> None:
        for elem in metric.d:
            if elem not in self.d:
                self.d[elem] = metric.d[elem]
            else:
                self.d[elem].merge(metric.d[elem])

    def report(self):
        return {key: metric.report() for key, metric in self.d.items()}


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


READ_POSITION_CUTOFFS = range(0, 5100, 100)
PER_DUPE_GROUP_MAX = 100
READ_PREFIX_LENGTH = 10
# Downsample reads when tracking top sequences to constrain mem usage
TOP_RAW_SEQUENCE_SAMPLE_RATE = 0.01
TOP_PROCESSED_SEQUENCE_SAMPLE_RATE = 0.1

_PREFIX_KEY = intern("prefixes")
_ARGS_KEY = intern("args")
_KWARGS_KEY = intern("kwargs")
_BARCODES = intern("barcodes")
_ALWAYS_ACTIVE = intern("always_active")


def _make_default_metrics() -> dict[
    str,
    tuple[
        type[Metric],
        dict[
            str,
            dict[str, str | bool] | list[str] | list[int | float] | list[range | list[int]],
        ],
    ],
]:
    # Common subsections of the default metrics get reused, reducing allocations
    # and making string comparisons in dict/set lookup faster.
    report_type = intern("report_type")
    references = "references"
    read_types = "read_types"
    barcode_quality_prefix = {
        _PREFIX_KEY: ["barcode_properties"]
    }  # type: Dict[str, Union[Dict[str, Union[str, bool]], List[str], List[Union[int, float]], List[Union[range, List[int]]]]]
    umi_prop_prefix = {
        _PREFIX_KEY: ["umi_properties"]
    }  # type: Dict[str, Union[Dict[str, Union[str, bool]], List[str], List[Union[int, float]], List[Union[range, List[int]]]]]

    empty_percent_metric = (PercentMetric, {})
    subsampled = {
        _KWARGS_KEY: {_ALWAYS_ACTIVE: True},
        _PREFIX_KEY: [references, "subsample_types", "subsample_depths"],
    }
    barcode_lib_dict_metric = (
        DictionaryMetric,
        {
            _KWARGS_KEY: {report_type: _BARCODES, _ALWAYS_ACTIVE: True},
            _PREFIX_KEY: ["library_types", references],
        },
    )
    refs_and_dups = [references, "dupe_types"]
    numpy_dup_metric = (
        NumpyHistogramMetric,
        {_ARGS_KEY: [PER_DUPE_GROUP_MAX], _PREFIX_KEY: refs_and_dups},
    )

    return {
        # Raw fastq metrics
        "total_reads": (CountMetric, {}),
        "total_read_pairs_per_library": (CountMetric, {_PREFIX_KEY: ["library_indices"]}),
        "total_read_pairs": (CountMetric, {}),
        "read_bases_with_q30_frac": empty_percent_metric,
        "read_N_bases_frac": empty_percent_metric,
        "top_read_prefixes": (
            SubsampledTopNMetric,
            {_ARGS_KEY: [cr_constants.TOP_N, TOP_RAW_SEQUENCE_SAMPLE_RATE]},
        ),
        "read2_bases_with_q30_frac": empty_percent_metric,
        "read2_N_bases_frac": empty_percent_metric,
        "top_read2_prefixes": (
            SubsampledTopNMetric,
            {_ARGS_KEY: [cr_constants.TOP_N, TOP_RAW_SEQUENCE_SAMPLE_RATE]},
        ),
        # Primer fastq metrics
        "read_with_perfect_primer_or_homopolymer_frac": empty_percent_metric,
        "perfect_homopolymer_frac": (PercentMetric, {_PREFIX_KEY: ["nucs"]}),
        "perfect_primer_frac": (PercentMetric, {_PREFIX_KEY: ["primer_names"]}),
        # Aligned bam metrics
        "good_bc_frac": empty_percent_metric,
        "barcodes_detected": (SetMetric, {}),
        "corrected_bc_frac": empty_percent_metric,
        "bc_bases_with_q30_frac": empty_percent_metric,
        "bc_N_bases_frac": empty_percent_metric,
        "sample_index_bases_with_q30_frac": empty_percent_metric,
        "sample_index_N_bases_frac": empty_percent_metric,
        "good_umi_frac": empty_percent_metric,
        "corrected_umi_frac": empty_percent_metric,
        "umi_bases_with_q30_frac": empty_percent_metric,
        "umi_N_bases_frac": empty_percent_metric,
        "failed_rescue_frac": empty_percent_metric,
        # Transcript-based metrics
        "read_distance_to_transcript_three_prime": (
            HistogramMetric,
            {_ARGS_KEY: [READ_POSITION_CUTOFFS]},
        ),
        "read_distance_to_transcript_five_prime": (
            HistogramMetric,
            {_ARGS_KEY: [READ_POSITION_CUTOFFS]},
        ),
        # Hierarchical (exclusive) sequence filters
        "umi_filter_frac": (PercentMetric, umi_prop_prefix),
        "barcode_filter_frac": (PercentMetric, barcode_quality_prefix),
        # Non-hierarchical (non-exclusive) sequence properties
        "umi_property_frac": (PercentMetric, umi_prop_prefix),
        "barcode_property_frac": (PercentMetric, barcode_quality_prefix),
        "median_insert_size": (MedianMetric, {_PREFIX_KEY: [references]}),
        "insert_size_histogram": (
            HistogramMetric,
            {_ARGS_KEY: [cr_constants.INSERT_SIZE_CUTOFFS], _PREFIX_KEY: [references]},
        ),
        "iqr_insert_size": (IQRMetric, {_PREFIX_KEY: [references]}),
        "reads_frac": (PercentMetric, {_PREFIX_KEY: [references, "regions", read_types]}),
        "unmapped_reads_frac": empty_percent_metric,
        "antisense_reads_frac": (PercentMetric, {_PREFIX_KEY: [references]}),
        "discordant_pairs_frac": (PercentMetric, {_PREFIX_KEY: [references]}),
        "dupe_reads_frac": (PercentMetric, {_PREFIX_KEY: refs_and_dups}),
        "dupe_reads_rate": (RateMetric, {_PREFIX_KEY: refs_and_dups}),
        "umi_hamming_distance_per_dupe_group_histogram": (
            DictionaryMetric,
            {_PREFIX_KEY: refs_and_dups},
        ),
        "reads_per_dupe_group_histogram": numpy_dup_metric,
        "umis_per_dupe_group_histogram": numpy_dup_metric,
        "reads_per_molecule_histogram": numpy_dup_metric,
        "top_raw_umis": (
            SubsampledTopNMetric,
            {_ARGS_KEY: [cr_constants.TOP_N, TOP_RAW_SEQUENCE_SAMPLE_RATE]},
        ),
        "top_processed_umis": (
            SubsampledTopNMetric,
            {
                _ARGS_KEY: [cr_constants.TOP_N, TOP_PROCESSED_SEQUENCE_SAMPLE_RATE],
                _PREFIX_KEY: [references, read_types],
            },
        ),
        "effective_umi_diversity": (
            EffectiveDiversityMetric,
            {_PREFIX_KEY: [references, read_types]},
        ),
        "low_support_umi_reads_frac": empty_percent_metric,
        "top_raw_barcodes": (
            SubsampledTopNMetric,
            {_ARGS_KEY: [cr_constants.TOP_N, TOP_RAW_SEQUENCE_SAMPLE_RATE]},
        ),
        "effective_barcode_diversity": (
            EffectiveDiversityMetric,
            {_PREFIX_KEY: [references, read_types]},
        ),
        "transcriptome_conf_mapped_barcoded_reads": barcode_lib_dict_metric,
        "transcriptome_conf_mapped_deduped_barcoded_reads": barcode_lib_dict_metric,
        ## Subsampled metrics
        "subsampled_duplication_frac": (PercentMetric, subsampled),
        # Subsampled matrix summary metrics
        "subsampled_filtered_bcs_median_unique_genes_detected": (CountMetric, subsampled),
        "subsampled_filtered_bcs_mean_unique_genes_detected": (CountMetric, subsampled),
        "subsampled_filtered_bcs_median_counts": (CountMetric, subsampled),
        "subsampled_filtered_bcs_mean_counts": (CountMetric, subsampled),
    }


METRICS = _make_default_metrics()

# UMI properties
LOW_MIN_QUAL_UMI_FILTER = "low_min_qual"
HAS_N_UMI_FILTER = "has_n"
HOMOPOLYMER_UMI_FILTER = "homopolymer"
# PRIMER_UMI_FILTER = "primer"
# POLYT_UMI_FILTER = "polyt_suffix"

# Barcode properties
# MISS_WHITELIST_BARCODE_FILTER = "miss_whitelist"
# HAS_N_BARCODE_FILTER = HAS_N_UMI_FILTER
# HOMOPOLYMER_BARCODE_FILTER = HOMOPOLYMER_UMI_FILTER
# LOW_MIN_QUAL_BARCODE_FILTER = LOW_MIN_QUAL_UMI_FILTER
BARCODE_MIN_QUAL_THRESHOLD: int = 10

# BARCODE_PROPERTIES = [
#     MISS_WHITELIST_BARCODE_FILTER,
#     LOW_MIN_QUAL_BARCODE_FILTER,
#     HAS_N_BARCODE_FILTER,
#     HOMOPOLYMER_BARCODE_FILTER,
# ]

UMI_POLYT_SUFFIX_LENGTH: int = 5


class Reporter:
    def __init__(
        self,
        umi_length: int | None = None,
        primers=None,
        reference_path: str | None = None,
        high_conf_mapq=None,
        chroms=None,
        barcode_whitelist: list[bytes] | None = None,
        barcode_summary=None,
        barcode_dist=None,
        gem_groups=None,
        metrics_dict: dict[
            str,
            tuple[
                type[Metric],
                dict[
                    str,
                    dict[str, str | bool] | list[str] | list[int | float] | list[range | list[int]],
                ],
            ],
        ] = METRICS,
        gene_index=None,
        umi_min_qual_threshold: int = 0,
        subsample_types=None,
        subsample_depths=None,
        genomes: list[str] | None = None,
        library_types=None,
        num_libraries: int | None = None,
    ):
        self.chroms = chroms
        self.umi_length = umi_length
        self.high_conf_mapq = high_conf_mapq
        self.barcode_whitelist = barcode_whitelist
        self.whitelist_set = set(barcode_whitelist) if barcode_whitelist else None
        self.barcode_summary = barcode_summary
        self.barcode_dist = barcode_dist
        self.gem_groups = gem_groups
        if num_libraries is None:
            self.library_indices = None  # type: Optional[List[int]]
        else:
            self.library_indices = list(range(num_libraries))

        self.gene_index = gene_index
        self.umi_min_qual_threshold = umi_min_qual_threshold
        random.seed(0)

        self.barcode_types = [cr_constants.ALL_BARCODES, cr_constants.FILTERED_BARCODES]
        self.molecule_types = [
            cr_constants.INSERT_MOLECULE_TYPE,
            cr_constants.FRAGMENT_MOLECULE_TYPE,
            cr_constants.CDNA_MOLECULE_TYPE,
            "cdna_candidate",
        ]
        self.read_types = [
            "all",
            "mapped",
            cr_constants.CONF_MAPPED_READ_TYPE,
            cr_constants.CONF_MAPPED_BC_READ_TYPE,
            cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
        ]
        self.dupe_types = [
            "cdna_pcr_uncorrected",
            cr_constants.CDNA_PCR_DUPE_TYPE,
            "si_pcr",
        ]
        self.umi_properties = [
            LOW_MIN_QUAL_UMI_FILTER,
            HAS_N_UMI_FILTER,
            HOMOPOLYMER_UMI_FILTER,
            "primer",
            "polyt_suffix",
        ]
        self.barcode_properties = [
            "miss_whitelist",
            LOW_MIN_QUAL_UMI_FILTER,
            HAS_N_UMI_FILTER,
            HOMOPOLYMER_UMI_FILTER,
        ]

        self.regions = [
            cr_constants.TRANSCRIPTOME_REGION,
            "genome",
            cr_constants.EXONIC_REGION,
            cr_constants.INTERGENIC_REGION,
            cr_constants.INTRONIC_REGION,
        ]

        assert not (reference_path is not None and genomes is not None)

        if reference_path is not None:
            self.genomes = get_reference_genomes(reference_path)
            self.references = self.genomes + [rna_library.MULTI_REFS_PREFIX]
            self.has_multiple_genomes = len(self.genomes) > 1
        elif genomes is not None:
            self.genomes = genomes
            self.references = self.genomes + [rna_library.MULTI_REFS_PREFIX]
            self.has_multiple_genomes = len(self.genomes) > 1
        else:
            self.genomes = []
            self.references = []
            self.has_multiple_genomes = False

        if library_types:
            self.library_types = [
                rna_library.get_library_type_metric_prefix(t) for t in library_types
            ]
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
                    "seq": seq,
                    "seq_rc": seq_rc,
                }
        else:
            self.primers = {}
            self.primer_names = []

        if self.umi_length is not None:
            self.umi_primers = {}
            self.poly_t = "T" * min(self.umi_length, UMI_POLYT_SUFFIX_LENGTH)
            for name, seq in self.primers.items():
                self.umi_primers[name] = {
                    "seq": seq["seq"][: self.umi_length],
                    "seq_rc": seq["seq_rc"][: self.umi_length],
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
        prefixes: list[list[str]] = []
        for base, (metric_cls, metric_dict) in self.metrics_dict.items():
            # Track prefixes for passing to the webshim
            prefixes.extend(metric_dict.get(_PREFIX_KEY, []))
            args = metric_dict.get(_ARGS_KEY, [])
            kwargs = metric_dict.get(_KWARGS_KEY, {})
            for key in self._get_metric_keys(base):
                self.metrics_data[key] = metric_cls(*args, **kwargs)

        self.prefixes = {}
        for prefix in prefixes:
            self.prefixes[prefix] = self.__dict__.get(prefix, [])

        self.reference_path = reference_path
        self.metadata = {}

    def _get_metric_keys(self, name: str) -> set[tuple[str]]:
        metric_cls, metric_dict = self.metrics_dict[name]
        prefixes: list[str] = metric_dict.get(_PREFIX_KEY, [])
        kwargs: dict[str, str | bool] = metric_dict.get(_KWARGS_KEY, {})

        always_active = kwargs.get(_ALWAYS_ACTIVE, False)

        parts = [[name]]
        for prefix in prefixes:
            prefix = getattr(self, prefix)

            if prefix:
                parts.append(prefix)

        # Check to make sure all specified metrics are present for metrics that are always active
        if always_active and len(parts) != len(prefixes) + 1:
            return set()

        # Return the set of keys
        keys = set(itertools.product(*parts))

        # Add bare keys
        keys.add((name,))

        return keys

    def _get_metric_attr(self, *args):
        return self.metrics_data[tuple(args)]

    def get_all_prefixes(self):
        return self.prefixes

    def store_pipeline_metadata(self, version):
        if self.metadata is None:
            self.metadata = {}
        self.metadata[cr_constants.CELLRANGER_VERSION_KEY] = version

    def store_reference_metadata(self, reference_path: str, ref_type: str, metric_prefix: str):
        """Stores reference metadata.

        Args:
            ref_type (str): e.g., 'Transcriptome'
            metric_prefix (str): e.g., 'vdj'
        """
        if self.metadata is None:
            self.metadata = {}

        with open(os.path.join(reference_path, cr_constants.REFERENCE_METADATA_FILE)) as f:
            ref_metadata = json.load(f)

        for key in [
            cr_constants.REFERENCE_FASTA_HASH_KEY,
            cr_constants.REFERENCE_GTF_HASH_KEY,
            cr_constants.REFERENCE_INPUT_FASTA_KEY,
            cr_constants.REFERENCE_INPUT_GTF_KEY,
            cr_constants.REFERENCE_MKREF_VERSION_KEY,
            cr_constants.REFERENCE_VERSION_KEY,
            cr_constants.REFERENCE_TYPE_KEY,  # New in CR 2.0.0
            cr_constants.REFERENCE_GENOMES_KEY,  # New in CR 2.0.0
        ]:
            value = ref_metadata.get(key, "")
            if value is None:
                value = ""

            # Backward compatibility with old reference metadata jsons that don't contain the type field
            if key == cr_constants.REFERENCE_TYPE_KEY and value == "":
                self.metadata[f"{metric_prefix}{cr_constants.REFERENCE_TYPE_KEY}"] = ref_type
                continue

            if np.isscalar(value):
                self.metadata[f"{metric_prefix}{key}"] = value
            elif key == cr_constants.REFERENCE_GENOMES_KEY:
                # Special case for genome key
                self.metadata[f"{metric_prefix}{key}"] = get_ref_name_from_genomes(value)
            else:
                self.metadata[f"{metric_prefix}{key}"] = ", ".join(str(x) for x in value)

    def store_chemistry_metadata(self, chemistry_def: dict):
        """Store the chemistry definition as metrics in the summary json."""
        if self.metadata is None:
            self.metadata = {}
        for key, value in chemistry_def.items():
            self.metadata["chemistry_%s" % key] = value

    def merge(self, reporter) -> None:
        for key, metric2 in reporter.metrics_data.items():
            metric1 = self.metrics_data.get(key)
            if isinstance(metric2, Metric):
                if metric1 and metric1.active and metric2.active:
                    metric1.merge(metric2)
                elif metric2.active:
                    self.metrics_data[key] = metric2

    def report(self, report_type):
        summary = {}
        for key, metric in self.metrics_data.items():
            if metric.active and metric.report_type == report_type:
                name = "_".join(
                    [str(x) for x in key[1:] + (key[0],)]
                )  # ensure keys are strings before output
                summary[name] = metric.report()

        if self.metadata is not None and report_type == DEFAULT_REPORT_TYPE:
            summary.update(self.metadata)

        return summary

    def to_json(self):
        """Return data as safe-for-json dictionary."""
        return tk_safe_json.json_sanitize(self.report(DEFAULT_REPORT_TYPE))

    def report_summary_json(self, filename):
        """Write to JSON file."""
        with open(filename, "w") as f:
            tk_safe_json.dump_numpy(self.to_json(), f, pretty=True)

    def report_barcodes_csv(self, filename):
        barcodes = sorted(self._get_metric_attr("barcodes_detected").d)
        with open(filename, "w") as f:
            for bc in barcodes:
                f.write(bc)
                f.write("\n")

    def save(self, filename):
        with open(filename, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(filename):
        with open(filename, "rb") as f:
            return pickle.load(f)


def merge_reporters(filenames):
    reporter = None
    for filename in filenames:
        tmp_reporter = Reporter.load(filename)
        if reporter is None:
            reporter = tmp_reporter
        else:
            reporter.merge(tmp_reporter)

    return reporter


def merge_jsons(in_filenames, out_filename, dicts: Iterable[dict] = []):
    """Merge a list of json files and optional dicts."""
    d = cr_utils.merge_jsons_as_dict(in_filenames)
    for data in dicts:
        d.update(data)

    with open(out_filename, "w") as f:
        tk_safe_json.dump_numpy(d, f, indent=4, sort_keys=True)
