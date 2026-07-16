#![deny(missing_docs)]

use crate::SequencingMetrics;
use crate::barcode_sort::ReadTrimmingMetrics;
use anyhow::Result;
use barcode::{BarcodeConstruct, BarcodeConstructMetric, BcSegSeq};
use cr_types::SequencingPlatform;
use cr_types::chemistry::ChemistryDef;
use cr_types::reference::feature_extraction::FeatureExtractor;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::rna_read::RnaRead;
use fastq_set::metric_utils::{ILLUMINA_QUAL_OFFSET, PatternCheck};
use fastq_set::read_pair::WhichRead;
use itertools::Itertools;
use json_report_derive::JsonReport;
use metric::{
    CountMetric, Histogram, JsonReporter, Metric, PercentMetric, SimpleHistogram, TxHashMap,
};
use metric_derive::Metric;
use regex::bytes::Regex;
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;
use std::convert::Into;
use std::sync::{Arc, LazyLock};

const HOMOPOLYMER_LENGTH: usize = 15;
const BARCODE_MIN_QUAL_THRESHOLD: u8 = 10;
const UMI_MIN_QUAL_THRESHOLD: u8 = 10;
const UMI_POLYT_SUFFIX_LENGTH: usize = 5;

#[allow(non_snake_case)]
#[derive(Debug, Clone, Default, Serialize, Deserialize, Metric, JsonReport)]
#[json_report(extend = "extra_reports")]
pub(super) struct MakeShardMetrics {
    pub bc_N_bases: PercentMetric,
    pub umi_N_bases: PercentMetric,
    pub read_N_bases: PercentMetric,
    pub read2_N_bases: Option<PercentMetric>,
    pub i1_N_bases: Option<PercentMetric>,
    pub i2_N_bases: Option<PercentMetric>,

    pub bc_bases_with_q30: PercentMetric,
    /// bc_bases_with_q30_in_gel_bead_frac and bc_bases_with_q30_in_probe_frac
    pub bc_bases_with_q30_in: BarcodeConstructMetric<PercentMetric>,
    pub umi_bases_with_q30: PercentMetric,
    pub read_bases_with_q30: PercentMetric,
    pub read2_bases_with_q30: Option<PercentMetric>,
    pub i1_bases_with_q30: Option<PercentMetric>,
    pub i2_bases_with_q30: Option<PercentMetric>,

    pub A_perfect_homopolymer: PercentMetric,
    pub C_perfect_homopolymer: PercentMetric,
    pub G_perfect_homopolymer: PercentMetric,
    pub T_perfect_homopolymer: PercentMetric,

    pub good_umi: PercentMetric,
    pub has_n_barcode_property: PercentMetric,
    pub has_n_umi_property: PercentMetric,
    pub homopolymer_barcode_property: PercentMetric,
    pub homopolymer_umi_property: PercentMetric,
    pub low_min_qual_barcode_property: PercentMetric,
    pub low_min_qual_umi_property: PercentMetric,
    pub polyt_suffix_umi_property: PercentMetric,
    pub miss_whitelist_barcode_property: PercentMetric,

    /// If there are a total of `r1` read pairs in the input fastqs, we process `r2` read pairs
    /// through `MAKE_SHARD`. `r2` is the minimum of these numbers:
    /// - `r1` read pairs subsampled at a rate specified in the top level mro or `SampleDef` (PD)
    /// - `initial_reads` specified at the top level mro
    ///
    /// Among the `r2` read pairs, we end up carrying `r3` reads forward from `MAKE_SHARD`, where
    /// `r3` is the number of reads among `r2` which are valid `RnaRead`s. The reads which are
    /// ignored are the ones which are too short. Typically, for customer data, `r1 = r2 = r3`
    ///
    /// The metric name is probably a misnomer, but by `sequenced_reads`, here we mean `r3`. This
    /// appears under the key `sequenced_reads_count` with the appropriate library type prefix
    /// in the metrics summary json
    ///
    /// There is a related metric, `total_read_pairs` in the metrics summary json which we have
    /// used historically for `r3`. `total_read_pairs` metric was originally computed in this stage,
    /// but has been moved to `ALIGN_AND_COUNT` except for VDJ which is why we have a
    /// `vdj_total_read_pairs` metric in this struct. `total_read_pairs` = `sequenced_reads_count`.
    pub sequenced_reads: CountMetric,

    /// Number of read pairs skipped because it cannot be a valid RnaRead, typically due to
    /// insufficient length. Equal to r2 - r3 (see comment above)
    pub unprocessed_read_pairs: i64,

    /// See the comment on `sequenced_reads`. Only reported for VDJ library types
    #[json_report(skip)]
    pub vdj_total_read_pairs: i64,

    /// Total read pairs per library
    #[json_report(skip)]
    pub total_read_pairs_per_library: TxHashMap<u16, i64>,

    /// Fraction of outstanding unknown feature barcode sequences
    #[json_report(block)]
    pub unknown_feature_bcs: Option<TxHashMap<String, f64>>,

    /// Metrics collected from trimming long reads.
    #[json_report(inline)]
    pub(crate) read_trimming_metrics: Option<ReadTrimmingMetrics>,

    /// Observed and inferred properties of the presumed sequencing platform.
    #[json_report(inline)]
    pub(super) platform: Option<SequencerPlatformMetrics>,
}

impl MakeShardMetrics {
    /// Emit total_read_pairs_per_library without the default _count suffix.
    fn extra_reports(&self) -> JsonReporter {
        let mut metrics: JsonReporter = self
            .total_read_pairs_per_library
            .iter()
            .map(|(library_id, count)| {
                (format!("{library_id}_total_read_pairs_per_library"), count)
            })
            .collect();

        if self.vdj_total_read_pairs > 0 {
            metrics.insert("total_read_pairs", self.vdj_total_read_pairs);
        }
        metrics
    }

    pub fn sequencing_metrics_for_fastq(&self, fastq_id: String) -> SequencingMetrics {
        // Only populate muli-part barcode metrics if both GEM and probe barcode are present
        let (q30_gem, q30_probe) = match self.bc_bases_with_q30_in.inner() {
            Some(BarcodeConstruct::GelBeadAndProbe(c)) => (Some(c.gel_bead), Some(c.probe)),
            _ => (None, None),
        };

        SequencingMetrics {
            fastq_id,
            number_of_reads: self.sequenced_reads.count() as usize,
            q30_barcode: self.bc_bases_with_q30,
            q30_gem_barcode: q30_gem,
            q30_probe_barcode: q30_probe,
            q30_umi: self.umi_bases_with_q30,
            q30_read1: self.read_bases_with_q30,
            unprocessed_reads: self.unprocessed_read_pairs as usize,
            q30_read2: self.read2_bases_with_q30,
            sequencing_platform: None,
            read_trimming_metrics: self.read_trimming_metrics.clone(),
        }
    }
    pub fn report_unknown_fbc(
        &mut self,
        library_id: u16,
        hist: &SimpleHistogram<String>,
        min_frac: f64,
    ) {
        // Report unknown feature barcode sequences that represent a non-trivial fraction of reads per library.
        if let Some(total_read_pairs) = self.total_read_pairs_per_library.get(&library_id) {
            for (seq, read_count) in hist {
                let frac_reads: f64 = read_count.count() as f64 / *total_read_pairs as f64;
                if frac_reads >= min_frac {
                    self.unknown_feature_bcs
                        .get_or_insert_default()
                        .insert(seq.clone(), frac_reads);
                }
            }
        } else {
            eprintln!("No total read pairs data available for library: {library_id}");
        }
    }
}
#[derive(Clone, Debug, Default, Serialize, Deserialize, JsonReport)]
pub(super) struct SequencerPlatformMetrics {
    /// The number of :-separated components in the fastq ID header.
    fastq_id_components: usize,
    /// The sequencing platform, if we were able to uniquely infer it from headers.
    inferred_platform: Option<SequencingPlatform>,
}

impl SequencerPlatformMetrics {
    pub fn from_fastq_id(id: &[u8]) -> Self {
        let fastq_id_components = id.iter().filter(|&&c| c == b':').count() + 1;

        static ID_PATTERNS: LazyLock<Vec<(SequencingPlatform, Regex)>> = LazyLock::new(|| {
            vec![
                (
                    SequencingPlatform::Illumina,
                    // <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
                    Regex::new(r"^[a-zA-Z0-9_]*:\d*:[a-zA-Z0-9]*:\d*:\d*:\d*:\d*(?: |$)").unwrap(),
                ),
                (
                    SequencingPlatform::Ultima,
                    // <instrument id>:<run id>:<wafer id>:<chuck>:<camera>:<ring>:<tile>:<x-pos>:<y-pos>:<segment-num>:<first-flow-signal>:<filtered>:RSQ
                    Regex::new(r"^[a-zA-Z0-9]*:\d*:[a-zA-Z0-9]*:\d*:\d*:\d*:\d*:\d*:\d*:\d*:\d*:[YN]:[0-9\.]*:[ACTG]*:[A-Z]*:\d*(?: |$)")
                        .unwrap(),
                ),
            ]
        });

        let inferred_platform = ID_PATTERNS
            .iter()
            .filter_map(|(platform, re)| re.is_match(id).then_some(*platform))
            .unique()
            .at_most_one()
            .ok() // if we matched more than one platform, return None
            .flatten();
        Self {
            fastq_id_components,
            inferred_platform,
        }
    }
}

impl Metric for SequencerPlatformMetrics {
    fn merge(&mut self, other: Self) {
        // Keep whichever FASTQ ID component count is higher.
        // TODO: this is arbitrary. It would be nice if we could report a value
        // per fastq, but telemetry does not support this currently.
        self.fastq_id_components = self.fastq_id_components.max(other.fastq_id_components);
        // We'll probably only ever have one of these...
        self.inferred_platform = self.inferred_platform.max(other.inferred_platform);
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MakeShardHistograms {
    pub(super) valid_bc_segment_counts: BarcodeConstruct<SimpleHistogram<BcSegSeq>>,
    pub(super) r1_lengths: SimpleHistogram<usize>,
    pub(super) unknown_feature_bcs: SimpleHistogram<String>,
}

impl MakeShardHistograms {
    pub fn new(barcode_construct: BarcodeConstruct<()>) -> Self {
        MakeShardHistograms {
            r1_lengths: SimpleHistogram::default(),
            valid_bc_segment_counts: barcode_construct.map(|()| SimpleHistogram::default()),
            unknown_feature_bcs: SimpleHistogram::default(),
        }
    }

    pub fn merge(&mut self, other: Self) {
        // NOTE: The functions new() and merge() look eerily similar to the
        // `Metric` trait. However` valid_bc_segment_counts` require a context
        // to define a `new()`. We expand the fields in the struct here explicitly
        // so that the compiler will remind you to update this function when
        // adding/removing fields from the struct.
        //
        // TODO: `Metric` can be defined for `Option<BarcodeConstruct<Metric>>`
        // in which case we can derive the trait
        let MakeShardHistograms {
            valid_bc_segment_counts,
            r1_lengths,
            unknown_feature_bcs,
        } = other;
        self.r1_lengths.merge(r1_lengths);
        self.unknown_feature_bcs.merge(unknown_feature_bcs);
        self.valid_bc_segment_counts
            .as_mut_ref()
            .zip(valid_bc_segment_counts)
            .iter()
            .for_each(|(this, that)| this.merge(that));
    }
    fn observe(&mut self, rna_read: &RnaRead) {
        if let Some(r1) = rna_read.seq(WhichRead::R1) {
            self.r1_lengths.observe(&r1.len());
        }
        rna_read
            .segmented_barcode
            .segments()
            .zip(self.valid_bc_segment_counts.as_mut_ref())
            .iter()
            .for_each(|(segment, hist)| {
                if segment.is_valid() {
                    hist.observe_owned(*segment.sequence());
                }
            });
    }
}

#[derive(Debug, Clone)]
pub struct MakeShardMetricsAccum {
    pub metrics: MakeShardMetrics,
    pub histograms: MakeShardHistograms,
    pub feature_counts: Vec<i64>,
    library_id: u16,
}

impl MakeShardMetricsAccum {
    pub fn new(feature_count: usize, library_id: u16, chemistry_def: &ChemistryDef) -> Self {
        let mut metrics = MakeShardMetrics::default();
        // Handle edge case where a library has no read counts (possibly because it wasn't sequenced
        // far enough to have a parseable read) and so would not be initialized with a count.
        metrics.total_read_pairs_per_library.insert(library_id, 0);
        Self {
            metrics,
            histograms: MakeShardHistograms::new(chemistry_def.barcode_construct().map(|_| ())),
            feature_counts: vec![0i64; feature_count],
            library_id,
        }
    }

    pub fn merge(&mut self, rhs: Self) {
        self.metrics.merge(rhs.metrics);
        self.histograms.merge(rhs.histograms);
        self.feature_counts
            .iter_mut()
            .zip_eq(rhs.feature_counts)
            .for_each(|(x, y)| *x += y);
        assert_eq!(self.library_id, rhs.library_id);
    }
}

pub struct MakeShardMetricProc {
    poly_a_pattern: PatternCheck,
    poly_c_pattern: PatternCheck,
    poly_g_pattern: PatternCheck,
    poly_t_pattern: PatternCheck,
    feature_extractor: Option<Arc<FeatureExtractor>>,
    pub feature_count: usize,
}

impl MakeShardMetricProc {
    pub fn new(reference: Option<FeatureReference>) -> Result<Self> {
        let (feature_count, feature_extractor) = if let Some(reference) = reference {
            (
                reference.feature_defs.len(),
                Some(Arc::new(FeatureExtractor::new(Arc::new(reference), None)?)),
            )
        } else {
            (0, None)
        };

        Ok(MakeShardMetricProc {
            poly_a_pattern: PatternCheck::new(&[b'A'; HOMOPOLYMER_LENGTH]),
            poly_c_pattern: PatternCheck::new(&[b'C'; HOMOPOLYMER_LENGTH]),
            poly_g_pattern: PatternCheck::new(&[b'G'; HOMOPOLYMER_LENGTH]),
            poly_t_pattern: PatternCheck::new(&[b'T'; HOMOPOLYMER_LENGTH]),
            feature_extractor,
            feature_count,
        })
    }

    pub fn visit_processed_read(&self, rna_read: &RnaRead, accum: &mut MakeShardMetricsAccum) {
        let homopolymer_metric = |pattern_check: &PatternCheck,
                                  read: &[u8],
                                  read2_option: Option<&[u8]>|
         -> PercentMetric {
            let mut exists = pattern_check.exists(read);
            if let Some(read2) = read2_option {
                exists = exists || pattern_check.exists(read2);
            }
            exists.into()
        };

        let read_metrics = MakeShardMetrics {
            bc_N_bases: frac_n_bases(rna_read.raw_bc_construct_seq().flat_iter()),
            umi_N_bases: frac_n_bases(rna_read.raw_umi_seq()),
            read_N_bases: frac_n_bases(rna_read.r1_seq()),
            read2_N_bases: rna_read.r2_seq().map(frac_n_bases),
            i1_N_bases: rna_read.raw_illumina_i1_seq().map(frac_n_bases),
            i2_N_bases: rna_read.raw_illumina_i2_seq().map(frac_n_bases),

            bc_bases_with_q30: frac_q30_bases(rna_read.raw_bc_construct_qual().flat_iter()),
            bc_bases_with_q30_in: rna_read.raw_bc_construct_qual().map(frac_q30_bases).into(),
            umi_bases_with_q30: frac_q30_bases(rna_read.raw_umi_qual()),
            read_bases_with_q30: frac_q30_bases(rna_read.r1_qual()),
            read2_bases_with_q30: rna_read.r2_qual().map(frac_q30_bases),
            i1_bases_with_q30: rna_read.raw_illumina_i1_qual().map(frac_q30_bases),
            i2_bases_with_q30: rna_read.raw_illumina_i2_qual().map(frac_q30_bases),

            A_perfect_homopolymer: homopolymer_metric(
                &self.poly_a_pattern,
                rna_read.r1_seq(),
                rna_read.r2_seq(),
            ),
            C_perfect_homopolymer: homopolymer_metric(
                &self.poly_c_pattern,
                rna_read.r1_seq(),
                rna_read.r2_seq(),
            ),
            G_perfect_homopolymer: homopolymer_metric(
                &self.poly_g_pattern,
                rna_read.r1_seq(),
                rna_read.r2_seq(),
            ),
            T_perfect_homopolymer: homopolymer_metric(
                &self.poly_t_pattern,
                rna_read.r1_seq(),
                rna_read.r2_seq(),
            ),

            good_umi: rna_read.umi().is_valid(&rna_read.raw_umi_qual()).into(),
            has_n_barcode_property: rna_read
                .raw_bc_construct_seq()
                .iter()
                .any(|s| s.has_n())
                .into(),
            has_n_umi_property: rna_read.umi().sseq().has_n().into(),
            homopolymer_barcode_property: rna_read
                .raw_bc_construct_seq()
                .iter()
                .any(|s| s.is_homopolymer())
                .into(),
            homopolymer_umi_property: rna_read.umi().sseq().is_homopolymer().into(),
            low_min_qual_barcode_property: (rna_read.barcode_min_qual()
                < BARCODE_MIN_QUAL_THRESHOLD)
                .into(),
            low_min_qual_umi_property: (rna_read.umi_min_qual() < UMI_MIN_QUAL_THRESHOLD).into(),
            polyt_suffix_umi_property: rna_read
                .umi()
                .sseq()
                .has_polyt_suffix(UMI_POLYT_SUFFIX_LENGTH)
                .into(),
            miss_whitelist_barcode_property: (!rna_read.barcode_is_valid()).into(),
            sequenced_reads: 1.into(),
            vdj_total_read_pairs: i64::from(rna_read.library_type.is_vdj()),
            // Optimization - manually merge this below to avoid allocating a hashmap for every read.
            total_read_pairs_per_library: TxHashMap::default(),
            read_trimming_metrics: None,
            unknown_feature_bcs: None,

            unprocessed_read_pairs: 0,
            platform: None,
        };
        accum.metrics.merge(read_metrics);
        *accum
            .metrics
            .total_read_pairs_per_library
            .get_mut(&accum.library_id)
            .unwrap() += 1;
        accum.histograms.observe(rna_read);

        // collect initial feature barcode counts
        if let Some(res) = self
            .feature_extractor
            .as_ref()
            .and_then(|e| e.match_read(rna_read))
        {
            if res.ids.len() == 1 {
                accum.feature_counts[res.ids[0].0] += 1;
            } else {
                let raw_feature_bc = String::from_utf8(res.barcode).unwrap();
                accum
                    .histograms
                    .unknown_feature_bcs
                    .observe_owned(raw_feature_bc);
            }
        }
    }
}

/// Compute the fraction of 'N' bases in the input sequence
/// as a PercentMetric
#[inline(always)]
fn frac_n_bases<C, D>(seq: D) -> PercentMetric
where
    C: Borrow<u8>,
    D: IntoIterator<Item = C>,
{
    let mut result = PercentMetric::default();
    for base in seq {
        result.increment(*base.borrow() == b'N');
    }
    result
}

/// Compute the fraction of Q30 bases in the input sequence
/// quality as a PercentMetric
/// Only bases with quality > 2 is counted in the denominator
fn frac_q30_bases<C, D>(qual: D) -> PercentMetric
where
    C: Borrow<u8>,
    D: IntoIterator<Item = C>,
{
    let mut num = 0;
    let mut den = 0;
    for q in qual {
        let q = *q.borrow();
        if q > 2 + ILLUMINA_QUAL_OFFSET {
            den += 1;
            if q >= 30 + ILLUMINA_QUAL_OFFSET {
                num += 1;
            }
        }
    }
    (num, den).into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barcode_sort::ReadTrimmingMetrics;
    #[test]
    fn test_sequencing_metrics_includes_trimming() {
        let m = MakeShardMetrics {
            read_trimming_metrics: Some(ReadTrimmingMetrics::default()),
            ..Default::default()
        };
        let seq = m.sequencing_metrics_for_fastq("test_fastq".to_string());
        assert!(seq.read_trimming_metrics.is_some());
    }

    fn check_fastq_id(
        id: &[u8],
        expected_components: usize,
        expected_platform: Option<SequencingPlatform>,
    ) {
        let metrics = SequencerPlatformMetrics::from_fastq_id(id);
        assert_eq!(
            metrics.fastq_id_components,
            expected_components,
            "Header {:?}: expected {} ID components, got {}",
            String::from_utf8_lossy(id),
            expected_components,
            metrics.fastq_id_components
        );
        assert_eq!(
            metrics.inferred_platform,
            expected_platform,
            "Header {:?}: expected platform {:?}, got {:?}",
            String::from_utf8_lossy(id),
            expected_platform,
            metrics.inferred_platform
        );
    }

    #[test]
    fn test_from_fastq_header() {
        // Illumina headers
        check_fastq_id(
            b"D000684:840:HFCNFBCXY:1:1101:4266:2693",
            7,
            Some(SequencingPlatform::Illumina),
        );

        // Production Ultima headers
        check_fastq_id(
            b"V143:442353:442353:442353:0:2:1:187:11:1:990:N:0.872:CTCAGCATA:NA:1286",
            16,
            Some(SequencingPlatform::Ultima),
        );

        // Unrecognized formats
        check_fastq_id(b"UNKNOWN_FORMAT_READ_12345", 1, None);
        check_fastq_id(b"", 1, None);
        check_fastq_id(b"INST:1:FLOWCELL:1:2", 5, None);
        check_fastq_id(b"read.123.something/1", 1, None);
        check_fastq_id(b"A:B:C", 3, None);
    }
}
