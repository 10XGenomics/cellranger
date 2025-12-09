#![deny(missing_docs)]

use crate::SequencingMetrics;
use crate::barcode_sort::ReadVisitor;
use anyhow::Result;
use barcode::{Barcode, BarcodeConstruct, BarcodeConstructMetric, BcSegSeq};
use cr_types::chemistry::ChemistryDef;
use cr_types::reference::feature_extraction::FeatureExtractor;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::rna_read::{RnaChunk, RnaRead};
use fastq_set::metric_utils::{ILLUMINA_QUAL_OFFSET, PatternCheck};
use fastq_set::read_pair::{ReadPair, ReadPart, WhichRead};
use json_report_derive::JsonReport;
use metric::{
    CountMetric, Histogram, JsonReporter, Metric, OrderedHistogram, PercentMetric, SimpleHistogram,
    TxHashMap,
};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;
use std::convert::Into;
use std::sync::Arc;

const HOMOPOLYMER_LENGTH: usize = 15;
const BARCODE_MIN_QUAL_THRESHOLD: u8 = 10;
const UMI_MIN_QUAL_THRESHOLD: u8 = 10;
const UMI_POLYT_SUFFIX_LENGTH: usize = 5;

#[allow(non_snake_case)]
#[derive(Default, Serialize, Deserialize, Metric, JsonReport)]
#[json_report(extend = "extra_reports")]
pub struct MakeShardMetrics {
    bc_N_bases: PercentMetric,
    umi_N_bases: PercentMetric,
    read_N_bases: PercentMetric,
    read2_N_bases: Option<PercentMetric>,
    i1_N_bases: Option<PercentMetric>,
    i2_N_bases: Option<PercentMetric>,

    bc_bases_with_q30: PercentMetric,
    /// bc_bases_with_q30_in_gel_bead_frac and bc_bases_with_q30_in_probe_frac
    bc_bases_with_q30_in: BarcodeConstructMetric<PercentMetric>,
    umi_bases_with_q30: PercentMetric,
    read_bases_with_q30: PercentMetric,
    read2_bases_with_q30: Option<PercentMetric>,
    i1_bases_with_q30: Option<PercentMetric>,
    i2_bases_with_q30: Option<PercentMetric>,

    A_perfect_homopolymer: PercentMetric,
    C_perfect_homopolymer: PercentMetric,
    G_perfect_homopolymer: PercentMetric,
    T_perfect_homopolymer: PercentMetric,

    good_umi: PercentMetric,
    has_n_barcode_property: PercentMetric,
    has_n_umi_property: PercentMetric,
    pub(super) homopolymer_barcode_property: PercentMetric,
    homopolymer_umi_property: PercentMetric,
    low_min_qual_barcode_property: PercentMetric,
    low_min_qual_umi_property: PercentMetric,
    polyt_suffix_umi_property: PercentMetric,
    miss_whitelist_barcode_property: PercentMetric,

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
    sequenced_reads: CountMetric,

    /// Number of read pairs skipped because it cannot be a valid RnaRead, typically due to
    /// insufficient length. Equal to r2 - r3 (see comment above)
    unprocessed_read_pairs: i64,

    /// See the comment on `sequenced_reads`. Only reported for VDJ library types
    #[json_report(skip)]
    vdj_total_read_pairs: i64,

    /// Total read pairs per library
    #[json_report(skip)]
    total_read_pairs_per_library: TxHashMap<u16, i64>,

    /// Fraction of outstanding unknown feature barcode sequences
    #[json_report(block)]
    unknown_feature_bcs: Option<TxHashMap<String, f64>>,
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

#[derive(Clone, Serialize, Deserialize)]
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
        if let Some(r1) = rna_read.read.get(WhichRead::R1, ReadPart::Seq) {
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

pub struct MakeShardVisitor<'a> {
    poly_a_pattern: PatternCheck,
    poly_c_pattern: PatternCheck,
    poly_g_pattern: PatternCheck,
    poly_t_pattern: PatternCheck,
    metrics: MakeShardMetrics,
    histograms: MakeShardHistograms,
    valid_bc_counts: OrderedHistogram<Barcode>,
    read_chunks: &'a [RnaChunk],
    pub(super) feature_counts: Vec<i64>,
    feature_extractor: Option<FeatureExtractor>,
}

impl<'a> MakeShardVisitor<'a> {
    pub fn new(
        reference: Option<FeatureReference>,
        read_chunks: &'a [RnaChunk],
        chunk_id: usize,
        chemistry_def: &ChemistryDef,
    ) -> Result<Self> {
        let (feature_counts, feature_extractor) = if let Some(reference) = reference {
            (
                vec![0i64; reference.feature_defs.len()],
                Some(FeatureExtractor::new(Arc::new(reference), None)?),
            )
        } else {
            (vec![], None)
        };
        // Handle edge case where a library has no read counts (possibly because it wasn't sequenced
        // far enough to have a parseable read) and so would not be initialized with a count.
        let mut metrics = MakeShardMetrics::default();
        metrics
            .total_read_pairs_per_library
            .insert(read_chunks[chunk_id].library_id(), 0);
        Ok(MakeShardVisitor {
            poly_a_pattern: PatternCheck::new(&[b'A'; HOMOPOLYMER_LENGTH]),
            poly_c_pattern: PatternCheck::new(&[b'C'; HOMOPOLYMER_LENGTH]),
            poly_g_pattern: PatternCheck::new(&[b'G'; HOMOPOLYMER_LENGTH]),
            poly_t_pattern: PatternCheck::new(&[b'T'; HOMOPOLYMER_LENGTH]),
            metrics,
            feature_counts,
            feature_extractor,
            read_chunks,
            histograms: MakeShardHistograms::new(chemistry_def.barcode_construct().map(|_| ())),
            valid_bc_counts: Default::default(),
        })
    }

    /// Return the associated metrics.
    pub fn metrics(&self) -> &MakeShardMetrics {
        &self.metrics
    }

    /// Return the various histograms.
    pub fn histograms(&self) -> &MakeShardHistograms {
        &self.histograms
    }

    /// Return the ordered valid barcode histogram.
    pub fn valid_bc_counts(&self) -> &OrderedHistogram<Barcode> {
        &self.valid_bc_counts
    }
}

impl ReadVisitor for MakeShardVisitor<'_> {
    type ReadType = RnaRead;

    fn visit_processed_read(&mut self, rna_read: &mut Self::ReadType) -> Result<()> {
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

            good_umi: rna_read.umi().is_valid().into(),
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
            total_read_pairs_per_library: std::iter::once((
                rna_read.read_chunk(self.read_chunks).library_id(),
                1,
            ))
            .collect(),
            unknown_feature_bcs: Default::default(),

            unprocessed_read_pairs: 0,
        };
        self.metrics.merge(read_metrics);
        self.histograms.observe(rna_read);
        if rna_read.barcode_is_valid() {
            self.valid_bc_counts.observe_owned(rna_read.barcode());
        }

        // collect initial feature barcode counts
        if let Some(res) = self
            .feature_extractor
            .as_ref()
            .and_then(|e| e.match_read(rna_read))
        {
            if res.ids.len() == 1 {
                self.feature_counts[res.ids[0].0] += 1;
            } else {
                let raw_feature_bc = String::from_utf8(res.barcode).unwrap();
                self.histograms
                    .unknown_feature_bcs
                    .observe_owned(raw_feature_bc);
            }
        }

        Ok(())
    }

    fn visit_unprocessed_read(&mut self, _: ReadPair, _: String) -> Result<()> {
        self.metrics.unprocessed_read_pairs += 1;
        Ok(())
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
