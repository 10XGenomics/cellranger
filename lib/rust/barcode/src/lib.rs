//! Crate for dealing with 10x barcodes.
//!
//! Contains tools for loading barcode whitelists, check barcodes against the whitelist
//! and correcting sequencing errors in barcodes.

use std::fmt;
use std::iter::{self, FromIterator};
use std::ops::Deref;
use std::str::FromStr;

pub mod binned;
pub mod corrector;
mod io_utils;
pub mod whitelist;

pub mod short_string;
use anyhow::{anyhow, bail, Result};
use arrayvec::ArrayVec;
pub use corrector::BarcodeCorrector;
use fastq_set::squality::SQualityGen;
use fastq_set::sseq::SSeqGen;
use itertools::Itertools;
use metric::{JsonReport, JsonReporter, Metric};
use serde::{Deserialize, Serialize};
pub use short_string::*;
pub use whitelist::{Whitelist, WhitelistSource, WhitelistSpec};
use BarcodeConstruct::{GelBeadAndProbe, GelBeadOnly, Segmented};

/* ---------------------------------------------------------------------------------------------- */

/// The maximum supported barcode length
pub const MAX_BARCODE_LENGTH: usize = 43;

/// The maximum length of each segment of the barcode we support. The
/// barcode segment on the gel bead determines the maximum length of 16
pub const MAX_BARCODE_SEGMENT_LENGTH: usize = 22;
pub const BARCODE_GG_SUFFIX_LENGTH: usize = 2; // Support "-1" to "-9" suffix.

/* ---------------------------------------------------------------------------------------------- */
// NOTE: DO NOT EDIT THESE DIRECTLY
// These types/constants are derived based on the types/constants defined above.

pub const MAX_BARCODE_AND_GEM_GROUP_LENGTH: usize = MAX_BARCODE_LENGTH + BARCODE_GG_SUFFIX_LENGTH;
/// Type for storing barcode sequence
pub type BcSeq = SSeqGen<MAX_BARCODE_LENGTH>;
/// Type for storing sequences of barcode segments
pub type BcSegSeq = SSeqGen<MAX_BARCODE_SEGMENT_LENGTH>;
/// Type for storing barcode quality
pub type BcQual = SQualityGen<MAX_BARCODE_LENGTH>;
/// Type for storing qualities of barcode segments
pub type BcSegQual = SQualityGen<MAX_BARCODE_SEGMENT_LENGTH>;

/* ---------------------------------------------------------------------------------------------- */
/// Represent a (possibly-corrected) 10x barcode sequence, and it's GEM group. If the barcode
/// contains multiple segments, the sequence will be the concatenation of sequences of individual
/// segments.
///
/// Note the natural sort order groups barcodes by GEM group, then whether they valid, then
/// by the barcode sequence.
#[derive(Serialize, Deserialize, Copy, Clone, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct Barcode {
    gem_group: u16,
    valid: bool,
    sequence: BcSeq,
}

impl Barcode {
    pub fn is_valid(self) -> bool {
        self.valid
    }
    pub fn from_bytes(seq_gg: &[u8]) -> Result<Self> {
        std::str::from_utf8(seq_gg)?.parse()
    }
    pub fn with_seq(gem_group: u16, sequence: BcSeq, valid: bool) -> Self {
        Barcode {
            gem_group,
            sequence,
            valid,
        }
    }
    pub fn new(gem_group: u16, sequence: &[u8], valid: bool) -> Self {
        Barcode {
            gem_group,
            sequence: BcSeq::from_bytes(sequence),
            valid,
        }
    }
    pub fn gem_group(self) -> u16 {
        self.gem_group
    }
    pub fn bcseq(&self) -> &BcSeq {
        &self.sequence
    }
    pub fn seq(&self) -> &[u8] {
        self.sequence.as_bytes()
    }

    /// ASCII string of corrected, GEM group appended form of
    /// barcode, suitable for use in BAM files (CB or BX tags)
    /// For example: "AGCCGATA-1"
    pub fn to_corrected_bytes(self) -> Vec<u8> {
        self.to_string().into_bytes()
    }
}

impl fmt::Display for Barcode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}", self.sequence, self.gem_group)
    }
}

impl FromStr for Barcode {
    type Err = anyhow::Error;

    /// Parse a barcode from its string representation "AACCGGTT-1".
    fn from_str(seq_gg: &str) -> Result<Barcode> {
        let Some((bc, gg_str)) = seq_gg.split_once('-') else {
            bail!("invalid barcode: '{}'", seq_gg)
        };
        let gg: u16 = gg_str.parse()?;
        Ok(Barcode::new(gg, bc.as_bytes(), true))
    }
}

impl From<&str> for Barcode {
    /// Parse a barcode from its string representation "AACCGGTT-1", and panic if it fails.
    fn from(seq_gg: &str) -> Self {
        seq_gg.parse().unwrap()
    }
}

/// Deserialize a Barcode from its string representation, "AACCGGTT-1".
#[derive(Deserialize, Serialize)]
#[serde(from = "&str")]
pub struct BarcodeFromString(pub Barcode);

impl Deref for BarcodeFromString {
    type Target = Barcode;

    fn deref(&self) -> &Barcode {
        &self.0
    }
}

impl From<&str> for BarcodeFromString {
    /// Convert a String to a BarcodeFromString.
    fn from(barcode_gg: &str) -> BarcodeFromString {
        BarcodeFromString(Barcode::from(barcode_gg))
    }
}

impl From<BarcodeFromString> for Barcode {
    /// Convert a BarcodeFromString to a Barcode.
    fn from(barcode: BarcodeFromString) -> Barcode {
        barcode.0
    }
}

/// Serialize and deserialize a Barcode, for use with `#[serde(with = "barcode_string")]`.
pub mod barcode_string {
    use crate::Barcode;
    use serde::de::Error;
    use serde::{Deserialize, Deserializer, Serializer};

    /// Deserialize a Barcode from its string representation, "AACCGGTT-1".
    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Barcode, D::Error> {
        <&str>::deserialize(deserializer)?
            .parse()
            .map_err(D::Error::custom)
    }

    /// Serialize a Barcode to its string representation, "AACCGGTT-1".
    pub fn serialize<S: Serializer>(barcode: &Barcode, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.collect_str(barcode)
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// Different states a barcode segment could be in.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub enum BarcodeSegmentState {
    /// The sequence is not checked against any whitelist
    NotChecked,
    /// The sequence exists in the whitelist
    ValidBeforeCorrection,
    /// The sequence exists in the whitelist after correction
    ValidAfterCorrection,
    /// The sequence does not exist in the whitelist
    Invalid,
}

impl BarcodeSegmentState {
    pub fn is_valid(self) -> bool {
        match self {
            BarcodeSegmentState::ValidBeforeCorrection
            | BarcodeSegmentState::ValidAfterCorrection => true,
            BarcodeSegmentState::NotChecked | BarcodeSegmentState::Invalid => false,
        }
    }
    /// Change the state given the evidence on whether the sequence is in the whitelist
    pub fn change(&mut self, is_in_whitelist: bool) {
        *self = match self {
            BarcodeSegmentState::NotChecked => {
                if is_in_whitelist {
                    BarcodeSegmentState::ValidBeforeCorrection
                } else {
                    BarcodeSegmentState::Invalid
                }
            }
            BarcodeSegmentState::Invalid => {
                if is_in_whitelist {
                    BarcodeSegmentState::ValidAfterCorrection
                } else {
                    BarcodeSegmentState::Invalid
                }
            }
            s => panic!("Cannot transition from state {s:?}"),
        }
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// A container to store two objects associated with a gel bead barcode and a probe barcode.
#[derive(Deserialize, Serialize, Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct GelBeadAndProbeConstruct<G> {
    pub gel_bead: G,
    pub probe: G,
}

impl<G> GelBeadAndProbeConstruct<G> {
    fn as_ref(&self) -> GelBeadAndProbeConstruct<&G> {
        GelBeadAndProbeConstruct {
            gel_bead: &self.gel_bead,
            probe: &self.probe,
        }
    }

    fn as_mut_ref(&mut self) -> GelBeadAndProbeConstruct<&mut G> {
        GelBeadAndProbeConstruct {
            gel_bead: &mut self.gel_bead,
            probe: &mut self.probe,
        }
    }

    pub fn map<F, K>(self, f: F) -> GelBeadAndProbeConstruct<K>
    where
        F: FnOnce(G) -> K + Copy,
    {
        GelBeadAndProbeConstruct {
            gel_bead: f(self.gel_bead),
            probe: f(self.probe),
        }
    }

    pub fn map_result<F, K>(self, f: F) -> Result<GelBeadAndProbeConstruct<K>>
    where
        F: FnOnce(G) -> Result<K> + Copy,
    {
        Ok(GelBeadAndProbeConstruct {
            gel_bead: f(self.gel_bead)?,
            probe: f(self.probe)?,
        })
    }

    pub fn zip<K>(self, other: GelBeadAndProbeConstruct<K>) -> GelBeadAndProbeConstruct<(G, K)> {
        GelBeadAndProbeConstruct {
            gel_bead: (self.gel_bead, other.gel_bead),
            probe: (self.probe, other.probe),
        }
    }
}

impl<G> IntoIterator for GelBeadAndProbeConstruct<G> {
    type Item = G;
    type IntoIter = iter::Chain<iter::Once<G>, iter::Once<G>>;

    fn into_iter(self) -> Self::IntoIter {
        iter::once(self.gel_bead).chain(iter::once(self.probe))
    }
}

impl<G: Metric> Metric for GelBeadAndProbeConstruct<G> {
    fn new() -> Self {
        Self {
            gel_bead: G::new(),
            probe: G::new(),
        }
    }

    fn merge(&mut self, other: Self) {
        self.gel_bead.merge(other.gel_bead);
        self.probe.merge(other.probe);
    }
}

impl<G: JsonReport> JsonReport for GelBeadAndProbeConstruct<G> {
    fn to_json_reporter(&self) -> JsonReporter {
        let mut g_reporter = self.gel_bead.to_json_reporter();
        g_reporter.add_prefix("gel_bead");
        let mut p_reporter = self.probe.to_json_reporter();
        p_reporter.add_prefix("probe");
        g_reporter + p_reporter
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// Container for storing information related to segments of barcode
///
/// NOTE: It would have been preferrable to use an ArrayVec, a stack allocated
/// small vector. But none of the top crates implement Copy for that type (arrayvec, smallvec).
/// `tinyvec` implements Copy but it required Default.
///
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct Segments<G> {
    pub segment1: G,
    pub segment2: G,
    pub segment3: Option<G>,
    pub segment4: Option<G>,
}

impl<G> Segments<G> {
    pub fn is_two_part(&self) -> bool {
        self.segment3.is_none() && self.segment4.is_none()
    }
    pub fn map_result<F, K>(self, f: F) -> Result<Segments<K>>
    where
        F: FnMut(G) -> Result<K> + Copy,
    {
        self.into_iter().map(f).collect()
    }
    pub fn map<F, K>(self, f: F) -> Segments<K>
    where
        F: FnMut(G) -> K + Copy,
    {
        self.into_iter().map(f).collect()
    }
    pub fn as_ref(&self) -> Segments<&G> {
        Segments {
            segment1: &self.segment1,
            segment2: &self.segment2,
            segment3: self.segment3.as_ref(),
            segment4: self.segment4.as_ref(),
        }
    }
    pub fn as_mut_ref(&mut self) -> Segments<&mut G> {
        Segments {
            segment1: &mut self.segment1,
            segment2: &mut self.segment2,
            segment3: self.segment3.as_mut(),
            segment4: self.segment4.as_mut(),
        }
    }
    pub fn zip<K>(self, other: Segments<K>) -> Segments<(G, K)> {
        self.into_iter().zip_eq(other.into_iter()).collect()
    }
    pub fn array_vec(self) -> ArrayVec<[G; 4]> {
        self.into()
    }
}

impl<G> From<Segments<G>> for ArrayVec<[G; 4]> {
    fn from(v: Segments<G>) -> ArrayVec<[G; 4]> {
        ArrayVec::from_iter(v.into_iter())
    }
}

impl<G> FromIterator<G> for Segments<G> {
    fn from_iter<I: IntoIterator<Item = G>>(it: I) -> Self {
        let mut it = it.into_iter();
        let segments = Segments {
            segment1: it.next().unwrap(),
            segment2: it.next().unwrap(),
            segment3: it.next(),
            segment4: it.next(),
        };
        assert!(
            it.next().is_none(),
            "Attempted to create Segments with an iterator yielding more than 4 elements"
        );
        segments
    }
}

impl<G> IntoIterator for Segments<G> {
    type Item = G;
    type IntoIter = std::iter::Flatten<std::array::IntoIter<Option<G>, 4>>;

    fn into_iter(self) -> Self::IntoIter {
        [
            Some(self.segment1),
            Some(self.segment2),
            self.segment3,
            self.segment4,
        ]
        .into_iter()
        .flatten()
    }
}

impl<M: Metric> Metric for Segments<M> {
    fn new() -> Self {
        Segments {
            segment1: M::new(),
            segment2: M::new(),
            segment3: Option::<M>::new(),
            segment4: Option::<M>::new(),
        }
    }

    fn merge(&mut self, other: Self) {
        self.segment1.merge(other.segment1);
        self.segment2.merge(other.segment2);
        self.segment3.merge(other.segment3);
        self.segment4.merge(other.segment4);
    }
}

impl<J: JsonReport> JsonReport for Segments<J> {
    /// Report a metric for each of the four barcode segments.
    fn to_json_reporter(&self) -> JsonReporter {
        let mut reporter1 = self.segment1.to_json_reporter();
        let mut reporter2 = self.segment2.to_json_reporter();
        let mut reporter3 = self.segment3.to_json_reporter();
        let mut reporter4 = self.segment4.to_json_reporter();
        reporter1.add_prefix("segment1");
        reporter2.add_prefix("segment2");
        reporter3.add_prefix("segment3");
        reporter4.add_prefix("segment4");
        reporter1 + reporter2 + reporter3 + reporter4
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// Defines the different flavors of 10x barcodes.
///
/// This is used as a container for storing information about barcode segments. For example
/// `BarcodeConstruct<Whitelist>` is the whitelist for all the segments in the barcode
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub enum BarcodeConstruct<G> {
    GelBeadOnly(G),
    GelBeadAndProbe(GelBeadAndProbeConstruct<G>),
    Segmented(Segments<G>),
}

impl<G> BarcodeConstruct<G> {
    pub fn new_gel_bead_and_probe(gel_bead: G, probe: G) -> Self {
        GelBeadAndProbe(GelBeadAndProbeConstruct { gel_bead, probe })
    }

    pub fn map_result<F, K>(self, mut f: F) -> Result<BarcodeConstruct<K>>
    where
        F: FnMut(G) -> Result<K> + Copy,
    {
        Ok(match self {
            GelBeadOnly(g) => GelBeadOnly(f(g)?),
            GelBeadAndProbe(x) => GelBeadAndProbe(x.map_result(f)?),
            Segmented(t) => Segmented(t.map_result(f)?),
        })
    }

    pub fn map_option<F, K>(self, f: F) -> Option<BarcodeConstruct<K>>
    where
        F: FnOnce(G) -> Option<K> + Copy,
    {
        self.map_result(|g| f(g).ok_or_else(|| anyhow!(""))).ok()
    }

    pub fn map<F, K>(self, mut f: F) -> BarcodeConstruct<K>
    where
        F: FnMut(G) -> K + Copy,
    {
        match self {
            GelBeadOnly(g) => GelBeadOnly(f(g)),
            GelBeadAndProbe(x) => GelBeadAndProbe(x.map(f)),
            Segmented(t) => Segmented(t.map(f)),
        }
    }

    pub fn any<F, K>(self, f: F) -> bool
    where
        F: FnMut(G) -> bool,
    {
        self.iter().any(f)
    }

    pub fn array_vec(self) -> ArrayVec<[G; 4]> {
        let mut xs = ArrayVec::<[G; 4]>::new();
        match self {
            GelBeadOnly(g) => xs.push(g),
            GelBeadAndProbe(x) => xs.extend(x),
            Segmented(x) => xs.extend(x),
        }
        xs
    }

    pub fn iter(self) -> impl Iterator<Item = G> {
        self.into_iter()
    }

    pub fn gel_bead(self) -> G {
        self.option_gel_bead().unwrap()
    }

    pub fn option_gel_bead(self) -> Option<G> {
        match self {
            GelBeadOnly(g) => Some(g),
            GelBeadAndProbe(x) => Some(x.gel_bead),
            Segmented(_) => None,
        }
    }

    pub fn probe(self) -> G {
        self.option_probe().unwrap()
    }

    pub fn option_probe(self) -> Option<G> {
        match self {
            GelBeadOnly(_) => None,
            GelBeadAndProbe(g) => Some(g.probe),
            Segmented(_) => None,
        }
    }

    pub fn option_segments(self) -> Option<Segments<G>> {
        match self {
            GelBeadOnly(_) => None,
            GelBeadAndProbe(_) => None,
            Segmented(s) => Some(s),
        }
    }

    pub fn segments(self) -> Segments<G> {
        self.option_segments().unwrap()
    }

    pub fn as_ref(&self) -> BarcodeConstruct<&G> {
        match self {
            GelBeadOnly(ref g) => GelBeadOnly(g),
            GelBeadAndProbe(ref x) => GelBeadAndProbe(x.as_ref()),
            Segmented(t) => Segmented(t.as_ref()),
        }
    }

    pub fn as_mut_ref(&mut self) -> BarcodeConstruct<&mut G> {
        match self {
            GelBeadOnly(ref mut g) => GelBeadOnly(g),
            GelBeadAndProbe(ref mut x) => GelBeadAndProbe(x.as_mut_ref()),
            Segmented(ref mut t) => Segmented(t.as_mut_ref()),
        }
    }

    pub fn zip<K>(self, other: BarcodeConstruct<K>) -> BarcodeConstruct<(G, K)> {
        match (self, other) {
            (GelBeadOnly(g1), GelBeadOnly(g2)) => GelBeadOnly((g1, g2)),
            (GelBeadAndProbe(x1), GelBeadAndProbe(x2)) => GelBeadAndProbe(x1.zip(x2)),
            (Segmented(t1), Segmented(t2)) => Segmented(t1.zip(t2)),
            _ => unreachable!("Cannot zip incompatible BarcodeConstruct"),
        }
    }
}

impl<G> IntoIterator for BarcodeConstruct<G> {
    type Item = G;
    type IntoIter = arrayvec::IntoIter<[G; 4]>;

    fn into_iter(self) -> Self::IntoIter {
        self.array_vec().into_iter()
    }
}

impl<G, T> BarcodeConstruct<G>
where
    G: IntoIterator<Item = T>,
{
    pub fn flat_iter(self) -> impl Iterator<Item = T> {
        self.iter().flat_map(IntoIterator::into_iter)
    }
}

/// A thin wrapper around Option<BarcodeConstruct> to implement traits like
/// Metric and JsonReport
#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct BarcodeConstructMetric<M> {
    inner: Option<BarcodeConstruct<M>>,
}

impl<M> BarcodeConstructMetric<M> {
    pub fn inner(self) -> Option<BarcodeConstruct<M>> {
        self.inner
    }
}

impl<M> From<BarcodeConstruct<M>> for BarcodeConstructMetric<M> {
    fn from(src: BarcodeConstruct<M>) -> Self {
        BarcodeConstructMetric { inner: Some(src) }
    }
}

impl<M> BarcodeConstructMetric<M> {
    pub fn map<F, K>(self, f: F) -> BarcodeConstructMetric<K>
    where
        F: FnMut(M) -> K + Copy,
    {
        BarcodeConstructMetric {
            inner: self.inner.map(|construct| construct.map(f)),
        }
    }
}

impl<M: Metric> Metric for BarcodeConstructMetric<M> {
    fn new() -> Self {
        BarcodeConstructMetric { inner: None }
    }

    fn merge(&mut self, other: Self) {
        if let Some(other_inner) = other.inner {
            match self.inner.as_mut() {
                Some(construct) => match (construct, other_inner) {
                    (GelBeadOnly(x), GelBeadOnly(other)) => x.merge(other),
                    (GelBeadAndProbe(x), GelBeadAndProbe(other)) => x.merge(other),
                    (Segmented(x), Segmented(other)) => x.merge(other),
                    (_, _) => unreachable!("Cannot merge incompatible BarcodeConstruct"),
                },
                None => {
                    self.inner = Some(other_inner);
                }
            }
        }
    }
}

impl<M: JsonReport> JsonReport for BarcodeConstructMetric<M> {
    fn to_json_reporter(&self) -> JsonReporter {
        match self.inner.as_ref() {
            Some(GelBeadOnly(_g)) => JsonReporter::new(),
            Some(GelBeadAndProbe(g)) => g.to_json_reporter(),
            Some(Segmented(g)) => g.to_json_reporter(),
            None => JsonReporter::new(),
        }
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// One segment of the entire barcode construct. A segment is a contiguous range of bases in the
/// read associated with a specific whitelist. Each segment has a state and a sequence.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct BarcodeSegment {
    pub state: BarcodeSegmentState,
    pub sequence: BcSegSeq,
}

impl BarcodeSegment {
    pub fn new(seq: &[u8], state: BarcodeSegmentState) -> Self {
        BarcodeSegment {
            sequence: BcSegSeq::from_bytes(seq),
            state,
        }
    }
    pub fn is_valid(self) -> bool {
        self.state.is_valid()
    }

    pub fn is_valid_after_correction(self) -> bool {
        matches!(self.state, BarcodeSegmentState::ValidAfterCorrection)
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// Collection of the individual segments that make up the Barcode.
///
/// IMPORTANT: We deliberately do not implement traits such as `PartialEq`, `Eq`, `PartialOrd`,
/// `Ord` and `Hash` for this struct to reduce errors when using this struct because a corrected
/// and an uncorrected barcode with the same sequence will not be equal. This struct is only
/// intended to be used in cases where you need to know about different segments of the barcodes,
/// which is typically only in the early part of the SLFE pipeline until barcode correction
#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct SegmentedBarcode {
    gem_group: u16,
    segments: BarcodeConstruct<BarcodeSegment>,
}

impl SegmentedBarcode {
    pub fn new(gem_group: u16, segments: BarcodeConstruct<BarcodeSegment>) -> SegmentedBarcode {
        SegmentedBarcode {
            gem_group,
            segments,
        }
    }

    pub fn gel_bead_only(
        gem_group: u16,
        sequence: &[u8],
        state: BarcodeSegmentState,
    ) -> SegmentedBarcode {
        SegmentedBarcode {
            gem_group,
            segments: BarcodeConstruct::GelBeadOnly(BarcodeSegment::new(sequence, state)),
        }
    }

    /// Does this represent a valid whitelist barcode
    pub fn is_valid(self) -> bool {
        self.segments
            .map(BarcodeSegment::is_valid)
            .iter()
            .all(|valid| valid)
    }

    /// Does it represent a valid whitelist barcode with correction
    /// applied to **at least** one of the barcode segments
    pub fn is_valid_and_at_least_one_segment_corrected(self) -> bool {
        self.is_valid()
            && self
                .segments
                .iter()
                .any(BarcodeSegment::is_valid_after_correction)
    }

    pub fn segments(self) -> BarcodeConstruct<BarcodeSegment> {
        self.segments
    }

    pub fn segments_valid(self) -> BarcodeConstruct<bool> {
        self.segments.map(BarcodeSegment::is_valid)
    }

    pub fn sequence(&self) -> BarcodeConstruct<&[u8]> {
        self.segments.as_ref().map(|seg| seg.sequence.seq())
    }

    pub fn seq_iter(&self) -> impl Iterator<Item = &u8> {
        self.sequence().flat_iter()
    }

    pub fn len(self) -> usize {
        self.segments.map(|s| s.sequence.len()).iter().sum()
    }
    pub fn is_empty(self) -> bool {
        !self.segments.iter().any(|s| !s.sequence.is_empty())
    }

    pub fn sseq(self) -> BarcodeConstruct<BcSegSeq> {
        self.segments.map(|s| s.sequence)
    }

    pub fn gem_group(self) -> u16 {
        self.gem_group
    }

    pub fn segments_mut(&mut self) -> BarcodeConstruct<&mut BarcodeSegment> {
        self.segments.as_mut_ref()
    }
}

impl fmt::Display for SegmentedBarcode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Barcode::from(*self))
    }
}

impl From<SegmentedBarcode> for Barcode {
    fn from(src: SegmentedBarcode) -> Self {
        Barcode {
            gem_group: src.gem_group,
            sequence: {
                let mut sequence = BcSeq::new();
                for s in src.sequence() {
                    sequence.push_unchecked(s);
                }
                sequence
            },
            valid: src.is_valid(),
        }
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// A trait for objects that carry a 10x barcode, allowing for querying the barcode,
/// and correcting the barcode.
pub trait HasBarcode {
    fn barcode(&self) -> Barcode;
    fn segmented_barcode(&self) -> SegmentedBarcode;
    fn segmented_barcode_mut(&mut self) -> &mut SegmentedBarcode;
    fn raw_bc_construct_seq(&self) -> BarcodeConstruct<BcSegSeq>;
    fn raw_bc_construct_qual(&self) -> BarcodeConstruct<BcSegQual>;
    fn raw_bc_seq(&self) -> BcSeq;
    fn raw_bc_qual(&self) -> BcQual;
    fn has_probe_barcode(&self) -> bool;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plain_bc_parse() {
        assert_eq!(
            "TCTCAGAATGAGCCCA-1".parse::<Barcode>().unwrap(),
            Barcode::new(1, b"TCTCAGAATGAGCCCA", true)
        );
        assert!("TCTCAGAATGAGCCCA".parse::<Barcode>().is_err());
    }

    #[test]
    fn test_plain_bc_display() {
        assert_eq!(
            Barcode::new(1, b"TCTCAGAATGAGCCCA", true).to_string(),
            "TCTCAGAATGAGCCCA-1"
        );
    }

    #[test]
    fn test_barcode_to_barcode() {
        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::GelBeadOnly(BarcodeSegment::new(
                    b"ACAGTCATGTCCAAAT",
                    BarcodeSegmentState::NotChecked,
                )),
            }),
            Barcode::new(1, b"ACAGTCATGTCCAAAT", false)
        );

        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::new_gel_bead_and_probe(
                    BarcodeSegment::new(b"ACAGTCATGTCCAAAT", BarcodeSegmentState::Invalid,),
                    BarcodeSegment::new(b"CTGCCACT", BarcodeSegmentState::ValidBeforeCorrection)
                ),
            }),
            Barcode::new(1, b"ACAGTCATGTCCAAATCTGCCACT", false)
        );

        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::new_gel_bead_and_probe(
                    BarcodeSegment::new(
                        b"ACAGTCATGTCCAAAT",
                        BarcodeSegmentState::ValidBeforeCorrection,
                    ),
                    BarcodeSegment::new(b"CTGCCA", BarcodeSegmentState::ValidAfterCorrection)
                ),
            }),
            Barcode::new(1, b"ACAGTCATGTCCAAATCTGCCA", true)
        );

        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::new_gel_bead_and_probe(
                    BarcodeSegment::new(b"ACAGTCATGTCCAAAT", BarcodeSegmentState::NotChecked,),
                    BarcodeSegment::new(b"CTGCCACT", BarcodeSegmentState::NotChecked),
                ),
            }),
            Barcode::new(1, b"ACAGTCATGTCCAAATCTGCCACT", false)
        );
    }

    #[test]
    fn test_segments_from_iter() {
        let x4 = [1, 2, 3, 4];
        let s4 = Segments::from_iter(x4);
        assert_eq!(
            s4,
            Segments {
                segment1: 1,
                segment2: 2,
                segment3: Some(3),
                segment4: Some(4)
            }
        );
        itertools::assert_equal(s4.into_iter(), x4);

        let x2 = [1, 2];
        let s2 = Segments::from_iter(x2);
        assert_eq!(
            s2,
            Segments {
                segment1: 1,
                segment2: 2,
                segment3: None,
                segment4: None
            }
        );
        itertools::assert_equal(s2.into_iter(), x2);
    }

    #[test]
    #[should_panic]
    fn test_segments_from_iter_panic() {
        let _ = Segments::from_iter([1, 2, 3, 4, 5]);
    }
}
