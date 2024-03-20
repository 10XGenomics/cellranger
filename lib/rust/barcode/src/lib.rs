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
use binned::{SquareBinIndex, SquareBinRowOrColumnIndex};
pub use corrector::BarcodeCorrector;
use fastq_set::squality::SQualityGen;
use fastq_set::sseq::SSeqGen;
use itertools::{zip_eq, Itertools};
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

/// Type for storing barcode content which could be a sequence or a spatial index.
/// The spatial index flavor is used in Visium HD where we have a square grid of spots.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub enum BarcodeContent {
    Sequence(BcSeq),
    SpatialIndex(SquareBinIndex),
}

impl fmt::Display for BarcodeContent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BarcodeContent::Sequence(seq) => write!(f, "{seq}"),
            BarcodeContent::SpatialIndex(index) => write!(f, "{index}"),
        }
    }
}

impl FromStr for BarcodeContent {
    type Err = anyhow::Error;

    /// Parse a barcode from its string representation "AACCGGTT-1".
    fn from_str(content_str: &str) -> Result<BarcodeContent> {
        if content_str.starts_with(binned::SQUARE_BIN_PREFIX) {
            Ok(BarcodeContent::SpatialIndex(content_str.parse()?))
        } else {
            Ok(BarcodeContent::Sequence(BcSeq::from_bytes(
                content_str.as_bytes(),
            )))
        }
    }
}

impl BarcodeContent {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        std::str::from_utf8(bytes)?.parse()
    }
    /// The spatial index of the barcode, if it has one. Panics for
    /// sequence barcodes
    pub fn spatial_index(&self) -> SquareBinIndex {
        match self {
            BarcodeContent::Sequence(_) => {
                panic!("Cannot get spatial index from sequence barcode")
            }
            BarcodeContent::SpatialIndex(index) => *index,
        }
    }
    /// The sequence of the barcode, if it has one. Panics for
    /// spatial index barcodes
    pub fn sequence(&self) -> &BcSeq {
        match self {
            BarcodeContent::Sequence(seq) => seq,
            BarcodeContent::SpatialIndex(_) => {
                panic!("Cannot get sequence from spatial barcode")
            }
        }
    }
}

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
    content: BarcodeContent,
}

impl Barcode {
    pub fn is_valid(self) -> bool {
        self.valid
    }
    pub fn from_bytes(seq_gg: &[u8]) -> Result<Self> {
        std::str::from_utf8(seq_gg)?.parse()
    }
    /// Create a barcode with the given sequence.
    pub fn with_seq(gem_group: u16, sequence: BcSeq, valid: bool) -> Self {
        Barcode {
            gem_group,
            content: BarcodeContent::Sequence(sequence),
            valid,
        }
    }
    /// Create a barcode with a spatial index. A barcode with a spatial index is a valid barcode
    pub fn with_square_bin_index(gem_group: u16, index: SquareBinIndex) -> Self {
        Barcode {
            gem_group,
            content: BarcodeContent::SpatialIndex(index),
            valid: true,
        }
    }

    pub fn with_content(gem_group: u16, content: BarcodeContent, valid: bool) -> Self {
        Barcode {
            gem_group,
            valid,
            content,
        }
    }
    pub fn gem_group(self) -> u16 {
        self.gem_group
    }
    /// The sequence of the barcode, if it has one. Panics for
    /// spatial index barcodes
    fn sequence(&self) -> &BcSeq {
        self.content.sequence()
    }
    /// The sequence of the barcode as byte slice, if it has one. Panics for
    /// spatial index barcodes
    pub fn sequence_bytes(&self) -> &[u8] {
        self.sequence().as_bytes()
    }
    pub fn content(&self) -> &BarcodeContent {
        &self.content
    }

    /// The spatial index of the barcode, if it has one. Panics otherwise
    pub fn spatial_index(&self) -> SquareBinIndex {
        match self.content {
            BarcodeContent::Sequence(_) => {
                panic!("Cannot get spatial index from sequence barcode")
            }
            BarcodeContent::SpatialIndex(index) => index,
        }
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
        write!(f, "{}-{}", self.content, self.gem_group)
    }
}

impl FromStr for Barcode {
    type Err = anyhow::Error;

    /// Parse a barcode from its string representation "AACCGGTT-1".
    fn from_str(seq_gg: &str) -> Result<Barcode> {
        let Some((content_str, gem_group_str)) = seq_gg.split_once('-') else {
            bail!("invalid barcode: '{seq_gg}'")
        };
        let gem_group: u16 = gem_group_str.parse()?;
        Ok(Barcode {
            gem_group,
            valid: true,
            content: content_str.parse()?,
        })
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
    fn merge(&mut self, other: Self) {
        self.gel_bead.merge(other.gel_bead);
        self.probe.merge(other.probe);
    }
}

impl<G: JsonReport> JsonReport for GelBeadAndProbeConstruct<G> {
    fn to_json_reporter(&self) -> JsonReporter {
        self.gel_bead.to_json_reporter().add_prefix("gel_bead")
            + self.probe.to_json_reporter().add_prefix("probe")
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
    pub fn map_result<F, K>(self, mut f: F) -> Result<Segments<K>>
    where
        F: FnMut(G) -> Result<K> + Copy,
    {
        Ok(Segments {
            segment1: f(self.segment1)?,
            segment2: f(self.segment2)?,
            segment3: self.segment3.map(f).transpose()?,
            segment4: self.segment4.map(f).transpose()?,
        })
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
        zip_eq(self, other).collect()
    }
    pub fn array_vec(self) -> ArrayVec<[G; 4]> {
        self.into()
    }
}

impl<G> From<Segments<G>> for ArrayVec<[G; 4]> {
    fn from(v: Segments<G>) -> ArrayVec<[G; 4]> {
        ArrayVec::from_iter(v)
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
        self.segment1.to_json_reporter().add_prefix("segment1")
            + self.segment2.to_json_reporter().add_prefix("segment2")
            + self.segment3.to_json_reporter().add_prefix("segment3")
            + self.segment4.to_json_reporter().add_prefix("segment4")
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

impl<G: IntoIterator<Item = T>, T> BarcodeConstruct<G> {
    pub fn flat_iter(self) -> impl Iterator<Item = T> {
        self.iter().flat_map(IntoIterator::into_iter)
    }
}

/// A thin wrapper around Option<BarcodeConstruct> to implement traits like
/// Metric and JsonReport
#[derive(Copy, Clone, Default, Serialize, Deserialize)]
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
            Some(GelBeadOnly(_g)) => JsonReporter::default(),
            Some(GelBeadAndProbe(g)) => g.to_json_reporter(),
            Some(Segmented(g)) => g.to_json_reporter(),
            None => JsonReporter::default(),
        }
    }
}

/// This enum represents the content of a barcode segment.
///
/// Each segment of the barcode can either be a sequence or a spatial index.
/// The spatial index flavor is used in Visium HD where we have a square grid of spots.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub enum BarcodeSegmentContent {
    Sequence(BcSegSeq),
    SpatialIndex(BcSegSeq, SquareBinRowOrColumnIndex),
}

impl From<BcSegSeq> for BarcodeSegmentContent {
    fn from(seq: BcSegSeq) -> BarcodeSegmentContent {
        BarcodeSegmentContent::Sequence(seq)
    }
}

impl BarcodeSegmentContent {
    fn spatial_index(&self) -> SquareBinRowOrColumnIndex {
        match self {
            BarcodeSegmentContent::Sequence(_) => {
                panic!("Cannot get spatial index from sequence barcode")
            }
            BarcodeSegmentContent::SpatialIndex(_, index) => *index,
        }
    }
}

/* ---------------------------------------------------------------------------------------------- */
/// One segment of the entire barcode construct. A segment is a contiguous range of bases in the
/// read associated with a specific whitelist. Each segment has a state and a sequence.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct BarcodeSegment {
    pub state: BarcodeSegmentState,
    pub content: BarcodeSegmentContent,
}

impl BarcodeSegment {
    pub fn new(content: BarcodeSegmentContent, state: BarcodeSegmentState) -> Self {
        BarcodeSegment { content, state }
    }
    pub fn with_sequence(seq: &[u8], state: BarcodeSegmentState) -> Self {
        BarcodeSegment {
            content: BarcodeSegmentContent::Sequence(BcSegSeq::from_bytes(seq)),
            state,
        }
    }
    pub fn with_sequence_unchecked(seq: &[u8]) -> Self {
        BarcodeSegment {
            content: BarcodeSegmentContent::Sequence(BcSegSeq::from_bytes_unchecked(seq)),
            state: BarcodeSegmentState::NotChecked,
        }
    }
    pub fn is_valid(self) -> bool {
        self.state.is_valid()
    }

    pub fn is_valid_after_correction(self) -> bool {
        matches!(self.state, BarcodeSegmentState::ValidAfterCorrection)
    }

    pub fn sequence(&self) -> &BcSegSeq {
        match &self.content {
            BarcodeSegmentContent::Sequence(seq) => seq,
            BarcodeSegmentContent::SpatialIndex(seq, _) => seq,
        }
    }

    pub fn sequence_mut(&mut self) -> &mut BcSegSeq {
        match &mut self.content {
            BarcodeSegmentContent::Sequence(seq) => seq,
            BarcodeSegmentContent::SpatialIndex(_, _) => {
                unreachable!("Cannot get mutable sequence from spatial barcode")
            }
        }
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
            segments: BarcodeConstruct::GelBeadOnly(BarcodeSegment::with_sequence(sequence, state)),
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
        self.segments.as_ref().map(|seg| seg.sequence().seq())
    }

    pub fn sseq(self) -> BarcodeConstruct<BcSegSeq> {
        self.segments.map(|s| *s.sequence())
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
        if src
            .segments
            .iter()
            .all(|seg| matches!(seg.content, BarcodeSegmentContent::SpatialIndex(_, _)))
        {
            assert!(
                src.segments.segments().is_two_part(),
                "Spatial barcode must have two segments"
            );
            let segments = src.segments.array_vec();
            Barcode {
                gem_group: src.gem_group,
                valid: src.is_valid(),
                content: BarcodeContent::SpatialIndex(SquareBinIndex {
                    // NOTE: The order of the segments is important here. The first segment
                    // is the X barcode (col) and the second segment is the Y barcode (row).
                    row: segments[1].content.spatial_index().index,
                    col: segments[0].content.spatial_index().index,
                    size_um: segments
                        .into_iter()
                        .map(|s| s.content.spatial_index().size_um)
                        .dedup()
                        .exactly_one()
                        .unwrap(),
                }),
            }
        } else {
            Barcode::with_seq(
                src.gem_group,
                {
                    let mut sequence = BcSeq::new();
                    for s in src.sequence() {
                        sequence.push_unchecked(s);
                    }
                    sequence
                },
                src.is_valid(),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn barcode(seq: &[u8], valid: bool) -> Barcode {
        Barcode::with_seq(1, BcSeq::from_bytes(seq), valid)
    }

    #[test]
    fn test_plain_bc_parse() {
        assert_eq!(
            "TCTCAGAATGAGCCCA-1".parse::<Barcode>().unwrap(),
            barcode(b"TCTCAGAATGAGCCCA", true)
        );
        assert!("TCTCAGAATGAGCCCA".parse::<Barcode>().is_err());
    }

    #[test]
    fn test_plain_bc_display() {
        assert_eq!(
            barcode(b"TCTCAGAATGAGCCCA", true).to_string(),
            "TCTCAGAATGAGCCCA-1"
        );
    }

    #[test]
    fn test_barcode_to_barcode() {
        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::GelBeadOnly(BarcodeSegment::with_sequence(
                    b"ACAGTCATGTCCAAAT",
                    BarcodeSegmentState::NotChecked,
                )),
            }),
            barcode(b"ACAGTCATGTCCAAAT", false)
        );

        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::new_gel_bead_and_probe(
                    BarcodeSegment::with_sequence(
                        b"ACAGTCATGTCCAAAT",
                        BarcodeSegmentState::Invalid,
                    ),
                    BarcodeSegment::with_sequence(
                        b"CTGCCACT",
                        BarcodeSegmentState::ValidBeforeCorrection
                    )
                ),
            }),
            barcode(b"ACAGTCATGTCCAAATCTGCCACT", false)
        );

        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::new_gel_bead_and_probe(
                    BarcodeSegment::with_sequence(
                        b"ACAGTCATGTCCAAAT",
                        BarcodeSegmentState::ValidBeforeCorrection,
                    ),
                    BarcodeSegment::with_sequence(
                        b"CTGCCA",
                        BarcodeSegmentState::ValidAfterCorrection
                    )
                ),
            }),
            barcode(b"ACAGTCATGTCCAAATCTGCCA", true)
        );

        assert_eq!(
            Barcode::from(SegmentedBarcode {
                gem_group: 1,
                segments: BarcodeConstruct::new_gel_bead_and_probe(
                    BarcodeSegment::with_sequence(
                        b"ACAGTCATGTCCAAAT",
                        BarcodeSegmentState::NotChecked,
                    ),
                    BarcodeSegment::with_sequence(b"CTGCCACT", BarcodeSegmentState::NotChecked),
                ),
            }),
            barcode(b"ACAGTCATGTCCAAATCTGCCACT", false)
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
        itertools::assert_equal(s4, x4);

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
        itertools::assert_equal(s2, x2);
    }

    #[test]
    #[should_panic]
    fn test_segments_from_iter_panic() {
        let _ = Segments::from_iter([1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_barcode_from_segmented_barcode() {
        let barcode = Barcode::from(SegmentedBarcode {
            gem_group: 1,
            segments: BarcodeConstruct::GelBeadOnly(BarcodeSegment::with_sequence(
                b"ACAGTCATGTCCAAAT",
                BarcodeSegmentState::Invalid,
            )),
        });
        assert_eq!(barcode.gem_group(), 1);
        assert_eq!(barcode.sequence_bytes(), b"ACAGTCATGTCCAAAT");
        assert!(!barcode.is_valid());

        let barcode = Barcode::from(SegmentedBarcode {
            gem_group: 1,
            segments: BarcodeConstruct::new_gel_bead_and_probe(
                BarcodeSegment::with_sequence(b"ACAGTCATGTCCAAAT", BarcodeSegmentState::Invalid),
                BarcodeSegment::with_sequence(
                    b"CTGCCACT",
                    BarcodeSegmentState::ValidBeforeCorrection,
                ),
            ),
        });
        assert_eq!(barcode.gem_group(), 1);
        assert_eq!(barcode.sequence_bytes(), b"ACAGTCATGTCCAAATCTGCCACT");
        assert!(!barcode.is_valid());

        let spatial_segment = |seq: &[u8], index: usize| {
            BarcodeSegmentContent::SpatialIndex(
                BcSegSeq::from_bytes(seq),
                SquareBinRowOrColumnIndex { index, size_um: 2 },
            )
        };

        let barcode = Barcode::from(SegmentedBarcode {
            gem_group: 1,
            segments: BarcodeConstruct::Segmented(Segments {
                segment1: BarcodeSegment {
                    state: BarcodeSegmentState::ValidBeforeCorrection,
                    content: spatial_segment(b"ACAGTCATGTCCAAAT", 2),
                },
                segment2: BarcodeSegment {
                    state: BarcodeSegmentState::ValidAfterCorrection,
                    content: spatial_segment(b"AAACCTGAGAACAACT", 4),
                },
                segment3: None,
                segment4: None,
            }),
        });
        assert_eq!(barcode.gem_group(), 1);
        assert_eq!(barcode.to_string(), "s_002um_00004_00002-1");
        assert!(barcode.is_valid());

        let barcode = Barcode::from(SegmentedBarcode {
            gem_group: 1,
            segments: BarcodeConstruct::Segmented(Segments {
                segment1: BarcodeSegment::with_sequence(
                    b"ACAGTCATGTCCAAAT",
                    BarcodeSegmentState::Invalid,
                ),
                segment2: BarcodeSegment {
                    state: BarcodeSegmentState::ValidAfterCorrection,
                    content: spatial_segment(b"AAACCTGAGAACAACT", 4),
                },
                segment3: None,
                segment4: None,
            }),
        });
        assert_eq!(barcode.gem_group(), 1);
        assert_eq!(barcode.to_string(), "ACAGTCATGTCCAAATAAACCTGAGAACAACT-1");
        assert!(!barcode.is_valid());
    }

    #[test]
    fn test_segments_map_result() {
        let segments = Segments {
            segment1: 1,
            segment2: 2,
            segment3: Some(3),
            segment4: Some(4),
        };
        assert!(segments.map_result(Ok).is_ok());
        assert!(segments.map_result::<_, i32>(|_| bail!("")).is_err());
    }
}
