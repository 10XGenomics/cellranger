//! Crate for dealing with UMI related types and functionalities.
#![deny(missing_docs)]

use fastq_set::squality::SQualityGen;
use fastq_set::sseq::{HammingIterOpt, SSeqGen, SSeqOneHammingIter};
use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter};

mod info;

pub use info::UmiInfo;

/// The maximum length of a UMI.
pub const MAX_UMI_LENGTH: usize = 16;

/// UMI sequence.
pub type UmiSeq = SSeqGen<MAX_UMI_LENGTH>;

/// UMI quality sequence.
pub type UmiQual = SQualityGen<MAX_UMI_LENGTH>;

/// A container for a read UMI sequence. Can hold upto 16 bases
#[derive(
    Clone, Copy, Debug, Default, Hash, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize,
)]
#[serde(transparent)]
pub struct Umi {
    sequence: UmiSeq,
}

impl From<UmiSeq> for Umi {
    fn from(src: UmiSeq) -> Self {
        Umi { sequence: src }
    }
}

impl Umi {
    /// Create a new UMI from a byte slice.
    pub fn new(sequence: &[u8]) -> Umi {
        Umi {
            sequence: UmiSeq::from_bytes(sequence),
        }
    }

    /// Create a new UMI from a byte slice without checking validity.
    pub fn new_unchecked(sequence: &[u8]) -> Umi {
        Umi {
            sequence: UmiSeq::from_bytes_unchecked(sequence),
        }
    }

    /// Return the UMI sequence as a byte slice.
    pub fn sequence(&self) -> &[u8] {
        self.sequence.seq()
    }

    /// Return the UMI sequence as a byte slice.
    pub fn seq(&self) -> &[u8] {
        self.sequence.seq()
    }

    /// Return the UMI sequence as a `UmiSeq`.
    pub fn sseq(self) -> UmiSeq {
        self.sequence
    }

    /// Return true if the UMI is valid.
    ///
    /// The qual scores are required to determine this.
    pub fn is_valid(self, qual: &UmiQual) -> bool {
        is_valid(&self.sequence, qual)
    }

    /// Create an iterator over UMIs one Hamming distance away from this one.
    pub fn one_hamming_iter(self, opt: HammingIterOpt) -> UmiOneHammingIter {
        UmiOneHammingIter {
            inner: self.sequence.one_hamming_iter(opt),
        }
    }
}

/// Return true if the provided sequence and qual are a valid UMI.
///
/// A valid UMI has no Ns, is not a homopolymer, and has no bases with a qual
/// score less than a minimum threshold.
pub fn is_valid(seq: &UmiSeq, qual: &UmiQual) -> bool {
    let has_n = seq.has_n();
    let is_homopolymer = seq.is_homopolymer();

    const UMI_MIN_QV: u8 = 10;
    let low_min_qual = qual.min_qual() < UMI_MIN_QV;
    !(has_n || is_homopolymer || low_min_qual)
}

/// An iterator over UMIs one Hamming distance away from a given UMI.
pub struct UmiOneHammingIter {
    inner: SSeqOneHammingIter<MAX_UMI_LENGTH>,
}

impl Iterator for UmiOneHammingIter {
    type Item = Umi;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|s| Umi { sequence: s })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl Display for Umi {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.sequence, f)
    }
}

/// Is a UMI consistent with splice junctions (Txomic) or not
/// This label is useful in intron mode to distinguish transcriptomic UMIs from
/// non-transcriptomic UMIs. The order of fields in this enum ensures that a
/// transcriptomic read (aligned with splice junctions) with the lexicographically
/// smallest qname is selected as the UMI representative of its dup cohort of reads
/// in intron mode. In regular mode this enum is always `Txomic` at the moment.
/// TODO: generalize this enum into a general bit flag where each bit will be
/// used to represent other molecule properties, like sense/anti-sense etc.
#[derive(Clone, Copy, Default, PartialOrd, Ord, PartialEq, Eq, Serialize, Deserialize, Debug)]
pub enum UmiType {
    /// Transcriptomic: exonic and consistent with splice junctions
    #[default]
    Txomic,
    /// Non-transcriptomic
    NonTxomic,
}

impl UmiType {
    /// Decode the UmiType enum into a u32 for storage in the molecule_info h5
    pub fn from_u32(value: &u32) -> UmiType {
        match value {
            0 => UmiType::NonTxomic,
            1 => UmiType::Txomic,
            v => unimplemented!(
                "Invalid value {} encountered while parsing the molecule
info. Valid values are 0 for non-transcriptomic and 1 for transcriptomic UMIs.",
                v
            ),
        }
    }

    /// Encode the UmiType enum for writing to the molecule_info h5
    pub fn to_u32(&self) -> u32 {
        match self {
            UmiType::NonTxomic => 0_u32,
            UmiType::Txomic => 1_u32,
        }
    }
}
