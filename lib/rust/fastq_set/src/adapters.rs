//! Module for storing adapter sequences.
#![deny(missing_docs)]
use crate::WhichEnd;
use crate::adapter_trimmer::{Adapter, AdapterLoc};
use martian::{AsMartianPrimaryType, MartianPrimaryType};
use serde::{Deserialize, Serialize};
use strum::{Display, EnumString};

/// TSO sequence used in Fiveprime chemistries.
pub const SPACER: &str = "TTTCTTATATGGG";

/// Reverse complement of the Fiveprime TSO sequence.
pub const SPACER_RC: &str = "CCCATATAAGAAA";

/// Non-poly(dT) reverse primer used in FivePrime chemistries post GEM-RT.
pub const NON_POLYDT: &str = "AAGCAGTGGTATCAACGCAGAGTAC";

/// Reverse complement of non-poly(dT) reverse primer used in FivePrime chemistries post GEM-RT.
pub const NON_POLYDT_RC: &str = "GTACTCTGCGTTGATACCACTGCTT";

/// TSO sequence used in Threeprime chemistries.
pub const RT_PRIMER: &str = "AAGCAGTGGTATCAACGCAGAGTACAT";

/// Reverse complement of the TSO sequence used in Threeprime chemistries.
pub const RT_PRIMER_RC: &str = "ATGTACTCTGCGTTGATACCACTGCTT";

/// Poly-A sequence used in the protocol.
pub const POLY_A: &str = "AAAAAAAAAAAAAAAAAAAA";

/// Poly-T sequence used in the protocol.
pub const POLY_T: &str = "TTTTTTTTTTTTTTTTTTTT";

/// Flex capture sequence, named pCS1, X12, or partial X22_6, as sequenced on R1.
pub const X12_CAPTURE_SEQ: &str = "TTGCTAGGACCG";

/// Reverse complement of the flex capture sequence.
pub const X12_CAPTURE_SEQ_RC: &str = "CGGTCCTAGCAA";

/// Full Illumina R1 adapter sequence.
pub const ILLUMINA_R1_RC: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

/// Illumina R1 adapter sequence included in gel-beads.
pub const ILLUMINA_R1_SHORT: &str = "CTACACGACGCTCTTCCGATCT";

/// Reverse complement of the Illumina R1 adapter sequence included in gel-beads.
pub const ILLUMINA_R1_SHORT_RC: &str = "AGATCGGAAGAGCGTCGTGTAG";

/// The last 12-bp of the Illumina R1 adapter sequence included in gel-beads.
pub const ILLUMINA_R1_BP12: &str = "CTCTTCCGATCT";

/// Reverse complement of the last 12-bp of the Illumina R1 adapter sequence included in gel-beads.
pub const ILLUMINA_R1_BP12_RC: &str = "AGATCGGAAGAG";

/// Illumina R2 adapter sequence.
pub const ILLUMINA_R2: &str = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

/// Reverse complement of the Illumina R2 adapter sequence.
pub const ILLUMINA_R2_RC: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

/// Reverse complement of the Illumina P5 adapter sequence.
pub const ILLUMINA_P5_RC: &str = "AGATCTCGGTGGTCGCCGTATCATT";

/// Reverse complement of the Illumina P7 adapter sequence.
pub const ILLUMINA_P7_RC: &str = "ATCTCGTATGCCGTCTTCTGCTTG";

/// Small RNA R2 adapter sequence.
pub const SMALL_RNA_R2: &str = "CCTTGGCACCCGAGAATTCCA";

/// Reverse complement of the Small RNA R2 adapter sequence.
pub const SMALL_RNA_R2_RC: &str = "TGGAATTCTCGGGTGCCAAGG";

/// The last 12-bp of the Small RNA R2 adapter sequence.
pub const SMALL_RNA_R2_BP12: &str = "CCGAGAATTCCA";

/// Reverse complement of the last 12-bp of the Small RNA R2 adapter sequence.
pub const SMALL_RNA_R2_BP12_RC: &str = "TGGAATTCTCGG";

/// The last 12-bp of the Nextera R2 adapter sequence.
pub const NEXTERA_R2_BP12: &str = "TATAAGAGACAG";

/// Reverse complement of the last 12-bp of the Nextera R2 adapter sequence.
pub const NEXTERA_R2_BP12_RC: &str = "CTGTCTCTTATA";

/// 5-bp truncated truseq R2 adapter.
pub const TRUNCATED_R2: &str = "GATCT";

/// Reverse complement of the 5-bp truncated truseq R2 adapter.
pub const TRUNCATED_R2_RC: &str = "AGATC";

/// Enumeration of adapter types used in chemistry definitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, EnumString, Display, Serialize, Deserialize)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum AdapterKind {
    /// Illumina R1 in gel-beads (forward orientation).
    R1,
    /// Illumina R1 in gel-beads(reverse complement).
    R1Rc,
    /// Illumina R1 last 12 bp (forward orientation).
    R1Bp12,
    /// Illumina R1 last 12 bp (reverse complement).
    R1Bp12Rc,
    /// Illumina R2 (forward orientation).
    R2,
    /// Illumina R2 (reverse complement).
    R2Rc,
    /// Small RNA R2 (forward orientation).
    SmallRnaR2,
    /// Small RNA R2 (reverse complement).
    SmallRnaR2Rc,
    /// Small RNA R2 last 12 bp (forward orientation).
    SmallRnaR2Bp12,
    /// Small RNA R2 last 12 bp (reverse complement).
    SmallRnaR2Bp12Rc,
    /// Nextera R2 last 12 bp (forward orientation).
    NexteraR2Bp12,
    /// Nextera R2 last 12 bp (reverse complement).
    NexteraR2Bp12Rc,
    /// Template Switch Oligo (TSO) (forward orientation).
    RtPrimer,
    /// Template Switch Oligo (TSO) (reverse complement).
    RtPrimerRc,
    /// Non-poly(dT) reverse primer (forward orientation).
    NonPolydt,
    /// Non-poly(dT) reverse primer (reverse complement).
    NonPolydtRc,
    /// Poly(dT) primer (forward orientation).
    PolyT,
    /// Poly(dT) primer (reverse complement).
    PolyA,
    /// Capture sequence used in 5' chemistiries (forward orientation).
    Spacer,
    /// Capture sequence used in 5' chemistiries (reverse complement).
    SpacerRc,
    /// Flex capture sequence (forward orientation).
    X12Capture,
    /// Flex capture sequence (reverse complement).
    X12CaptureRc,
    /// 5-bp truncated truseq R2 (forward orientation).
    TruncatedR2,
    /// 5-bp truncated truseq R2 reverse-complement (reverse complement).
    TruncatedR2Rc,
    /// No adapter.
    End,
}

impl AsMartianPrimaryType for AdapterKind {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl AdapterKind {
    /// Build a 5' Adapter instance from this adapter kind.
    pub fn build_5p(self, forward: bool) -> Adapter {
        if forward {
            Adapter::new(
                self.to_string(),
                WhichEnd::FivePrime,
                AdapterLoc::Anywhere,
                self.sequence(),
            )
        } else {
            Adapter::new(
                self.get_revcomp_adapter().to_string(),
                WhichEnd::ThreePrime,
                AdapterLoc::Anywhere,
                self.get_revcomp_adapter().sequence(),
            )
        }
    }
    /// Build a 3' Adapter instance from this adapter kind.
    pub fn build_3p(self, forward: bool) -> Adapter {
        if forward {
            Adapter::new(
                self.to_string(),
                WhichEnd::ThreePrime,
                AdapterLoc::Anywhere,
                self.sequence(),
            )
        } else {
            Adapter::new(
                self.get_revcomp_adapter().to_string(),
                WhichEnd::FivePrime,
                AdapterLoc::Anywhere,
                self.get_revcomp_adapter().sequence(),
            )
        }
    }

    /// Returns the corresponding adapter sequence.
    pub fn sequence(&self) -> &str {
        match self {
            AdapterKind::R1 => ILLUMINA_R1_SHORT,
            AdapterKind::R1Rc => ILLUMINA_R1_SHORT_RC,
            AdapterKind::R1Bp12 => ILLUMINA_R1_BP12,
            AdapterKind::R1Bp12Rc => ILLUMINA_R1_BP12_RC,
            AdapterKind::R2 => ILLUMINA_R2,
            AdapterKind::R2Rc => ILLUMINA_R2_RC,
            AdapterKind::SmallRnaR2 => SMALL_RNA_R2,
            AdapterKind::SmallRnaR2Rc => SMALL_RNA_R2_RC,
            AdapterKind::SmallRnaR2Bp12 => SMALL_RNA_R2_BP12,
            AdapterKind::SmallRnaR2Bp12Rc => SMALL_RNA_R2_BP12_RC,
            AdapterKind::NexteraR2Bp12 => NEXTERA_R2_BP12,
            AdapterKind::NexteraR2Bp12Rc => NEXTERA_R2_BP12_RC,
            AdapterKind::RtPrimer => RT_PRIMER,
            AdapterKind::RtPrimerRc => RT_PRIMER_RC,
            AdapterKind::NonPolydt => NON_POLYDT,
            AdapterKind::NonPolydtRc => NON_POLYDT_RC,
            AdapterKind::PolyT => POLY_T,
            AdapterKind::PolyA => POLY_A,
            AdapterKind::Spacer => SPACER,
            AdapterKind::SpacerRc => SPACER_RC,
            AdapterKind::X12Capture => X12_CAPTURE_SEQ,
            AdapterKind::X12CaptureRc => X12_CAPTURE_SEQ_RC,
            AdapterKind::TruncatedR2 => TRUNCATED_R2,
            AdapterKind::TruncatedR2Rc => TRUNCATED_R2_RC,
            AdapterKind::End => "",
        }
    }

    /// Returns the corresponding reverse complement adapter kind.
    pub fn get_revcomp_adapter(&self) -> Self {
        match self {
            AdapterKind::R1 => AdapterKind::R1Rc,
            AdapterKind::R1Rc => AdapterKind::R1,
            AdapterKind::R1Bp12 => AdapterKind::R1Bp12Rc,
            AdapterKind::R1Bp12Rc => AdapterKind::R1Bp12,
            AdapterKind::R2 => AdapterKind::R2Rc,
            AdapterKind::R2Rc => AdapterKind::R2,
            AdapterKind::SmallRnaR2 => AdapterKind::SmallRnaR2Rc,
            AdapterKind::SmallRnaR2Rc => AdapterKind::SmallRnaR2,
            AdapterKind::SmallRnaR2Bp12 => AdapterKind::SmallRnaR2Bp12Rc,
            AdapterKind::SmallRnaR2Bp12Rc => AdapterKind::SmallRnaR2Bp12,
            AdapterKind::NexteraR2Bp12 => AdapterKind::NexteraR2Bp12Rc,
            AdapterKind::NexteraR2Bp12Rc => AdapterKind::NexteraR2Bp12,
            AdapterKind::RtPrimer => AdapterKind::RtPrimerRc,
            AdapterKind::RtPrimerRc => AdapterKind::RtPrimer,
            AdapterKind::NonPolydt => AdapterKind::NonPolydtRc,
            AdapterKind::NonPolydtRc => AdapterKind::NonPolydt,
            AdapterKind::PolyT => AdapterKind::PolyA,
            AdapterKind::PolyA => AdapterKind::PolyT,
            AdapterKind::Spacer => AdapterKind::SpacerRc,
            AdapterKind::SpacerRc => AdapterKind::Spacer,
            AdapterKind::X12Capture => AdapterKind::X12CaptureRc,
            AdapterKind::X12CaptureRc => AdapterKind::X12Capture,
            AdapterKind::TruncatedR2 => AdapterKind::TruncatedR2Rc,
            AdapterKind::TruncatedR2Rc => AdapterKind::TruncatedR2,
            AdapterKind::End => AdapterKind::End,
        }
    }
}
