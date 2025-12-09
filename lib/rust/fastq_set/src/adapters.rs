//! Module for storing adapter sequences.
#![deny(missing_docs)]

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
pub const X12_CAPTURE_SEQ: &[u8] = b"TTGCTAGGACCG";

/// Full Illumina R1 adapter sequence.
pub const ILLUMINA_R1_RC: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

/// Illumina R1 adapter sequence included in gel-beads.
pub const ILLUMINA_R1_SHORT: &str = "CTACACGACGCTCTTCCGATCT";

/// Reverse complement of the Illumina R1 adapter sequence included in gel-beads.
pub const ILLUMINA_R1_SHORT_RC: &str = "AGATCGGAAGAGCGTCGTGTAG";

/// Illumina R2 adapter sequence.
pub const ILLUMINA_R2: &str = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

/// Reverse complement of the Illumina R2 adapter sequence.
pub const ILLUMINA_R2_RC: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

/// Reverse complement of the Illumina P5 adapter sequence.
pub const ILLUMINA_P5_RC: &str = "AGATCTCGGTGGTCGCCGTATCATT";

/// Reverse complement of the Illumina P7 adapter sequence.
pub const ILLUMINA_P7_RC: &str = "ATCTCGTATGCCGTCTTCTGCTTG";
