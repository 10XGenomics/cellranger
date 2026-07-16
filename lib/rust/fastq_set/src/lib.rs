// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
//! Tools for working with groups FASTQ of files.
//! Major functionality includes:
//! * Find groups FASTQs (R1/R2/I1/I2) following Illumina filename conventions
//! * Parsing flowcell information from Illumina FASTQ headers
//! * High-speed FASTQ I/O (via the `fastq` crate), with careful validation of
//!   FASTQ correctness and good error message.
//! * Containers for FASTQ read-pairs (along with index reads), providing access
//!   to 'technical' read components like cell barcode and UMI sequences.
//! * Flexible read trimming inspired by `cutadapt`
#![deny(missing_docs)]

pub mod adapter_trimmer;
pub mod adapters;
mod array;
mod background_iterator;
pub mod filenames;
mod illumina_header_info;
pub mod length_trim_iter;
pub mod long_read_adapter_trimmer;
pub mod metric_utils;
pub mod read_pair;
pub mod read_pair_iter;
pub mod sample_index_map;
pub mod squality;
pub mod sseq;

pub use crate::background_iterator::{
    BackgroundIterator, BackgroundIteratorExt, BackgroundIteratorScopedExt,
};
use crate::read_pair::ReadPair;
pub use crate::read_pair_iter::FastqReadPairIter;
pub use crate::sseq::SSeq;
pub use fastq::{OwnedRecord, Record};
pub use read_pair::WhichRead;
use serde::{Deserialize, Serialize};

/// See this URL for more information on an error code.
pub const ERROR_CODE_INFO: &str = "For more information: https://10xgen.com/pipeline-error-codes";

/// The result of processing a read-pair.
pub enum ProcessResult<T> {
    /// A processed read
    Processed(T),
    /// An unprocessable read
    Unprocessed {
        /// The unprocessed read
        read: ReadPair,
        /// The reason the read was unprocessed
        reason: String,
    },
}

/// Which end of a transcript reads come from
#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq, Eq)]
pub enum WhichEnd {
    /// The 3' end of a transcript
    #[serde(rename = "three_prime")]
    ThreePrime,
    /// The 5' end of a transcript
    #[serde(rename = "five_prime")]
    FivePrime,
}
