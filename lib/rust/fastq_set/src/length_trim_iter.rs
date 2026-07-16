// Copyright (c) 10x Genomics, Inc. All rights reserved.

//! Iterator adapter for post-hoc trimming of read pair lengths.
//!
//! This is the preferred mechanism for applying Illumina read-length trimming
//! because it operates as a final step after all read-source adapters have been
//! composed (e.g. the long-read adapter trimmer).  This ensures that trimming
//! is applied uniformly regardless of whether reads come from standard FASTQ
//! parsing or from pseudo-paired reads produced by adapter splitting.

use crate::read_pair::{ReadPair, ReadPairStorage, TrimLengths, trim_read_pair_into};
use anyhow::Result;
use bytes::BytesMut;

/// Iterator adapter that trims reads in a stream of [`ReadPair`]s according to
/// [`TrimLengths`].
///
/// Constructed via [`trim_read_lengths`].
pub struct ReadPairLengthTrimIter<I> {
    inner: I,
    params: TrimLengths,
    /// Persistent buffer passed to [`ReadPair::trimmed_into`] on each trimmed read.
    /// Avoids a per-read `BytesMut` allocation; behavior depends on [`storage`](Self::storage).
    buffer: BytesMut,
    storage: ReadPairStorage,
}

impl<I: Iterator<Item = Result<ReadPair>>> Iterator for ReadPairLengthTrimIter<I> {
    type Item = Result<ReadPair>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|r| {
            r.map(|rp| trim_read_pair_into(rp, self.params, &mut self.buffer, self.storage))
        })
    }
}

impl<I> ReadPairLengthTrimIter<I> {
    /// Set the [`ReadPairStorage`] strategy used when trimming is applied.
    ///
    /// - [`ReadPairStorage::SharedBuffer`]: trimmed bytes are zero-copy slice-views
    ///   into a shared backing allocation that grows amortized (O(log n) allocations
    ///   for n reads).  Best throughput when all reads in the batch are consumed
    ///   together (the typical pipeline case).
    /// - [`ReadPairStorage::PerReadAllocation`] *(default)*: the internal buffer is
    ///   reused as a staging area; each output [`ReadPair`] owns its bytes
    ///   independently.  Preferable when reads may be held or dropped at different
    ///   times.
    ///
    /// Has no effect on reads that pass through unchanged (no-op trim params).
    pub fn storage(mut self, storage: ReadPairStorage) -> Self {
        self.storage = storage;
        self
    }
}

/// Wrap an iterator of [`ReadPair`] results to trim read lengths per `params`.
///
/// When `params` is a no-op, the underlying [`ReadPair::trimmed_into`] short-circuits
/// with no allocation.  When limits are set, each [`ReadPair`] whose stored length
/// exceeds the limit is rebuilt with only the trimmed bytes physically stored —
/// preserving correct serialization semantics.
///
/// The iterator holds a persistent [`BytesMut`] buffer and reuses it across reads,
/// avoiding per-read allocation overhead.  Use [`.storage()`](ReadPairLengthTrimIter::storage)
/// to select between [`ReadPairStorage::SharedBuffer`] (amortized, zero-copy) and
/// [`ReadPairStorage::PerReadAllocation`] (independent per-read ownership).
///
/// # Example
/// ```no_run
/// use fastq_set::length_trim_iter::trim_read_lengths;
/// use fastq_set::read_pair::TrimLengths;
/// # use fastq_set::read_pair_iter::FastqReadPairIter;
/// # let iter: FastqReadPairIter = todo!();
/// // Trim R1 to 90 bp, leave R2 at full length.
/// let trimmed = trim_read_lengths(iter, TrimLengths { r1: Some(90), r2: None });
/// ```
pub fn trim_read_lengths<I>(iter: I, params: TrimLengths) -> ReadPairLengthTrimIter<I> {
    ReadPairLengthTrimIter {
        inner: iter,
        params,
        buffer: BytesMut::new(),
        storage: ReadPairStorage::default(),
    }
}
