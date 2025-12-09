// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

//! Container for the FASTQ data from a single sequencing 'cluster',
//! including the primary 'R1' and 'R2' and index 'I1' and 'I2' reads.

use crate::WhichEnd;
use anyhow::{Result, bail, ensure};
use bytes::{Bytes, BytesMut};
use fastq::Record;
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::{fmt, ops};

// ▪ fastq crate imposes an upper limit of 68*1024 bytes for a fastq record (id & seq & qual)
// ▪ bit-packed RpRange storage imposes a u15 limit of 32,767 on offsets and lengths
const MAX_READ_HEAD: u16 = 767;
const MAX_READ_LEN: u16 = 32_000;

/// Pointers into a buffer that identify the positions of lines from a FASTQ record
/// header exists at buf[start .. head], seq exists at buf[head .. seq], etc.
#[derive(Deserialize, Serialize, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct ReadOffset {
    exists: bool,
    start: u16,
    head: u16,
    seq: u16,
    qual: u16,
}

impl ReadOffset {
    fn seq_len(&self) -> Option<usize> {
        if self.exists {
            Some((self.seq - self.head) as usize)
        } else {
            None
        }
    }
}

/// The possible reads from a Illumina cluster. R1 and R2 are the two
/// 'primary' reads, I1 and I2 are the two 'index' samples. I1 is
/// often referred to as the 'sample index read', or I7.  I2 contains
/// the 10x barcode sequence in some 10x assays.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum WhichRead {
    R1 = 0,
    R2 = 1,
    I1 = 2,
    I2 = 3,
}

macro_rules! whichread_from {
    ($type:ident) => {
        impl From<$type> for WhichRead {
            fn from(i: $type) -> Self {
                match i {
                    0 => WhichRead::R1,
                    1 => WhichRead::R2,
                    2 => WhichRead::I1,
                    3 => WhichRead::I2,
                    _ => panic!(
                        "Values other than 0,1,2,3 cannot be converted to WhichRead. Got {}",
                        i
                    ),
                }
            }
        }
    };
}
whichread_from!(usize);
whichread_from!(u32);

impl WhichRead {
    pub const fn read_types() -> [WhichRead; 4] {
        [WhichRead::R1, WhichRead::R2, WhichRead::I1, WhichRead::I2]
    }
}

impl fmt::Display for WhichRead {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            WhichRead::R1 => "read1",
            WhichRead::R2 => "read2",
            WhichRead::I1 => "index1",
            WhichRead::I2 => "index2",
        })
    }
}

impl std::str::FromStr for WhichRead {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        Ok(match s {
            "R1" => WhichRead::R1,
            "R2" => WhichRead::R2,
            "I1" => WhichRead::I1,
            "I2" => WhichRead::I2,
            _ => bail!("could not parse WhichRead from '{s}'"),
        })
    }
}

/// Components of a FASTQ record.
#[derive(Debug, Copy, Clone)]
pub enum ReadPart {
    Header,
    Seq,
    Qual,
}

/// Compact representation of selected read and an interval in that read.
/// Supports offsets and lengths up to 32K.
/// Internally it is stored as a `u32` with the following bit layout
/// ```text
/// +-----------+--------------+-----------+
/// | WhichRead | Start Offset | Length    |
/// | (2 bits)  | (15 bits)    | (15 bits) |
/// +-----------+--------------+-----------+
/// ```
/// Length is optional, with `None` indicating everything until the end of the read.
///
/// # Example
/// ```rust
/// extern crate fastq_set;
/// extern crate fastq;
/// use fastq_set::read_pair::{RpRange, WhichRead, ReadPart, ReadPair};
/// use fastq_set::WhichEnd;
/// use fastq::OwnedRecord;
/// let read1 = OwnedRecord {
///     head: b"some_name".to_vec(),
///     seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
///     qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
///     sep: None,
/// };
/// let input = [Some(read1), None, None, None];
///
/// // WARNING: DO NOT USE THIS FUNCTION IF YOU ARE STREAMING FASTQ DATA
/// // AND WANT TO CREATE ReadPair STRUCTS.
/// // USE `fastq_set::read_pair_iter::ReadPairIter` INSTEAD.
/// let read_pair = ReadPair::new(input);
/// // Let's say the read1 is of the form BC(16)-UMI(10)-Insert. This example will
/// // setup different RpRanges to represent these ranges.
///
/// let barcode_range = RpRange::new(WhichRead::R1, 0, Some(16));
/// assert_eq!(read_pair.get_range(barcode_range, ReadPart::Seq).unwrap(), b"GTCGCACTGATCTGGG".to_vec().as_slice());
///
/// let umi_range = RpRange::new(WhichRead::R1, 16, Some(10));
/// assert_eq!(read_pair.get_range(umi_range, ReadPart::Seq).unwrap(), b"TTAGGCGCGG".to_vec().as_slice());
///
/// let mut r1_range = RpRange::new(WhichRead::R1, 26, None); // None => everything beyond offset
/// assert_eq!(read_pair.get_range(r1_range, ReadPart::Seq).unwrap(), b"AGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec().as_slice());
///
/// // Let's say you want to trim first 5 bases in r1
/// r1_range.trim(WhichEnd::FivePrime, 5);
/// assert_eq!(r1_range.offset(), 31);
/// assert_eq!(r1_range.len(), None);
///
/// // Let's say you only want to consider first 50 bases of the read1.seq
/// // and update the r1_range accordingly
/// let useful_range = RpRange::new(WhichRead::R1, 0, Some(50));
/// r1_range.intersect(useful_range);
/// assert_eq!(r1_range.offset(), 31);
/// assert_eq!(r1_range.len(), Some(19));
///
/// ```
#[derive(Serialize, Deserialize, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct RpRange {
    val: u32,
}

impl fmt::Debug for RpRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RpRange {{ read: {:?}, offset: {}, len: {:?}, val: {:#b} }}",
            self.read(),
            self.offset(),
            self.len(),
            self.val
        )
    }
}

#[allow(clippy::len_without_is_empty)]
impl RpRange {
    /// Create a `RpRange` that represent the interval [`offset`, `offset + len`) in
    /// the `read` component of a ReadPair.
    ///
    /// # Args
    /// * `read` - Specify `WhichRead`
    /// * `offset` - Start of the interval. Must be less than 2^15 (=32,768)
    /// * `len` - Optional length that determines the end of the interval.
    ///   A value `None` indicates everything from `offset` until the end of the `read`.
    ///   Must be less than 2^15 (=32,768)
    ///
    /// # Panics
    /// * If `offset` or `len` is >= `2^15`
    ///
    /// # Example
    /// ```rust
    /// use fastq_set::read_pair::RpRange;
    /// use fastq_set::read_pair::WhichRead;
    /// let range = RpRange::new(WhichRead::R1, 10, Some(50));
    /// assert!(range.read() == WhichRead::R1);
    /// assert!(range.offset() == 10);
    /// assert!(range.len() == Some(50));
    /// ```
    ///
    /// # Tests
    /// * `test_rprange_invalid_offset()` - Test that this function panics with
    ///   an offset that is too large
    /// * `test_rprange_invalid_len()` - Test that this function panics with
    ///   a length that is too large
    /// * `prop_test_rprange_representation()` - Test that arbitrary construction of RpRange
    ///   stores the values correctly.
    pub fn new(read: WhichRead, offset: usize, len: Option<usize>) -> RpRange {
        assert!(offset < (1 << 15));
        let len_bits = match len {
            Some(v) => {
                assert!(v < 0x7FFF);
                v
            }
            None => 0x7FFF,
        };

        let val = ((read as u32) << 30) | ((offset as u32) << 15) | len_bits as u32;
        RpRange { val }
    }

    #[inline]
    /// Retrieive the read from the internal representation
    pub fn read(self) -> WhichRead {
        let k = self.val >> 30;
        WhichRead::from(k)
    }

    #[inline]
    /// Retrieive the offset from the internal representation
    pub fn offset(self) -> usize {
        ((self.val >> 15) & 0x7FFF) as usize
    }

    #[inline]
    /// Retrieive the (optional) length from the internal representation
    pub fn len(self) -> Option<usize> {
        let len_bits = self.val & 0x7FFF;
        if len_bits == 0x7FFF {
            None
        } else {
            Some(len_bits as usize)
        }
    }

    // Slice the input to the range [offset, offset + len).
    fn slice<'a>(&self, input: &'a [u8]) -> Option<&'a [u8]> {
        let o = self.offset();
        match self.len() {
            Some(l) if o + l <= input.len() => Some(&input[o..o + l]),
            None if o <= input.len() => Some(&input[o..]),
            _ => None,
        }
    }

    // Set the length
    // # Tests
    // * `prop_test_rprange_setter()`
    fn set_len(&mut self, len: usize) {
        assert!(len < 0x7FFF);
        self.val = (self.val & (!0x7FFFu32)) | len as u32;
    }

    // Set the offset
    // # Tests
    // * `prop_test_rprange_setter()`
    fn set_offset(&mut self, offset: usize) {
        assert!(offset < (1 << 15));
        self.val = (self.val & (!(0x7FFFu32 << 15))) | ((offset as u32) << 15);
    }

    /// Intersect the `RpRange` with another.
    ///
    /// # Input
    /// * `other` - Intersect `self` with this `RpRange`
    ///
    /// # Panics
    /// * If `self` and `other` does not point to the same read
    ///
    /// # Example
    /// ```rust
    /// use fastq_set::read_pair::{RpRange, WhichRead};
    /// let read_seq = b"AAGCAGGGGCGGGCAAATCCAGCCGTTACCTTACACGCCCCACTGGGAAG";
    /// let full_range = RpRange::new(WhichRead::R1, 0, Some(read_seq.len()));
    /// let mut r1_range = RpRange::new(WhichRead::R1, 20, None);
    /// r1_range.intersect(full_range);
    /// assert_eq!(r1_range.read(), WhichRead::R1);
    /// assert_eq!(r1_range.offset(), 20);
    /// assert_eq!(r1_range.len(), Some(read_seq.len()-20));
    /// ```
    ///
    /// # Tests
    /// * `test_rprange_intersect_panic()` - Make sure that this function panics
    ///   if the reads do not match
    /// * `test_rprange_intersect_both_open()` - Test a case when both lengths are not set
    /// * `test_rprange_intersect_self_open()` - Test a case when only self length is set
    /// * `test_rprange_intersect_other_open()` - Test a case when only other length is set
    /// * `test_rprange_intersect_both_closed()` - Test a case when both lengths are set
    pub fn intersect(&mut self, other: RpRange) {
        use std::cmp::{max, min};
        assert!(
            self.read() == other.read(),
            "Self = {self:?} and other = {other:?} does not have the same read"
        );

        let self_offset = self.offset();
        let other_offset = other.offset();
        let new_end = match (self.len(), other.len()) {
            (Some(self_len), Some(other_len)) => {
                Some(min(self_offset + self_len, other_offset + other_len))
            }
            (Some(self_len), None) => Some(self_offset + self_len),
            (None, Some(other_len)) => Some(other_offset + other_len),
            (None, None) => None,
        };

        let new_offset = max(self_offset, other_offset);
        match new_end {
            Some(end) => {
                let final_offset = min(end, new_offset);
                self.set_offset(final_offset);
                self.set_len(end - final_offset);
            }
            None => {
                self.set_offset(new_offset);
            }
        }
    }

    /// Shrink an `RpRange` using the specifed local range. This is useful
    /// in adapter trimming, for example when you find a 3' adapter at some
    /// position `x` from the start of this `RpRange` and want to only retain
    /// the range `[0..x)` within this `RpRange`
    ///
    /// # Input
    /// * `shrink_range` - `Range<usize>` to shrink the `RpRange` to
    ///
    /// # Panics
    /// * If the length is set and the `shrink_range` is not within `0..self.len().unwrap()`
    /// * If shrink_range.start > shrink_range.end
    ///
    ///
    /// # Example
    /// ```rust
    /// use fastq_set::read_pair::{RpRange, WhichRead};
    /// // 100 bases in R1 starting from base 10
    /// let mut rp_range = RpRange::new(WhichRead::R1, 10, Some(100));
    /// let shrink_range = 20..60;
    /// rp_range.shrink(&shrink_range);
    /// assert!(rp_range.read() == WhichRead::R1);
    /// assert!(rp_range.offset() == 30); // 20 + 10
    /// assert!(rp_range.len() == Some(40)); // 60-20
    /// ```
    ///
    /// # Tests
    /// * `test_shrink_invalid_range_1()`: Test for panic if shrink range start is > length
    /// * `test_shrink_invalid_range_2()`: Test for panic if shrink range end is > length
    /// * `test_shrink_invalid_range_3()`: Test for panic if shrink range start > end
    /// * `test_rprange_trivial_shrink()`: Test shrink to an empty range.
    /// * `prop_test_rprange_shrink()`: Test shrink for arbitrary values of
    ///   `RpRange` and valid `shrink_range`
    pub fn shrink(&mut self, shrink_range: &ops::Range<usize>) {
        assert!(
            shrink_range.start <= shrink_range.end,
            "RpRange shrink() expects a valid range with start<=end. Received {shrink_range:?}"
        );

        if let Some(len) = self.len() {
            assert!(
                shrink_range.end <= len,
                "Attempting to shrink more than the current length. shrink_range = {shrink_range:?}, RpRange = {self:?}",
            );
        }

        let new_offset = self.offset() + shrink_range.start;
        let new_len = shrink_range.end - shrink_range.start;
        self.set_offset(new_offset);
        self.set_len(new_len);
    }

    /// Trim the `RpRange` by specifying the `end` and the `amount`
    ///
    /// # Inputs
    /// * `end` - whether to trim the `FivePrime` end or the `ThreePrime` end
    /// * `amount` - how many bases to trim.
    ///
    /// # Panics
    /// * To trim the `ThreePrime` end, the length needs to be known. Panics otherwise
    /// * If length is known, panics if the amount to trim is more than the length
    ///
    /// # Example
    /// ```rust
    /// use fastq_set::read_pair::{RpRange, WhichRead};
    /// use fastq_set::WhichEnd;
    /// let mut rp_range1 = RpRange::new(WhichRead::R1, 40, None);
    /// // Trim 10 bases in the 5' end
    /// rp_range1.trim(WhichEnd::FivePrime, 10);
    /// assert!(rp_range1.read() == WhichRead::R1); // Unchanged
    /// assert!(rp_range1.offset() == 50); // 40 + 10
    /// assert!(rp_range1.len() == None); // Still until the end
    ///
    /// let mut rp_range2 = RpRange::new(WhichRead::R1, 40, Some(110));
    /// // Trim 20 bases in the 3' end
    /// // Trimming from 3' end needs the length to be set. In this
    /// // case, the range defined is [40, 150)
    /// rp_range2.trim(WhichEnd::ThreePrime, 20);
    /// assert!(rp_range2.read() == WhichRead::R1); // Unchanged
    /// assert!(rp_range2.offset() == 40); // Unchanged
    /// assert!(rp_range2.len() == Some(90)); // 110-20
    /// ```
    ///
    /// # Tests
    /// * `test_rprange_trim_without_len()` - Make sure 3' trimming without length panics
    /// * `prop_test_rprange_trim()` - Proptest to make sure the trim results are correct
    pub fn trim(&mut self, end: WhichEnd, amount: usize) {
        if amount == 0 {
            // Trivial case
            return;
        }

        // Panic if we know the length and are asked to trim more than the length
        if let Some(l) = self.len() {
            assert!(
                amount <= l,
                "Attempt to trim more than the length of RpRange"
            );
        }

        match end {
            WhichEnd::ThreePrime => {
                match self.len() {
                    Some(len) => self.set_len(len - amount), // Won't underflow because of the assert above
                    None => {
                        panic!("ThreePrime trim is only possible for RpRange with known length!")
                    }
                }
            }
            WhichEnd::FivePrime => {
                let new_offset = self.offset() + amount;
                self.set_offset(new_offset);
                if let Some(l) = self.len() {
                    self.set_len(l - amount); // Won't underflow because of the assert above
                }
            }
        }
    }
}

/// Storage patterns for a read pair. There are two
/// options which is a compromise between performance
/// and memory usage.
#[derive(Clone, Copy, Default)]
pub enum ReadPairStorage {
    /// Multiple `ReadPair` objects will be backed slices into
    /// the same buffer. This reductes the allocation overhead.
    /// However, a buffer is only dropped when all the read
    /// pairs associated with it is dropped. This could lead to
    /// higher than expected memory usage when you filter the reads
    /// and store only a subset of reads in memory, because the memory
    /// associated with other reads are not necessarily freed.
    SharedBuffer,
    /// Perform an allocation per read, so that each read has its
    /// own memory which will be freed when the read is dropped. This
    /// has a deterministic memory usage particularly when you are
    /// filtering reads. However this will have inferior performance
    /// compared to `ReadPairStorage::SharedBuffer` due to higher
    /// [de]allocation overhead
    #[default]
    PerReadAllocation,
}

/// Helper struct used during construction of a ReadPair. The data for the ReadPair is
/// accumulated in the buffer bytes::BytesMut. When all the data has been added, call
/// `freeze()` to convert this into an immutable `ReadPair` object. Multiple `ReadPair` objects
/// will be backed slices into the same buffer, which reduces allocation overhead.
pub(super) struct MutReadPair<'a> {
    offsets: [ReadOffset; 4],
    data: &'a mut BytesMut,
    storage: ReadPairStorage,
}

impl<'a> MutReadPair<'a> {
    pub(super) fn empty(buffer: &mut BytesMut) -> MutReadPair<'_> {
        let offsets = [ReadOffset::default(); 4];
        MutReadPair {
            offsets,
            data: buffer,
            storage: ReadPairStorage::default(),
        }
    }

    pub fn new<R: Record>(buffer: &'a mut BytesMut, rr: &[Option<R>; 4]) -> MutReadPair<'a> {
        let mut rp = MutReadPair::empty(buffer);
        for (rec, which) in rr.iter().zip(WhichRead::read_types()) {
            if let Some(rec) = rec {
                rp.push_read(rec, which);
            }
            // default ReadOffsets is exists = false
        }
        rp
    }

    pub(super) fn storage(mut self, storage: ReadPairStorage) -> Self {
        self.storage = storage;
        self
    }

    pub(super) fn push_read<R: Record>(&mut self, rec: &R, which: WhichRead) {
        assert!(!self.offsets[which as usize].exists);
        let start = self.data.len() as u16;
        self.data.extend_from_slice(rec.head());
        let head = self.data.len() as u16;
        assert!(
            head - start <= MAX_READ_HEAD,
            "Invalid fastq record with id:{}. The fastq record id must be {} characters or less.",
            std::str::from_utf8(rec.head()).unwrap(),
            MAX_READ_HEAD
        );
        self.data.extend_from_slice(rec.seq());
        let seq = self.data.len() as u16;
        assert!(
            seq - head <= MAX_READ_LEN,
            "Invalid fastq record with id:{}. The fastq record sequence must be {} bases or less.",
            std::str::from_utf8(rec.head()).unwrap(),
            MAX_READ_LEN
        );
        self.data.extend_from_slice(rec.qual());
        let qual = self.data.len() as u16;
        assert!(
            qual - seq == seq - head,
            "Invalid fastq record with id:{}. The fastq record sequence and quality length do not match!",
            std::str::from_utf8(rec.head()).unwrap()
        );
        let read_offset = ReadOffset {
            exists: true,
            start,
            head,
            seq,
            qual,
        };
        self.offsets[which as usize] = read_offset;
    }

    pub fn freeze(self) -> ReadPair {
        ReadPair {
            offsets: self.offsets,
            data: match self.storage {
                ReadPairStorage::SharedBuffer => self.data.split().freeze(),
                ReadPairStorage::PerReadAllocation => {
                    Bytes::from(self.data.split().freeze().to_vec()) // Allocate a vector and then make Bytes
                }
            },
        }
    }

    #[inline]
    /// Get a ReadPart `part` from a read `which` in this cluster
    pub fn get(&self, which: WhichRead, part: ReadPart) -> Option<&[u8]> {
        if self.offsets[which as usize].exists {
            let w = self.offsets[which as usize];
            match part {
                ReadPart::Header => Some(&self.data[w.start as usize..w.head as usize]),
                ReadPart::Seq => Some(&self.data[w.head as usize..w.seq as usize]),
                ReadPart::Qual => Some(&self.data[w.seq as usize..w.qual as usize]),
            }
        } else {
            None
        }
    }
}

/// Container for all read data from a single Illumina cluster. Faithfully represents
/// the FASTQ data from all available reads, if available.
/// Generally should be created by a `ReadPairIter`.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct ReadPair {
    offsets: [ReadOffset; 4],

    // Single vector with all the raw FASTQ data.
    // Use with = "serde_bytes" to get much faster perf
    data: Bytes,
}

impl ReadPair {
    #[inline]
    /// Get a ReadPart `part` from a read `which` in this cluster
    pub fn get(&self, which: WhichRead, part: ReadPart) -> Option<&[u8]> {
        if self.offsets[which as usize].exists {
            let w = self.offsets[which as usize];
            match part {
                ReadPart::Header => Some(&self.data[w.start as usize..w.head as usize]),
                ReadPart::Seq => Some(&self.data[w.head as usize..w.seq as usize]),
                ReadPart::Qual => Some(&self.data[w.seq as usize..w.qual as usize]),
            }
        } else {
            None
        }
    }

    #[inline]
    /// Get the range in `RpRange`, return the chosen `part` (sequence or qvs).
    pub fn get_range(&self, rp_range: RpRange, part: ReadPart) -> Option<&[u8]> {
        let read = self.get(rp_range.read(), part);
        read.and_then(|r| rp_range.slice(r))
    }

    pub fn check_range(&self, range: &RpRange, region_name: &str) -> Result<()> {
        let req_len = range.offset() + range.len().unwrap_or(0);
        let Some(read) = self.get(range.read(), ReadPart::Seq) else {
            bail!(
                "{region_name} is missing from FASTQ. Read {} is not present.",
                range.read()
            );
        };
        ensure!(
            read.len() >= req_len,
            "{region_name} is expected in positions {}-{} in Read {}, but read is {} bp long.",
            range.offset(),
            req_len,
            range.read(),
            read.len()
        );
        Ok(())
    }

    /// Read length of the selected read.
    pub fn len(&self, which: WhichRead) -> Option<usize> {
        self.offsets[which as usize].seq_len()
    }

    pub fn is_empty(&self) -> bool {
        !self
            .offsets
            .iter()
            .any(|o| matches!(o.seq_len(), Some(x) if x != 0))
    }

    /// Write read selected by `which` in FASTQ format to `writer`.
    /// This method will silently do nothing if the selected read doesn't exist.
    pub fn write_fastq<W: Write>(&self, which: WhichRead, writer: &mut W) -> Result<()> {
        if self.offsets[which as usize].exists {
            let head = self.get(which, ReadPart::Header).unwrap();
            writer.write_all(b"@")?;
            writer.write_all(head)?;
            writer.write_all(b"\n")?;

            let seq = self.get(which, ReadPart::Seq).unwrap();
            writer.write_all(seq)?;
            writer.write_all(b"\n+\n")?;

            let qual = self.get(which, ReadPart::Qual).unwrap();
            writer.write_all(qual)?;
            writer.write_all(b"\n")?;
        }

        Ok(())
    }

    /// WARNING: DO NOT USE THIS FUNCTION IF YOU ARE STREAMING FASTQ DATA
    /// This function is intended for testing and illustration purposes
    /// only. Use `ReadPairIter` if you are iterating over a fastq.
    pub fn new<R: Record>(rr: [Option<R>; 4]) -> ReadPair {
        let mut buffer = BytesMut::with_capacity(4096);
        MutReadPair::new(&mut buffer, &rr).freeze()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastq::OwnedRecord;
    use proptest::arbitrary::any;
    use proptest::prelude::Just;
    use proptest::proptest;
    use proptest::strategy::Strategy;
    use std::cmp::{max, min};

    const MAX_RPRANGE_ENTRY: usize = 1usize << 15;

    fn generate_fastq(head_len: usize, seq_len: usize) -> OwnedRecord {
        OwnedRecord {
            head: "H".repeat(head_len).into(),
            seq: b"A".repeat(seq_len),
            sep: None,
            qual: b"I".repeat(seq_len),
        }
    }

    #[test]
    fn test_rprange_slice() {
        let data = &[1u8, 2, 3, 4, 5];
        let r1 = RpRange::new(WhichRead::R1, 0, Some(5));
        assert_eq!(r1.slice(data), Some(&data[..]));
        let r2 = RpRange::new(WhichRead::R1, 0, Some(6));
        assert_eq!(r2.slice(data), None);
        let r3 = RpRange::new(WhichRead::R1, 0, None);
        assert_eq!(r3.slice(data), Some(&data[..]));
        let r4 = RpRange::new(WhichRead::R1, 5, None);
        assert_eq!(r4.slice(data), Some(&[][..]));
        let r5 = RpRange::new(WhichRead::R1, 6, None);
        assert_eq!(r5.slice(data), None);
    }

    #[test]
    #[should_panic]
    fn test_rprange_invalid_offset() {
        // Offset too large
        let _ = RpRange::new(WhichRead::R1, MAX_RPRANGE_ENTRY, None);
    }

    #[test]
    #[should_panic]
    fn test_rprange_invalid_len() {
        // Length too large
        let _ = RpRange::new(WhichRead::R1, 10, Some(MAX_RPRANGE_ENTRY));
    }

    #[test]
    #[should_panic]
    fn test_rprange_trim_without_len() {
        // 3' Trimming is not allowed if length is unknown
        let mut r = RpRange::new(WhichRead::R1, 10, None);
        r.trim(WhichEnd::ThreePrime, 5);
    }

    #[test]
    fn test_rprange_trivial_shrink() {
        // Trimming is not allowed if length is unknows
        let mut r = RpRange::new(WhichRead::R1, 0, Some(1));
        let shrink_range = 0..0;
        r.shrink(&shrink_range);
        assert_eq!(r.read(), WhichRead::R1);
        assert_eq!(r.offset(), 0);
        assert_eq!(r.len(), Some(0));
    }

    #[test]
    #[should_panic]
    fn test_shrink_invalid_range_1() {
        let mut r = RpRange::new(WhichRead::R1, 10, Some(20));
        r.shrink(&(40..50));
    }

    #[test]
    #[should_panic]
    fn test_shrink_invalid_range_2() {
        let mut r = RpRange::new(WhichRead::R1, 10, Some(20));
        r.shrink(&(10..50));
    }

    #[test]
    #[should_panic]
    fn test_shrink_invalid_range_3() {
        let mut r = RpRange::new(WhichRead::R1, 10, Some(20));
        #[allow(clippy::reversed_empty_ranges)]
        r.shrink(&(10..5));
    }

    proptest! {
        #[test]
        fn prop_test_rprange_representation(
            read in (0..4usize).prop_map(WhichRead::from),
            offset in 0..MAX_RPRANGE_ENTRY,
            len in any::<usize>()
        ) {
            // Make sure out internal compact representation is valid for random inputs
            let len = if len < MAX_RPRANGE_ENTRY { Some(len) } else { None };
            let rprange = RpRange::new(read, offset, len);
            assert_eq!(rprange.read(), read);
            assert_eq!(rprange.offset(), offset);
            assert_eq!(rprange.len(), len);
        }
    }

    proptest! {
        #[test]
        fn prop_test_rprange_setter(
            read in (0..4usize).prop_map(WhichRead::from),
            offset in 0..MAX_RPRANGE_ENTRY,
            len in any::<usize>(),
            new_offset in 0..MAX_RPRANGE_ENTRY,
            new_len in 0..MAX_RPRANGE_ENTRY-1
        ) {
            // Make sure we set length and offset correctly using the setter functions
            // for random inputs
            let len = if len < MAX_RPRANGE_ENTRY { Some(len) } else { None };
            let mut rprange = RpRange::new(read, offset, len);
            rprange.set_offset(new_offset);
            rprange.set_len(new_len);
            assert_eq!(rprange.read(), read);
            assert_eq!(rprange.offset(), new_offset);
            assert_eq!(rprange.len(), Some(new_len));
        }
    }

    proptest! {
        #[test]
        fn prop_test_rprange_trim(
            val in any::<u32>(),
            end in any::<bool>(),
            amount in any::<usize>()
        ) {
            // Make sure we trim the range correctly
            let mut rprange = RpRange { val };
            let old_offset = rprange.offset();
            let old_read = rprange.read();
            if let Some(l) = rprange.len() {
                let amount = amount % max(l, 1);
                if amount + old_offset < MAX_RPRANGE_ENTRY {
                    if end {
                        rprange.trim(WhichEnd::ThreePrime, amount);
                        assert_eq!(rprange.offset(), old_offset);
                    } else {
                        rprange.trim(WhichEnd::FivePrime, amount);
                        assert_eq!(rprange.offset(), old_offset + amount);
                    }
                    assert_eq!(rprange.read(), old_read);
                    assert_eq!(rprange.len(), Some(l-amount));
                }
            } else if amount < (MAX_RPRANGE_ENTRY - old_offset) {
                rprange.trim(WhichEnd::FivePrime, amount);
                assert_eq!(rprange.read(), old_read);
                assert_eq!(rprange.offset(), old_offset + amount);
                assert_eq!(rprange.len(), None);
            }

        }
    }

    proptest! {
        #[test]
        fn prop_test_rprange_shrink(
            val in any::<u32>(),
            x in any::<usize>(),
            y in any::<usize>()
        ) {
            // Make sure we shrink the range correctly
            let mut rprange = RpRange { val };

            let max_val = rprange.len().unwrap_or(MAX_RPRANGE_ENTRY);
            let x = min(x, max_val);
            let y = min(y, max_val);
            let shrink_range = min(x, y)..max(x, y);

            let old_offset = rprange.offset();
            let old_read = rprange.read();

            let expected_offset = old_offset + shrink_range.start;
            let expected_len = Some(shrink_range.end - shrink_range.start);
            if expected_offset < MAX_RPRANGE_ENTRY {
                rprange.shrink(&shrink_range);
                assert_eq!(rprange.read(), old_read);
                assert_eq!(rprange.offset(), expected_offset);
                assert_eq!(rprange.len(), expected_len);
            }

        }
    }

    proptest! {
        #[test]
        fn prop_test_readpair_get(
            ref head in proptest::collection::vec(any::<u8>(), 0usize..500usize),
            (seq, qual) in proptest::collection::vec(any::<u8>(), 0usize..1000usize).prop_flat_map(|seq| {
                let len = seq.len();
                (Just(seq), proptest::collection::vec(any::<u8>(), len))
            }),
            pos in 0..4usize
        ) {
            let mut buffer = BytesMut::with_capacity(4096);
            let owned = OwnedRecord {
                head: head.clone(),
                seq: seq.clone(),
                qual: qual.clone(),
                sep: None,
            };
            let mut input = [None, None, None, None];
            input[pos] = Some(owned);
            let read_pair = MutReadPair::new(&mut buffer, &input).freeze();
            let read = WhichRead::from(pos);
            assert_eq!(read_pair.get(read, ReadPart::Header), Some(head.as_slice()));
            assert_eq!(read_pair.get(read, ReadPart::Qual), Some(qual.as_slice()));
            assert_eq!(read_pair.get(read, ReadPart::Seq), Some(seq.as_slice()));
        }
    }

    #[test]
    #[should_panic]
    fn test_rprange_intersect_panic() {
        let mut rp_range = RpRange::new(WhichRead::R1, 0, None);
        rp_range.intersect(RpRange::new(WhichRead::R2, 0, None));
    }

    #[test]
    fn test_rprange_intersect_both_open() {
        let mut rp_range = RpRange::new(WhichRead::R1, 50, None);
        rp_range.intersect(RpRange::new(WhichRead::R1, 10, None));
        assert_eq!(rp_range.read(), WhichRead::R1);
        assert_eq!(rp_range.offset(), 50);
        assert_eq!(rp_range.len(), None);
    }

    #[test]
    fn test_rprange_intersect_self_open() {
        let mut rp_range = RpRange::new(WhichRead::R1, 30, None);
        rp_range.intersect(RpRange::new(WhichRead::R1, 0, Some(100)));
        assert_eq!(rp_range.read(), WhichRead::R1);
        assert_eq!(rp_range.offset(), 30);
        assert_eq!(rp_range.len(), Some(70));
    }

    #[test]
    fn test_rprange_intersect_other_open() {
        let mut rp_range = RpRange::new(WhichRead::R1, 30, Some(120));
        rp_range.intersect(RpRange::new(WhichRead::R1, 65, None));
        assert_eq!(rp_range.read(), WhichRead::R1);
        assert_eq!(rp_range.offset(), 65);
        assert_eq!(rp_range.len(), Some(85));
    }

    #[test]
    fn test_rprange_intersect_both_closed() {
        let mut rp_range = RpRange::new(WhichRead::R1, 40, Some(110));
        rp_range.intersect(RpRange::new(WhichRead::R1, 30, Some(70)));
        assert_eq!(rp_range.read(), WhichRead::R1);
        assert_eq!(rp_range.offset(), 40);
        assert_eq!(rp_range.len(), Some(60));
    }

    #[test]
    fn test_readpair_max_size() {
        let r2 = generate_fastq(MAX_READ_HEAD.into(), MAX_READ_LEN.into());
        let read_pair = ReadPair::new([None, Some(r2), None, None]);

        let expected_r2_offsets = ReadOffset {
            exists: true,
            start: 0,
            head: MAX_READ_HEAD,
            seq: MAX_READ_HEAD + MAX_READ_LEN,
            qual: MAX_READ_HEAD + MAX_READ_LEN * 2,
        };
        assert_eq!(
            read_pair.offsets[WhichRead::R2 as usize],
            expected_r2_offsets
        );
    }

    #[test]
    #[should_panic]
    fn test_readpair_id_toolong() {
        let r2 = generate_fastq((MAX_READ_HEAD + 1) as usize, 100);
        ReadPair::new([None, Some(r2), None, None]);
    }

    #[test]
    #[should_panic]
    fn test_readpair_seq_toolong() {
        let r2 = generate_fastq(100, (MAX_READ_LEN + 1) as usize);
        ReadPair::new([None, Some(r2), None, None]);
    }
}
