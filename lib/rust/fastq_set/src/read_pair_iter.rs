//! Read a set of FASTQs, convert into an Iterator over ReadPairs.
// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

use crate::ERROR_CODE_INFO;
use crate::read_pair::{MutReadPair, ReadHeader, ReadPair, ReadPairStorage, ReadPart, WhichRead};
use anyhow::{Context, Result, anyhow, bail, ensure};
use bytes::{BufMut, BytesMut};
use fastq::{self, OwnedRecord, Record, RecordRefIter};
use itertools::{Itertools, zip_eq};
use martian_derive::MartianType;
use rand::SeedableRng;
use rand::distr::{Distribution, Uniform};
use rand::rngs::SmallRng;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, ErrorKind, Read, Seek};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

const GZ_BUF_SIZE: usize = 1 << 16;

/// Format a FASTQ error message.
fn format_error(message: &str, path: &Path, line: usize) -> anyhow::Error {
    anyhow!("{message}: {path:?} line {line}")
}

trait FileIoError<T> {
    fn open_err(self, path: &Path) -> Result<T>;
    fn fastq_err(self, path: &Path, line: usize) -> Result<T>;
}

impl<T> FileIoError<T> for std::io::Result<T> {
    #[cold]
    fn open_err(self, path: &Path) -> Result<T> {
        self.map_err(|err| {
            anyhow!("TXRNGR10010: Error opening FASTQ file: {path:?}: {err}. {ERROR_CODE_INFO}")
        })
    }

    #[cold]
    fn fastq_err(self, path: &Path, line: usize) -> Result<T> {
        self.map_err(|err| {
            if err.kind() == ErrorKind::InvalidData {
                anyhow!("TXRNGR10012: {err}: {path:?} line {line}. {ERROR_CODE_INFO}")
            } else {
                anyhow!(
                    "TXRNGR10011: IO error in FASTQ file: {path:?} line {line}: {err}. \
                     {ERROR_CODE_INFO}"
                )
            }
        })
    }
}

/// A set of corresponding FASTQ representing the different read components from a set of flowcell 'clusters'
/// All reads are optional except for R1. For an interleaved R1/R2 file, set the filename in the `r1` field,
/// and set `r1_interleaved = true`.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, MartianType)]
pub struct FastqSet {
    pub r1: String,
    pub r2: Option<String>,
    pub i1: Option<String>,
    pub i2: Option<String>,
    pub r1_interleaved: bool,
}

/// Detect whether a gzip file uses BGZF block structure.
///
/// Reads up to 4 gzip members and checks each for the BGZF extra subfield
/// (SI1=0x42, SI2=0x43, SLEN=2). Returns true only if all members inspected
/// carry the BGZF tag; returns false on the first non-BGZF member or if the
/// file has no extra field.
fn is_bgzf(path: &Path) -> std::io::Result<bool> {
    use std::io::SeekFrom;

    let mut file = File::open(path)?;
    const MAX_MEMBERS: usize = 4;

    for i in 0..MAX_MEMBERS {
        let member_start = file.stream_position()?;

        let mut hdr = [0u8; 10];
        match file.read_exact(&mut hdr) {
            Ok(()) => {}
            Err(e) if i > 0 && e.kind() == ErrorKind::UnexpectedEof => return Ok(true),
            Err(e) => return Err(e),
        }
        if hdr[0] != 0x1f || hdr[1] != 0x8b {
            return Ok(false);
        }
        let flags = hdr[3];
        if flags & 4 == 0 {
            // No FEXTRA field → cannot be BGZF
            return Ok(false);
        }
        let mut xlen_buf = [0u8; 2];
        file.read_exact(&mut xlen_buf)?;
        let xlen = u16::from_le_bytes(xlen_buf) as usize;
        let mut extra = vec![0u8; xlen];
        file.read_exact(&mut extra)?;

        // Walk subfields looking for SI1=B(0x42) SI2=C(0x43) SLEN=2
        let mut bsize: Option<usize> = None;
        let mut i = 0;
        while i + 4 <= extra.len() {
            let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            if extra[i] == 0x42 && extra[i + 1] == 0x43 && slen == 2 && i + 6 <= extra.len() {
                bsize = Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize + 1);
                break;
            }
            i += 4 + slen;
        }
        let Some(bs) = bsize else {
            return Ok(false);
        };

        // Skip to the next member using the absolute block size from the member start.
        // BGZF BSIZE covers the entire block including header, deflate data, and trailer.
        file.seek(SeekFrom::Start(member_start + bs as u64))?;
    }
    Ok(true)
}

impl FastqSet {
    /// Returns true if any FASTQ file in this set appears to be BGZF-compressed.
    /// Used by stage split logic to decide whether to request extra decompression threads.
    pub fn is_bgzf(&self) -> Result<bool> {
        for path in [
            Some(&self.r1),
            self.r2.as_ref(),
            self.i1.as_ref(),
            self.i2.as_ref(),
        ]
        .into_iter()
        .flatten()
        {
            let path = Path::new(path);
            if is_bgzf(path).fastq_err(path, 0)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Open an iterator over the reads in this set of fastq files.
    ///
    /// When `allowed_bgzf_threads` is Some(N) and the input files are BGZF-compressed,
    /// a multithreaded bgzf reader is used instead of flate2. Use N
    /// workers for decompression.
    pub fn read_iter(
        &self,
        header_handling: ReadHeader,
        allowed_bgzf_threads: Option<usize>,
    ) -> Result<FastqReadPairIter> {
        FastqReadPairIter::new(
            Some(self.r1.clone().into()),
            self.r2.clone().map(Into::into),
            self.i1.clone().map(Into::into),
            self.i2.clone().map(Into::into),
            self.r1_interleaved,
            header_handling,
            None,
            allowed_bgzf_threads,
        )
    }

    /// Open an iterator over the reads in this set of fastq files.
    /// Use a custom buffer size for parsing reads.
    pub fn read_iter_with_buf_size(
        &self,
        header_handling: ReadHeader,
        buf_size: usize,
        allowed_bgzf_threads: Option<usize>,
    ) -> Result<FastqReadPairIter> {
        FastqReadPairIter::new(
            Some(self.r1.clone().into()),
            self.r2.clone().map(Into::into),
            self.i1.clone().map(Into::into),
            self.i2.clone().map(Into::into),
            self.r1_interleaved,
            header_handling,
            Some(buf_size),
            allowed_bgzf_threads,
        )
    }
}

const BUF_SIZE: usize = 4096 * 4;

/// Track the number of reads read from a FASTQ file.
/// Provides a helper method for getting an accurate line number for a specific
/// part of the read.
#[derive(Default, Clone, Copy)]
pub(crate) struct RecordCount(usize);

impl RecordCount {
    /// Return a line number, starting from 1, pointing to the specified part
    /// of the read.
    ///
    /// Note: since this type stores a count of "the number of records read",
    /// the appropriate line number for the current state will be one less than
    /// the current count - if we've read one read, then any reporting about
    /// the last read read should be for "record 0".
    pub(crate) fn line_number(&self, which: Option<ReadPart>) -> usize {
        let record_index = self.0.checked_sub(1).unwrap_or_default();
        let header_line = (record_index * 4) + 1;
        header_line
            + match which {
                None => 0,
                Some(ReadPart::Seq) => 1,
                Some(ReadPart::Qual) => 3,
            }
    }

    /// Increment the record count.
    fn increment(&mut self) {
        self.0 += 1;
    }
}

struct ReadPairIterSet {
    pub iters: [Option<RecordRefIter<Box<dyn BufRead + Send>>>; 4],
    pub paths: [Option<PathBuf>; 4],
    pub records_read: [RecordCount; 4],
}

/// Read sequencing data from a parallel set of FASTQ files.
/// Illumina sequencers typically emit a parallel set of FASTQ files, with one file
/// for each read component taken by the sequencer. Up to 4 reads are possible (R1, R2, I1, and I2).
/// The reader supports any combination of R1/R2/I1/I2 read files,
/// as well as an interleaved R1/R2 file. Supports plain or gzipped FASTQ files, which
/// will be detected based on the filename extension.
pub struct FastqReadPairIter {
    // Each input file can interleave up to 2 -- declare those here
    r1_interleaved: bool,
    buffer: BytesMut,
    rand: SmallRng,
    uniform: Uniform<f64>,
    subsample_rate: f64,
    storage: ReadPairStorage,
    iters: Box<ReadPairIterSet>,
    is_single_ended: bool,
    /// Keep or discard read IDs; read IDs will still be checked for concurrence.
    header_handling: ReadHeader,
}

impl FastqReadPairIter {
    /// Open a FASTQ file that is uncompressed, gzipped compressed, or lz4 compressed.
    /// The extension of the file is ignored & the filetype is determined by looking
    /// for magic bytes at the of the file.
    ///
    /// If bgzf is detected and allowed_bgzf_threads is Some(N), use N worker
    /// threads for decompression.
    fn open_fastq_from_file(
        p: &Path,
        mut file: File,
        allowed_bgzf_threads: Option<usize>,
    ) -> Result<Box<dyn BufRead + Send>> {
        let mut buf = [0u8; 4];
        file.read_exact(&mut buf[..]).fastq_err(p, 0)?;
        file.rewind().fastq_err(p, 0)?;

        if buf[0..2] == [0x1F, 0x8B] {
            if let Some(n_workers) = allowed_bgzf_threads
                && is_bgzf(p).fastq_err(p, 0)?
            {
                use noodles_bgzf::io::MultithreadedReader;
                let workers = NonZeroUsize::new(n_workers).unwrap();
                let reader = MultithreadedReader::with_worker_count(workers, file);
                Ok(Box::new(BufReader::with_capacity(GZ_BUF_SIZE, reader)))
            } else {
                let gz = flate2::read::MultiGzDecoder::new(file);
                Ok(Box::new(BufReader::with_capacity(GZ_BUF_SIZE, gz)))
            }
        } else if buf[0..4] == [0x04, 0x22, 0x4D, 0x18] {
            let lz = lz4::Decoder::new(file).fastq_err(p, 0)?;
            let buf_reader = BufReader::with_capacity(GZ_BUF_SIZE, lz);
            Ok(Box::new(buf_reader))
        } else if buf[0] == b'@' {
            let buf_reader = BufReader::with_capacity(32 * 1024, file);
            Ok(Box::new(buf_reader))
        } else {
            bail!(
                "FASTQ file does not appear to be valid. Input FASTQ file must be gzip or lz4 \
                 compressed, or must begin with the '@' symbol: {p:?}",
            )
        }
    }

    /// Open a (possibly gzip or lz4 compressed) FASTQ file & read some records to confirm the format looks good.
    fn open_fastq_confirm_fmt(
        p: &Path,
        allowed_bgzf_threads: Option<usize>,
    ) -> Result<Box<dyn BufRead + Send>> {
        let mut file = File::open(p).open_err(p)?;
        // Clone the file descriptor so we can read the file header.
        // The alternative would be to open the file at that path a second time,
        // but that would be less performant (because it wouldn't share i/o
        // cache as effectively, and because it would involve more round trips
        // to the underlying fileystem) and also technically incorrect
        // since the file at that path could have changed.
        //
        // Note: always use allowed_bgzf_threads=None for the validation pass. dup()-based clones
        // share the same open file description (and therefore the same file offset) as the
        // original. If the validation reader spawns bgzf worker threads, they would race
        // against the real reader's workers over the shared offset, producing corrupted reads.
        // The single-threaded gzip impl is safe to use with a dup clone.
        let reader = Self::open_fastq_from_file(p, file.try_clone().open_err(p)?, None)?;
        let parser = fastq::Parser::new(reader);

        // make sure we can successfully read some records
        // try and give a useful message if we can't
        let mut iter = parser.ref_iter();

        for rec in 0..10 {
            iter.advance().fastq_err(p, rec * 4)?;
            let rec = iter.get();
            if rec.is_none() {
                break;
            }
        }

        // re-open file so we re-read the initial records
        file.rewind().fastq_err(p, 0)?;
        Self::open_fastq_from_file(p, file, allowed_bgzf_threads)
    }

    fn is_single_ended<P: AsRef<Path>>(
        r1: &Option<P>,
        r2: &Option<P>,
        r1_interleaved: bool,
    ) -> Result<bool> {
        if !r1_interleaved {
            // We only allow parsing single ended RA files
            Ok(r2.is_none())
        } else if let Some(p) = r1.as_ref().map(P::as_ref) {
            let rdr = Self::open_fastq_confirm_fmt(p, None)?;
            let parser = fastq::Parser::new(rdr);
            let mut r1_iter = parser.ref_iter();
            r1_iter
                .advance()
                .with_context(|| format_error("Reading first read", p, 0))?;
            let first_read_header = r1_iter.get().and_then(|x| {
                x.head()
                    .split(|&x| matches!(x, b' ' | b'/'))
                    .next()
                    .map(Vec::from)
            });
            r1_iter
                .advance()
                .with_context(|| format_error("Reading second read", p, 4))?;
            let second_read_header = r1_iter.get().and_then(|x| {
                x.head()
                    .split(|&x| matches!(x, b' ' | b'/'))
                    .next()
                    .map(Vec::from)
            });
            Ok(first_read_header != second_read_header)
        } else {
            Ok(false)
        }
    }

    /// Open a `ReadPairIter` given of FASTQ files.
    /// For interleaved R1/R2 files, set `r2 = None`, and set
    /// `r1_interleaved = true`.
    #[expect(clippy::too_many_arguments)]
    pub(crate) fn new(
        r1: Option<PathBuf>,
        r2: Option<PathBuf>,
        i1: Option<PathBuf>,
        i2: Option<PathBuf>,
        r1_interleaved: bool,
        header_handling: ReadHeader,
        buffer_size: Option<usize>,
        allowed_bgzf_threads: Option<usize>,
    ) -> Result<FastqReadPairIter> {
        let is_single_ended = Self::is_single_ended(&r1, &r2, r1_interleaved)?;
        let paths = [r1, r2, i1, i2];

        let mut iters = [None, None, None, None];

        for (idx, r) in paths.iter().enumerate() {
            if let Some(p) = r {
                let rdr = Self::open_fastq_confirm_fmt(p, allowed_bgzf_threads)?;
                let parser = if let Some(buf_size) = buffer_size {
                    fastq::Parser::with_buf_size(rdr, buf_size)
                } else {
                    fastq::Parser::new(rdr)
                };
                iters[idx] = Some(parser.ref_iter());
            }
        }

        let buffer = BytesMut::with_capacity(BUF_SIZE);

        Ok(FastqReadPairIter {
            r1_interleaved,
            buffer,
            rand: SmallRng::seed_from_u64(0),
            uniform: Uniform::new(0.0, 1.0).unwrap(),
            subsample_rate: 1.0,
            storage: ReadPairStorage::default(),
            iters: Box::new(ReadPairIterSet {
                paths,
                iters,
                records_read: [RecordCount::default(); 4],
            }),
            is_single_ended,
            header_handling,
        })
    }

    pub fn storage(mut self, storage: ReadPairStorage) -> Self {
        self.storage = storage;
        self
    }

    pub fn seed(mut self, seed: u64) -> Self {
        self.rand = SmallRng::seed_from_u64(seed);
        self
    }

    pub fn get_is_single_ended(&self) -> bool {
        self.is_single_ended
    }

    pub fn subsample_rate(mut self, subsample_rate: f64) -> Self {
        self.subsample_rate = subsample_rate;
        self
    }

    fn get_next(&mut self) -> Result<Option<ReadPair>> {
        // Recycle the buffer if it's almost full.
        if self.buffer.remaining_mut() < 512 {
            self.buffer = BytesMut::with_capacity(BUF_SIZE);
        }
        assert!(self.buffer.is_empty());

        // need these local reference to avoid borrow checker problem
        let iters = &mut *self.iters;
        let paths = &iters.paths;
        let rec_num = &mut iters.records_read;
        let iters = &mut iters.iters;

        let mut rp = MutReadPair::empty(&mut self.buffer).storage(self.storage);

        loop {
            // If we're subsampling, decide whether or not we'll be skipping
            // this record.
            let sample = self.subsample_rate >= 1.0
                || self.uniform.sample(&mut self.rand) < self.subsample_rate;

            // Track which reader was the first to finish.
            let mut iter_ended = [false; 4];

            for (idx, iter_opt) in iters.iter_mut().enumerate() {
                let Some(ref mut iter) = *iter_opt else {
                    continue;
                };
                rec_num[idx].increment();
                iter.advance()
                    .fastq_err(paths[idx].as_ref().unwrap(), rec_num[idx].line_number(None))?;

                let current_read_record = {
                    let record = iter.get();

                    // Check for non-ACGTN characters
                    if let Some(ref rec) = record {
                        ensure!(
                            Record::validate_dnan(rec),
                            format_error(
                                "FASTQ contains sequence base with character other than [ACGTN].",
                                paths[idx].as_ref().unwrap(),
                                rec_num[idx].line_number(Some(ReadPart::Seq)),
                            )
                        );
                        if sample {
                            let which = WhichRead::read_types()[idx];
                            rp.push_read(rec, which)
                                .map_err(|err| err.augment_id_mismatch(rec_num, paths))?;
                        }
                    } else {
                        // track which reader finished
                        iter_ended[idx] = true;
                    }

                    record
                };

                // If R1 is interleaved, read another entry
                // and store it as R2
                if idx == 0 && self.r1_interleaved && !iter_ended[idx] {
                    if !self.is_single_ended {
                        rec_num[idx].increment();
                        iter.advance().fastq_err(
                            paths[idx].as_ref().unwrap(),
                            rec_num[idx].line_number(None),
                        )?;
                        let record = iter.get();

                        // Check for non-ACGTN characters
                        if let Some(rec) = &record {
                            ensure!(
                                Record::validate_dnan(rec),
                                format_error(
                                    "FASTQ contains sequence base with character other than [ACGTN].",
                                    paths[idx].as_ref().unwrap(),
                                    rec_num[idx].line_number(Some(ReadPart::Seq)),
                                )
                            );

                            if sample {
                                let which = WhichRead::read_types()[idx + 1];
                                rp.push_read(rec, which)
                                    .map_err(|err| err.augment_id_mismatch(rec_num, paths))?;
                            }
                        } else {
                            // We should only hit this if the FASTQ has an odd number of records.
                            bail!(format_error(
                                "Input FASTQ file was input as interleaved R1 and R2, \
                                     but contains an odd number of records",
                                paths[idx].as_ref().unwrap(),
                                rec_num[idx].line_number(None),
                            ));
                        }
                    } else if let Some(current_read_record) = current_read_record
                        && sample
                    {
                        let fake_r2_read = OwnedRecord {
                            head: current_read_record.head().to_owned(),
                            sep: None,
                            seq: vec![],
                            qual: vec![],
                        };
                        rp.push_read(&fake_r2_read, WhichRead::R2)
                            .map_err(|err| err.augment_id_mismatch(rec_num, paths))?;
                    }
                }
            }

            // At least one of the readers got to the end -- make sure they all did.
            if let Some((ended_index, _)) = iter_ended.iter().find_position(|&&x| x) {
                // Are there any incomplete iterators?
                let all_complete =
                    zip_eq(iter_ended, paths).all(|(ended, path)| ended || path.is_none());
                ensure!(
                    all_complete,
                    format_error(
                        "Input FASTQ file ended prematurely",
                        paths[ended_index].as_ref().unwrap(),
                        rec_num[ended_index].line_number(None),
                    )
                );
                return Ok(None);
            }

            if sample {
                return Ok(Some(rp.freeze(self.header_handling)));
            }
        }
    }
}

impl Iterator for FastqReadPairIter {
    type Item = Result<ReadPair>;

    /// Iterate over ReadPair objects
    fn next(&mut self) -> Option<Result<ReadPair>> {
        self.get_next().transpose()
    }
}

#[cfg(test)]
mod test_read_pair_iter {
    use super::*;
    use file_diff::diff_files;
    use std::io::Write;

    fn read_pair_iter(
        r1: Option<&str>,
        r2: Option<&str>,
        i1: Option<&str>,
        i2: Option<&str>,
        interleaved: bool,
    ) -> Result<FastqReadPairIter> {
        FastqReadPairIter::new(
            r1.map(PathBuf::from),
            r2.map(PathBuf::from),
            i1.map(PathBuf::from),
            i2.map(PathBuf::from),
            interleaved,
            ReadHeader::KeepId,
            None,
            None,
        )
    }

    // Verify that we can parse and write to the identical FASTQ.
    // Verify that we can parse and write to the identical FASTQ. Note that this
    // test works because we have removed the non-ID read metadata from good-RA.fastq
    #[test]
    fn test_round_trip() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();
        let res: Result<Vec<ReadPair>> = it.collect();
        let res = res.unwrap();

        {
            {
                let mut output = File::create("tests/fastq_round_trip.fastq").unwrap();
                for rec in res {
                    rec.write_fastq(WhichRead::R1, &mut output).unwrap();
                    rec.write_fastq(WhichRead::R2, &mut output).unwrap();
                }
                output.flush().unwrap();
            }

            let mut output = File::open("tests/fastq_round_trip.fastq").unwrap();
            let mut input = File::open("tests/read_pair_iter/good-RA.fastq").unwrap();
            assert!(diff_files(&mut input, &mut output));
        }

        std::fs::remove_file("tests/fastq_round_trip.fastq").unwrap();
    }

    #[test]
    fn test_correct() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        assert!(!it.get_is_single_ended());

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_ok());
        assert_eq!(res.unwrap().len(), 8);
    }

    #[test]
    fn test_single_ended_correct() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA-single-end.fastq"),
            None,
            None,
            None,
            true,
        )
        .unwrap();

        assert!(it.get_is_single_ended());

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_ok());
        println!("{:#?}", &res);
        assert_eq!(res.unwrap().len(), 8);
    }

    #[test]
    fn test_mgi() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/slash1.fastq"),
            Some("tests/read_pair_iter/slash2.fastq"),
            None,
            None,
            false,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_ok());
        assert_eq!(res.unwrap().len(), 6);
    }

    #[test]
    fn test_csi_1376() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/csi-1376-R1.fastq"),
            Some("tests/read_pair_iter/csi-1376-R2.fastq"),
            None,
            None,
            false,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        println!("res: {res:?}");
        assert!(res.is_ok());
        assert_eq!(res.unwrap().len(), 3);
    }

    #[test]
    fn test_not_gzipped() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/not-gzipped-RA.fastq.gz"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_ok());
    }

    #[test]
    fn test_gzipped() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-gzipped-RA.fastq.gz"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_ok());
    }

    #[test]
    fn test_lz4() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-lz4-RA.fastq.lz4"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_ok());
    }

    #[test]
    fn test_missing_pair() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/short-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_missing_single_end() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/imbalanced-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_short_i1() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/short-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_bad_char_i1() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/bad-char-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_short_i2() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/short-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_mismatched_header() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/bad-header-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_mismatched_fastq_error() {
        let it = read_pair_iter(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/bad-format-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        );

        // this error gets caught in the opening pre-check
        assert!(it.is_err());
    }

    #[cfg(target_os = "linux")]
    fn get_rss_pages() -> u64 {
        std::fs::read_to_string("/proc/self/statm")
            .unwrap()
            .split_ascii_whitespace()
            .nth(1)
            .unwrap()
            .parse::<u64>()
            .unwrap()
    }

    #[cfg(target_os = "linux")]
    fn test_mem_single(every: usize, storage: ReadPairStorage) -> Result<u64> {
        let iter = read_pair_iter(
            Some("tests/read_pair_iter/vdj_micro_50k.fastq"),
            None,
            None,
            None,
            true,
        )
        .unwrap()
        .storage(storage)
        .step_by(every);
        let rss_before = get_rss_pages();
        let rp: Vec<ReadPair> = iter.collect::<Result<_, _>>()?;
        let rss_after = get_rss_pages();
        let elements = rp.len();
        drop(rp);
        let rss_after_drop = get_rss_pages();

        // BUG: This test is incorrect.  Process RSS usage will generally not go
        // down after a deallocation, because most malloc implementations
        // will hang on to recently-freed memory for a little while in case
        // it can be reused for subsequent allications.  So rss_used, measured
        // this way, will almost always show 0.
        let rss_used = rss_after.saturating_sub(rss_after_drop);

        let ps = rustix::param::page_size() as u64;
        println!(
            "{:<10} {:<10} {:<12} {:<12} {:<12} {:<12}",
            every,
            elements,
            ps * rss_before / 1024,
            ps * rss_after / 1024,
            ps * rss_used / 1024,
            ps * rss_after_drop / 1024
        );
        Ok(ps * rss_used)
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn test_mem_usage() {
        println!(
            "{:10} {:10} {:12} {:12} {:12} {:12}",
            "Every", "Elements", "Before(kB)", "After(kB)", "Diff(kB)", "Drop(kB)"
        );
        let every = 4;
        let used_rss = test_mem_single(every, ReadPairStorage::PerReadAllocation).unwrap();
        // let used_rss = test_mem_single(every, ReadPairStorage::SharedBuffer).unwrap(); // This fails

        // 50k lines = 6250 reads,
        // With approx 1kB per read the total memory should be ~ 6.25MB / every
        let max_used_rss = 7 * 1024 * 1024 / (every as u64);
        assert!(
            used_rss <= max_used_rss,
            "Used {used_rss} kB whereas max is {max_used_rss}"
        );
    }

    #[test]
    fn test_bgzf_detection_positive() {
        use noodles_bgzf as bgzf;
        use std::io::Write as _;
        use tempfile::NamedTempFile;

        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_path_buf();
        {
            let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
            writer.write_all(b"ACGT").unwrap();
        }
        let _ = tmp; // keep alive so the path isn't deleted
        assert!(is_bgzf(&path).unwrap());
    }

    #[test]
    fn test_bgzf_detection_negative() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write as _;
        use tempfile::NamedTempFile;

        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_path_buf();
        {
            let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
            encoder.write_all(b"ACGT").unwrap();
            encoder.finish().unwrap();
        }
        let _ = tmp; // keep alive so the path isn't deleted
        assert!(!is_bgzf(&path).unwrap());
    }
}
