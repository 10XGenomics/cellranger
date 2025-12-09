#![deny(missing_docs)]

use anyhow::Result;
use barcode::{Barcode, BarcodeConstruct, BarcodeConstructMetric, Whitelist};
use cr_bam::constants::{ALN_BC_DISK_CHUNK_SZ, ALN_BC_ITEM_BUFFER_SZ};
use cr_types::mempool::{MemPool, PoolMember};
use cr_types::rna_read::{RnaChunk, RnaProcessor, RnaRead};
use fastq_set::read_pair::{ReadPair, ReadPairStorage};
use fastq_set::{FastqProcessor, ProcessResult};
use metric::{CountMetric, Metric};
use shardio::{BufHandler, Compressor, SortAndWriteHandler, SortKey};
use std::borrow::Cow;
use std::path::PathBuf;
use std::sync::mpsc::{Receiver, Sender, channel};

/// The maximum number of reads we will pre-buffer in the fastq reader.
const BG_FASTA_READAHEAD: usize = 1024;

/// The largest size of each buffer sent to a shard writer.
/// shardio splits the requested item buffer size into two ping-pong buffers,
/// so we match the previous implementation's buffer size by dividing by 2.
const MAX_READ_BUF_SIZE: usize = ALN_BC_ITEM_BUFFER_SZ / 2;

/// Use a much smaller buffer for invalid reads, to ensure we are promptly
/// writing them into shards and releasing shared backing read memory.
/// This will result in larger overhead when reading the invalid shards, but
/// this should not be an issue for the vast majority of analyses.
/// The ratio of the invalid and valid read buffers will be dynamically adjusted
/// to balance high-quality sharding vs. the need to promptly release shared
/// memory backing the reads.
const INITIAL_INVALID_READ_BUF_SIZE: usize = MAX_READ_BUF_SIZE / 16;

pub trait ReadVisitor {
    type ReadType;
    fn visit_processed_read(&mut self, read: &mut Self::ReadType) -> Result<()>;
    fn visit_unprocessed_read(&mut self, read: ReadPair, reason: String) -> Result<()>;
}

pub fn execute_barcode_sort_with_visitor<V: ReadVisitor<ReadType = RnaRead>>(
    chunk: RnaChunk,
    valid_path: PathBuf,
    invalid_path: PathBuf,
    whitelist: BarcodeConstruct<Whitelist>,
    max_iters: Option<usize>,
    visitor: &mut V,
) -> Result<BarcodeSortMetrics> {
    let mut metrics = BarcodeSortMetrics::default();
    let mut nreads: usize = 0;
    let (send_err, recv_err) = channel();

    // Spawn shard writer threads.
    std::thread::scope(|s| {
        let (send_valid, recv_valid) = channel();
        let (send_invalid, recv_invalid) = channel();

        let mut read_buffers = ReadBuffers::new(send_valid, send_invalid);

        s.spawn(|| {
            if let Err(err) = process_read_buffers(valid_path, recv_valid) {
                let _ = send_err.send(err);
            }
        });
        s.spawn(|| {
            if let Err(err) = process_read_buffers(invalid_path, recv_invalid) {
                let _ = send_err.send(err);
            }
        });

        let processor = RnaProcessor::new(chunk, whitelist);

        for read_result in processor
            .iter_background_with_storage(BG_FASTA_READAHEAD, ReadPairStorage::SharedBuffer)
            .unwrap()
            .take(max_iters.unwrap_or(usize::MAX))
        {
            nreads += 1;

            match read_result? {
                ProcessResult::Processed(mut read) => {
                    visitor.visit_processed_read(&mut read)?;
                    metrics.valid_items_in.merge(
                        read.segmented_barcode
                            .segments_valid()
                            .map(CountMetric::from)
                            .into(),
                    );

                    if (if read.barcode_is_valid() {
                        metrics.valid_items += 1;
                        read_buffers.put_valid(read)
                    } else {
                        metrics.invalid_items += 1;
                        read_buffers.put_invalid(read)
                    })
                    .is_err()
                    {
                        // one of the processors hung up, stop sending and wait
                        // for error
                        break;
                    }
                }
                ProcessResult::Unprocessed { read, reason } => {
                    visitor.visit_unprocessed_read(read, reason)?;
                }
            }
        }
        anyhow::Ok(())
    })?;

    println!("Processed {nreads} reads");

    // Check for processor error.
    if let Ok(err) = recv_err.try_recv() {
        return Err(err);
    }

    Ok(metrics)
}

/// Sort by barcode.
pub struct BarcodeOrder;
impl SortKey<RnaRead> for BarcodeOrder {
    type Key = Barcode;
    fn sort_key(v: &RnaRead) -> Cow<'_, Barcode> {
        // Note: we do not use shardio's default sorting mechanism for
        // MAKE_SHARD because sorting RnaReads using merge sort results in too
        // much memcpy. Instead we use a cached key sort.
        Cow::Owned(v.barcode())
    }
}

#[derive(Copy, Clone, Default)]
pub struct BarcodeSortMetrics {
    pub valid_items: i64,
    pub valid_items_in: BarcodeConstructMetric<CountMetric>,
    pub invalid_items: i64,
}

struct ReadBuffer {
    pool: MemPool<Vec<RnaRead>>,
    buf: Option<PoolMember<Vec<RnaRead>>>,
    /// The maximum number of items to buffer before sending.
    ///
    /// This can be changed dynamically, but buffer allocations will only
    /// grow. Capacity will not be reclaimed if shrunk.
    pub size: usize,
    send: Sender<PoolMember<Vec<RnaRead>>>,
}

impl ReadBuffer {
    pub fn new(size: usize, send: Sender<PoolMember<Vec<RnaRead>>>) -> Self {
        Self {
            // Use two read buffers, to allow pre-buffering reads while the other
            // buffer is being processed by the writer.
            pool: MemPool::new(2, || Vec::with_capacity(size)),
            send,
            buf: None,
            size,
        }
    }

    /// Put a read in this buffer.
    /// Send the buffer if full.
    /// Fetch a new buffer if necessary.
    ///
    /// Return Err if a send operation failed due to the processor hanging up.
    pub fn put(&mut self, read: RnaRead) -> Result<()> {
        let mut buf = self.buf.take().unwrap_or_else(|| {
            let mut buf = self.pool.fetch();
            buf.clear();
            buf
        });
        buf.push(read);
        if buf.len() >= self.size {
            self.send.send(buf)?;
        } else {
            self.buf = Some(buf);
        }
        Ok(())
    }
}

impl Drop for ReadBuffer {
    /// Send any remaining reads when dropped.
    fn drop(&mut self) {
        if let Some(buf) = self.buf.take() {
            self.send.send(buf).unwrap();
        }
    }
}

fn process_read_buffers(path: PathBuf, recv: Receiver<PoolMember<Vec<RnaRead>>>) -> Result<()> {
    let mut handler: SortAndWriteHandler<RnaRead, BarcodeOrder> =
        SortAndWriteHandler::new(ALN_BC_DISK_CHUNK_SZ, Compressor::Lz4, &path)?;
    for mut buf in &recv {
        // RnaRead structs are quite large - sort_by_cached_key generates a vec
        // of just the keys+indices and sorts that, then performs the final
        // re-arrangement of the read data itself in an O(N) operation. This
        // results in drastically lower CPU time wasted in memcpy. This also
        // avoids recomputing the barcode for each comparision operation, as in
        // the default shardio implementation.
        buf.sort_by_cached_key(RnaRead::barcode);
        handler.process_buf(&mut buf)?;
    }
    handler.write_index()?;
    Ok(())
}

/// Own the read buffers and balance their sizes based on read quality.
struct ReadBuffers {
    valid: ReadBuffer,
    invalid: ReadBuffer,
    frac: ValidFrac,
}

impl ReadBuffers {
    pub fn new(
        send_valid: Sender<PoolMember<Vec<RnaRead>>>,
        send_invalid: Sender<PoolMember<Vec<RnaRead>>>,
    ) -> Self {
        Self {
            valid: ReadBuffer::new(MAX_READ_BUF_SIZE, send_valid),
            invalid: ReadBuffer::new(INITIAL_INVALID_READ_BUF_SIZE, send_invalid),
            frac: ValidFrac::new(MAX_READ_BUF_SIZE / 4),
        }
    }

    fn update_buf_sizes(&mut self, valid: bool) {
        let new_buffer_sizing = self.frac.observe(valid);
        if new_buffer_sizing {
            let (valid_buf_scale, invalid_buf_scale) = self.frac.valid_to_invalid_read_buf_ratio();
            self.valid.size = MAX_READ_BUF_SIZE / valid_buf_scale;
            self.invalid.size = MAX_READ_BUF_SIZE / invalid_buf_scale;
            println!(
                "Adjusting buffer ratios: valid: {valid_buf_scale}, invalid: {invalid_buf_scale}.  Valid rate: {:.2}",
                self.frac.prev_valid as f64 / self.frac.chunk_size as f64
            );
        }
    }

    pub fn put_valid(&mut self, read: RnaRead) -> Result<()> {
        self.update_buf_sizes(true);
        self.valid.put(read)
    }

    pub fn put_invalid(&mut self, read: RnaRead) -> Result<()> {
        self.update_buf_sizes(false);
        self.invalid.put(read)
    }
}

/// Track a crude chunked average of the fraction of valid reads.
struct ValidFrac {
    valid: usize,
    total: usize,
    /// The number of reads to observe before updating the valid frac.
    chunk_size: usize,
    /// The number of reads in the last chunk that were valid.
    ///
    /// Initialized as the chunk size, so we start out assuming that all reads
    /// are valid.
    prev_valid: usize,
}

impl ValidFrac {
    fn new(chunk_size: usize) -> Self {
        Self {
            valid: 0,
            total: 0,
            chunk_size,
            prev_valid: chunk_size,
        }
    }

    /// Observe a count.
    ///
    /// Return true if we filled a chunk and updated our estimate of the rate.
    fn observe(&mut self, valid: bool) -> bool {
        self.total += 1;
        if valid {
            self.valid += 1;
        }
        // We've observed an entire chunk; save the valid count and reset.
        if self.total == self.chunk_size {
            self.prev_valid = self.valid;
            self.total = 0;
            self.valid = 0;
            return true;
        }
        false
    }

    /// Return a suggested ratio for the size of the valid vs invalid read
    /// buffer.
    ///
    /// Return (valid_buffer_reduction_factor, invalid_buffer_reduction_factor).
    ///
    /// This should allow us to promptly release shared memory for the case
    /// where we have few invalid reads, but avoid creating many small, poorly-
    /// collated invalid read shards for the case where we have many invalid
    /// reads (and vice-versa for very poor sequencing data).
    ///
    /// We constrain these buffer scaling factors to be
    ///  - no larger 1
    ///  - no smaller than 1/16
    ///  - scaled by powers of 2
    ///
    /// These values were chosen fairly arbitrarily and have not been empirically tuned.
    fn valid_to_invalid_read_buf_ratio(&self) -> (usize, usize) {
        let prev_invalid = self.chunk_size - self.prev_valid;

        let invalid_ratio = if self.prev_valid == 0 {
            f64::MAX
        } else {
            prev_invalid as f64 / self.prev_valid as f64
        };

        // Use the largest buffer size for whichever track has the larger
        // fraction of reads; shrink the other track's buffers proportionate to
        // that.
        if invalid_ratio < 1.0 {
            // More valid than invalid reads.
            (1, buffer_scaling_factor(invalid_ratio))
        } else {
            // More invalid than valid reads.
            (buffer_scaling_factor(1.0 / invalid_ratio), 1)
        }
    }
}

fn buffer_scaling_factor(ratio: f64) -> usize {
    // We need to take the inverse/swap the scaling if we have more than half of
    // whatever read category this is.
    assert!(ratio <= 1.0);
    if ratio > 1. / 2. {
        1
    } else if ratio > 1. / 4. {
        2
    } else if ratio > 1. / 8. {
        4
    } else if ratio > 1. / 16. {
        8
    } else {
        16
    }
}
