#![deny(missing_docs)]

use crate::make_shard_metrics::{MakeShardMetricProc, MakeShardMetricsAccum};
use anyhow::Result;
use barcode::{Barcode, BarcodeConstructMetric};
use cr_bam::constants::{ALN_BC_DISK_CHUNK_SZ, ALN_BC_ITEM_BUFFER_SZ};
use cr_types::chemistry::ChemistryDef;
use cr_types::rna_read::{RnaChunk, RnaProcessor, RnaRead};
use fastq_set::length_trim_iter::trim_read_lengths;
use fastq_set::long_read_adapter_trimmer::{InsertDefinition, LongReadTrimmer, TrimmedReadVisitor};
use fastq_set::read_pair::{ReadHeader, ReadPair, ReadPairStorage};
use fastq_set::{BackgroundIteratorScopedExt, ProcessResult};
use itertools::{Either, Itertools, process_results};
use json_report_derive::JsonReport;
use metric::{CountMetric, Histogram, Metric, PercentMetric, SimpleHistogram};
use metric_derive::Metric;
use pariter::IteratorExt;
use serde::{Deserialize, Serialize};
use shardio::{Compressor, SortKey, process_sorted_bufs};
use std::borrow::Cow;
use std::path::PathBuf;
use std::sync::mpsc::{Receiver, SyncSender, channel, sync_channel};
use std::sync::{Arc, Mutex};

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

/// How many worker threads to dedicate to computing read metrics and doing
/// serialization/compression? Note: there's no intrinsic reason these must be
/// the same, it just happens to be the case that the workload per read is
/// roughly equal between these two tasks. Serialization/compression is probably
/// a bit heavier.
const NUM_WORKERS: usize = 2;

/// Number of downstream processor workers to use if we are doing bgzf decompression.
/// This was tuned to ensure that the read producer is never stalled.
const NUM_WORKERS_BGZF: usize = NUM_WORKERS * 3;

/// Use this many threads with bgzf decompression itself. This gives about 2.7 GB/s
/// throughput into fastq-rs's ~2 GB/s throughput, ensuring that decompression
/// is not the bottleneck.
const NUM_DECOMPRESSION_WORKERS_BGZF: usize = 4;

pub(super) struct BarcodeSortResults {
    pub metrics: BarcodeSortMetrics,
    pub valid_bc_counts: SimpleHistogram<Barcode>,
    pub accum: MakeShardMetricsAccum,
}

/// Number of elements to accumulate in a buffer before sending on down the pipeline.
/// The main purpose of batching here is to reduce inter-thread communication overhead;
/// the value just needs to be large enough so that we suppress the number of inter-
/// thread sync events, as well as limit allocator pressure by not making too
/// many allocations.
const BATCH_SIZE: usize = 1024;

#[expect(clippy::too_many_arguments)]
pub fn execute_barcode_sort_with_visitor(
    chunk: RnaChunk,
    chem: ChemistryDef,
    valid_path: PathBuf,
    invalid_path: PathBuf,
    max_iters: Option<usize>,
    proc: &MakeShardMetricProc,
    header_handling: ReadHeader,
    use_multithread_bgzf: bool,
) -> Result<BarcodeSortResults> {
    let mut valid_bc_counts = SimpleHistogram::default();
    let mut metrics = BarcodeSortMetrics::default();

    let mut nreads: usize = 0;
    let (send_err, recv_err) = channel();

    let processor = &RnaProcessor::new(chunk.clone(), chem.clone())?;

    let metric_shard_sink: Arc<Mutex<Vec<MakeShardMetricsAccum>>> = Arc::default();
    // This is cloned once per worker in this state, so each worker gets an
    // initialized copy of the metrics accumulator.
    let mut metric_shard = MetricShard {
        state: Some(MakeShardMetricsAccum::new(
            proc.feature_count,
            chunk.library_id(),
            &chem,
        )),
        sink: metric_shard_sink.clone(),
    };

    let mut read_trimming_metrics = ReadTrimmingMetrics::default();

    // Use a scope to manage all worker threads.
    std::thread::scope(|s| {
        let (send_valid, recv_valid) = sync_channel(1);
        let (send_invalid, recv_invalid) = sync_channel(1);

        let mut read_buffers = ReadBuffers::new(send_valid, send_invalid);

        let shard_workers = if use_multithread_bgzf {
            NUM_WORKERS_BGZF
        } else {
            NUM_WORKERS
        };
        let send_err = &send_err;
        s.spawn(move || {
            if let Err(err) = process_read_buffers(valid_path, recv_valid, shard_workers) {
                let _ = send_err.send(err);
            }
        });
        s.spawn(move || {
            if let Err(err) = process_read_buffers(invalid_path, recv_invalid, shard_workers) {
                let _ = send_err.send(err);
            }
        });

        let iter = chunk
            .read_iter(
                header_handling,
                use_multithread_bgzf.then_some(NUM_DECOMPRESSION_WORKERS_BGZF),
            )?
            .storage(ReadPairStorage::SharedBuffer)
            .background_scoped(s, BATCH_SIZE)
            .flatten();

        let iter = match chunk.adapter_trimmer_config.as_ref() {
            None => Either::Left(iter),
            Some(config) => Either::Right(pseudo_paired_from_adapters_iter(
                iter,
                config.clone().into_trimmer(),
                header_handling,
                &mut read_trimming_metrics,
            )),
        };

        // Apply read-length trimming as the final step so it works uniformly for
        // all read sources, including pseudo-paired reads from the adapter trimmer.
        let trim = chunk.trim_lengths();
        let iter = if trim.is_noop() {
            Either::Left(iter)
        } else {
            Either::Right(trim_read_lengths(iter, trim))
        };

        let iter = iter.take(max_iters.unwrap_or(usize::MAX));

        process_results(iter, |reads| {
            let batched = reads.batching(|i| {
                let mut batch = Vec::with_capacity(BATCH_SIZE);
                batch.extend(i.take(BATCH_SIZE));
                nreads += batch.len();
                (!batch.is_empty()).then_some(batch)
            });

            let iter = batched.parallel_map_scoped_custom(
                s,
                |o| {
                    o.threads(if use_multithread_bgzf {
                        NUM_WORKERS_BGZF
                    } else {
                        NUM_WORKERS
                    })
                },
                move |batch| {
                    let metric_accum = metric_shard.state.as_mut().unwrap();
                    batch
                        .into_iter()
                        .filter_map(|r| match processor.process_read(r) {
                            ProcessResult::Processed(read) => {
                                proc.visit_processed_read(&read, metric_accum);
                                Some(read)
                            }
                            ProcessResult::Unprocessed { .. } => {
                                metric_accum.metrics.unprocessed_read_pairs += 1;
                                None
                            }
                        })
                        .collect::<Vec<_>>()
                },
            );

            for proc_read in iter.flatten() {
                metrics.valid_items_in.merge(
                    proc_read
                        .segmented_barcode
                        .segments_valid()
                        .map(CountMetric::from)
                        .into(),
                );
                if (if proc_read.barcode_is_valid() {
                    valid_bc_counts.observe_owned(proc_read.barcode());
                    metrics.valid_items += 1;
                    read_buffers.put_valid(proc_read)
                } else {
                    metrics.invalid_items += 1;
                    read_buffers.put_invalid(proc_read)
                })
                .is_err()
                {
                    // one of the processors hung up, stop sending and wait
                    // for error
                    break;
                }
            }
        })?;

        anyhow::Ok(())
    })?;

    // Check for processor error.
    drop(send_err);
    if let Some(err) = recv_err.into_iter().next() {
        return Err(err);
    }

    println!("Processed {nreads} reads");

    // Fold up all the metric accumulators. Note that this must be outside of
    // the scope to ensure that all scoped workers are guaranteed to have joined,
    // thus why we can get the Vec out of the Arc below. We could refactor the
    // scope to use a regular & reference instead of Arc to statically prevent
    // this from being used incorrectly.
    let mut accum = Arc::try_unwrap(metric_shard_sink)
        .expect("all worker threads should have terminated")
        .into_inner()
        .unwrap()
        .into_iter()
        .reduce(|mut x, y| {
            x.merge(y);
            x
        })
        .unwrap();
    if read_trimming_metrics.total_raw_reads.count() > 0 {
        accum.metrics.read_trimming_metrics = Some(read_trimming_metrics);
    }

    Ok(BarcodeSortResults {
        metrics,
        valid_bc_counts,
        accum,
    })
}

/// We need one metrics accumulator per worker - but pariter doesn't
/// have something like fold, so we need a way to get this back out
/// of the worker thread at the end of processing. We make use of drop
/// semantics so that each metric accumulator hands itself back to
/// the original calling context when processing is complete.
#[derive(Clone)]
struct MetricShard {
    // Needs to be Option so we can take out the value on Drop.
    state: Option<MakeShardMetricsAccum>,
    sink: Arc<Mutex<Vec<MakeShardMetricsAccum>>>,
}
impl Drop for MetricShard {
    fn drop(&mut self) {
        let Some(s) = self.state.take() else {
            return;
        };
        let Ok(mut guard) = self.sink.lock() else {
            // If the mutex is poisoned, someone panicked, we're going to crash.
            // Just skip this and drop the metrics.
            return;
        };
        guard.push(s);
    }
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

/// These metrics don't escape the stage, but are used for deriving
#[derive(Copy, Clone, Default)]
pub(super) struct BarcodeSortMetrics {
    pub valid_items: i64,
    pub valid_items_in: BarcodeConstructMetric<CountMetric>,
    pub invalid_items: i64,
}

struct ReadBuffer {
    buf: Vec<RnaRead>,
    /// The maximum number of items to buffer before sending.
    ///
    /// This can be changed dynamically.
    pub size: usize,
    send: SyncSender<Vec<RnaRead>>,
}

impl ReadBuffer {
    pub fn new(size: usize, send: SyncSender<Vec<RnaRead>>) -> Self {
        Self {
            send,
            buf: Vec::with_capacity(size),
            size,
        }
    }

    /// Put a read in this buffer.
    /// Send the buffer if full.
    /// Fetch a new buffer if necessary.
    ///
    /// Return Err if a send operation failed due to the processor hanging up.
    pub fn put(&mut self, read: RnaRead) -> Result<()> {
        self.buf.push(read);
        if self.buf.len() >= self.size {
            let buf = std::mem::replace(&mut self.buf, Vec::with_capacity(self.size));
            self.send.send(buf)?;
        }
        Ok(())
    }
}

impl Drop for ReadBuffer {
    /// Send any remaining reads when dropped.
    fn drop(&mut self) {
        if !self.buf.is_empty() {
            self.send.send(std::mem::take(&mut self.buf)).unwrap();
        }
    }
}

fn process_read_buffers(
    path: PathBuf,
    recv: Receiver<Vec<RnaRead>>,
    n_workers: usize,
) -> Result<()> {
    let sorted_buf_iter = recv.into_iter().map(|mut buf| {
        // RnaRead structs are quite large - sort_by_cached_key generates a vec
        // of just the keys+indices and sorts that, then performs the final
        // re-arrangement of the read data itself in an O(N) operation. This
        // results in drastically lower CPU time wasted in memcpy. This also
        // avoids recomputing the barcode for each comparision operation, as in
        // the default shardio implementation.
        buf.sort_by_cached_key(RnaRead::barcode);
        buf
    });

    process_sorted_bufs::<RnaRead, BarcodeOrder>(
        sorted_buf_iter,
        ALN_BC_DISK_CHUNK_SZ,
        n_workers,
        Compressor::Lz4,
        &path,
    )
}

/// Own the read buffers and balance their sizes based on read quality.
struct ReadBuffers {
    valid: ReadBuffer,
    invalid: ReadBuffer,
    frac: ValidFrac,
}

impl ReadBuffers {
    pub fn new(
        send_valid: SyncSender<Vec<RnaRead>>,
        send_invalid: SyncSender<Vec<RnaRead>>,
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

/// Number of threads in the read pre-processor pool.
// TODO: tune this value
const N_ADAPTER_TRIMMER_THREADS: usize = 1;

/// Maximum number of reads to buffer during read pre-processing.
// TODO: tune this value
const ADAPTER_TRIMMER_BUFFER_SIZE: usize = 32;

pub fn pseudo_paired_from_adapters_iter<'a, I: Iterator<Item = Result<ReadPair>> + 'a>(
    iter: I,
    trimmer: LongReadTrimmer,
    header_handling: ReadHeader,
    metrics: &'a mut ReadTrimmingMetrics,
) -> impl Iterator<Item = Result<ReadPair>> + 'a {
    iter.parallel_map_custom(
        |o| {
            o.threads(N_ADAPTER_TRIMMER_THREADS)
                .buffer_size(ADAPTER_TRIMMER_BUFFER_SIZE)
        },
        move |read_result| {
            let read: ReadPair = read_result?;
            let mut per_read_metrics = PerReadMetrics::default();
            let trimmed = trimmer.trim_read(read, header_handling, &mut per_read_metrics)?;
            Ok((trimmed, per_read_metrics))
        },
    )
    .map_ok(|(trimmed_read, per_read_metrics)| {
        metrics.update(per_read_metrics);
        trimmed_read
    })
    .filter_map_ok(|opt_read| opt_read)
}

/// Classification of a read after adapter-trimming boundary detection.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TrimmedReadCategory {
    /// Read contains a canonical insert with both terminal adapters.
    InsertCanonical,
    /// Read contains an insert whose final terminal adapter is missing.
    InsertTruncated,
    /// Read matches a truncated insert, but truncation is disallowed.
    RejectedTruncated,
    /// Read does not match canonical or truncated insert topology.
    Unrecognized,
}

/// Add necessary per read metrics here.
#[derive(Default)]
pub struct PerReadMetrics {
    insert_category: Option<TrimmedReadCategory>,
    insert_size: Option<usize>,
}

impl TrimmedReadVisitor for PerReadMetrics {
    fn visit_trimmed_read(&mut self, config: Option<&InsertDefinition>, is_truncated: bool) {
        self.insert_category = Some(match (&config, is_truncated) {
            (Some(_), false) => TrimmedReadCategory::InsertCanonical,
            (Some(_), true) => TrimmedReadCategory::InsertTruncated,
            (None, true) => TrimmedReadCategory::RejectedTruncated,
            (None, false) => TrimmedReadCategory::Unrecognized,
        });
        self.insert_size = config.map(|c| c.get_retain_range().len());
    }
}

#[derive(Debug, Clone, Default, Deserialize, Serialize, Metric, JsonReport)]
#[json_report(extend = "extra_reports")]
pub struct ReadTrimmingMetrics {
    total_raw_reads: CountMetric,
    reads_with_inserts: CountMetric,
    insert_full_length: PercentMetric,
    insert_truncated: PercentMetric,
    rejected_truncated: CountMetric,
    unrecognized: CountMetric,
    insert_size_distribution: SimpleHistogram<i32>,
}

impl ReadTrimmingMetrics {
    /// Computed fields added to the JSON reporter: median insert size.
    fn extra_reports(&self) -> metric::JsonReporter {
        let mut reporter = metric::JsonReporter::default();
        if let Some(median) = self.insert_size_median() {
            reporter.insert("insert_size_median", median);
        }
        reporter
    }

    /// Returns the median insert size in bases, or `None` if no insert sizes were observed.
    pub fn insert_size_median(&self) -> Option<i64> {
        Histogram::percentiles(&self.insert_size_distribution, [50.0_f64])
            .and_then(|p| p.into_iter().next())
            .map(|v| v as i64)
    }

    pub fn update(&mut self, other: PerReadMetrics) {
        self.total_raw_reads.increment();
        if let Some(category) = other.insert_category {
            if category == TrimmedReadCategory::InsertCanonical
                || category == TrimmedReadCategory::InsertTruncated
            {
                self.reads_with_inserts.increment();
                self.insert_full_length
                    .increment(category == TrimmedReadCategory::InsertCanonical);
                self.insert_truncated
                    .increment(category == TrimmedReadCategory::InsertTruncated);
            }
            if category == TrimmedReadCategory::RejectedTruncated {
                self.rejected_truncated.increment();
            }
            if category == TrimmedReadCategory::Unrecognized {
                self.unrecognized.increment();
            }
        }
        if let Some(insert_size) = other.insert_size {
            Histogram::observe(&mut self.insert_size_distribution, &(insert_size as i32));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::chemistry::adapter::layouts;
    use fastq_set::OwnedRecord;
    use fastq_set::long_read_adapter_trimmer::{
        DEFAULT_MAX_R2_LENGTH, ReadSplittingStrategy, TrimmerConfig,
    };
    use fastq_set::read_pair::{ReadHeader, ReadPair};

    const RTL_GEX: TrimmerConfig = TrimmerConfig {
        layout: layouts::FLEX_GEX,
        splitting_strategy: ReadSplittingStrategy::AdapterLocated {},
        allow_truncated_constructs: false,
        is_directional: true,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    const UNIVERSAL_GEX: TrimmerConfig = TrimmerConfig {
        layout: layouts::THREE_PRIME_GEX,
        splitting_strategy: ReadSplittingStrategy::FixedOffset { offset: 4 },
        allow_truncated_constructs: true,
        is_directional: true,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    fn clean_read(read: &[u8]) -> Vec<u8> {
        read.iter().filter(|&&b| b != b'|').copied().collect()
    }

    fn make_read_pair(seq: Vec<u8>) -> ReadPair {
        let qual = vec![b'I'; seq.len()];
        let record = OwnedRecord {
            head: b"test".to_vec(),
            seq,
            qual,
            sep: None,
        };
        ReadPair::new([Some(record), None, None, None], ReadHeader::KeepId)
    }

    fn check_per_read_metrics(
        config: &TrimmerConfig,
        seq: &[u8],
        expected_category: Option<TrimmedReadCategory>,
        expected_insert_size: Option<usize>,
    ) {
        let trimmer = config.clone().into_trimmer();
        let rp = make_read_pair(clean_read(seq));
        let mut metrics = PerReadMetrics::default();
        trimmer
            .trim_read(rp, ReadHeader::KeepId, &mut metrics)
            .unwrap();
        assert_eq!(metrics.insert_category, expected_category);
        assert_eq!(metrics.insert_size, expected_insert_size);
    }

    #[test]
    fn test_rtl_gex_canonical_full_length() {
        // (r1 22 bases)(GGGG)(x12 12 bases)(spacer 13 bases)(ATGCA)(small_rna_r2_rc 12 bases)
        check_per_read_metrics(
            &RTL_GEX,
            b"CTACACGACGCTCTTCCGATCT|GGGG|TTGCTAGGACCG|TTTCTTATATGGG|ATGCA|TGGAATTCTCGG",
            Some(TrimmedReadCategory::InsertCanonical),
            Some(34),
        );
    }

    #[test]
    fn test_rtl_gex_truncated_rejected() {
        // (r1 22 bases)(GGGG)(x12 12 bases)(spacer 13 bases)(ATGCA) — missing terminal 3p adapter
        check_per_read_metrics(
            &RTL_GEX,
            b"CTACACGACGCTCTTCCGATCT|GGGG|TTGCTAGGACCG|TTTCTTATATGGG|ATGCA",
            Some(TrimmedReadCategory::RejectedTruncated),
            None,
        );
    }

    #[test]
    fn test_rtl_gex_missing_x12_unrecognized() {
        // (r1 22 bases)(GGGG)(spacer 13 bases)(ATGCA)(small_rna_r2_rc 12 bases) — missing x12 adapter
        check_per_read_metrics(
            &RTL_GEX,
            b"CTACACGACGCTCTTCCGATCT|GGGG|TTTCTTATATGGG|ATGCA|TGGAATTCTCGG",
            Some(TrimmedReadCategory::Unrecognized),
            None,
        );
    }

    #[test]
    fn test_universal_gex_canonical_full_length() {
        // (r1_bp12 12 bases)(4 bases)(poly_t 20 bases)(4 bases)(rt_primer_rc 28 bases)
        check_per_read_metrics(
            &UNIVERSAL_GEX,
            b"CTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|GATC|ATGTACTCTGCGTTGATACCACTGCTT",
            Some(TrimmedReadCategory::InsertCanonical),
            Some(28),
        );
    }

    #[test]
    fn test_universal_gex_truncated() {
        // (r1_bp12 12 bases)(4 bases)(poly_t 20 bases)(4 bases) — missing terminal 3p adapter, truncated allowed
        check_per_read_metrics(
            &UNIVERSAL_GEX,
            b"CTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|GATC",
            Some(TrimmedReadCategory::InsertTruncated),
            Some(28),
        );
    }

    #[test]
    fn test_universal_gex_unrecognized() {
        // No 5' terminal adapter, so no canonical/truncated topology can be formed.
        check_per_read_metrics(
            &UNIVERSAL_GEX,
            b"TTCC|TTTTTTTTTTTTTTTTTTTT|GATC|ATGTACTCTGCGTTGATACCACTGCTT",
            Some(TrimmedReadCategory::Unrecognized),
            None,
        );
    }
}
