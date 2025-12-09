//! Martian stage BARCODE_CORRECTION
#![deny(missing_docs)]

use crate::barcode_correction_metrics::{
    BarcodeCorrectionMetricsFormat, BarcodeCorrectionVisitor, BarcodeDiversityMetrics,
};
use crate::barcode_sort::{BarcodeOrder, ReadVisitor};
use crate::stages::make_correction_map::CorrectionMapFormat;
use crate::types::ReadShardFile;
use anyhow::Result;
use barcode::{
    Barcode, BarcodeConstruct, BarcodeCorrector, BarcodeSegment, BarcodeSegmentState, BcSegQual,
    SegmentedBarcode, Segments,
};
use barcode_extensions::select_barcode_corrector;
use cr_bam::constants::{ALN_BC_DISK_CHUNK_SZ, ALN_BC_ITEM_BUFFER_SZ, ALN_BC_SEND_BUFFER_SZ};
use cr_types::chemistry::{BarcodeExtraction, ChemistryDefs};
use cr_types::rna_read::{RnaRead, find_with_one_mismatch};
use cr_types::types::BcSegmentCountFormat;
use cr_types::{
    LibraryType, MergeSortedCounts, MetricsFile, PerLibrarySortedBarcodeCounts,
    SortedBarcodeCountFile, write_sorted_count_file_from_histogram,
};
use fastq_set::WhichRead;
use fastq_set::adapters::X12_CAPTURE_SEQ;
use fastq_set::read_pair::{ReadPair, ReadPart, RpRange};
use itertools::{Itertools, izip};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::{Histogram, JsonReport, Metric, OrderedHistogram, TxHashMap};
use serde::{Deserialize, Serialize};
use shardio::{ShardWriter, UnsortedShardReader};
use stats::effective_diversity;
use std::ops::RangeInclusive;
use std::path::PathBuf;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct BarcodeCorrectionStageInputs {
    pub invalid_uncorrected: Vec<ReadShardFile>,
    pub chemistry_defs: ChemistryDefs,
    pub barcode_segment_counts: BcSegmentCountFormat,
    pub barcode_counts: PerLibrarySortedBarcodeCounts,
    pub valid_read_metrics: BarcodeCorrectionMetricsFormat,
    pub correction_map: Option<CorrectionMapFormat>,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct BarcodeCorrectionChunkInputs {
    index: usize,
    count: usize,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct BarcodeCorrectionChunkOutputs {
    valid_shard: ReadShardFile,
    invalid_shard: ReadShardFile,
    corrected_barcode_counts: PerLibrarySortedBarcodeCounts,
    invalid_barcode_counts: SortedBarcodeCountFile,
    chunk_summary: BarcodeCorrectionMetricsFormat,
}

#[derive(Serialize, MartianStruct)]
pub struct BarcodeCorrectionStageOutputs {
    pub valid_corrected: Vec<ReadShardFile>,
    pub invalid: Vec<ReadShardFile>,
    pub summary: MetricsFile,
    pub corrected_barcode_counts: PerLibrarySortedBarcodeCounts,
    pub invalid_barcode_counts: SortedBarcodeCountFile,
}

pub const MAX_BC_CORRECT_CHUNKS_PER_GG: usize = 100;
pub const MIN_BC_CORRECT_READ_PAIRS_PER_CHUNK: usize = 500_000;

/// Return the position of a perfect match of a prefix of `query` to a suffix of `target`.
/// The match has no minimum length. One matching character is sufficient.
fn find_suffix_overlap_match(target: &[u8], query: &[u8]) -> Option<usize> {
    (1..=query.len()).rev().find_map(|query_len| {
        target
            .ends_with(&query[..query_len])
            .then_some(target.len() - query_len)
    })
}

/// Return whether the molecule is unusual for Flex chemistry.
/// Return true if either the capture sequence is found in an unusual position
/// or the read contains two complete capture sequences
/// or the read contains a partial capture sequence at its end.
fn is_unusual_flex_molecule(
    r1_seq: &[u8],
    probe_barcode_offset_range: &RangeInclusive<usize>,
) -> bool {
    const MINIMUM_SUFFIX_MATCH: usize = 4;
    if find_suffix_overlap_match(r1_seq, X12_CAPTURE_SEQ)
        .is_some_and(|pos| r1_seq.len() - pos >= MINIMUM_SUFFIX_MATCH)
    {
        return true;
    }

    let Some(x12_pos) = find_with_one_mismatch(r1_seq, X12_CAPTURE_SEQ) else {
        return false;
    };

    let barcode_pos = x12_pos + X12_CAPTURE_SEQ.len();
    !probe_barcode_offset_range.contains(&barcode_pos)
        || find_with_one_mismatch(&r1_seq[barcode_pos..], X12_CAPTURE_SEQ).is_some()
}

/// Martian stage BARCODE_CORRECTION
/// Correct sequencing errors in barcodes, up to one mismatch.
pub struct BarcodeCorrection;

/// function to correct barcode in read
pub fn correct_barcode_in_read(
    rna_read: &mut RnaRead,
    extractor: Option<&BarcodeExtraction>,
    corrector_and_length_range: BarcodeConstruct<&(
        BarcodeCorrector,
        Option<std::ops::Range<usize>>,
    )>,
) {
    #[allow(clippy::enum_glob_use)]
    use BarcodeSegmentState::*;

    match extractor {
        None => {
            let bc_qual = rna_read.raw_bc_construct_qual();
            let corrector = corrector_and_length_range.map(|(corrector, _range)| corrector);

            for (segment, qual, corrector) in izip!(
                rna_read.segmented_barcode.segments_mut(),
                bc_qual,
                corrector
            ) {
                if !segment.is_valid() {
                    corrector.correct_barcode(segment, Some(qual));
                }
            }
        }

        Some(&BarcodeExtraction::VariableMultiplexingBarcode {
            min_offset,
            max_offset,
        }) => {
            // We should only ever have GelBeadAndProbe constructs in this case.
            let (corrector_gelbead, corrector_probe) = corrector_and_length_range
                .map(|(corrector, _range)| corrector)
                .gel_bead_and_probe();

            // Correct the gel-bead barcode.
            if !rna_read.segmented_barcode.segments().gel_bead().is_valid() {
                let gelbead_qual = rna_read.raw_bc_construct_qual().gel_bead();
                corrector_gelbead.correct_barcode(
                    rna_read.segmented_barcode.segments_mut().gel_bead(),
                    Some(gelbead_qual),
                );
            }

            // Nothing to do for the probe barcode if it is already valid.
            if rna_read.segmented_barcode.segments().probe().is_valid() {
                return;
            }

            // Correct the probe barcode.
            let default_bc_range = rna_read.bc_range().probe();
            let default_bc_range_offset = default_bc_range.offset() as i64;
            let probe_barcode_offset_range = usize::try_from(default_bc_range_offset + min_offset)
                .unwrap()
                ..=usize::try_from(default_bc_range_offset + max_offset).unwrap();

            if default_bc_range.read() == WhichRead::R1
                && is_unusual_flex_molecule(
                    rna_read.raw_illumina_read1_seq(),
                    &probe_barcode_offset_range,
                )
            {
                // Discard the read because the molecule has unusual structure.
                rna_read.segmented_barcode.segments_mut().probe().state = InvalidStructure;
                return;
            }

            // Correct the probe barcode, taking variable position into account.
            let corrections = probe_barcode_offset_range
                .filter_map(|offset| {
                    let (range, corrected_seg, dist) = try_correct(
                        rna_read.readpair(),
                        RpRange::new(default_bc_range.read(), offset, default_bc_range.len()),
                        corrector_probe,
                    )?;
                    assert_eq!(dist, 1);
                    Some((range, corrected_seg))
                })
                .coalesce(|x, y| {
                    // Merge adjacent positions that correct to the same barcode sequence,
                    // preferring to keep the substitution correction over the indel correction.
                    if x.1.content != y.1.content {
                        return Err((x, y));
                    }
                    Ok(match (x.1.state, y.1.state) {
                        (InvalidAmbiguous, InvalidAmbiguous)
                        | (ValidAfterCorrection, ValidAfterCorrection | InvalidIndel)
                        | (InvalidIndel, InvalidIndel) => x,
                        (InvalidIndel, ValidAfterCorrection) => y,
                        (
                            NotChecked
                            | ValidBeforeCorrection
                            | Invalid
                            | InvalidAmbiguous
                            | InvalidStructure,
                            _,
                        )
                        | (
                            _,
                            NotChecked
                            | ValidBeforeCorrection
                            | Invalid
                            | InvalidAmbiguous
                            | InvalidStructure,
                        ) => {
                            unreachable!("Unexpected barcode state: {x:?} {y:?}")
                        }
                    })
                });

            let (corrected_range, corrected_seg) = match corrections.at_most_one() {
                // No matches
                Ok(None) => return,
                // One match
                Ok(Some(correction)) => correction,
                // Multiple matches
                Err(mut iter) => {
                    // Return the position of the first correction, but flag it as ambiguous.
                    let (corrected_range, corrected_seg) = iter.next().unwrap();
                    (
                        corrected_range,
                        BarcodeSegment {
                            content: corrected_seg.content,
                            state: InvalidAmbiguous,
                        },
                    )
                }
            };

            // Set the position of the raw barcode, whether or not it's a valid correction.
            rna_read.bc_range = BarcodeConstruct::new_gel_bead_and_probe(
                rna_read.bc_range.gel_bead(),
                corrected_range,
            );

            match corrected_seg.state {
                // One match with one substitution
                ValidAfterCorrection => {
                    *rna_read.segmented_barcode.segments_mut().probe() = corrected_seg;
                }
                // One match with one indel or multiple matches
                InvalidIndel | InvalidAmbiguous => {
                    let raw_seq = rna_read
                        .readpair()
                        .get_range(corrected_range, ReadPart::Seq)
                        .unwrap();
                    // Set the barcode sequence to the raw uncorrected sequence at this position.
                    *rna_read.segmented_barcode.segments_mut().probe() =
                        BarcodeSegment::with_sequence(raw_seq, corrected_seg.state);
                }
                // Invalid state
                NotChecked | ValidBeforeCorrection | Invalid | InvalidStructure => {
                    unreachable!("Unexpected barcode state: {:?}", corrected_seg.state)
                }
            }
        }

        Some(BarcodeExtraction::JointBc1Bc2 {
            min_offset,
            max_offset,
        }) => {
            let valid = rna_read.segmented_barcode.segments_valid().array_vec();

            let corrector_and_length = corrector_and_length_range.segments().array_vec();
            let bc_range = rna_read.bc_range();
            let bc_range_vec = bc_range.array_vec();
            let which_read = bc_range_vec[0].read();
            let read_len = rna_read.readpair().len(which_read).unwrap();

            const SEARCH_PADDING: usize = 1;
            let search_range = |length_range: std::ops::Range<usize>| -> std::ops::Range<usize> {
                length_range.start.saturating_sub(SEARCH_PADDING)
                    ..(length_range.end + SEARCH_PADDING).min(read_len)
            };

            assert_eq!(valid.len(), 2);
            match (valid[0], valid[1]) {
                (true, true) => {}
                (true, false) => {
                    // bc1 is valid. Look for bc2 starting from the end of bc1
                    let bc1_end = bc_range_vec[0].offset() + bc_range_vec[0].len().unwrap();
                    let (corrector, lengths) = corrector_and_length[1];

                    let length_range = lengths.clone().unwrap();

                    if let Some((range, segment, _)) = search_range(length_range)
                        .filter_map(|len2| {
                            // See if we can find a corrected bc2 with this length
                            try_correct(
                                rna_read.readpair(),
                                RpRange::new(which_read, bc1_end, Some(len2)),
                                corrector,
                            )
                        })
                        // Among the options pick the one with the least correction
                        .min_by_key(|(_, _, dist)| *dist)
                    {
                        *rna_read.bc_range.as_mut_ref().segments().segment2 = range;
                        *rna_read
                            .segmented_barcode
                            .segments_mut()
                            .segments()
                            .segment2 = segment;
                    };
                }
                (false, true) => {
                    // bc2 is valid. Search for bc1 that ends right before the start of bc2.
                    let bc2_start = bc_range_vec[1].offset();
                    let (corrector, lengths) = corrector_and_length[0];
                    let max_bc1_length = lengths.clone().unwrap().end + SEARCH_PADDING - 1;

                    if let Some((range, segment, _)) = (*min_offset..=*max_offset)
                        .filter_map(|offset| {
                            try_correct(
                                rna_read.readpair(),
                                RpRange::new(
                                    which_read,
                                    offset,
                                    Some((bc2_start - offset).min(max_bc1_length)),
                                ),
                                corrector,
                            )
                        })
                        .min_by_key(|(_, _, dist)| *dist)
                    {
                        *rna_read.bc_range.as_mut_ref().segments().segment1 = range;
                        *rna_read
                            .segmented_barcode
                            .segments_mut()
                            .segments()
                            .segment1 = segment;
                    }
                }
                (false, false) => {
                    const LARGE_DISTANCE: u16 = 1_000;
                    let corrector = corrector_and_length_range.map(|(corrector, _range)| corrector);
                    if let Some(corrected) = (*min_offset..=*max_offset)
                        .flat_map(|offset| {
                            // For the given bc1 start position (offset), find the possible ranges for bc1 and bc2
                            // that could give us an approximate match. To find this, consider all combinations of
                            // the lengths of both bc1 and bc2 with a padding of 1.
                            search_range(corrector_and_length[0].1.clone().unwrap())
                                .cartesian_product(search_range(
                                    corrector_and_length[1].1.clone().unwrap(),
                                ))
                                .filter_map(move |(len1, len2)| {
                                    // Skip cases where the end of bc2 goes past the end of read
                                    (offset + len1 + len2 <= read_len).then_some(
                                        BarcodeConstruct::Segmented(Segments::from_iter([
                                            RpRange::new(which_read, offset, Some(len1)),
                                            RpRange::new(which_read, offset + len1, Some(len2)),
                                        ])),
                                    )
                                })
                        })
                        .map(|bc_range| {
                            bc_range.zip(corrector).map(|(range, segment_corrector)| {
                                try_correct(rna_read.readpair(), range, segment_corrector)
                            })
                        })
                        .min_by_key(|value| {
                            value
                                .iter()
                                .map(|v| v.as_ref().map_or(LARGE_DISTANCE, |(_, _, dist)| *dist))
                                .sum::<u16>()
                        })
                    {
                        let bc_range = bc_range
                            .zip(corrected)
                            .map(|(default, new)| new.map_or(default, |x| x.0));
                        let barcode = SegmentedBarcode::new(
                            rna_read.segmented_barcode.gem_group(),
                            rna_read
                                .segmented_barcode
                                .segments()
                                .zip(corrected)
                                .map(|(default, new)| new.map_or(default, |x| x.1)),
                        );
                        rna_read.bc_range = bc_range;
                        rna_read.segmented_barcode = barcode;
                    }
                }
            }
        }
    }
}

fn try_correct(
    read: &ReadPair,
    range: RpRange,
    corrector: &BarcodeCorrector,
) -> Option<(RpRange, BarcodeSegment, u16)> {
    read.get_range(range, ReadPart::Seq).and_then(|seq| {
        let bc_qual =
            BcSegQual::from_bytes_unchecked(read.get_range(range, ReadPart::Qual).unwrap());
        let mut bc = BarcodeSegment::with_sequence(seq, BarcodeSegmentState::Invalid);
        corrector
            .correct_barcode(&mut bc, Some(bc_qual))
            .map(|dist| (range, bc, dist))
    })
}

#[make_mro(volatile = strict)]
impl MartianStage for BarcodeCorrection {
    type StageInputs = BarcodeCorrectionStageInputs;
    type StageOutputs = BarcodeCorrectionStageOutputs;
    type ChunkInputs = BarcodeCorrectionChunkInputs;
    type ChunkOutputs = BarcodeCorrectionChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // Figure out how many chunks to create.
        let n = UnsortedShardReader::<RnaRead, BarcodeOrder>::len(&args.invalid_uncorrected)?;
        let num_chunks =
            (n / MIN_BC_CORRECT_READ_PAIRS_PER_CHUNK).clamp(1, MAX_BC_CORRECT_CHUNKS_PER_GG);
        let chunk_size = n.div_ceil(num_chunks);
        Ok((0..num_chunks)
            .map(|n| {
                (
                    BarcodeCorrectionChunkInputs {
                        index: n,
                        count: chunk_size,
                    },
                    Resource::with_mem_gb(5),
                )
            })
            .collect::<StageDef<_>>()
            .join_resource(Resource::with_mem_gb(4)))
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let valid_shard: ReadShardFile = rover.make_path("valid_shard");
        let invalid_shard: ReadShardFile = rover.make_path("invalid_shard");

        let mut reader =
            UnsortedShardReader::<RnaRead, BarcodeOrder>::open_set(&args.invalid_uncorrected);

        // Seek to the start of this chunk.
        reader.skip_lazy(chunk_args.index * chunk_args.count)?;

        let mut valid_writer: ShardWriter<RnaRead, BarcodeOrder> = ShardWriter::with_compressor(
            &valid_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            ALN_BC_ITEM_BUFFER_SZ,
            shardio::Compressor::Lz4,
        )?;
        let mut valid_sender = valid_writer.get_sender();

        let mut invalid_writer: ShardWriter<RnaRead, BarcodeOrder> = ShardWriter::with_compressor(
            &invalid_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            ALN_BC_ITEM_BUFFER_SZ,
            shardio::Compressor::Lz4,
        )?;
        let mut invalid_sender = invalid_writer.get_sender();

        let mut bc_counts_corrected: TxHashMap<LibraryType, OrderedHistogram<Barcode>> =
            Default::default();
        let mut bc_counts_invalid = OrderedHistogram::default();

        let bc_segment_counts = args.barcode_segment_counts.read()?;

        let mut per_lib_correction_map = args
            .correction_map
            .map(|f| f.read())
            .transpose()?
            .unwrap_or_default();

        let extractors_and_correctors: TxHashMap<_, _> = bc_segment_counts
            .into_iter()
            .map(|(library_type, counts)| {
                let chemistry_def = &args.chemistry_defs[&library_type];
                let extractor = chemistry_def.barcode_extraction();
                let whitelist_spec = chemistry_def.barcode_whitelist_spec();
                let whitelist = chemistry_def.barcode_whitelist()?;

                anyhow::Ok((
                    library_type,
                    (
                        extractor,
                        select_barcode_corrector(
                            whitelist_spec.zip(whitelist).zip(counts),
                            extractor,
                            per_lib_correction_map.remove(&library_type),
                        ),
                    ),
                ))
            })
            .try_collect()?;

        let mut visitor = BarcodeCorrectionVisitor::new();
        for rna_read in reader.take(chunk_args.count) {
            let mut rna_read = rna_read?;
            let (extractor, corrector) = &extractors_and_correctors[&rna_read.library_type];
            correct_barcode_in_read(&mut rna_read, *extractor, corrector.as_ref());

            visitor.visit_processed_read(&mut rna_read)?;
            let bc = rna_read.barcode();
            if bc.is_valid() {
                bc_counts_corrected
                    .entry(rna_read.library_type)
                    .or_default()
                    .observe(&bc);
                valid_sender.send(rna_read)?;
            } else {
                bc_counts_invalid.observe_owned(bc);
                invalid_sender.send(rna_read)?;
            }
        }
        valid_sender.finished()?;
        valid_writer.finish()?;
        invalid_sender.finished()?;
        invalid_writer.finish()?;

        let chunk_summary: BarcodeCorrectionMetricsFormat = rover.make_path("chunk_summary");
        chunk_summary.write(&visitor.metrics)?;

        let corrected_barcode_counts =
            PerLibrarySortedBarcodeCounts::write_histograms(&bc_counts_corrected, |lib_type| {
                rover.make_path(format!("{lib_type}_sorted_barcode_counts"))
            })?;

        let invalid_barcode_counts = write_sorted_count_file_from_histogram(
            &bc_counts_invalid,
            &rover.make_path::<PathBuf>("invalid_barcode_counts"),
        )?;

        Ok(BarcodeCorrectionChunkOutputs {
            valid_shard,
            invalid_shard,
            corrected_barcode_counts,
            invalid_barcode_counts,
            chunk_summary,
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        // Merge the original valid and corrected barcode counts.
        let corrected_barcode_counts = PerLibrarySortedBarcodeCounts::merge(
            chunk_outs
                .iter()
                .map(|c| c.corrected_barcode_counts.clone())
                .chain(std::iter::once(args.barcode_counts)),
            |lib_type| rover.make_path(format!("{lib_type}_corrected_barcode_counts")),
        )?;

        let mut metrics = args.valid_read_metrics.read()?;
        for chunk_out in &chunk_outs {
            metrics.merge(chunk_out.chunk_summary.read()?);
        }

        // VDJ prefix is lowercase for these metrics.
        let metrics: TxHashMap<_, _> = metrics
            .0
            .into_iter()
            .map(|(lib_type, v)| {
                if lib_type.is_vdj() {
                    (Some("vdj"), v)
                } else {
                    (lib_type.as_metric_prefix_static(), v)
                }
            })
            .collect();

        // FIXME: streaming this collection from disk twice isn't ideal, but its
        // awkward to try to avoid.
        let barcode_diversity_metrics: TxHashMap<_, _> = corrected_barcode_counts
            .iter_each_lib()?
            .into_iter()
            .map(|(lib_type, bc_counts)| {
                let mut barcodes_detected: usize = 0;

                // Count the number of elements in the iterator while computing
                // effective diversity.
                let effective_barcode_diversity = bc_counts.process_results(|iter| {
                    effective_diversity(iter.map(|(_bc, count)| {
                        barcodes_detected += 1;
                        count as f64
                    }))
                })?;

                let key = if lib_type.is_vdj() {
                    None
                } else {
                    lib_type.as_metric_prefix_static()
                };
                anyhow::Ok((
                    key,
                    BarcodeDiversityMetrics {
                        barcodes_detected,
                        effective_barcode_diversity,
                    },
                ))
            })
            .try_collect()?;

        let summary = MetricsFile::from_reporter(
            &rover,
            "summary",
            &(metrics.to_json_reporter() + barcode_diversity_metrics.to_json_reporter()),
        )?;

        let invalid_barcode_counts = merge_invalid_barcode_counts(&chunk_outs, &rover)?;

        let (valid_corrected, invalid) = chunk_outs
            .into_iter()
            .map(|x| (x.valid_shard, x.invalid_shard))
            .unzip();
        Ok(BarcodeCorrectionStageOutputs {
            valid_corrected,
            invalid,
            summary,
            corrected_barcode_counts,
            invalid_barcode_counts,
        })
    }
}

/// Merge all invalid barcode counts from all chunks.
///
/// The barcode counts must have been written in sorted order for this merge
/// to work properly.
fn merge_invalid_barcode_counts(
    chunks: &[BarcodeCorrectionChunkOutputs],
    rover: &MartianRover,
) -> Result<SortedBarcodeCountFile> {
    let f: SortedBarcodeCountFile = rover.make_path("invalid_barcode_counts");
    let mut writer = f.lazy_writer()?;
    // Merge all invalid barcode counts.
    // These were written in sorted order, so we can group their merged iterators
    // without extra allocations.
    for result in MergeSortedCounts::from_readers(
        chunks
            .iter()
            .map(|c| c.invalid_barcode_counts.lazy_reader())
            .try_collect()?,
    ) {
        writer.write_item(&result?)?;
    }
    writer.finish()?;
    Ok(f)
}
