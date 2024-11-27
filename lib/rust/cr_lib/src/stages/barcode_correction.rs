//! Martian stage BARCODE_CORRECTION

use crate::barcode_correction_metrics::{
    BarcodeCorrectionMetricsFormat, BarcodeCorrectionVisitor, BarcodeDiversityMetrics,
};
use crate::barcode_sort::{BarcodeOrder, ReadVisitor};
use crate::stages::make_correction_map::CorrectionMapFormat;
use crate::types::ReadShardFile;
use anyhow::Result;
use barcode::{
    Barcode, BarcodeConstruct, BarcodeCorrector, BarcodeSegment, BarcodeSegmentState, BcSegQual,
    SegmentedBarcode, Segments, Whitelist,
};
use barcode_extensions::select_barcode_corrector;
use cr_bam::constants::{ALN_BC_DISK_CHUNK_SZ, ALN_BC_ITEM_BUFFER_SZ, ALN_BC_SEND_BUFFER_SZ};
use cr_types::chemistry::{BarcodeExtraction, ChemistryDefs};
use cr_types::rna_read::RnaRead;
use cr_types::types::{
    BcCountDataType, BcCountFormat, BcSegmentCountFormat, GemWell, TotalBcCountFormat,
};
use cr_types::MetricsFile;
use fastq_set::read_pair::{ReadPair, ReadPart, RpRange};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{self, JsonReport, JsonReporter, Metric, SimpleHistogram, TxHashMap};
use serde::{Deserialize, Serialize};
use shardio::{ShardWriter, UnsortedShardReader};

#[derive(Deserialize, Clone, MartianStruct)]
pub struct BarcodeCorrectionStageInputs {
    pub gem_well: GemWell,
    pub invalid_uncorrected: Vec<ReadShardFile>,
    pub chemistry_defs: ChemistryDefs,
    pub barcode_segment_counts: BcSegmentCountFormat,
    pub barcode_counts: BcCountFormat,
    pub valid_read_metrics: BarcodeCorrectionMetricsFormat,
    pub correction_map: Option<CorrectionMapFormat>,
    // Counts for all barcodes over this amount will be
    // be reported in the total_barcode_counts output.
    pub min_reads_to_report_bc: i64,
    pub barcodes_under_tissue: Option<JsonFile<Vec<String>>>,
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
    corrected_barcode_counts: BcCountFormat,
    total_barcode_counts: TotalBcCountFormat,
    chunk_summary: BarcodeCorrectionMetricsFormat,
}

#[derive(Serialize, MartianStruct)]
pub struct BarcodeCorrectionStageOutputs {
    pub valid_corrected: Vec<ReadShardFile>,
    pub invalid: Vec<ReadShardFile>,
    pub summary: MetricsFile,
    pub corrected_barcode_counts: BcCountFormat,
    pub total_barcode_counts: TotalBcCountFormat,
}

pub const MAX_BC_CORRECT_CHUNKS_PER_GG: usize = 100;
pub const MIN_BC_CORRECT_READ_PAIRS_PER_CHUNK: usize = 500_000;

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
    // This is needed to get around the borrow checker
    let bc_qual = rna_read.raw_bc_construct_qual();
    let corrector = corrector_and_length_range.map(|(corrector, _range)| corrector);

    match extractor {
        Some(BarcodeExtraction::Independent) | None => {
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
        Some(BarcodeExtraction::JointBc1Bc2 {
            min_offset,
            max_offset,
        }) => {
            fn try_correct(
                read: &ReadPair,
                range: RpRange,
                corrector: &BarcodeCorrector,
                bc_qual: BcSegQual,
            ) -> Option<(RpRange, BarcodeSegment, u16)> {
                read.get_range(range, ReadPart::Seq).and_then(|seq| {
                    let mut bc = BarcodeSegment::with_sequence(seq, BarcodeSegmentState::Invalid);
                    corrector
                        .correct_barcode(&mut bc, Some(bc_qual))
                        .map(|dist| (range, bc, dist))
                })
            }
            let valid = rna_read.segmented_barcode.segments_valid().array_vec();

            let corrector_and_length = corrector_and_length_range.segments().array_vec();
            let bc_range = rna_read.bc_range();
            let bc_range_vec = bc_range.array_vec();
            let bc_qual_vec = bc_qual.array_vec();
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
                                bc_qual_vec[1],
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
                                bc_qual_vec[0],
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
                            bc_range.zip(bc_qual).zip(corrector).map(
                                |((range, qual), segment_corrector)| {
                                    try_correct(rna_read.readpair(), range, segment_corrector, qual)
                                },
                            )
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
            .join_resource(Resource::with_mem_gb(6)))
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

        let mut valid_writer: ShardWriter<RnaRead, BarcodeOrder> = ShardWriter::new(
            &valid_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            ALN_BC_ITEM_BUFFER_SZ,
        )?;
        let mut valid_sender = valid_writer.get_sender();

        let mut invalid_writer: ShardWriter<RnaRead, BarcodeOrder> = ShardWriter::new(
            &invalid_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            ALN_BC_ITEM_BUFFER_SZ,
        )?;
        let mut invalid_sender = invalid_writer.get_sender();

        let mut bc_counts_corrected: BcCountDataType = TxHashMap::default();
        let mut bc_counts_total: SimpleHistogram<Barcode> = SimpleHistogram::default();

        let bc_segment_counts = args.barcode_segment_counts.read()?;

        let mut per_lib_correction_map = args
            .correction_map
            .map(|f| f.read())
            .transpose()?
            .unwrap_or_default();

        let extractors_and_correctors: TxHashMap<_, _> = bc_segment_counts
            .into_iter()
            .map(|(library_type, v)| {
                let chemistry_def = &args.chemistry_defs[&library_type];
                let extractor = chemistry_def.barcode_extraction();
                let whitelist = Whitelist::construct(chemistry_def.barcode_whitelist())?;

                anyhow::Ok((
                    library_type,
                    (
                        extractor,
                        select_barcode_corrector(
                            whitelist.zip(v),
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
                invalid_sender.send(rna_read)?;
            }
            bc_counts_total.observe_owned(bc);
        }
        valid_sender.finished()?;
        valid_writer.finish()?;
        invalid_sender.finished()?;
        invalid_writer.finish()?;

        let chunk_summary: BarcodeCorrectionMetricsFormat = rover.make_path("chunk_summary");
        chunk_summary.write(&visitor.metrics)?;

        let corrected_barcode_counts_file: BcCountFormat =
            rover.make_path("corrected_barcode_counts");
        corrected_barcode_counts_file.write(&bc_counts_corrected)?;

        let total_barcode_counts_file: TotalBcCountFormat = rover.make_path("total_barcode_counts");

        bc_counts_total.retain(|_, v| args.min_reads_to_report_bc <= v.count());
        total_barcode_counts_file.write(&bc_counts_total)?;

        Ok(BarcodeCorrectionChunkOutputs {
            valid_shard,
            invalid_shard,
            corrected_barcode_counts: corrected_barcode_counts_file,
            total_barcode_counts: total_barcode_counts_file,
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
        let bc_counts_raw = args.barcode_counts.read()?;

        let bc_counts_total = {
            let mut result = SimpleHistogram::default();
            for co in &chunk_outs {
                result += co.total_barcode_counts.read()?;
            }
            for counts in bc_counts_raw.values() {
                let mut retained = SimpleHistogram::default();
                counts
                    .distribution()
                    .iter()
                    .filter(|(_, v)| args.min_reads_to_report_bc <= v.count())
                    .for_each(|(k, v)| retained.observe_by(k, v.count()));
                result += retained;
            }
            result
        };

        let total_barcode_counts: TotalBcCountFormat = rover.make_path("total_barcode_counts");
        total_barcode_counts.write(&bc_counts_total)?;

        let mut bc_counts_corrected = bc_counts_raw;
        for co in &chunk_outs {
            bc_counts_corrected.merge(co.corrected_barcode_counts.read()?);
        }

        let corrected_barcode_counts: BcCountFormat = rover.make_path("correct_barcode_counts");
        corrected_barcode_counts.write(&bc_counts_corrected)?;

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
        let barcode_diversity_metrics: TxHashMap<_, _> = bc_counts_corrected
            .iter()
            .map(|(lib_type, bc_counts)| {
                (
                    if lib_type.is_vdj() {
                        None
                    } else {
                        lib_type.as_metric_prefix_static()
                    },
                    BarcodeDiversityMetrics {
                        barcodes_detected: bc_counts.distribution().len(),
                        effective_barcode_diversity: bc_counts.effective_diversity(),
                    },
                )
            })
            .collect();

        let total_barcodes_detected: usize = bc_counts_corrected
            .into_values()
            .flatten()
            .map(|(k, _)| k)
            .chain(
                args.barcodes_under_tissue
                    .as_ref()
                    .map(martian_filetypes::FileTypeRead::read)
                    .transpose()?
                    .unwrap_or_default()
                    .into_iter()
                    .map(|x| x.parse().unwrap()),
            )
            .unique()
            .count();

        let summary = MetricsFile::from_reporter(
            &rover,
            "summary",
            &(metrics.to_json_reporter()
                + barcode_diversity_metrics.to_json_reporter()
                + JsonReporter::from(("total_barcodes_detected", total_barcodes_detected))),
        )?;

        let (valid_corrected, invalid) = chunk_outs
            .into_iter()
            .map(|x| (x.valid_shard, x.invalid_shard))
            .unzip();
        Ok(BarcodeCorrectionStageOutputs {
            valid_corrected,
            invalid,
            summary,
            corrected_barcode_counts,
            total_barcode_counts,
        })
    }
}
