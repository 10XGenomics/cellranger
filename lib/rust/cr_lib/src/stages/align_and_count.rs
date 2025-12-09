//! Martian stage ALIGN_AND_COUNT
#![deny(missing_docs)]

use super::write_minimap_index::MinimapIndex;
use crate::BcUmiInfoShardFile;
use crate::align_and_count_metrics::StageVisitor;
use crate::align_metrics::{BarcodeMetrics, LibFeatThenBarcodeOrder};
use crate::aligner::{Align, Aligner, AlignmentOptions, BarcodeSummary, MAX_ANNOTATIONS_IN_MEM};
use crate::barcode_sort::BarcodeOrder;
use crate::minimap2::MinimapReference;
use crate::parquet_file::{ParquetFile, ParquetWriter, PerReadGapAlignRow};
#[cfg(feature = "tenx_internal")]
use crate::stages::internal::V1PatternFixParams;
#[cfg(feature = "tenx_source_available")]
use crate::stages::stubs::V1PatternFixParams;
use crate::types::{
    BarcodeMetricsShardFile, FeatureReferenceFormat, ReadShardFile, ReadSpillFormat,
};
use anyhow::{Error, Result, bail, ensure};
use barcode::Barcode;
use cr_bam::bam::{BamPosSort, is_unmapped_to_any};
use cr_bam::constants::{
    ALN_BC_DISK_CHUNK_SZ, ALN_BC_GIB, ALN_BC_ITEM_BUFFER_SZ, ALN_BC_SEND_BUFFER_SZ,
};
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::probe_set::{
    MappedGap, MappedGapAlignmentInfo, ProbeSetReference, ProbeSetReferenceMetadata,
};
use cr_types::reference::feature_checker::compute_feature_dist;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::rna_read::{RnaChunk, RnaRead};
use cr_types::spill_vec::SpillVec;
use cr_types::types::{BarcodeSetFormat, BcUmiInfo, FeatureBarcodeCount, ProbeBarcodeCount};
use cr_types::{
    AlignShardFile, AlignerParam, BarcodeThenFeatureOrder, CountShardFile, FeatureCountFormat,
    PerLibrarySortedBarcodeCounts, SamHeaderFile, SortedBarcodeCountFile,
};
use fastq_set::WhichEnd;
use itertools::{Itertools, zip_eq};
use log::warn;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::TxHashSet;
use minimap2::Preset;
use orbit::{StarReference, StarSettings};
use par_proc::{MAX_ITEMS_IN_MEM, Proc, group_by_processor};
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::{HeaderView, Record};
use serde::{Deserialize, Serialize};
use shardio::{Range, ShardReader, ShardSender, ShardWriter};
use std::borrow::Borrow;
use std::cmp::{Reverse, max};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use transcriptome::Bed12Format;
use tx_annotation::read::{
    AnnotationFiles, AnnotationInfo, ReadAnnotations, ReadAnnotationsFormat, ReadAnnotator,
    RecordAnnotation,
};
use tx_annotation::visitor::AnnotatedReadVisitor;

/// Default value for the stage parameter transcriptome_min_score, STAR parameter --outFilterScoreMin.
pub const DEFAULT_TRANSCRIPTOME_MIN_SCORE: usize = 30;

const RNA_READ_SZ: usize = 1024;
const RNA_READ_LONG_SZ: usize = 1024 * 32;
const ANN_SZ: usize = RNA_READ_SZ * 3; // RnaRead + BamRecords
const ANN_LONG_SZ: usize = RNA_READ_LONG_SZ * 3; // RnaRead + BamRecords
const BAM_BUFFER_SZ: usize = 500_000;
const NUM_CHUNK_THREADS: usize = 4;

/// maximum numbers of reads to sample for VISIT_ANNOTATED_READS_PD
/// we don't need to look at all reads for those metrics.
pub const MAX_READ_ANN_SAMPLE: usize = 50_000_000;

martian_filetype!(BarcodeSummaryFile, "bsf");

/// ALIGN_AND_COUNT stage metrics
#[derive(Serialize)]
pub struct AlignAndCountMetrics {
    /// The aligner used to align reads to the transcriptome reference.
    alignment_aligner: AlignerParam,
    /// The minimum MAPQ threshold for a confidently-mapped read.
    alignment_high_conf_mapq: i64,
}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct StageInputs {
    pub read_chunks: Vec<RnaChunk>,
    pub reference_path: Option<PathBuf>,

    pub read_shards: ReadShards,

    pub feature_counts: FeatureCountFormat,
    pub feature_reference_binary: FeatureReferenceFormat,

    /// The target panel CSV file.
    pub target_set: Option<TargetSetFile>,

    /// The minimap genome index.
    pub minimap_index: Option<MinimapIndex>,
    pub tx_bed: Option<Bed12Format>,

    pub chemistry_defs: ChemistryDefs,

    pub aligner: Option<AlignerParam>,
    pub include_exons: bool,
    pub include_introns: bool,
    pub is_pd: bool,
    pub no_bam: bool,

    /// Minimum alignment score to align to the transcriptome.
    /// Set the --outFilterScoreMin parameter of STAR.
    /// Alignment will be output only if its score is higher than or equal to this value.
    pub transcriptome_min_score: Option<usize>,

    /// Minimum alignment score to trim polyA sequence.
    /// Set to None to disable polyA trimming.
    pub trim_polya_min_score: Option<usize>,

    /// Minimum alignment score to trim template switch oligo (TSO) sequence.
    /// Set to None to disable TSO trimming.
    pub trim_tso_min_score: Option<usize>,

    /// Hidden parameter: Command line parameters of STAR
    pub star_parameters: Option<String>,

    /// All barcode counts that were valid after correction.
    pub corrected_barcode_counts: PerLibrarySortedBarcodeCounts,
    /// All barcode counts that were uncorrectable.
    pub invalid_barcode_counts: SortedBarcodeCountFile,

    /// Optionally specify a set of barcodes in a file. If supplied,
    /// only the barcodes listed in the file are processed through
    /// the aligner.
    /// TODO: Currently, we read all the reads from the shardio disk
    /// and then apply the filter prior to alignment. It would be more
    /// efficient to instead read only the list of barcodes we are
    /// interested in from disk. This requires a bit more invasive
    /// set of changes and is deferred to the future
    pub barcode_subset: Option<BarcodeSetFormat>,

    pub v1_pattern_fix_params: Option<V1PatternFixParams>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ChunkInputs {
    #[mro_type = "map"]
    range: Range<Barcode>,
    read_ann_subsample_rate: f32,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ChunkOutputs {
    counts_bc_order_shard: CountShardFile,
    probe_barcode_counts_shard: Option<CountShardFile>,
    bc_umi_info_shard: BcUmiInfoShardFile,
    pos_sorted_mapped_shard: AlignShardFile,
    pos_sorted_unmapped_shard: AlignShardFile,
    bam_header: Option<SamHeaderFile>,
    barcode_summary_shard: BinaryFormat<BarcodeSummaryFile, Vec<BarcodeSummary>>,
    annotation_files: Option<AnnotationFiles>,
    metrics_shard: BarcodeMetricsShardFile,
    no_tx_alignments: bool,
    per_read_gap_align: Vec<ParquetFile>,
}

#[derive(Serialize, MartianStruct)]
pub struct StageOutputs {
    /// feature x barcode counts in barcode order
    pub counts_bc_order: Vec<CountShardFile>,

    /// probe x barcode counts in barcode order
    pub probe_barcode_counts: Option<Vec<CountShardFile>>,

    /// UMI/molecule information for each barcode
    /// This includes read count per UMI and used
    /// to construct the molecule_info.h5 file.
    pub bc_umi_info: Vec<BcUmiInfoShardFile>,

    /// BAM records sorted in position-sorted order
    // in shardio that can be used to construct a BAM file
    pub pos_sorted: Vec<AlignShardFile>,

    /// Header of the BAM file to create, as a text file
    pub bam_header: Option<SamHeaderFile>,

    /// Barcode summary info, for use in filter barcodes
    pub barcode_summary: CsvFile<BarcodeSummary>,

    /// Annotation files for consumption by a PD stage
    pub annotation_files: Option<AnnotationFiles>,

    /// parquet files containing gap align information
    pub per_read_gap_align: Vec<ParquetFile>,

    /// Metrics computed per barcode sorted by (library type, barcode)
    pub per_barcode_metrics: Vec<BarcodeMetricsShardFile>,

    /// Metrics summary JSON.
    pub summary: JsonFile<AlignAndCountMetrics>,

    /// Used to disable useless passes over the asf files if
    /// they don't contain information from the star aligner
    pub no_tx_alignments: bool,
}

const MAX_ALIGN_CHUNKS_PER_LIB: usize = 100;
const READS_PER_CHUNK: usize = 15_000_000;
const READS_PER_CHUNK_LONG: usize = 3_750_000;

/// Martian stage ALIGN_AND_COUNT
/// Align the reads to the transcriptome and count reads and molecules per feature and barcode.
/// This stage receives all reads in corrected-barcode order and performs
/// alignment, UMI correction & deduplication, and gene counting for each
/// barcode. Gene counts and alignment records are output in shardio formats
/// for use downstream.
pub struct AlignAndCount;

/// A list of shard files that contains all the reads that need to be processed:
///  - `valid_reads`: A list of files containing`RnaRead`s with a valid barcode
///  - `corrected_reads`: A list of files containing `RnaRead`s which originally
///    had an invalid barcode, but we were able to correct it to a barcode in
///    our whitelist
/// - `invalid_reads`: a list of files containing `RnaRead`s with an invalid barcode
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ReadShards {
    /// Reads with valid barcodes
    pub valid_reads: Vec<ReadShardFile>,
    /// Reads with corrected barcodes
    pub corrected_reads: Vec<ReadShardFile>,
    /// Reads with invalid barcodes
    pub invalid_reads: Vec<ReadShardFile>,
}

impl ReadShards {
    /// Return a shard reader over all the reads.
    pub fn reader(&self) -> Result<ShardReader<RnaRead, BarcodeOrder>> {
        let mut shards = Vec::new();
        shards.extend(self.valid_reads.clone());
        shards.extend(self.corrected_reads.clone());
        shards.extend(self.invalid_reads.clone());
        let reader: ShardReader<RnaRead, BarcodeOrder> = ShardReader::open_set(&shards)?;
        Ok(reader)
    }
}

// Refer to either all the barcodes or a subset of barcodes
enum BarcodeSet {
    All,
    SubSet(TxHashSet<Barcode>),
}

impl BarcodeSet {
    fn new(barcode_subset: Option<&BarcodeSetFormat>) -> Result<Self> {
        Ok(match barcode_subset {
            Some(f) => BarcodeSet::SubSet(f.read()?),
            None => BarcodeSet::All,
        })
    }
    fn contains(&self, barcode: Barcode) -> bool {
        match self {
            BarcodeSet::All => true,
            BarcodeSet::SubSet(s) => s.contains(&barcode),
        }
    }
}

enum TxReference {
    Star(StarReference),
    Minimap2(MinimapReference),
}

impl TxReference {
    fn header_view(&self) -> HeaderView {
        match self {
            TxReference::Star(r) => r.header_view().to_owned(),
            TxReference::Minimap2(r) => HeaderView::from_header(&r.header),
        }
    }

    fn get_aligner(&self) -> Align {
        match self {
            TxReference::Star(r) => Align::Star(r.get_aligner()),
            TxReference::Minimap2(r) => Align::Minimap2(r.get_aligner()),
        }
    }
}

pub struct AlignThreadProcessor<V>
where
    V: AnnotatedReadVisitor,
{
    args: StageInputs,
    aligner: Aligner,
    bc_umi_info_sender: ShardSender<BcUmiInfo, BcUmiInfo>,
    bc_counts_sender: ShardSender<FeatureBarcodeCount, BarcodeThenFeatureOrder>,
    bc_probe_counts_sender: Option<ShardSender<ProbeBarcodeCount>>,
    per_read_gap_align_writer: Option<ParquetWriter<PerReadGapAlignRow>>,
    pos_mapped_reads_sender: ShardSender<Record, BamPosSort>,
    pos_unmapped_reads_sender: ShardSender<Record, BamPosSort>,
    visitor: V,
    // set of barcodes to align/annotate
    barcode_set: BarcodeSet,
    barcodes_to_subsample: TxHashSet<Barcode>,
    barcode_subsample_rate: Option<f64>,
}

impl<V> AlignThreadProcessor<V>
where
    V: AnnotatedReadVisitor,
{
    fn handle_gap_align_annotation(
        per_read_gap_align_writer: &mut Option<ParquetWriter<PerReadGapAlignRow>>,
        ann: &ReadAnnotations,
    ) -> Result<()> {
        if let Some(writer) = per_read_gap_align_writer
            && let RecordAnnotation::Probe(_, mapped_probe) = &ann.primary
            && let Some(gap_info) = mapped_probe.gap_info()
        {
            let read = &ann.read;
            assert!(
                !read.r2_exists(),
                "Paired end reads not expected in gap align!"
            );
            let read_gap_sequence =
                String::from_utf8(gap_info.gap_seq(read.r1_seq()).to_vec()).unwrap();
            let gap = mapped_probe.gap();
            let mgi = gap.as_ref().and_then(|x| {
                if x.gap_within_max_error() {
                    Some(MappedGapAlignmentInfo {
                        gap_seq: x.gap_seq.clone(),
                        expected_gap_seq: x.expected_gap_seq.clone(),
                        alignment_operations: x.get_alignment().operations.clone(),
                    })
                } else {
                    None
                }
            });

            writer.push(PerReadGapAlignRow {
                barcode: read.barcode().to_string(),
                is_valid_barcode: ann.read.barcode_is_valid(),
                umi: read.umi().to_string(),
                probe_type: gap_info.probe_type.to_string(),
                left_probe_id: mapped_probe.lhs_probe().unwrap().probe_id.to_string(),
                right_probe_id: mapped_probe.rhs_probe().unwrap().probe_id.to_string(),
                read_gap_sequence: read_gap_sequence.clone(),
                read_gap_len: read_gap_sequence.len(),
                expected_gap_sequence: gap.as_ref().map(MappedGap::get_expected_gap_seq),
                gap_levenshtein_distance: gap.as_ref().map(MappedGap::get_gap_levenshtein_distance),
                gap_within_max_error: gap.as_ref().map(MappedGap::gap_within_max_error),
                gap_exactly_matches_expected: gap.map(|x| read_gap_sequence == x.expected_gap_seq),
                num_matches: mgi.as_ref().map(MappedGapAlignmentInfo::get_num_matches),
                num_mismatches: mgi.as_ref().map(MappedGapAlignmentInfo::get_num_mismatches),
                num_deletions: mgi.as_ref().map(MappedGapAlignmentInfo::get_num_deletions),
                num_insertions: mgi.as_ref().map(MappedGapAlignmentInfo::get_num_insertions),
                starts_with_deletion: mgi
                    .as_ref()
                    .map(MappedGapAlignmentInfo::starts_with_deletion),
                starts_with_insertion: mgi
                    .as_ref()
                    .map(MappedGapAlignmentInfo::starts_with_insertion),
                ends_with_deletion: mgi.as_ref().map(MappedGapAlignmentInfo::ends_with_deletion),
                ends_with_insertion: mgi
                    .as_ref()
                    .map(MappedGapAlignmentInfo::ends_with_insertion),
            })?;
        }
        Ok(())
    }
}

impl<V> Proc for AlignThreadProcessor<V>
where
    V: AnnotatedReadVisitor,
{
    type Item = (Barcode, SpillVec<RnaRead, ReadSpillFormat>);
    type Err = Error;

    fn process(&mut self, item: Self::Item) -> Result<()> {
        let (bc, read_vec) = item;

        if !self.barcode_set.contains(bc) {
            return Ok(());
        }

        let barcode_subsample_rate = if self.barcodes_to_subsample.contains(&bc) {
            self.barcode_subsample_rate.unwrap()
        } else {
            1.0
        };
        let chemistry_name = self.args.chemistry_defs.primary().name;
        let alignment_options = AlignmentOptions {
            trim_polya_min_score: self.args.trim_polya_min_score,
            trim_tso_min_score: self.args.trim_tso_min_score,
            is_pd: self.args.is_pd,
        };
        let annotations_iter = self.aligner.process_barcode_se(
            chemistry_name,
            &alignment_options,
            self.args.no_bam,
            bc,
            read_vec.iter()?,
            barcode_subsample_rate,
        )?;
        let mut umi_counts = Vec::new();

        for ann in annotations_iter {
            let ann = ann?;
            if let Some(count) = ann.umi_count() {
                umi_counts.push(count);
            }

            self.visitor.visit_read_annotation(&ann);
            Self::handle_gap_align_annotation(&mut self.per_read_gap_align_writer, &ann)?;
            for r in ann.records() {
                if is_unmapped_to_any(&r) {
                    self.pos_unmapped_reads_sender.send(r)?;
                } else {
                    self.pos_mapped_reads_sender.send(r)?;
                }
            }
        }

        if bc.is_valid() {
            // Required for ordering in molecule info
            umi_counts.sort();
            let bc_umi_info = BcUmiInfo {
                barcode: bc,
                umi_counts,
            };

            // Generate individual feature counts & send
            for c in bc_umi_info.feature_counts() {
                self.bc_counts_sender.send(c)?;
            }

            // Generate individual probe counts & send
            if let Some(bc_probe_counts_sender_unwrapped) = self.bc_probe_counts_sender.as_mut() {
                for c in bc_umi_info.probe_counts() {
                    bc_probe_counts_sender_unwrapped.send(c)?;
                }
            }

            self.bc_umi_info_sender.send(bc_umi_info)?;
        }

        Ok(())
    }
}

fn star_settings(
    reference_path: &Path,
    transcriptome_min_score: Option<usize>,
    star_parameters: Option<&str>,
) -> Result<StarSettings> {
    let mut settings = StarSettings::new(reference_path.join("star").to_str().unwrap());
    let parameters = star_parameters.unwrap_or("");
    if parameters.contains("--outFilterScoreMin") {
        warn!("using --outFilterScoreMin from star_parameters");
    } else {
        settings = settings.arg(&format!(
            "--outFilterScoreMin={}",
            transcriptome_min_score.unwrap_or(DEFAULT_TRANSCRIPTOME_MIN_SCORE)
        ));
    }
    for arg in parameters.split_ascii_whitespace() {
        settings = settings.arg(arg);
    }
    Ok(settings)
}

/// Choose an aligner based on chemistry and input aligner argument.
/// If target_set=null or target set CSV is not a probe set return STAR as the default aligner.
/// If inferred chemistry is RTL return Hurtle as the default aligner.
fn choose_aligner(
    arg_aligner: Option<AlignerParam>,
    is_rtl: Option<bool>,
    target_set: Option<&TargetSetFile>,
) -> Result<AlignerParam> {
    use AlignerParam::{Hurtle, Minimap2, Star};
    match target_set {
        None => match arg_aligner {
            None => Ok(Star),
            Some(Hurtle) => bail!("aligner=hurtle is incompatible with target_set=null"),
            Some(aligner) => Ok(aligner),
        },
        Some(target_set) => {
            let is_probe_set =
                ProbeSetReferenceMetadata::load_from(target_set)?.is_probe_set_metadata();

            match (is_rtl.unwrap_or(is_probe_set), arg_aligner) {
                (false, None) => Ok(Star),
                (false, Some(Hurtle)) => bail!(
                    "aligner=hurtle is incompatible with target_set={}",
                    target_set.display(),
                ),
                (false, Some(aligner)) => Ok(aligner),
                (true, None) => Ok(Hurtle),
                (true, Some(Minimap2)) => bail!("aligner=minimap2 is incompatible with rtl"),
                (true, Some(aligner)) => Ok(aligner),
            }
        }
    }
}

/// Choose an aligner based on the chemistry and target_set CSV file.
fn choose_aligner_from_args(args: &StageInputs) -> Result<AlignerParam> {
    choose_aligner(
        args.aligner,
        args.chemistry_defs.is_rtl(),
        args.target_set.as_ref(),
    )
}

/// When not aligning the reads with STAR, we assign the tid as the lhs probe index and the
/// pos is equal to the rhs index, this gives us a dummy "reference" to use with this.
fn make_dummy_header_for_probes(n_probes: usize) -> HeaderView {
    let mut header = Header::new();
    for i in 0..n_probes {
        let mut header_rec = HeaderRecord::new(b"SQ");
        header_rec.push_tag(b"SN", format!("PROBE{i}"));
        header_rec.push_tag(b"LN", n_probes + 1);
        header.push_record(&header_rec);
    }
    HeaderView::from_header(&header)
}

fn write_output_bam_header<T: Borrow<HeaderView>>(
    bam_header_data: T,
    rover: &MartianRover,
) -> Result<Option<SamHeaderFile>> {
    let bam_header: SamHeaderFile = rover.make_path("bam_header");
    let mut f = File::create(&bam_header)?;
    f.write_all(bam_header_data.borrow().as_bytes())?;
    Ok(Some(bam_header))
}

#[make_mro(mem_gb = 9, volatile = strict)]
impl MartianStage for AlignAndCount {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;
    type ChunkInputs = ChunkInputs;
    type ChunkOutputs = ChunkOutputs;

    /// Create upto `MAX_ALIGN_CHUNKS_PER_LIB` chunks per physical library with at least
    /// `READS_PER_CHUNK` read pairs per chunk
    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let aligner = choose_aligner_from_args(&args)?;
        let ref_gb = if args.no_bam && aligner == AlignerParam::Hurtle {
            // For probe assays, STAR is needed only to produce alignments for the BAM file.
            // The memory used by the Hurtle reference is non-zero, but it's less than the
            // wiggle room added to mem_gb below.
            0.0
        } else if let Some(reference_path) = &args.reference_path {
            let settings = star_settings(
                reference_path,
                args.transcriptome_min_score,
                args.star_parameters.as_deref(),
            )?;
            settings.est_mem()? as f64 / ALN_BC_GIB
        } else {
            0.0
        };

        let make_resource =
            |range: &Range<Barcode>,
             top_n_counts: &TopCounts,
             reader: &ShardReader<RnaRead, BarcodeOrder>| {
                let max_reads = top_n_counts.sum_with_limit(MAX_ITEMS_IN_MEM);
                let max_anns = top_n_counts.sum_with_limit(MAX_ANNOTATIONS_IN_MEM);
                // if we presume ~1KB per RnaRead for short-reads
                let rna_read_size = match aligner {
                    AlignerParam::Minimap2 => RNA_READ_LONG_SZ,
                    _ => RNA_READ_SZ,
                };
                let reads_gb = (max_reads * rna_read_size) as f64 / ALN_BC_GIB;
                // if we presume ~3KB per ReadAnnotations for short-reads
                let ann_size = match aligner {
                    AlignerParam::Minimap2 => ANN_LONG_SZ,
                    _ => ANN_SZ,
                };
                let anns_gb = ((max_anns + BAM_BUFFER_SZ) * ann_size) as f64 / ALN_BC_GIB;
                // account for maximum overhead of shardio reading (ShardIter and a Record)
                let shardio_gb = ((reader.shard_iter_size() + ann_size)
                    * reader.est_len_range(range)
                    / ALN_BC_DISK_CHUNK_SZ) as f64
                    / ALN_BC_GIB;
                // blow up in really extreme cases, this has been tuned against a
                // 13.04GB reference memory reservation
                let mem_gb = (reads_gb + anns_gb + shardio_gb).min(24.0);
                let mem_gb = (2.0 + mem_gb + ref_gb).ceil() as isize;
                Resource::new()
                    .threads(NUM_CHUNK_THREADS as isize)
                    .mem_gb(mem_gb)
                    .vmem_gb(mem_gb + 6)
            };

        // Figure out how many chunks to create.
        let reader = args.read_shards.reader()?;

        let n = reader.len(); // This is the total number of read pairs across all the shard files in this gem well

        // DETECT_CHEMISTRY should generally prevent us from hitting this, but kill the pipeline
        ensure!(
            n > 0,
            "Zero read-pairs in the input data were detected as valid for chemistry {}.\n\
             Check the --chemistry argument to Cell Ranger",
            args.chemistry_defs.name(),
        );

        // Create upto `MAX_ALIGN_CHUNKS_PER_LIB` chunks, with at least `READS_PER_CHUNK` read pairs per chunk
        let reads_per_chunk = match aligner {
            AlignerParam::Minimap2 => READS_PER_CHUNK_LONG,
            _ => READS_PER_CHUNK,
        };
        let num_chunks = (n / reads_per_chunk).clamp(1, MAX_ALIGN_CHUNKS_PER_LIB);

        // Stream the invalid counts from disk, chained with the valid counts.
        let chunks = args
            .invalid_barcode_counts
            .lazy_reader()?
            .chain(args.corrected_barcode_counts.iter_counts()?)
            .process_results(move |counts_iter| make_chunks(n, num_chunks, counts_iter))?;

        // Pick a sampling rate for AnnotatedReads output. This output is only usef by PD stages.
        let read_ann_subsample_rate = (MAX_READ_ANN_SAMPLE as f32) / (max(n, 1) as f32);
        let read_ann_subsample_rate = if read_ann_subsample_rate > 1.0 {
            1.0
        } else {
            read_ann_subsample_rate
        };

        let stage_def: StageDef<_> = chunks
            .into_iter()
            .map(|chunk| {
                // Create a chunk that will process all the reads that fall within this barcode range.
                (
                    ChunkInputs {
                        range: chunk.range,
                        read_ann_subsample_rate,
                    },
                    make_resource(&chunk.range, &chunk.top_n_counts, &reader),
                )
            })
            .collect();

        Ok(stage_def.join_resource(Resource::with_mem_gb(1)))
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        // Disable polyA and TSO trimming for 5' gene expression assay.
        let args = if args.chemistry_defs.endedness() == Some(WhichEnd::FivePrime) {
            Self::StageInputs {
                trim_polya_min_score: None,
                trim_tso_min_score: None,
                ..args
            }
        } else {
            args
        };

        let aligner_param = choose_aligner_from_args(&args)?;

        let probe_set_reference = if aligner_param == AlignerParam::Hurtle {
            Some(ProbeSetReference::from_path(
                args.target_set.as_ref().unwrap(),
                args.reference_path.as_deref(),
                args.transcriptome_min_score
                    .unwrap_or(DEFAULT_TRANSCRIPTOME_MIN_SCORE),
            )?)
        } else {
            None
        };

        let tx_reference = match (aligner_param, args.no_bam) {
            // STAR is needed only to create the BAM file when using Hurtle.
            (AlignerParam::Hurtle, true) => None,
            (AlignerParam::Star, _) | (AlignerParam::Hurtle, false) => {
                // load up the aligner and annotator
                // do this first because Orbit loads the reference in the background
                // with mmap, so we want to start this ASAP.
                if let Some(reference_path) = &args.reference_path {
                    let star = StarReference::load(star_settings(
                        reference_path,
                        args.transcriptome_min_score,
                        args.star_parameters.as_deref(),
                    )?)?;
                    Some(TxReference::Star(star))
                } else {
                    None
                }
            }
            (AlignerParam::Minimap2, _) => {
                let mmap = MinimapReference::load(
                    Preset::Splice,
                    args.chemistry_defs.primary().strandedness,
                    args.minimap_index.clone().unwrap().as_ref(),
                    args.tx_bed.clone().unwrap().as_ref(),
                );
                Some(TxReference::Minimap2(mmap))
            }
        };

        let reader = args.read_shards.reader()?;

        // Output files
        let counts_bc_order_shard: CountShardFile = rover.make_path("bc_sort");
        let probe_barcode_counts_shard: Option<CountShardFile> =
            if aligner_param == AlignerParam::Hurtle {
                Some(rover.make_path("bc_probe_sort"))
            } else {
                None
            };
        let bc_umi_info_shard: BcUmiInfoShardFile = rover.make_path("bc_umi_info");
        let pos_sorted_mapped_shard: AlignShardFile = rover.make_path("pos_sorted_mapped");
        let pos_sorted_unmapped_shard: AlignShardFile = rover.make_path("pos_sorted_unmapped");
        let metrics_shard: BarcodeMetricsShardFile = rover.make_path("metrics_shard");

        // Count data ordered by barcode.
        let mut bc_counts: ShardWriter<FeatureBarcodeCount, BarcodeThenFeatureOrder> =
            ShardWriter::with_compressor(
                &counts_bc_order_shard,
                ALN_BC_SEND_BUFFER_SZ,
                ALN_BC_DISK_CHUNK_SZ,
                ALN_BC_ITEM_BUFFER_SZ,
                shardio::Compressor::Lz4,
            )?;

        let mut bc_probe_counts: Option<ShardWriter<ProbeBarcodeCount>> =
            probe_barcode_counts_shard.as_ref().map(|x| {
                ShardWriter::with_compressor(
                    x,
                    ALN_BC_SEND_BUFFER_SZ,
                    ALN_BC_DISK_CHUNK_SZ,
                    ALN_BC_ITEM_BUFFER_SZ,
                    shardio::Compressor::Lz4,
                )
                .unwrap()
            });

        // umi info data
        let mut bc_umi_info: ShardWriter<BcUmiInfo, BcUmiInfo> = ShardWriter::with_compressor(
            &bc_umi_info_shard,
            1,
            256,
            2048,
            shardio::Compressor::Lz4,
        )?;

        // Position sorted mapped read data
        let mut pos_mapped_reads: ShardWriter<Record, BamPosSort> = ShardWriter::with_compressor(
            &pos_sorted_mapped_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            BAM_BUFFER_SZ,
            shardio::Compressor::Lz4,
        )?;

        // Position sorted unmapped read data; we write a separate shard file
        // for the unmapped reads and then group them with the mapped reads in
        // a single collection of files to enable transparent downstream optimizations
        // when we only need to read the mapped or unmapped reads subset.
        let mut pos_unmapped_reads: ShardWriter<Record, BamPosSort> = ShardWriter::with_compressor(
            &pos_sorted_unmapped_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            BAM_BUFFER_SZ,
            shardio::Compressor::Lz4,
        )?;

        // Per barcode metrics
        let mut metrics_writer: ShardWriter<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardWriter::with_compressor(
                &metrics_shard,
                1024,
                8192,
                100_000,
                shardio::Compressor::Lz4,
            )?;

        // Get the BAM header & write it out.
        let bam_header = if let Some(ref tx_reference) = tx_reference {
            write_output_bam_header(tx_reference.header_view(), &rover)?
        } else if aligner_param == AlignerParam::Hurtle && (args.no_bam || args.is_pd) {
            write_output_bam_header(
                make_dummy_header_for_probes(
                    probe_set_reference.as_ref().unwrap().number_of_probes(),
                ),
                &rover,
            )?
        } else {
            None
        };

        let feature_reference = args.feature_reference_binary.read()?;
        let feature_dist = compute_feature_dist(args.feature_counts.read()?, &feature_reference)?;
        let target_genes = feature_reference.target_genes();

        let annotator = if let Some(reference_path) = &args.reference_path {
            Some(ReadAnnotator::new(
                reference_path,
                args.chemistry_defs.primary(),
                args.include_exons,
                args.include_introns,
                aligner_param,
            )?)
        } else {
            None
        };

        let n_threads = rover.get_threads().max(1);

        let per_read_gap_align_files = match &probe_set_reference {
            Some(probe_set_ref) if probe_set_ref.has_gap_probes => (0..n_threads)
                .map(|i| {
                    Some(rover.make_path::<ParquetFile>(format!("per_read_gap_alignments_{i}")))
                })
                .collect(),
            _ => vec![None; n_threads],
        };
        let per_read_gap_align_writers: Vec<_> = per_read_gap_align_files
            .iter()
            .map(|f| {
                f.as_ref()
                    .map(|f| f.writer(PerReadGapAlignRow::ROW_GROUP_SIZE))
                    .transpose()
            })
            .try_collect()?;

        let no_tx_alignments = tx_reference.is_none();
        let aligner = Aligner::new(
            tx_reference.map(|x| x.get_aligner()),
            feature_reference,
            feature_dist,
            annotator,
            args.read_chunks.clone(),
            rover.files_path().into(),
            probe_set_reference,
        )?;

        let ann_files: Vec<ReadAnnotationsFormat> = (0..n_threads)
            .map(|i| rover.make_path(format!("read_annotations_{i}")))
            .collect();

        let (barcodes_to_subsample, barcode_subsample_rate) = match &args.v1_pattern_fix_params {
            Some(v1_pattern_fix_params) => v1_pattern_fix_params.barcode_subsampling()?,
            None => (TxHashSet::default(), None),
        };

        // Each thread that handles reads uses one clone() of AlignThreadProcessor struct
        // It contains:
        // - An owned StarAligner that has it's own scratch data for aligning reads
        // - An Arc<StarReference> which is the shared reference data structure used by STAR
        // - Shared FeatureReference and ReadAnnotator structs used to
        let processors = zip_eq(&ann_files, per_read_gap_align_writers)
            .enumerate()
            .map(|(idx, (ann_file, per_read_gap_align_writer))| {
                anyhow::Ok(AlignThreadProcessor {
                    args: args.clone(),
                    aligner: aligner.clone(),
                    bc_umi_info_sender: bc_umi_info.get_sender(),
                    bc_counts_sender: bc_counts.get_sender(),
                    pos_mapped_reads_sender: pos_mapped_reads.get_sender(),
                    pos_unmapped_reads_sender: pos_unmapped_reads.get_sender(),
                    bc_probe_counts_sender: bc_probe_counts.as_mut().map(|x| x.get_sender()),
                    per_read_gap_align_writer,
                    visitor: if args.is_pd {
                        let rng = SmallRng::seed_from_u64(idx as u64);
                        StageVisitor::with_ann_writer_sample(
                            metrics_writer.get_sender(),
                            target_genes.clone(),
                            ann_file.lazy_writer()?,
                            chunk_args.read_ann_subsample_rate,
                            rng,
                        )
                    } else {
                        StageVisitor::new(metrics_writer.get_sender(), target_genes.clone())
                    },
                    barcode_set: BarcodeSet::new(args.barcode_subset.as_ref())?,
                    barcodes_to_subsample: barcodes_to_subsample.clone(),
                    barcode_subsample_rate,
                })
            })
            .try_collect()?;

        // get the range of barcode we're going to analyze
        let read_iter = reader.iter_range(&chunk_args.range)?;

        // Process barcode groups, using one thread for each processor
        // Return the processors when complete.
        // See impl Proc for AlignThreadProcessor above for the per-barcode inner loop.
        let processors =
            group_by_processor(read_iter, processors, RnaRead::barcode, rover.files_path())?;

        // Pull together the lightweight bc summary info.
        let mut bc_summaries = Vec::new();
        // Total number of reads written to the annotation files
        let mut num_reads_annotation_files = 0;
        for mut p in processors {
            bc_summaries.extend(p.visitor.barcode_summaries());
            num_reads_annotation_files += p.visitor.ann_writer_num_reads();
            p.visitor.finish()?;

            // shut down ShardSenders
            p.bc_umi_info_sender.finished()?;
            p.bc_counts_sender.finished()?;
            if let Some(x) = p.bc_probe_counts_sender.as_mut() {
                x.finished().unwrap();
            };
            p.pos_mapped_reads_sender.finished()?;
            p.pos_unmapped_reads_sender.finished()?;
        }

        // shutdown shardio writers
        bc_umi_info.finish()?;
        bc_counts.finish()?;
        bc_probe_counts.as_mut().map(|x| x.finish().unwrap());
        pos_mapped_reads.finish()?;
        pos_unmapped_reads.finish()?;
        metrics_writer.finish()?;

        let barcode_summary_shard: BinaryFormat<_, _> = rover.make_path("barcode_summary");
        barcode_summary_shard.write(&bc_summaries)?;

        Ok(ChunkOutputs {
            counts_bc_order_shard,
            probe_barcode_counts_shard,
            bc_umi_info_shard,
            pos_sorted_mapped_shard,
            pos_sorted_unmapped_shard,
            bam_header,
            barcode_summary_shard,
            annotation_files: if args.is_pd {
                Some(AnnotationFiles {
                    num_reads: num_reads_annotation_files,
                    files: ann_files,
                })
            } else {
                None
            },
            per_read_gap_align: per_read_gap_align_files.into_iter().flatten().collect(),
            metrics_shard,
            no_tx_alignments,
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let bam_header = if let Some(ref chunk_bam_header) = chunk_outs[0].bam_header {
            let bam_header = rover.make_path("bam_header");
            std::fs::copy(chunk_bam_header, &bam_header)?;
            Some(bam_header)
        } else {
            None
        };

        let barcode_summary: CsvFile<BarcodeSummary> = rover.make_path("barcode_summary");
        let mut writer = barcode_summary.lazy_writer()?;
        for chunk in &chunk_outs {
            for record in chunk
                .barcode_summary_shard
                .read()?
                .into_iter()
                // Threading causes some flutter in the order.
                .sorted_unstable()
            {
                writer.write_item(&record)?;
            }
        }
        writer.finish()?;

        let summary: JsonFile<AlignAndCountMetrics> = rover.make_path("summary");
        let aligner = choose_aligner_from_args(&args)?;
        summary.write(&AlignAndCountMetrics {
            alignment_aligner: aligner,
            alignment_high_conf_mapq: aligner.high_conf_mapq() as i64,
        })?;

        Ok(StageOutputs {
            counts_bc_order: chunk_outs
                .iter()
                .map(|v| v.counts_bc_order_shard.clone())
                .collect(),
            probe_barcode_counts: {
                let x: Vec<_> = chunk_outs
                    .iter()
                    .filter_map(|v| v.probe_barcode_counts_shard.clone())
                    .collect();
                if x.is_empty() { None } else { Some(x) }
            },
            bc_umi_info: chunk_outs
                .iter()
                .map(|v| v.bc_umi_info_shard.clone())
                .collect(),
            pos_sorted: chunk_outs
                .iter()
                .flat_map(|v| {
                    [
                        v.pos_sorted_mapped_shard.clone(),
                        v.pos_sorted_unmapped_shard.clone(),
                    ]
                })
                .collect(),
            bam_header,
            barcode_summary,
            annotation_files: if args.is_pd {
                let mut annotation_files = AnnotationFiles {
                    num_reads: 0,
                    files: Vec::new(),
                };
                for chunk_out in &chunk_outs {
                    // Will be present if `is_pd` is true
                    let chunk_ann_files = chunk_out.annotation_files.as_ref().unwrap();
                    annotation_files.num_reads += chunk_ann_files.num_reads;
                    annotation_files
                        .files
                        .extend(chunk_ann_files.files.iter().cloned());
                }
                Some(annotation_files)
            } else {
                None
            },
            per_read_gap_align: chunk_outs
                .iter()
                .flat_map(|x| x.per_read_gap_align.clone())
                .collect(),
            per_barcode_metrics: chunk_outs.iter().map(|x| x.metrics_shard.clone()).collect(),
            summary,
            no_tx_alignments: chunk_outs[0].no_tx_alignments,
        })
    }
}

/// Chunk the barcode space into num_chunks that approximately evenly divide up all reads.
/// Track the top barcodes in each chunk for use in memory estimation.
fn make_chunks(
    num_reads: usize,
    num_chunks: usize,
    bc_count_iter: impl Iterator<Item = (Barcode, usize)>,
) -> Vec<ChunkData> {
    // Augment the count iterator to maintain a running sum.
    // We use this running sum to chunk the iterator.
    let bc_count_iter = bc_count_iter.scan(0, |sum, (bc, count)| {
        *sum += count;
        Some((bc, count, *sum))
    });

    let approx_n_per_chunk = num_reads / num_chunks;

    // Derive a chunk key for each barcode.
    let chunker =
        bc_count_iter.chunk_by(|(_, _, sum)| (sum / approx_n_per_chunk).min(num_chunks - 1));

    let mut chunks: Vec<_> = chunker
        .into_iter()
        .map(|(_, counts_iter)| {
            let mut counts_iter = counts_iter.map(|(bc, count, _sum)| (bc, count)).peekable();
            let start_bc = counts_iter.peek().unwrap().0;
            let top_n_counts = counts_iter.fold(
                TopCounts::new(NUM_CHUNK_THREADS),
                |mut top_n, (_, count)| {
                    top_n.observe(count);
                    top_n
                },
            );
            ChunkData {
                range: Range {
                    start: Some(start_bc),
                    end: None, // we will populate end later
                },
                top_n_counts,
            }
        })
        .collect();
    // Populate range ends from next start.
    let mut chunk_iter = chunks.iter_mut().peekable();
    while let Some(chunk) = chunk_iter.next() {
        chunk.range.end = chunk_iter
            .peek()
            .and_then(|next_chunk| next_chunk.range.start);
    }
    chunks
}

/// Track top N counts.
struct TopCounts(Vec<Reverse<usize>>);

impl TopCounts {
    pub fn new(n: usize) -> Self {
        Self(vec![Reverse(0); n])
    }

    pub fn observe(&mut self, count: usize) {
        self.0.push(Reverse(count));
        self.0.sort_unstable();
        self.0.pop();
    }

    /// Sum all counts, clipping each value at limit.
    pub fn sum_with_limit(&self, limit: usize) -> usize {
        self.0.iter().map(|c| c.0.min(limit)).sum()
    }
}

struct ChunkData {
    /// The barcode range to read for the chunk.
    pub range: Range<Barcode>,
    /// The counts of the top N most populous barcodes for the chunk.
    pub top_n_counts: TopCounts,
}

#[cfg(test)]
mod test {
    use super::*;
    use cr_bam::bam::BamPosSort;
    use rust_htslib::bam::{HeaderView, Record};
    use shardio::SortKey;

    #[test]
    fn test_top_n() {
        let counts = [0, 5, 3, 12, 0, 1, 4, 5];
        let mut top_n = TopCounts::new(3);
        for c in counts {
            top_n.observe(c);
        }
        assert_eq!(
            top_n.0.iter().map(|c| c.0).collect::<Vec<_>>(),
            vec![12, 5, 5]
        );
        assert_eq!(top_n.sum_with_limit(4), 12);
        assert_eq!(top_n.sum_with_limit(10), 20);
    }

    #[test]
    fn test_choose_aligner() -> Result<()> {
        use AlignerParam::{Hurtle, Star};
        use cr_types::chemistry::ChemistryName::{
            MFRP_RNA, SpatialThreePrimeV1, ThreePrimeV3PolyA,
        };

        let hybcap_file = "test/target_panels/Immunology_targeting_hybrid.csv".into();
        let hybcap = Some(&hybcap_file);
        let rtl_file = "test/target_panels/Immunology_targeting_templated_ligation.csv".into();
        let rtl = Some(&rtl_file);

        // Test aligner and target_set=null.
        assert_eq!(choose_aligner(None, None, None)?, Star);
        assert!(choose_aligner(Some(Hurtle), None, None).is_err());
        assert_eq!(choose_aligner(Some(Star), None, None)?, Star);

        // Test aligner and target_panel_file_format.
        assert_eq!(choose_aligner(None, None, hybcap)?, Star);
        assert!(choose_aligner(Some(Hurtle), None, hybcap).is_err());
        assert_eq!(choose_aligner(Some(Star), None, hybcap)?, Star);

        // Test aligner and probe_set_file_format.
        assert_eq!(choose_aligner(None, None, rtl)?, Hurtle);
        assert_eq!(choose_aligner(Some(Hurtle), None, rtl)?, Hurtle);
        assert_eq!(choose_aligner(Some(Star), None, rtl)?, Star);

        // Test chemistry and target_panel_file_format.
        let spatial_is_rtl = SpatialThreePrimeV1.is_rtl();
        assert_eq!(
            choose_aligner(None, ThreePrimeV3PolyA.is_rtl(), hybcap)?,
            Star
        );
        assert!(choose_aligner(Some(Hurtle), None, hybcap).is_err());
        assert_eq!(choose_aligner(None, spatial_is_rtl, hybcap)?, Star);

        // Test chemistry and probe_set_file_format.
        assert_eq!(choose_aligner(None, ThreePrimeV3PolyA.is_rtl(), rtl)?, Star);
        assert_eq!(choose_aligner(None, MFRP_RNA.is_rtl(), rtl)?, Hurtle);
        assert_eq!(choose_aligner(None, spatial_is_rtl, rtl)?, Hurtle);

        Ok(())
    }

    #[test]
    fn test_bam_pos_sort_key() {
        /*  The pyO3 code assumes that the sort_key function creates certain values for
        unmapped reads, we verify that here. If this test fails, you need to update the
        unmappedreadsiterator to be compatible with the assumptions about the sort order key
        that this test is verifying. */
        let fake_header = "@SQ\tSN:CHROMOSOME_I\tLN:15072423\n".to_string();
        let header = HeaderView::from_bytes(&fake_header.into_bytes());
        let low_read = "FAKEQNAME\t4\t*\t0\t0\t*\t*\t0\t4\tACGT\t<<<<"
            .to_string()
            .into_bytes();
        // same as above but with higher position at 2^31 -1
        let high_read = "FAKEQNAME\t4\t*\t2147483647\t0\t*\t*\t0\t4\tACGT\t<<<<"
            .to_string()
            .into_bytes();
        let low_record = Record::from_sam(&header, &low_read).unwrap();
        let high_record = Record::from_sam(&header, &high_read).unwrap();
        let low_key = <BamPosSort as SortKey<Record>>::sort_key(&low_record);
        let high_key = <BamPosSort as SortKey<Record>>::sort_key(&high_record);
        assert_eq!(low_key.target_id, 4294967295);
        assert_eq!(high_key.target_id, 4294967295);
        assert_eq!(low_key.position, -1);
        assert_eq!(high_key.position, 2147483646);
    }
}
