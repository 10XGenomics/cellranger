//! Martian stage ALIGN_AND_COUNT

use crate::align_and_count_metrics::StageVisitor;
use crate::align_metrics::{BarcodeMetrics, LibFeatThenBarcodeOrder};
use crate::aligner::{Aligner, BarcodeSummary, MAX_ANNOTATIONS_IN_MEM};
use crate::barcode_sort::BarcodeOrder;
#[cfg(feature = "tenx_internal")]
use crate::stages::internal::get_barcode_subsampling;
#[cfg(feature = "tenx_source_available")]
use crate::stages::stubs::get_barcode_subsampling;
use crate::types::{
    BarcodeMetricsShardFile, FeatureReferenceFormat, ReadShardFile, ReadSpillFormat,
};
use crate::{AlignShardFile, BcUmiInfoShardFile};
use anyhow::{ensure, Result};
use barcode::Barcode;
use cr_bam::bam::BamPosSort;
use cr_bam::constants::{
    ALN_BC_DISK_CHUNK_SZ, ALN_BC_GIB, ALN_BC_ITEM_BUFFER_SZ, ALN_BC_SEND_BUFFER_SZ,
};
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::probe_set::{ProbeSetReference, ProbeSetReferenceMetadata};
use cr_types::reference::feature_checker::compute_feature_dist;
use cr_types::reference::feature_reference::TargetSetFile;
use cr_types::rna_read::{RnaChunk, RnaRead, HIGH_CONF_MAPQ};
use cr_types::spill_vec::SpillVec;
use cr_types::types::{
    BarcodeSetFormat, BcUmiInfo, FeatureBarcodeCount, GemWell, ProbeBarcodeCount,
};
use cr_types::{
    AlignerParam, BarcodeThenFeatureOrder, CountShardFile, FeatureCountFormat, TotalBcCountFormat,
};
use fastq_set::WhichEnd;
use itertools::Itertools;
use log::warn;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO};
use metric::{CountMetric, TxHashSet};
use orbit::{StarReference, StarSettings};
use par_proc::par_proc::{group_by_processor, Proc, MAX_ITEMS_IN_MEM};
use parameters_toml::star_parameters;
use rand::SeedableRng;
// rand::StdRng does not have a stable RNG across versions, so we explicitly choose one
use rand_chacha::ChaCha20Rng;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::{HeaderView, Record};
use serde::{Deserialize, Serialize};
use shardio::{Range, ShardReader, ShardSender, ShardWriter, SHARD_ITER_SZ as SHARD_SZ};
use std::borrow::Borrow;
use std::cmp::{max, Reverse};
use std::collections::{BTreeMap, BinaryHeap};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use tx_annotation::read::{AnnotationFiles, AnnotationInfo, ReadAnnotationsFormat, ReadAnnotator};
use tx_annotation::visitor::AnnotatedReadVisitor;

/// Default value for the stage parameter transcriptome_min_score, STAR parameter --outFilterScoreMin.
pub const DEFAULT_TRANSCRIPTOME_MIN_SCORE: usize = 30;

const RNA_READ_SZ: usize = 1024;
const ANN_SZ: usize = 3072; // RnaRead + BamRecords
const BAM_BUFFER_SZ: usize = 500_000;
const NUM_CHUNK_THREADS: usize = 4;
const PRESUMED_BARCODE_COUNT: usize = 2_000;

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
    pub gem_well: GemWell,
    pub read_chunks: Vec<RnaChunk>,
    pub reference_path: PathBuf,

    pub read_shards: ReadShards,

    pub feature_counts: FeatureCountFormat,
    pub feature_reference: FeatureReferenceFormat,

    /// The target panel CSV file.
    pub target_set: Option<TargetSetFile>,

    pub chemistry_defs: ChemistryDefs,

    pub aligner: Option<AlignerParam>,
    pub include_exons: bool,
    pub include_introns: bool,
    pub is_pd: bool,
    pub no_bam: bool,

    /// If this is Some(r), we filter out (UMI, genes) pairs
    /// with **less than** r reads and not include them in UMI counts.
    /// Useful in targeting
    pub targeted_umi_min_read_count: Option<u64>,

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

    pub total_barcode_counts: TotalBcCountFormat,

    /// Optionally specify a set of barcodes in a file. If supplied,
    /// only the barcodes listed in the file are processed through
    /// the aligner.
    /// TODO: Currently, we read all the reads from the shardio disk
    /// and then apply the filter prior to alignment. It would be more
    /// efficient to instead read only the list of barcodes we are
    /// interested in from disk. This requires a bit more invasive
    /// set of changes and is deferred to the future
    pub barcode_subset: Option<BarcodeSetFormat>,

    #[allow(unused)]
    pub chevron_correction_factor: Option<f64>,
    #[allow(unused)]
    pub chevron_affected_barcodes: Option<JsonFile<Vec<String>>>,
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
    pos_sorted_shard: AlignShardFile,
    bam_header: Option<PathBuf>,
    barcode_summary_shard: BinaryFormat<BarcodeSummaryFile, Vec<BarcodeSummary>>,
    annotation_files: Option<AnnotationFiles>,
    metrics_shard: BarcodeMetricsShardFile,
    no_star_alignments: bool,
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
    pub bam_header: Option<PathBuf>,

    /// Barcode summary info, for use in filter barcodes
    pub barcode_summary: CsvFile<BarcodeSummary>,

    /// Annotation files for consumption by a PD stage
    pub annotation_files: Option<AnnotationFiles>,

    /// Metrics computed per barcode sorted by (library type, barcode)
    pub per_barcode_metrics: Vec<BarcodeMetricsShardFile>,

    /// Metrics summary JSON.
    pub summary: JsonFile<AlignAndCountMetrics>,

    /// Used to disable useless passes over the asf files if
    /// they don't contain information from the star aligner
    pub no_star_alignments: bool,
}

pub const MAX_ALIGN_CHUNKS_PER_LIB: usize = 100;
pub const READS_PER_CHUNK: usize = 15_000_000;

/// Align the reads to the transcriptome and count reads and molecules per feature and barcode.
/// This stage receives all reads in corrected-barcode order and performs
/// alignment, UMI correction & deduplication, and gene counting for each
/// barcode. Gene counts and alignment records are output in shardio formats
/// for use downstream.
pub struct AlignAndCount;

/// A list of shard files that contains all the reads that need to be processed:
///  - `valid_reads`: A list of iles containing`RnaRead`s with a valid barcode
///  - `corrected_reads`: A list of files containing `RnaRead`s which originally had an invalid
///                       barcode, but we were able to correct it to a barcode in our whitelist
/// - `invalid_reads`: a list of files containing `RnaRead`s with an invalid barcode
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ReadShards {
    pub valid_reads: Vec<ReadShardFile>,
    pub corrected_reads: Vec<ReadShardFile>,
    pub invalid_reads: Vec<ReadShardFile>,
}

impl ReadShards {
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

pub struct AlignThreadProcessor<V>
where
    V: AnnotatedReadVisitor,
{
    args: StageInputs,
    aligner: Aligner,
    bc_umi_info_sender: ShardSender<BcUmiInfo, BcUmiInfo>,
    bc_counts_sender: ShardSender<FeatureBarcodeCount, BarcodeThenFeatureOrder>,
    bc_probe_counts_sender: Option<ShardSender<ProbeBarcodeCount>>,
    pos_reads_sender: ShardSender<Record, BamPosSort>,
    visitor: V,
    // set of barcodes to align/annotate
    barcode_set: BarcodeSet,
    barcodes_to_subsample: TxHashSet<Barcode>,
    barcode_subsample_rate: Option<f64>,
}

impl<V> Proc for AlignThreadProcessor<V>
where
    V: AnnotatedReadVisitor,
{
    type Item = (Barcode, SpillVec<RnaRead, ReadSpillFormat>);
    type Err = anyhow::Error;

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

        let annotations_iter = self.aligner.process_barcode_se(
            &self.args,
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
            for r in ann.records() {
                self.pos_reads_sender.send(r)?;
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
) -> Result<StarSettings> {
    let mut settings = StarSettings::new(reference_path.join("star").to_str().unwrap());
    let parameters = star_parameters()?;
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

fn max_items_in_range(
    counts: &BTreeMap<Barcode, CountMetric>,
    range: &Range<Barcode>,
    limit: usize,
) -> usize {
    let counts = match (range.start, range.end) {
        (Some(ref s), Some(ref e)) => counts.range(s..e),
        (Some(ref s), None) => counts.range(s..),
        (None, Some(ref e)) => counts.range(..e),
        (None, None) => counts.range(..),
    };
    counts
        .fold(
            std::iter::repeat(Reverse(PRESUMED_BARCODE_COUNT))
                .take(NUM_CHUNK_THREADS)
                .collect::<BinaryHeap<_>>(),
            |mut acc, (_, v)| {
                acc.push(Reverse((v.count() as usize).min(limit)));
                while acc.len() > NUM_CHUNK_THREADS {
                    let _ = acc.pop();
                }
                acc
            },
        )
        .into_iter()
        .map(|Reverse(v)| v)
        .sum()
}

/// Choose an aligner based on the chemistry and target_set CSV file.
/// Return STAR if target_set is null.
/// Return the aligner specified by the aligner argument if specified.
/// Return Hurtle if the chemistry is RTL and STAR otherwise.
/// Return Hurtle if the target_set CSV file is a probe set and STAR otherwise.
fn choose_aligner(
    arg_aligner: Option<AlignerParam>,
    is_rtl: Option<bool>,
    target_set: Option<&TargetSetFile>,
) -> Result<AlignerParam> {
    use AlignerParam::{Hurtle, Star};

    let Some(target_set) = target_set else {
        ensure!(
            arg_aligner != Some(Hurtle),
            "aligner=hurtle is incompatible with target_set=null"
        );
        return Ok(Star);
    };

    let is_probe_set = ProbeSetReferenceMetadata::load_from(target_set)?.is_probe_set_metadata();
    let aligner = arg_aligner.unwrap_or(if is_rtl.unwrap_or(is_probe_set) {
        Hurtle
    } else {
        Star
    });

    if aligner == Hurtle {
        ensure!(
            is_probe_set,
            "aligner=hurtle is incompatible with target_set={}",
            target_set.display()
        );
    }

    Ok(aligner)
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
        header_rec.push_tag(b"SN", &format!("PROBE{i}"));
        header_rec.push_tag(b"LN", n_probes + 1);
        header.push_record(&header_rec);
    }
    HeaderView::from_header(&header)
}

fn write_output_bam_header<T: Borrow<HeaderView>>(
    bam_header_data: T,
    rover: &MartianRover,
) -> Result<Option<PathBuf>> {
    let bam_header = rover.make_path("bam_header");
    let mut f = File::create(&bam_header)?;
    f.write_all(bam_header_data.borrow().as_bytes())?;
    Ok(Some(bam_header))
}

#[make_mro(mem_gb = 4, volatile = strict)]
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
        let ref_gb = if args.no_bam && choose_aligner_from_args(&args)? == AlignerParam::Hurtle {
            // For probe assays, STAR is needed only to produce alignments for the BAM file.
            // The memory used by the Hurtle reference is non-zero, but it's less than the
            // wiggle room added to mem_gb below.
            0.0
        } else {
            let settings = star_settings(&args.reference_path, args.transcriptome_min_score)?;
            settings.est_mem()? as f64 / ALN_BC_GIB
        };

        let barcode_counts: BTreeMap<Barcode, CountMetric> =
            args.total_barcode_counts.read()?.into_iter().collect();
        let make_resource =
            |range: &Range<Barcode>, reader: &ShardReader<RnaRead, BarcodeOrder>| -> Result<_> {
                let max_reads = max_items_in_range(&barcode_counts, range, MAX_ITEMS_IN_MEM);
                let max_anns = max_items_in_range(&barcode_counts, range, MAX_ANNOTATIONS_IN_MEM);
                // if we presume ~1KB per RnaRead
                let reads_gb = (max_reads * RNA_READ_SZ) as f64 / ALN_BC_GIB;
                // if we presume ~3KB per ReadAnnotations
                let anns_gb = ((max_anns + BAM_BUFFER_SZ) * ANN_SZ) as f64 / ALN_BC_GIB;
                // account for maximum overhead of shardio reading (ShardIter and a Record)
                let shardio_gb = ((SHARD_SZ + ANN_SZ) * reader.est_len_range(range)
                    / ALN_BC_DISK_CHUNK_SZ) as f64
                    / ALN_BC_GIB;
                // blow up in really extreme cases, this has been tuned against a
                // 13.04GB reference memory reservation
                let mem_gb = (reads_gb + anns_gb + shardio_gb).min(17.96);
                let mem_gb = (mem_gb + ref_gb + 1.5).ceil() as isize;
                Ok(Resource::new()
                    .threads(NUM_CHUNK_THREADS as isize)
                    .mem_gb(mem_gb))
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
        let num_chunks = (n / READS_PER_CHUNK).clamp(1, MAX_ALIGN_CHUNKS_PER_LIB);

        // The shard files are barcode ordered. We are asking the shard reader to chunk the
        // barcode space into `num_chunks` ranges that have roughly the same number of reads.
        // The output we get is a list of barcode ranges.
        let chunks = reader.make_chunks(num_chunks, &Range::all());

        // Pick a sampling rate for AnnotatedReads output. This output is only usef by PD stages.
        let read_ann_subsample_rate = (MAX_READ_ANN_SAMPLE as f32) / (max(n, 1) as f32);
        let read_ann_subsample_rate = if read_ann_subsample_rate > 1.0 {
            1.0
        } else {
            read_ann_subsample_rate
        };

        let stage_def: StageDef<_> = chunks
            .into_iter()
            .map(|range| {
                // Create a chunk that will process all the reads that fall within this barcode range.
                anyhow::Ok((
                    ChunkInputs {
                        range,
                        read_ann_subsample_rate,
                    },
                    make_resource(&range, &reader)?,
                ))
            })
            .try_collect()?;

        // PD uses more memory to produce annotation files.
        Ok(stage_def.join_resource(Resource::with_mem_gb(if args.is_pd { 16 } else { 6 })))
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

        let probe_set_reference = if choose_aligner_from_args(&args)? == AlignerParam::Hurtle {
            Some(ProbeSetReference::from_path(
                args.target_set.as_ref().unwrap(),
                &args.reference_path,
                args.transcriptome_min_score
                    .unwrap_or(DEFAULT_TRANSCRIPTOME_MIN_SCORE),
            )?)
        } else {
            None
        };
        let is_targeted_rtl = probe_set_reference.is_some();

        let star_reference = if is_targeted_rtl && args.no_bam {
            // For probe assays, STAR is needed only to produce alignments for the BAM file.
            None
        } else {
            // load up the aligner and annotator
            // do this first because Orbit loads the reference in the background
            // with mmap, so we want to start this ASAP.
            Some(StarReference::load(star_settings(
                &args.reference_path,
                args.transcriptome_min_score,
            )?)?)
        };

        let reader = args.read_shards.reader()?;

        // Output files
        let counts_bc_order_shard: CountShardFile = rover.make_path("bc_sort");
        let probe_barcode_counts_shard: Option<CountShardFile> = if is_targeted_rtl {
            Some(rover.make_path("bc_probe_sort"))
        } else {
            None
        };
        let bc_umi_info_shard: BcUmiInfoShardFile = rover.make_path("bc_umi_info");
        let pos_sorted_shard: AlignShardFile = rover.make_path("pos_sorted");
        let metrics_shard: BarcodeMetricsShardFile = rover.make_path("metrics_shard");

        // Count data ordered by barcode.
        let mut bc_counts: ShardWriter<FeatureBarcodeCount, BarcodeThenFeatureOrder> =
            ShardWriter::new(
                &counts_bc_order_shard,
                ALN_BC_SEND_BUFFER_SZ,
                ALN_BC_DISK_CHUNK_SZ,
                ALN_BC_ITEM_BUFFER_SZ,
            )?;

        let mut bc_probe_counts: Option<ShardWriter<ProbeBarcodeCount>> =
            probe_barcode_counts_shard.as_ref().map(|x| {
                ShardWriter::new(
                    x,
                    ALN_BC_SEND_BUFFER_SZ,
                    ALN_BC_DISK_CHUNK_SZ,
                    ALN_BC_ITEM_BUFFER_SZ,
                )
                .unwrap()
            });

        // umi info data
        let mut bc_umi_info: ShardWriter<BcUmiInfo, BcUmiInfo> =
            ShardWriter::new(&bc_umi_info_shard, 1, 256, 2048)?;

        // Position sorted read data
        let mut pos_reads: ShardWriter<Record, BamPosSort> = ShardWriter::new(
            &pos_sorted_shard,
            ALN_BC_SEND_BUFFER_SZ,
            ALN_BC_DISK_CHUNK_SZ,
            BAM_BUFFER_SZ,
        )?;

        // Per barcode metrics
        let mut metrics_writer: ShardWriter<BarcodeMetrics, LibFeatThenBarcodeOrder> =
            ShardWriter::new(&metrics_shard, 1024, 8192, 100_000)?;

        // Get the BAM header & write it out.
        let bam_header = if let Some(ref star_reference) = star_reference {
            write_output_bam_header(star_reference.header_view(), &rover)?
        } else if is_targeted_rtl && args.is_pd {
            write_output_bam_header(
                make_dummy_header_for_probes(
                    probe_set_reference.as_ref().unwrap().number_of_probes(),
                ),
                &rover,
            )?
        } else {
            None
        };

        let feature_reference = args.feature_reference.read()?;
        let feature_dist = compute_feature_dist(args.feature_counts.read()?, &feature_reference)?;
        let target_genes = feature_reference.target_genes();

        let annotator = ReadAnnotator::new(
            &args.reference_path,
            args.chemistry_defs.primary(),
            args.include_exons,
            args.include_introns,
        )?;

        let mut subsample_rate = chunk_args.read_ann_subsample_rate;

        // TODO: Don't subsample the the ReadAnnotation output
        // if we're using RTL mode. In RTL mode the ReadAnnotations are used
        // by WRITE_PROBE_INFO_PD to create a CSV with lhs-rhs info, that
        // shouldn't be subsampled. Consider moving the probe_info.csv creation
        // into this stage to avoid needing to write the complete ReadAnnotations
        // in this case.
        if probe_set_reference.is_some() {
            subsample_rate = 1.0;
        }
        let no_star_alignments = star_reference.is_none();
        let aligner = Aligner::new(
            star_reference.map(|x| x.get_aligner()),
            feature_reference,
            feature_dist,
            annotator,
            args.read_chunks.clone(),
            rover.files_path().into(),
            args.targeted_umi_min_read_count,
            probe_set_reference,
        )?;

        let n_threads = rover.get_threads().max(1);
        let ann_files: Vec<ReadAnnotationsFormat> = (0..n_threads)
            .map(|i| rover.make_path(format!("read_annotations_{i}")))
            .collect();

        let (barcodes_to_subsample, barcode_subsample_rate) = get_barcode_subsampling(&args)?;

        // Each thread that handles reads uses one clone() of AlignThreadProcessor struct
        // It contains:
        // - An owned StarAligner that has it's own scratch data for aligning reads
        // - An Arc<StarReference> which is the shared reference data structure used by STAR
        // - Shared FeatureReference and ReadAnnotator structs used to
        let is_pd = args.is_pd;
        let processors = ann_files
            .iter()
            .enumerate()
            .map(|(idx, f)| {
                anyhow::Ok(AlignThreadProcessor {
                    args: args.clone(),
                    aligner: aligner.clone(),
                    bc_umi_info_sender: bc_umi_info.get_sender(),
                    bc_counts_sender: bc_counts.get_sender(),
                    pos_reads_sender: pos_reads.get_sender(),
                    bc_probe_counts_sender: bc_probe_counts.as_mut().map(|x| x.get_sender()),
                    visitor: if is_pd {
                        let rng = ChaCha20Rng::seed_from_u64(idx as u64);
                        StageVisitor::with_ann_writer_sample(
                            metrics_writer.get_sender(),
                            target_genes.clone(),
                            f.lazy_writer()?,
                            subsample_rate,
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
            p.pos_reads_sender.finished()?;
        }

        // shutdown shardio writers
        bc_umi_info.finish()?;
        bc_counts.finish()?;
        bc_probe_counts.as_mut().map(|x| x.finish().unwrap());
        pos_reads.finish()?;
        metrics_writer.finish()?;

        let barcode_summary_shard: BinaryFormat<_, _> = rover.make_path("barcode_summary");
        barcode_summary_shard.write(&bc_summaries)?;

        Ok(ChunkOutputs {
            counts_bc_order_shard,
            probe_barcode_counts_shard,
            bc_umi_info_shard,
            pos_sorted_shard,
            bam_header,
            barcode_summary_shard,
            annotation_files: if is_pd {
                Some(AnnotationFiles {
                    num_reads: num_reads_annotation_files,
                    files: ann_files,
                })
            } else {
                None
            },
            metrics_shard,
            no_star_alignments,
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

        let barcode_summary_filename: CsvFile<_> = rover.make_path("barcode_summary");
        let barcode_summary = {
            let mut barcode_summary: Vec<_> = chunk_outs
                .iter()
                .map(|chunk| chunk.barcode_summary_shard.read())
                .flatten_ok()
                .try_collect()?;
            // Threading causes some flutter in the order.
            barcode_summary.sort();
            barcode_summary
        };
        barcode_summary_filename.write(&barcode_summary)?;

        let summary: JsonFile<AlignAndCountMetrics> = rover.make_path("summary");
        summary.write(&AlignAndCountMetrics {
            alignment_aligner: choose_aligner_from_args(&args)?,
            alignment_high_conf_mapq: HIGH_CONF_MAPQ as i64,
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
                if x.is_empty() {
                    None
                } else {
                    Some(x)
                }
            },
            bc_umi_info: chunk_outs
                .iter()
                .map(|v| v.bc_umi_info_shard.clone())
                .collect(),
            pos_sorted: chunk_outs
                .iter()
                .map(|v| v.pos_sorted_shard.clone())
                .collect(),
            bam_header,
            barcode_summary: barcode_summary_filename,
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
            per_barcode_metrics: chunk_outs.iter().map(|x| x.metrics_shard.clone()).collect(),
            summary,
            no_star_alignments: chunk_outs[0].no_star_alignments,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use cr_bam::bam::BamPosSort;
    use rust_htslib::bam::{HeaderView, Record};
    use shardio::SortKey;

    #[test]
    fn test_choose_aligner() -> Result<()> {
        use cr_types::chemistry::ChemistryName::{SpatialThreePrimeV1, ThreePrimeV3, MFRP_RNA};
        use AlignerParam::{Hurtle, Star};

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
        assert_eq!(choose_aligner(None, ThreePrimeV3.is_rtl(), hybcap)?, Star);
        assert!(choose_aligner(None, MFRP_RNA.is_rtl(), hybcap).is_err());
        assert_eq!(choose_aligner(None, spatial_is_rtl, hybcap)?, Star);

        // Test chemistry and probe_set_file_format.
        assert_eq!(choose_aligner(None, ThreePrimeV3.is_rtl(), rtl)?, Star);
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
        assert_eq!(low_key.0, 4294967295);
        assert_eq!(high_key.0, 4294967295);
        assert_eq!(low_key.1, -1);
        assert_eq!(high_key.1, 2147483646);
    }
}
