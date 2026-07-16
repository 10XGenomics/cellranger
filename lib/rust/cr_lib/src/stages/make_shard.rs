//! Martian stage MAKE_SHARD
#![deny(missing_docs)]

use crate::SequencingMetricsFormat;
use crate::barcode_correction_metrics::{BarcodeCorrectionMetrics, BarcodeCorrectionMetricsFormat};
use crate::barcode_sort::{BarcodeSortResults, execute_barcode_sort_with_visitor};
use crate::make_shard_metrics::{
    MakeShardHistograms, MakeShardMetricProc, MakeShardMetrics, SequencerPlatformMetrics,
};
use crate::types::{ReadPrefixCountFile, ReadShardFile, UmiCountFile};
use anyhow::{Result, bail, ensure};
use barcode::{BarcodeConstruct, BcSegSeq};
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::reference::feature_reference::{
    FeatureConfig, FeatureReference, FeatureReferenceFile, FeatureType, TargetGeneIndicesFile,
};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::ReferenceInfo;
use cr_types::rna_read::RnaChunk;
use cr_types::types::{BcSegmentCountFormat, FeatureBarcodeType, LibraryType, VdjChainType};
use cr_types::{
    ERROR_CODE_INFO, FeatureCountFormat, FeatureReferenceFormat, MetricsFile,
    PerLibrarySortedBarcodeCounts, SortedBarcodeCountFile,
    write_sorted_count_file_from_unordered_histogram,
};
use fastq_set::read_pair::{ReadHeader, WhichRead};
use itertools::{Itertools, zip_eq};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{Histogram, JsonReport, JsonReporter, Metric, SimpleHistogram, TxHashMap, TxHashSet};
use multi::config::multiconst::{FEATURE, GENE_EXPRESSION, VDJ};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::PathBuf;
use std::str::FromStr;

// Report up to 50 unknown feature barcode sequences which occupies at least 0.1% library reads
const MAX_UNKNOWN_FBC_TO_REPORT: usize = 50;
const MIN_FRAC_READS_UNKNOWN_FBC: f64 = 0.001;

martian_filetype!(MakeShardMetricsFile, "msm");

martian_filetype!(MakeShardHistogramsFile, "msh");

#[derive(Deserialize, Clone, MartianStruct)]
pub struct MakeShardStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub read_chunks: JsonFile<Vec<RnaChunk>>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub subsample_rate: Option<f64>,
    pub initial_reads: Option<usize>,
    pub reference_info: Option<ReferenceInfo>,
    pub feature_reference: Option<FeatureReferenceFile>,
    pub target_features: Option<TargetGeneIndicesFile>,
    pub target_set: Option<TargetSetFile>,
    pub target_set_name: Option<String>,
    pub feature_config: Option<FeatureConfig>,
    pub umi_min_read_length: Option<usize>,
    /// If true, do not include read IDs in the output data.
    pub skip_read_id: bool,
    /// VDJ receptor string (e.g. "VDJ-T", "VDJ-B"). When set, chain-type-specific
    /// platform metric keys are emitted alongside the generic "VDJ_" prefix keys.
    pub receptor: Option<String>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeShardChunkInputs {
    chunk_id: usize,
    feature_reference_binary: Option<FeatureReferenceFormat>,
    use_multithread_bgzf: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeShardChunkOutputs {
    valid_shard: ReadShardFile,
    invalid_shard: ReadShardFile,
    feature_counts: FeatureCountFormat,
    read_prefix_counts: ReadPrefixCountFile,
    umi_counts: UmiCountFile,
    chunk_summary: BinaryFormat<MakeShardMetricsFile, MakeShardMetrics>,
    chunk_hist: BinaryFormat<MakeShardHistogramsFile, MakeShardHistograms>,
    valid_bc_counts: SortedBarcodeCountFile,
    bc_correct_summary: BarcodeCorrectionMetricsFormat,
    total_read_pairs: i64,
}

#[derive(Clone, Serialize, MartianStruct)]
pub struct MakeShardStageOutputs {
    pub valid: Vec<ReadShardFile>,
    pub invalid: Vec<ReadShardFile>,
    pub barcode_counts: PerLibrarySortedBarcodeCounts,
    pub barcode_segment_counts: BcSegmentCountFormat,
    pub feature_counts: FeatureCountFormat,
    pub summary: MetricsFile,
    // this is for VDJ, because exclusive can collapse or choose VDJ
    pub total_read_pairs: i64,
    pub feature_reference_binary: Option<FeatureReferenceFormat>,
    /// The barcode correction stage only iterates through reads with invalid barcodes.
    /// This means that the valid read counts need to be supplied by this stage in order to
    /// compute metrics such as "corrected_bc_frac", "good_bc_frac".
    pub bc_correct_summary: BarcodeCorrectionMetricsFormat,
    pub sequencing_metrics: SequencingMetricsFormat,
}

/// Read the FASTQ files and extract the barcode and UMI sequences.
pub struct MakeShard;

#[make_mro(mem_gb = 2, volatile = strict)]
impl MartianStage for MakeShard {
    type StageInputs = MakeShardStageInputs;
    type StageOutputs = MakeShardStageOutputs;
    type ChunkInputs = MakeShardChunkInputs;
    type ChunkOutputs = MakeShardChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let reference_path = args
            .reference_info
            .as_ref()
            .and_then(|x| x.get_reference_path());

        let feature_reference = if reference_path.is_some()
            || args.target_set.is_some()
            || args.feature_reference.is_some()
        {
            let feature_reference: FeatureReferenceFormat = rover.make_path("feature_reference");
            feature_reference.write(&FeatureReference::from_paths(
                reference_path,
                args.feature_reference.as_ref(),
                args.target_set_name.as_deref(),
                args.target_set.as_ref(),
                args.target_features.as_ref(),
                args.feature_config.as_ref(),
            )?)?;
            Some(feature_reference)
        } else {
            None
        };

        let read_chunks = args.read_chunks.read()?;
        Ok(read_chunks
            .iter()
            .enumerate()
            .map(|(chunk_id, chunk)| -> Result<_> {
                let use_multithread_bgzf = chunk.fastqs.is_bgzf()?;
                Ok((
                    MakeShardChunkInputs {
                        chunk_id,
                        feature_reference_binary: feature_reference.clone(),
                        use_multithread_bgzf,
                    },
                    Resource::with_mem_gb(9).threads(if use_multithread_bgzf { 8 } else { 4 }),
                ))
            })
            .collect::<Result<StageDef<_>>>()?
            .join_resource(Resource::with_mem_gb(13)))
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let id = chunk_args.chunk_id;
        let read_chunks = args.read_chunks.read()?;

        // Apply overrides to the chunk and chemistry.
        let (rna_chunk, chemistry) = {
            let mut chunk = read_chunks[id].clone();
            let mut chemistry = args.chemistry_defs[&chunk.library_type].clone();
            if let Some(subsample_rate) = args.subsample_rate {
                chunk.set_subsample_rate(subsample_rate);
            }
            if let Some(r1_length) = args.r1_length {
                chunk.set_illumina_r1_trim_length(r1_length);
            }
            if let Some(r2_length) = args.r2_length {
                chunk.set_illumina_r2_trim_length(r2_length);
            }
            if let Some(umi_min_read_length) = args.umi_min_read_length {
                for umi in &mut chemistry.umi {
                    if umi.length > umi_min_read_length {
                        umi.min_length = Some(umi_min_read_length);
                    }
                }
            }
            (chunk, chemistry)
        };

        let proc = MakeShardMetricProc::new(
            chunk_args
                .feature_reference_binary
                .as_ref()
                .map(FileTypeRead::read)
                .transpose()?,
        )?;

        // Output filenames
        let valid_shard: ReadShardFile = rover.make_path("valid");
        let invalid_shard: ReadShardFile = rover.make_path("invalid");
        let feature_counts: FeatureCountFormat = rover.make_path("feature_counts");
        let read_prefix_counts: ReadPrefixCountFile = rover.make_path("read_prefix_counts");
        let umi_counts: UmiCountFile = rover.make_path("umi_counts");
        let chunk_summary: BinaryFormat<_, _> = rover.make_path("chunk_summary");
        let chunk_hist: BinaryFormat<_, _> = rover.make_path("chunk_hist");

        let initial_reads = args.initial_reads.map(|i| {
            let frac = i as f64 / read_chunks.len() as f64;
            let lwr = (frac * id as f64).round() as usize;
            let upr = (frac * (1 + id) as f64).round() as usize;
            upr - lwr
        });

        let BarcodeSortResults {
            metrics,
            valid_bc_counts,
            mut accum,
        } = execute_barcode_sort_with_visitor(
            rna_chunk.clone(),
            chemistry,
            valid_shard.to_path_buf(),
            invalid_shard.to_path_buf(),
            initial_reads,
            &proc,
            if args.skip_read_id {
                ReadHeader::SkipId
            } else {
                ReadHeader::KeepId
            },
            chunk_args.use_multithread_bgzf,
        )?;

        accum.metrics.platform = Some(infer_sequencing_platform(&rna_chunk)?);

        feature_counts.write(&accum.feature_counts)?;
        chunk_summary.write(&accum.metrics)?;
        chunk_hist.write(&accum.histograms)?;

        let valid_bc_counts = write_sorted_count_file_from_unordered_histogram(
            &valid_bc_counts,
            &rover.make_path::<PathBuf>("valid_bc_counts"),
        )?;

        // For the per-barcode segment metrics good_bc_in_*_frac
        // - MAKE_SHARD counts the reads with a valid barcode, that is, all segments are valid.
        // - BARCODE_CORRECTION counts the reads with an invalid barcode, that is, any segment is invalid.
        // Counting which individual segments are valid/invalid here in MAKE_SHARD would double count the reads with any invalid segment.
        let num_valid_barcodes = metrics.valid_items;
        let num_valid_barcode_segments = metrics.valid_items_in.map(|_x| num_valid_barcodes.into());
        let bc_correct_summary: BarcodeCorrectionMetricsFormat =
            rover.make_path("bc_correct_summary");
        bc_correct_summary.write(&BarcodeCorrectionMetrics::with_valid(
            rna_chunk.library_type(),
            num_valid_barcodes,
            num_valid_barcode_segments,
        ))?;

        let total_read_pairs = metrics.valid_items + metrics.invalid_items;

        ensure!(
            total_read_pairs > 0,
            "There were no reads to process. Please check whether the read lengths in the input \
             fastqs satisfy the minumum read length requirements for the chemistry."
        );

        Ok(MakeShardChunkOutputs {
            total_read_pairs,
            valid_shard,
            invalid_shard,
            feature_counts,
            read_prefix_counts,
            umi_counts,
            chunk_summary,
            chunk_hist,
            valid_bc_counts,
            bc_correct_summary,
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let feature_reference_binary = chunk_defs
            .into_iter()
            .next()
            .and_then(|x| x.feature_reference_binary);

        let read_chunks = args.read_chunks.read()?;
        let feature_counts_file: FeatureCountFormat = rover.make_path("feature_counts");
        let mut feature_counts = Vec::new();
        for co in &chunk_outs {
            feature_counts.merge(co.feature_counts.read()?);
        }
        feature_counts_file.write(&feature_counts)?;
        let feature_metrics = if let Some(feature_reference_binary) = &feature_reference_binary {
            Some(calculate_feature_metrics(
                &feature_counts,
                &feature_reference_binary.read()?,
            ))
        } else {
            None
        };
        drop(feature_counts);

        let bc_correct_metrics: BarcodeCorrectionMetrics = chunk_outs
            .iter()
            .map(|x| x.bc_correct_summary.read())
            .sum::<Result<_>>()?;
        let bc_correct_summary: BarcodeCorrectionMetricsFormat =
            rover.make_path("bc_correct_summary");
        bc_correct_summary.write(&bc_correct_metrics)?;
        drop(bc_correct_metrics);

        let total_read_pairs = chunk_outs.iter().map(|co| co.total_read_pairs).sum();

        let mut per_lib_per_fastq_id_metrics: TxHashMap<_, TxHashMap<_, MakeShardMetrics>> =
            TxHashMap::default();
        let mut per_lib_hist = TxHashMap::default();

        for (read_chunk, chunk_out) in zip_eq(&read_chunks, &chunk_outs) {
            let library_type = read_chunk.library_type();
            let barcode_construct = args.chemistry_defs[&library_type]
                .barcode_construct()
                .map(|_| ());

            let metric = chunk_out.chunk_summary.read()?;
            per_lib_per_fastq_id_metrics
                .entry(library_type)
                .or_default()
                .entry(read_chunk.fastq_id.as_ref())
                .or_default()
                .merge(metric);

            let hist = chunk_out.chunk_hist.read()?;
            per_lib_hist
                .entry(library_type)
                .or_insert_with(|| MakeShardHistograms::new(barcode_construct))
                .merge(hist);
        }

        let lib_types: Vec<_> = per_lib_hist.keys().copied().collect();

        let PerLibHistograms {
            valid_bc_segment_counts,
            per_lib_r1_lengths,
            unknown_feature_bcs,
        } = decompose_per_lib_histograms(per_lib_hist);

        let top_unknown_fbc_seqs: TxHashMap<_, _> = unknown_feature_bcs
            .into_iter()
            .map(|(k, v)| (k, v.top_n(MAX_UNKNOWN_FBC_TO_REPORT)))
            .collect();

        let mut per_lib_seq_metrics = TxHashMap::default();
        let mut per_lib_metrics = TxHashMap::default();
        for (library_id, (library_type, per_fastq_id_metrics)) in
            (0_u16..).zip(per_lib_per_fastq_id_metrics.into_iter())
        {
            let mut lib_metrics = MakeShardMetrics::default();
            let mut seq_metrics_rows = Vec::new();
            for (fastq_id, metrics) in per_fastq_id_metrics {
                seq_metrics_rows.push(
                    metrics.sequencing_metrics_for_fastq(fastq_id.cloned().unwrap_or_default()),
                );
                lib_metrics.merge(metrics);
                if let Some(feature_bc_hist) = top_unknown_fbc_seqs.get(&library_type) {
                    lib_metrics.report_unknown_fbc(
                        library_id,
                        feature_bc_hist,
                        MIN_FRAC_READS_UNKNOWN_FBC,
                    );
                }
            }
            per_lib_metrics.insert(library_type, lib_metrics);
            per_lib_seq_metrics.insert(library_type, seq_metrics_rows);
        }

        let sequencing_metrics_file: SequencingMetricsFormat =
            rover.make_path("sequencing_metrics");
        sequencing_metrics_file.write(&per_lib_seq_metrics)?;

        let barcode_segment_counts = rover
            .make_path::<BcSegmentCountFormat>("barcode_segment_counts")
            .with_content(&valid_bc_segment_counts)?;
        drop(valid_bc_segment_counts);

        // BARCODE_CORRECTION may exceed available memory when a single barcode
        // has a very large number of reads. Print a warning message.
        for (library, metrics) in &per_lib_metrics {
            let frac = metrics
                .homopolymer_barcode_property
                .fraction()
                .unwrap_or_default();
            if frac >= 0.1 {
                eprintln!("Warning: {library}: homopolymer_barcode_property_frac is high: {frac}");
            }
        }

        // Also if r1_lengths are different, abort the run and complain
        let mut r1_lengths = SimpleHistogram::default();
        for this_lib_r1_len in per_lib_r1_lengths.into_values() {
            r1_lengths.merge(this_lib_r1_len);
        }
        if r1_lengths.len() > 1 {
            let r1_min_observed = r1_lengths.min_key().unwrap();
            let r1_max_observed = r1_lengths.max_key().unwrap();
            let umi_max_length = read_chunks
                .iter()
                .flat_map(|chunk| {
                    args.chemistry_defs[&chunk.library_type]
                        .umi
                        .iter()
                        .filter(|x| x.read_type == WhichRead::R1)
                        .map(|x| x.offset + x.min_length.unwrap_or(x.length))
                })
                .max()
                .unwrap_or(0);

            ensure!(
                r1_min_observed >= umi_max_length,
                "TXRNGR10009: We detected a mixture of different R1 lengths ([{min}-{max}]), which \
                 breaks assumptions in how UMIs are tabulated and corrected. To process these \
                 data, you will need to truncate to the shortest observed R1 length by providing \
                 {min} to the --r1-length argument if you are running count/vdj, via the \
                 --rna-r1-length parameter if you are running Multiome, or via the r1-length \
                 parameter in each of the following tables of your multi config CSV if you are \
                 running multi: [{tbl}]. {ERROR_CODE_INFO}",
                min = r1_min_observed,
                max = r1_max_observed,
                tbl = lib_types
                    .iter()
                    .map(|&key| match key {
                        LibraryType::Atac => unreachable!(),
                        LibraryType::FeatureBarcodes(_) => FEATURE,
                        LibraryType::Gex => GENE_EXPRESSION,
                        LibraryType::Vdj(_) => VDJ,
                    })
                    .format("], ["),
            );
        }

        let mut per_lib_reporter = per_lib_metrics.to_json_reporter();

        // HERE BE DRAGONS
        //   there are oddities to the way metrics are prefixed and reported,
        //   these oddities are implemented here.
        for library_id in read_chunks.iter().map(RnaChunk::library_id) {
            // VDJ is expected to be disjoint from GEX/FB and so these metrics are expected unprefixed
            if let Some(value) =
                per_lib_reporter.remove(&format!("VDJ_{library_id}_total_read_pairs_per_library"))
            {
                per_lib_reporter
                    .insert(format!("{library_id}_total_read_pairs_per_library"), value);
            }
        }
        if let Some(receptor) = &args.receptor
            && let Ok(LibraryType::Vdj(chain_type)) = LibraryType::from_str(receptor)
            && chain_type != VdjChainType::Auto
        {
            let chain_prefix = chain_type.as_str().replace('-', "_");
            for metric in ["inferred_platform", "fastq_id_components"] {
                let vdj_key = format!("VDJ_{metric}");
                let chain_key = format!("{chain_prefix}_{metric}");
                if let Some(value) = per_lib_reporter.get(&vdj_key).cloned() {
                    per_lib_reporter.insert(chain_key, value);
                }
            }
        }
        // SO SAY THE DRAGONS

        let summary = JsonReporter::from(("cellranger_version", rover.pipelines_version()))
            + per_lib_reporter
            + feature_metrics.unwrap_or_default()
            + args.chemistry_defs.primary().to_json_reporter();

        Ok(MakeShardStageOutputs {
            valid: chunk_outs.iter().map(|x| x.valid_shard.clone()).collect(),
            invalid: chunk_outs.iter().map(|x| x.invalid_shard.clone()).collect(),
            // Merge the valid barcode counts.
            barcode_counts: PerLibrarySortedBarcodeCounts::merge_from_iter(
                chunk_outs
                    .into_iter()
                    .zip(&read_chunks)
                    .map(|(chunk_out, read_chunk)| {
                        (read_chunk.library_type, chunk_out.valid_bc_counts)
                    }),
                |lib_type| rover.make_path(format!("{lib_type}_valid_barcode_counts")),
            )?,
            barcode_segment_counts,
            feature_counts: feature_counts_file,
            summary: MetricsFile::from_reporter(&rover, "summary", &summary)?,
            total_read_pairs,
            feature_reference_binary,
            bc_correct_summary,
            sequencing_metrics: sequencing_metrics_file,
        })
    }
}

#[derive(Default)]
struct PerLibHistograms {
    valid_bc_segment_counts: TxHashMap<LibraryType, BarcodeConstruct<SimpleHistogram<BcSegSeq>>>,
    per_lib_r1_lengths: TxHashMap<LibraryType, SimpleHistogram<usize>>,
    unknown_feature_bcs: TxHashMap<LibraryType, SimpleHistogram<String>>,
}

/// Decompose the per-library histograms into individual per-library maps by transferring ownership.
fn decompose_per_lib_histograms(
    per_lib_histos: TxHashMap<LibraryType, MakeShardHistograms>,
) -> PerLibHistograms {
    let (valid_bc_segment_counts, per_lib_r1_lengths, unknown_feature_bcs) = per_lib_histos
        .into_iter()
        .map(|(lib_type, x)| {
            (
                (lib_type, x.valid_bc_segment_counts),
                (lib_type, x.r1_lengths),
                (lib_type, x.unknown_feature_bcs),
            )
        })
        .multiunzip();
    PerLibHistograms {
        valid_bc_segment_counts,
        per_lib_r1_lengths,
        unknown_feature_bcs,
    }
}

/// Read the hashtag barcode sequences.
fn read_hashtag_sequences() -> TxHashSet<String> {
    // id,name,read,pattern,sequence,feature_type
    const SEQUENCE_COLUMN_INDEX: usize = 4;

    [
        include_str!("data/PTG_HashMax.csv"),
        include_str!("data/TotalSeq_Hashtag.csv"),
    ]
    .into_iter()
    .flat_map(|csv_str| {
        let mut reader = csv::Reader::from_reader(csv_str.as_bytes());
        assert_eq!(
            &reader.headers().unwrap()[SEQUENCE_COLUMN_INDEX],
            "sequence"
        );
        reader
            .records()
            .map(|record: csv::Result<csv::StringRecord>| {
                record.unwrap()[SEQUENCE_COLUMN_INDEX].to_string()
            })
            .collect::<Vec<_>>()
    })
    .collect()
}

/// Calculate metrics on feature counts.
/// `detected_hashtag_feature_count` is the number of features with known
/// commercial antibody hashtag barcode sequences and that are detected to be
/// present. A hashtag oligo whose proportion of reads of all hashtag oligos
/// is at least 0.1% is considered to be detected.
/// Return no metric for no hashtag feature types in the feature reference.
/// Return null for no hashtag feature sequences in the feature reference.
/// Return 0 for no hashtag feature reads detected in the data.
fn calculate_feature_metrics(
    feature_counts: &[i64],
    feature_reference: &FeatureReference,
) -> JsonReporter {
    use FeatureBarcodeType::{Antibody, Custom, Multiplexing};
    const DETECTED_HASHTAG_THRESHOLD: f64 = 0.001;

    // Consider only feature types that could potentially be antibody hashtags.
    let is_possible_hashtag_feature_type = |feature: &FeatureType| {
        matches!(
            feature,
            FeatureType::Barcode(Antibody | Custom | Multiplexing)
        )
    };

    if !feature_reference
        .feature_defs
        .iter()
        .any(|feature| is_possible_hashtag_feature_type(&feature.feature_type))
    {
        return JsonReporter::default();
    }

    let detected_hashtag_feature_count = 'detected_hashtag_feature_count: {
        let hashtag_sequences = read_hashtag_sequences();
        let hashtag_counts: Vec<i64> = zip_eq(&feature_reference.feature_defs, feature_counts)
            .filter(|(feature, _count)| {
                is_possible_hashtag_feature_type(&feature.feature_type)
                    && hashtag_sequences.contains(&feature.sequence)
            })
            .map(|(_feature, &count)| count)
            .collect();
        if hashtag_counts.is_empty() {
            break 'detected_hashtag_feature_count Value::Null;
        }

        let hashtag_count_sum: i64 = hashtag_counts.iter().sum();
        if hashtag_count_sum == 0 {
            break 'detected_hashtag_feature_count Value::from(0);
        }

        let detected_hashtag_feature_count = hashtag_counts
            .into_iter()
            .filter(|&count| count as f64 / hashtag_count_sum as f64 >= DETECTED_HASHTAG_THRESHOLD)
            .count();
        assert_ne!(detected_hashtag_feature_count, 0);
        Value::from(detected_hashtag_feature_count)
    };

    [
        (
            "detected_hashtag_feature_count",
            detected_hashtag_feature_count,
        ),
        (
            "hashtag_reference_version",
            Value::from(HASHTAG_REFERENCE_VERSION),
        ),
    ]
    .into_iter()
    .collect()
}

/// Attempt to infer the sequencing platform from the first read ID.
fn infer_sequencing_platform(rna_chunk: &RnaChunk) -> Result<SequencerPlatformMetrics> {
    let Some(read) = rna_chunk.read_iter(ReadHeader::KeepId, None)?.next() else {
        bail!("could not infer sequencing platform, no reads in RNA chunk");
    };
    Ok(SequencerPlatformMetrics::from_fastq_id(
        read?.id().unwrap_or_default(),
    ))
}

/// The hashtag reference version.
/// Update this version when the hashtag reference sequences are changed.
const HASHTAG_REFERENCE_VERSION: i64 = 1;

#[cfg(test)]
mod tests {
    use super::*;
    use metric::TxHasher;
    use std::hash::{Hash, Hasher};

    #[test]
    /// Ensure that the hashtag sequences are not changed unintentionally.
    /// Update `HASHTAG_REFERENCE_VERSION` when the hashtag reference sequences are changed.
    fn test_read_hashtag_sequences() {
        let sequences = read_hashtag_sequences();
        assert_eq!(sequences.len(), 80);
        let hash_value = sequences
            .iter()
            .sorted()
            .fold(TxHasher::hasher(), |mut hasher, seq| {
                seq.hash(&mut hasher);
                hasher
            })
            .finish();
        assert_eq!(hash_value, 0x52ce28113d0e57b8);
    }
}
