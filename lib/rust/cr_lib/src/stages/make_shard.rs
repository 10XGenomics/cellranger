//! Martian stage MAKE_SHARD

use crate::barcode_correction_metrics::{BarcodeCorrectionMetrics, BarcodeCorrectionMetricsFormat};
use crate::barcode_sort::execute_barcode_sort_with_visitor;
use crate::make_shard_metrics::{MakeShardHistograms, MakeShardMetrics, MakeShardVisitor};
use crate::types::{FeatureReferenceFormat, ReadPrefixCountFile, ReadShardFile, UmiCountFile};
use crate::SequencingMetricsFormat;
use anyhow::{ensure, Result};
use barcode::Whitelist;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::reference::feature_reference::{
    FeatureConfig, FeatureReference, FeatureReferenceFile, FeatureType, TargetGeneIndicesFile,
};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::rna_read::RnaChunk;
use cr_types::types::{
    BcCountFormat, BcSegmentCountFormat, FeatureBarcodeType, GemWell, LibraryType,
};
use cr_types::{FeatureCountFormat, MetricsFile};
use fastq_set::read_pair::WhichRead;
use itertools::{zip_eq, Itertools};
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{JsonReport, JsonReporter, Metric, SimpleHistogram, TxHashMap, TxHashSet};
use multi::config::multiconst;
use parameters_toml::umi_min_read_length;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

// Report up to 50 unknown feature barcode sequences which occupies at least 0.1% library reads
const MAX_UNKNOWN_FBC_TO_REPORT: usize = 50;
const MIN_FRAC_READS_UNKNOWN_FBC: f64 = 0.001;

martian_filetype!(MakeShardMetricsFile, "msm");

martian_filetype!(MakeShardHistogramsFile, "msh");

#[derive(Deserialize, Clone, MartianStruct)]
pub struct MakeShardStageInputs {
    pub chemistry_defs: ChemistryDefs,
    // Make shard works on data from a single gem well
    pub gem_well: GemWell,
    pub read_chunks: Vec<RnaChunk>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub subsample_rate: Option<f64>,
    pub initial_read_pairs: Option<usize>,
    pub reference_path: Option<PathBuf>,
    pub feature_reference_path: Option<FeatureReferenceFile>,
    pub target_features: Option<TargetGeneIndicesFile>,
    pub target_set: Option<TargetSetFile>,
    pub target_set_name: Option<String>,
    pub feature_config: Option<FeatureConfig>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeShardChunkInputs {
    chunk_id: usize,
    feature_reference: Option<FeatureReferenceFormat>,
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
    bc_correct_summary: BarcodeCorrectionMetricsFormat,
    total_read_pairs: i64,
}

#[derive(Clone, Serialize, MartianStruct)]
pub struct MakeShardStageOutputs {
    pub valid: Vec<ReadShardFile>,
    pub invalid: Vec<ReadShardFile>,
    pub barcode_counts: BcCountFormat,
    pub barcode_segment_counts: BcSegmentCountFormat,
    pub feature_counts: FeatureCountFormat,
    pub summary: MetricsFile,
    // this is for VDJ, because exclusive can collapse or choose VDJ
    pub total_read_pairs: i64,
    pub feature_reference: Option<FeatureReferenceFormat>,
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
        // Ensure that we got data only from the specified gem well
        ensure!(
            args.read_chunks
                .iter()
                .all(|chunk| GemWell(chunk.gem_group()) == args.gem_well),
            "You have provided data from multiple gem wells, which is not supported in \
             count, vdj and multi pipelines. To perform analysis on the combined data, run \
             data from individual GEM wells using one of count/vdj/multi pipelines and \
             aggregate the outputs using the aggr pipeline."
        );

        // TODO: rather than special-casing antibody, perhaps only
        // generate features for input library types?
        let antibody_only = {
            let mut r = TxHashSet::default();
            r.insert(FeatureType::Barcode(FeatureBarcodeType::Antibody));
            r
        };
        let use_feature_types = if args.read_chunks.iter().all(|chunk| {
            chunk
                .library_type()
                .is_fb_type(FeatureBarcodeType::Antibody)
        }) {
            Some(antibody_only)
        } else {
            None
        };

        let feature_reference = if args.reference_path.is_some()
            || args.target_set.is_some()
            || args.feature_reference_path.is_some()
        {
            let feature_reference: FeatureReferenceFormat = rover.make_path("feature_reference");
            feature_reference.write(&FeatureReference::from_paths(
                args.reference_path.as_deref(),
                args.feature_reference_path.as_ref(),
                use_feature_types,
                args.target_set_name.as_deref(),
                args.target_set.as_ref(),
                args.target_features.as_ref(),
                args.feature_config.as_ref(),
            )?)?;
            Some(feature_reference)
        } else {
            None
        };

        Ok((0..args.read_chunks.len())
            .map(|chunk_id| {
                (
                    MakeShardChunkInputs {
                        chunk_id,
                        feature_reference: feature_reference.clone(),
                    },
                    Resource::with_mem_gb(4).threads(4),
                )
            })
            .collect::<StageDef<_>>()
            .join_resource(Resource::with_mem_gb(9)))
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let id = chunk_args.chunk_id;

        let rna_chunk = {
            let mut chunk = args.read_chunks[id].clone();
            if let Some(subsample_rate) = args.subsample_rate {
                chunk.set_subsample_rate(subsample_rate);
            }
            if let Some(r1_length) = args.r1_length {
                chunk.set_illumina_r1_trim_length(r1_length);
            }
            if let Some(r2_length) = args.r2_length {
                chunk.set_illumina_r2_trim_length(r2_length);
            }
            if let Some(umi_min_read_length) = *umi_min_read_length()? {
                for umi in &mut chunk.chemistry.umi {
                    if umi.length > umi_min_read_length {
                        umi.min_length = Some(umi_min_read_length);
                    }
                }
            }
            chunk
        };

        let mut visitor = MakeShardVisitor::new(
            chunk_args
                .feature_reference
                .as_ref()
                .map(martian_filetypes::FileTypeRead::read)
                .transpose()?,
            &args.read_chunks,
            id,
            &rna_chunk.chemistry,
        )?;

        // Output filenames
        let valid_shard: ReadShardFile = rover.make_path("valid");
        let invalid_shard: ReadShardFile = rover.make_path("invalid");
        let feature_counts: FeatureCountFormat = rover.make_path("feature_counts");
        let read_prefix_counts: ReadPrefixCountFile = rover.make_path("read_prefix_counts");
        let umi_counts: UmiCountFile = rover.make_path("umi_counts");
        let chunk_summary: BinaryFormat<_, _> = rover.make_path("chunk_summary");
        let chunk_hist: BinaryFormat<_, _> = rover.make_path("chunk_hist");

        let initial_read_pairs = args.initial_read_pairs.map(|i| {
            let frac = i as f64 / args.read_chunks.len() as f64;
            let lwr = (frac * id as f64).round() as usize;
            let upr = (frac * (1 + id) as f64).round() as usize;
            upr - lwr
        });

        let metrics = execute_barcode_sort_with_visitor(
            rna_chunk.clone(),
            valid_shard.to_path_buf(),
            invalid_shard.to_path_buf(),
            Whitelist::construct(rna_chunk.chemistry.barcode_whitelist())?,
            initial_read_pairs,
            &mut visitor,
        )?;
        feature_counts.write(&visitor.feature_counts)?;
        chunk_summary.write(visitor.metrics())?;
        chunk_hist.write(visitor.histograms())?;

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
        let feature_counts_file: FeatureCountFormat = rover.make_path("feature_counts");

        let mut feature_counts = Vec::new();
        for co in &chunk_outs {
            feature_counts.merge(co.feature_counts.read()?);
        }
        feature_counts_file.write(&feature_counts)?;

        let bc_correct_metrics: BarcodeCorrectionMetrics = chunk_outs
            .iter()
            .map(|x| x.bc_correct_summary.read())
            .sum::<Result<_>>()?;
        let bc_correct_summary: BarcodeCorrectionMetricsFormat =
            rover.make_path("bc_correct_summary");
        bc_correct_summary.write(&bc_correct_metrics)?;

        let total_read_pairs = chunk_outs.iter().map(|co| co.total_read_pairs).sum();

        let mut per_lib_per_fastq_id_metrics: TxHashMap<_, TxHashMap<_, MakeShardMetrics>> =
            TxHashMap::default();
        let mut per_lib_hist = TxHashMap::default();

        for (read_chunk, chunk_out) in zip_eq(&args.read_chunks, &chunk_outs) {
            let library_type = read_chunk.library_type();
            let barcode_construct = read_chunk.chemistry.barcode_construct().map(|_| ());

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

        let top_unknown_fbc_seqs: Vec<(_, SimpleHistogram<String>)> = per_lib_hist
            .iter()
            .map(|(&k, v)| {
                (
                    k,
                    v.unknown_feature_bcs
                        .top_n(MAX_UNKNOWN_FBC_TO_REPORT)
                        .clone(),
                )
            })
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
                if let Some((_, feature_bc_hist)) = top_unknown_fbc_seqs
                    .iter()
                    .find(|(lib_type, _)| *lib_type == library_type)
                {
                    lib_metrics.report_unknown_fbc(
                        library_id,
                        feature_bc_hist.clone(),
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

        let barcode_counts: BcCountFormat = rover.make_path("barcode_counts");
        barcode_counts.write(
            &per_lib_hist
                .iter()
                .map(|(&k, v)| (k, v.valid_bc_counts.clone()))
                .collect(),
        )?;

        let barcode_segment_counts: BcSegmentCountFormat =
            rover.make_path("barcode_segment_counts");
        barcode_segment_counts.write(
            &per_lib_hist
                .iter()
                .map(|(&k, v)| (k, v.valid_bc_segment_counts.clone()))
                .collect(),
        )?;

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
        for histograms in per_lib_hist.values() {
            r1_lengths.merge(histograms.r1_lengths.clone());
        }
        if r1_lengths.distribution().len() > 1 {
            let r1_min_length = r1_lengths.min_key().unwrap();
            let r1_max_length = r1_lengths.max_key().unwrap();
            let umi_max_length = args
                .read_chunks
                .iter()
                .flat_map(|chunk| {
                    chunk
                        .chemistry
                        .umi
                        .iter()
                        .filter(|x| x.read_type == WhichRead::R1)
                        .map(|x| x.offset + x.length)
                })
                .max()
                .unwrap_or(0);

            ensure!(
                r1_min_length >= umi_max_length,
                "We detected a mixture of different R1 lengths ([{min}-{max}]), which \
                 breaks assumptions in how UMIs are tabulated and corrected. To process these \
                 data, you will need to truncate to the shortest observed R1 length by providing \
                 {min} to the --r1-length argument if are running count/vdj, or via the r1-length \
                 parameter in each of the following tables of your multi config CSV if you are \
                 running multi: [{tbl}]",
                min = r1_min_length,
                max = r1_max_length,
                tbl = per_lib_hist
                    .keys()
                    .map(|&key| match key {
                        LibraryType::Gex => multiconst::GENE_EXPRESSION,
                        LibraryType::Vdj(_) => multiconst::VDJ,
                        LibraryType::FeatureBarcodes(_) => multiconst::FEATURE,
                        LibraryType::Atac => unreachable!(),
                    })
                    .format("], ["),
            );
        }

        let mut per_lib_reporter = per_lib_metrics.to_json_reporter();

        // HERE BE DRAGONS
        //   there are oddities to the way metrics are prefixed and reported,
        //   these oddities are implemented here.
        for library_id in args.read_chunks.iter().map(RnaChunk::library_id) {
            // VDJ is expected to be disjoint from GEX/FB and so these metrics are expected unprefixed
            if let Some(value) =
                per_lib_reporter.remove(&format!("VDJ_{library_id}_total_read_pairs_per_library"))
            {
                per_lib_reporter
                    .insert(format!("{library_id}_total_read_pairs_per_library"), value);
            }
        }
        // SO SAY THE DRAGONS

        let summary = JsonReporter::from(("cellranger_version", rover.pipelines_version()))
            + per_lib_reporter
            + args.chemistry_defs.primary().to_json_reporter()
            + JsonReporter::from(("chemistry_defs", args.chemistry_defs.to_value()));

        let feature_reference = chunk_defs.first().and_then(|x| x.feature_reference.clone());

        Ok(MakeShardStageOutputs {
            valid: chunk_outs.iter().map(|x| x.valid_shard.clone()).collect(),
            invalid: chunk_outs.iter().map(|x| x.invalid_shard.clone()).collect(),
            barcode_counts,
            barcode_segment_counts,
            feature_counts: feature_counts_file,
            summary: MetricsFile::from_reporter(&rover, "summary", &summary)?,
            total_read_pairs,
            feature_reference,
            bc_correct_summary,
            sequencing_metrics: sequencing_metrics_file,
        })
    }
}
