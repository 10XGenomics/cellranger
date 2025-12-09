//! Martian stage RUST_BRIDGE
#![deny(missing_docs)]

use crate::barcode_sort::BarcodeOrder;
use crate::types::ReadShardFile;
use anyhow::Result;
use barcode::Barcode;
use cr_types::chemistry::ChemistryDefs;
use cr_types::rna_read::RnaRead;
use cr_types::types::{GemWell, LibraryType};
use cr_types::{MetricsFile, PerLibrarySortedBarcodeCounts};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::bin_file::BincodeFile;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::lz4_file::Lz4;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::{Histogram, JsonReporter, SimpleHistogram, TxHashMap, TxHasher};
use serde::{Deserialize, Serialize};
use shardio::{Range, ShardReader};
use stats::ReservoirSampler;
use std::fs::File;
use vdj_types::get_max_read_pairs_per_barcode;

/// This will determine the number of chunks in the stages
/// after barcode correction in VDJ. At 5k read pairs per cell,
/// this means a minimum of 100 cells per chunk.
pub const MIN_READ_PAIRS_PER_CHUNK: usize = 500_000;

/// Maximum number of chunks to create within a stage.
pub const MAX_CHUNKS: usize = 250;

/// JSON file containing mapping from barcode to count.
pub type VdjBarcodeCountsFile = JsonFile<TxHashMap<String, usize>>;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct RustBridgeStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub gem_well: GemWell,
    pub valid_uncorrected: Vec<ReadShardFile>,
    pub valid_corrected: Vec<ReadShardFile>,
    pub corrected_barcode_counts: PerLibrarySortedBarcodeCounts,
    pub max_reads_per_barcode: Option<usize>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct RustBridgeChunkInputs {
    #[mro_type = "map"]
    range: Range<Barcode>,
    valid_shards: Vec<ReadShardFile>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct RustBridgeChunkOutputs {
    chunk_bc_sorted_rna_reads: Lz4<BincodeFile<Vec<RnaRead>>>,
    barcodes_shard: JsonFile<()>,
    n50s_shard: BincodeFile<Vec<u32>>,
    processed_read_pairs: u64,
}

#[derive(Serialize, MartianStruct)]
pub struct RustBridgeStageOutputs {
    pub bc_sorted_rna_reads: Vec<Lz4<BincodeFile<Vec<RnaRead>>>>,
    pub gem_groups: Vec<u16>,
    pub barcodes: Vec<JsonFile<()>>,
    pub corrected_barcode_counts_json: VdjBarcodeCountsFile,
    pub summary: MetricsFile,
    pub n50_n50_rpu: u32,
    pub processed_read_pairs: u64,
}

/// Martian stage RUST_BRIDGE
pub struct RustBridge;

#[make_mro(mem_gb = 4)]
impl MartianStage for RustBridge {
    type StageInputs = RustBridgeStageInputs;
    type StageOutputs = RustBridgeStageOutputs;
    type ChunkInputs = RustBridgeChunkInputs;
    type ChunkOutputs = RustBridgeChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let shards: Vec<_> = args
            .valid_uncorrected
            .into_iter()
            .chain(args.valid_corrected)
            .collect();

        // Figure out how many chunks to create.
        let reader: ShardReader<RnaRead, BarcodeOrder> = ShardReader::open_set(&shards)?;
        let n = reader.len();
        let num_chunks = (n / MIN_READ_PAIRS_PER_CHUNK).clamp(1, MAX_CHUNKS);
        Ok(reader
            .make_chunks(num_chunks, &Range::all())
            .into_iter()
            .map(|range| RustBridgeChunkInputs {
                range,
                valid_shards: shards.clone(),
            })
            .collect())
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let reader: ShardReader<RnaRead, BarcodeOrder> =
            ShardReader::open_set(&chunk_args.valid_shards)?;

        let mut processed_read_pairs: u64 = 0;

        let mut barcodes = Vec::new();

        let chunk_bc_sorted_rna_reads: Lz4<_> = rover.make_path("chunk_bc_sorted_rna_reads");

        let barcodes_shard: JsonFile<_> = rover.make_path("barcode_shard");
        let mut read_lazy_writer = chunk_bc_sorted_rna_reads.lazy_writer()?;
        let mut n50s = Vec::new();

        // CELLRANGER-7889 "VDJ" is hardcoded in the MRO chemistry map.
        let chemistry_def = &args.chemistry_defs[&LibraryType::VdjAuto];
        let max_read_pairs_per_barcode = get_max_read_pairs_per_barcode(
            chemistry_def.is_paired_end(),
            args.max_reads_per_barcode,
        );

        for (barcode, barcode_read_iter) in &reader
            .iter_range(&chunk_args.range)?
            .chunk_by(|read| read.as_ref().ok().map(RnaRead::barcode))
        {
            let subsampled_reads = ReservoirSampler::sample_from_iter(
                barcode_read_iter,
                max_read_pairs_per_barcode,
                TxHasher::hash(barcode), // seed for subsampling
            );
            let mut umi_read_counts = SimpleHistogram::default();
            processed_read_pairs += subsampled_reads.len() as u64;
            for rna_read in subsampled_reads {
                let rna_read = rna_read?;
                read_lazy_writer.write_item(&rna_read)?;
                umi_read_counts.observe_owned(rna_read.umi());
            }
            if let Some(n50) = stats::n50(umi_read_counts.raw_counts().filter(|&&count| count > 1))
            {
                n50s.push(n50 as u32);
            }
            barcodes.push(barcode.unwrap());
        }
        read_lazy_writer.finish()?;
        serde_json::to_writer_pretty(File::create(&barcodes_shard)?, &barcodes)?;

        let n50s_shard: BincodeFile<_> = rover.make_path("n50s_shard");
        n50s_shard.write(&n50s)?;

        Ok(RustBridgeChunkOutputs {
            chunk_bc_sorted_rna_reads,
            barcodes_shard,
            n50s_shard,
            processed_read_pairs,
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        // FIXME: downstream stages should not depend on this JSON count file,
        // and should instead be updated to read the sorted counts format
        // directly.
        let vdj_lib_type = args
            .corrected_barcode_counts
            .library_types()
            .filter(LibraryType::is_vdj)
            .exactly_one()
            .ok()
            .unwrap();
        // Older code expects a JSON object mapping string barcode to numeric counts.
        // We have to collect this in order to write it, unfortunately.
        let counts_map: TxHashMap<_, _> = args
            .corrected_barcode_counts
            .iter_for_library(vdj_lib_type)?
            .process_results(|counts_iter| {
                counts_iter
                    .map(|(bc, count)| (bc.to_string(), count))
                    .collect()
            })?;

        let corrected_barcode_counts_json = rover
            .make_path::<VdjBarcodeCountsFile>("corrected_barcode_counts")
            .with_content(&counts_map)?;

        let mut n50s = Vec::new();
        for n50s_shard in chunk_outs.iter().map(|x| &x.n50s_shard) {
            let this_n50s: Vec<u32> = n50s_shard.read()?;
            n50s.extend(this_n50s);
        }
        let n50_n50_rpu = stats::n50(&n50s).unwrap_or(0);

        let mut metrics = JsonReporter::default();
        metrics.insert("n50_n50_rpu", n50_n50_rpu);

        Ok(RustBridgeStageOutputs {
            bc_sorted_rna_reads: chunk_outs
                .iter()
                .map(|x| x.chunk_bc_sorted_rna_reads.clone())
                .collect(),
            gem_groups: vec![args.gem_well.inner(); chunk_defs.len()],
            barcodes: chunk_outs
                .iter()
                .map(|x| x.barcodes_shard.clone())
                .collect(),
            corrected_barcode_counts_json,
            summary: MetricsFile::from_reporter(&rover, "summary", &metrics)?,
            n50_n50_rpu,
            processed_read_pairs: chunk_outs.iter().map(|x| x.processed_read_pairs).sum(),
        })
    }
}
