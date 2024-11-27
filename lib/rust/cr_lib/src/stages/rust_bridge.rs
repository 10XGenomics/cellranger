//! Martian stage RUST_BRIDGE

use crate::barcode_sort::BarcodeOrder;
use crate::types::ReadShardFile;
use anyhow::Result;
use barcode::Barcode;
use cr_types::chemistry::ChemistryDefs;
use cr_types::rna_read::RnaRead;
use cr_types::types::{GemWell, LibraryType};
use cr_types::{BcCountFormat, MetricsFile};
use itertools::{GroupBy, Itertools};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::bin_file::BincodeFile;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::lz4_file::Lz4;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::{JsonReport, JsonReporter, SimpleHistogram};
use parameters_toml::vdj_max_reads_per_barcode;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use shardio::{Range, ShardReader, SortKey};
use stats::ReservoirSampler;
use std::fs::File;

/// This will determine the number of chunks in the stages
/// after barcode correction in VDJ. At 5k read pairs per cell,
/// this means a minimum of 100 cells per chunk.
pub const MIN_READ_PAIRS_PER_CHUNK: usize = 500_000;

/// Maximum number of chunks to create within a stage.
pub const MAX_CHUNKS: usize = 250;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct RustBridgeStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub gem_well: GemWell,
    pub valid_uncorrected: Vec<ReadShardFile>,
    pub valid_corrected: Vec<ReadShardFile>,
    pub raw_barcode_counts: BcCountFormat,
    pub corrected_barcode_counts: BcCountFormat,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct RustBridgeChunkInputs {
    #[mro_type = "map"]
    range: shardio::Range<Barcode>,
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
    pub raw_barcode_counts_json: JsonFile<()>,
    pub corrected_barcode_counts_json: JsonFile<SimpleHistogram<String>>,
    pub summary: MetricsFile,
    pub n50_n50_rpu: u32,
    pub processed_read_pairs: u64,
}

fn iter_range_grouped<'a, T, S>(
    reader: &'a ShardReader<T, S>,
    range: &Range<<S as SortKey<T>>::Key>,
) -> Result<
    GroupBy<
        Option<<S as SortKey<T>>::Key>,
        impl Iterator<Item = Result<T>> + 'a,
        impl FnMut(&Result<T>) -> Option<<S as SortKey<T>>::Key> + 'a,
    >,
    anyhow::Error,
>
where
    T: 'a + DeserializeOwned,
    <S as SortKey<T>>::Key: Clone + Ord + DeserializeOwned,
    S: 'a + SortKey<T>,
{
    Ok(reader.iter_range(range)?.group_by(|read| {
        read.as_ref()
            .ok()
            .map(|r| <S as SortKey<T>>::sort_key(r).into_owned())
    }))
}

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
        let max_read_pairs_per_barcode = if chemistry_def.is_paired_end() {
            *vdj_max_reads_per_barcode()? / 2
        } else {
            *vdj_max_reads_per_barcode()?
        };

        for (barcode, barcode_read_iter) in &iter_range_grouped(&reader, &chunk_args.range)? {
            let subsampled_reads = ReservoirSampler::sample_from_iter(
                barcode_read_iter,
                max_read_pairs_per_barcode,
                fxhash::hash64(&barcode), // seed for subsampling
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
        let raw_bc_counts: JsonFile<()> = rover.make_path("raw_barcode_counts");
        let corrected_bc_counts: JsonFile<_> = rover.make_path("corrected_barcode_counts");

        args.raw_barcode_counts
            .read()?
            .into_iter()
            .filter_map(|(k, v)| match k {
                LibraryType::Vdj(_) => Some(v),
                _ => None,
            })
            .exactly_one()
            .unwrap()
            .report(&raw_bc_counts)?;

        args.corrected_barcode_counts
            .read()?
            .into_iter()
            .filter_map(|(k, v)| match k {
                LibraryType::Vdj(_) => Some(v),
                _ => None,
            })
            .exactly_one()
            .unwrap()
            .report(&corrected_bc_counts.as_ref())?;

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
            raw_barcode_counts_json: raw_bc_counts,
            corrected_barcode_counts_json: corrected_bc_counts,
            summary: MetricsFile::from_reporter(&rover, "summary", &metrics)?,
            n50_n50_rpu,
            processed_read_pairs: chunk_outs.iter().map(|x| x.processed_read_pairs).sum(),
        })
    }
}

/*
// Functions to create non-SLFE fastq files with augmented header. This format was used in
// CR 3.1 and earlier.
fn create_owned_records(read: &RnaRead) -> (OwnedRecord, Option<OwnedRecord>) {
    let head = augmented_fastq_header(read);
    let seq = read.r1_seq().to_vec();
    let qual = read.r1_qual().to_vec();

    let r1_record = OwnedRecord {
        head: head.clone(),
        seq,
        sep: None,
        qual,
    };

    let r2_record = if read.r2_exists() {
        let seq = read.r2_seq().unwrap().to_vec();
        let qual = read.r2_qual().unwrap().to_vec();
        Some(OwnedRecord {
            head,
            seq,
            sep: None,
            qual,
        })
    } else {
        None
    };

    (r1_record, r2_record)
}


fn augmented_fastq_header(read: &RnaRead) -> Vec<u8> {
    const SAMPLE_INDEX_TAG: &'static [u8] = b"BC";
    const SAMPLE_INDEX_QUAL_TAG: &'static [u8] = b"QT";
    const RAW_BARCODE_TAG: &'static [u8] = b"CR";
    const RAW_BARCODE_QUAL_TAG: &'static [u8] = b"CY";
    const RAW_UMI_TAG: &'static [u8] = b"UR";
    const RAW_UMI_QUAL_TAG: &'static [u8] = b"UY";
    const PROCESSED_BARCODE_TAG: &'static [u8] = b"CB";
    const SEP: &'static [u8] = b"|||";

    let header = read.header().to_vec();
    let mut iter = header.split(|&m| m == b' ');
    let mut header_augment: Vec<u8> = Vec::new();
    header_augment.extend(iter.next().unwrap());

    header_augment.extend(SEP);
    header_augment.extend(SAMPLE_INDEX_TAG);
    header_augment.extend(SEP);
    header_augment.extend(&[]);
    header_augment.extend(SEP);
    header_augment.extend(SAMPLE_INDEX_QUAL_TAG);
    header_augment.extend(SEP);
    header_augment.extend(&[]);
    header_augment.extend(SEP);
    header_augment.extend(RAW_BARCODE_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_bc_seq().as_bytes());
    header_augment.extend(SEP);
    header_augment.extend(RAW_BARCODE_QUAL_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_bc_qual().as_bytes());
    header_augment.extend(SEP);
    header_augment.extend(RAW_UMI_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_umi_seq());
    header_augment.extend(SEP);
    header_augment.extend(RAW_UMI_QUAL_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_umi_qual());
    header_augment.extend(SEP);
    header_augment.extend(PROCESSED_BARCODE_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.barcode().to_string().as_bytes());
    header_augment.push(b' ');
    header_augment.extend(iter.collect::<Vec<_>>().join(&b' '));

    header_augment
}

*/
