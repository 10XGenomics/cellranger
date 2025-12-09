//! Martian stage COUNT_ALLELES
//! Count alleles of each cell using VarTrix
#![deny(missing_docs)]

use crate::BamFile;
use anyhow::{Context, Ok, Result};
use bio::io::fasta;
use cr_types::filtered_barcodes::FilteredBarcodesCsv;
use cr_types::{BarcodeIndex, GeneticDemuxParams, Mtx, SamHeaderFile};
use itertools::Itertools;
use martian::prelude::{MartianRover, MartianStage};
use martian::{Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::bin_file::BincodeFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use rust_htslib::bam::IndexedReader as BamReader;
use rust_htslib::bcf::Reader as BcfReader;
use serde::{Deserialize, Serialize};
use sprs::TriMat;
use sprs::io::write_matrix_market;
use std::path::PathBuf;
use vartrix::{AlleleType, Snp, SnpResults, VartrixSettings, VcfChunk, evaluate_snp, snp_records};

type SnpBinFile = BincodeFile<Vec<Snp>>;
type SnpResultBinFile = BincodeFile<Vec<SnpResults>>;

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct CountAllelesStageInputs {
    /// Path to the reference genome.
    pub reference_path: PathBuf,
    /// Position sorted BAM file.
    pub bam_file: BamFile,
    /// BAM file header used to map contig names to contig IDs.
    pub bam_header: SamHeaderFile,
    /// Path to the filtered barcodes CSV file.
    pub filtered_barcodes: FilteredBarcodesCsv,
    /// Genetic demultiplexing parameters.
    pub genetic_demux_params: Option<GeneticDemuxParams>,
}

/// The Martian stage outputs.
#[derive(Serialize, MartianStruct)]
pub struct CountAllelesStageOutputs {
    /// Path to the ref.mtx file.
    pub ref_mtx: Option<Mtx>,
    /// Path to the alt.mtx file.
    pub alt_mtx: Option<Mtx>,
}

/// COUNT_ALLELES chunk inputs
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ChunkInputs {
    /// Serialized SNP records.
    pub chunk_file: BincodeFile<VcfChunk>,
    /// Max SNP index in all chunks.
    pub max_snp_idx: usize,
}

/// COUNT_ALLELES chunk outputs
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ChunkOutputs {
    /// Scores for the SNP records.
    pub results: SnpResultBinFile,
    /// Larget SNP index
    pub max_snp_idx: usize,
}

/// Martian stage COUNT_ALLELES
/// Count alleles of each cell using VarTrix
pub struct CountAlleles;

#[make_mro(volatile = strict)]
impl MartianStage for CountAlleles {
    type StageInputs = CountAllelesStageInputs;
    type StageOutputs = CountAllelesStageOutputs;
    type ChunkInputs = ChunkInputs;
    type ChunkOutputs = ChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let Some(genetic_demux_params) = args.genetic_demux_params else {
            println!("Genetic demultiplexing parameters are not provided");
            return Ok(StageDef::new());
        };

        let mut vcf_reader = BcfReader::from_path(&genetic_demux_params.candidate_snps)?;
        let contig_info = args.bam_header.contig_info_by_name();
        let snps = snp_records(&mut vcf_reader, &contig_info)?;
        let (chunk_files, max_snp_idx) =
            snps.process_results(|iter| distribute_records(iter, 10000, 64, rover))??;

        let stage_def: StageDef<_> = chunk_files
            .into_iter()
            .map(|chunk_file| {
                (
                    ChunkInputs {
                        chunk_file,
                        max_snp_idx,
                    },
                    Resource::with_mem_gb(2).threads(1),
                )
            })
            .collect();

        Ok(stage_def.join_resource(Resource::with_mem_gb(4)))
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let settings = VartrixSettings::default();
        let mut fasta_reader =
            fasta::IndexedReader::from_file(&args.reference_path.join("fasta/genome.fa"))?;
        let cb_index = BarcodeIndex::from_sorted(
            args.filtered_barcodes
                .lazy_reader()?
                .process_results(|iter| iter.map(|x| x.barcode).sorted())?
                .dedup(),
        );
        let mut bam_reader = BamReader::from_path(args.bam_file)?;
        let contig_info_by_tid = args.bam_header.contig_info_by_tid();

        let snp_reader = chunk_args.chunk_file.lazy_reader()?;
        let outfile: SnpResultBinFile = rover.make_path("chunk_scores.bin");
        let mut writer = outfile.lazy_writer()?;

        for snp in snp_reader {
            let snp = snp?;
            let snp_results = evaluate_snp(
                &snp,
                &mut bam_reader,
                &mut fasta_reader,
                &contig_info_by_tid,
                &cb_index,
                &settings,
            )?;
            writer.write_item(&snp_results)?;
        }

        writer.finish()?;

        Ok(Self::ChunkOutputs {
            results: outfile,
            max_snp_idx: chunk_args.max_snp_idx,
        })
    }

    /// Join the chunk outputs and process their results.
    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        if chunk_outs.is_empty() {
            return Ok(Self::StageOutputs {
                ref_mtx: None,
                alt_mtx: None,
            });
        }

        let ref_mtx: Mtx = rover.make_path("ref");
        let alt_mtx: Mtx = rover.make_path("alt");
        let num_cols = chunk_outs[0].max_snp_idx + 1;
        let num_rows = args
            .filtered_barcodes
            .lazy_reader()?
            .process_results(|iter| iter.map(|x| x.barcode).sorted())?
            .dedup()
            .count();
        let mut ref_matrix = TriMat::new((num_cols, num_rows));
        let mut alt_matrix = TriMat::new((num_cols, num_rows));
        for chunk_outs in chunk_outs {
            for results in chunk_outs.results.lazy_reader()? {
                let results = results?;
                for calls in results.cell_barcode_calls {
                    ref_matrix.add_triplet(
                        results.rec_idx,
                        calls.cb_idx,
                        calls.allele_counts[&AlleleType::Ref],
                    );
                    alt_matrix.add_triplet(
                        results.rec_idx,
                        calls.cb_idx,
                        calls.allele_counts[&AlleleType::Alt],
                    );
                }
            }
        }
        write_matrix_market(&ref_mtx, &ref_matrix)
            .with_context(|| format!("failed to write ref matrix to {ref_mtx:?}"))?;
        write_matrix_market(&alt_mtx, &alt_matrix)
            .with_context(|| format!("failed to write alt matrix to {alt_mtx:?}"))?;

        Ok(Self::StageOutputs {
            ref_mtx: Some(ref_mtx),
            alt_mtx: Some(alt_mtx),
        })
    }
}

/// Distribute SNPs into chunks for Martian processing.
/// The functions write `min_chunk_size` SNPs to a chunk at a time.
/// It then moves to the next chunk (and creating if needed).
/// It goes back to the first chunk if it has created `max_num_chunks` chunks.
/// Rover is used to create the chunk files.
/// Returns a vector of chunk files and the maximum Snp.idx value
fn distribute_records(
    records: impl Iterator<Item = Snp>,
    min_chunk_size: usize,
    max_num_chunks: usize,
    rover: MartianRover,
) -> Result<(Vec<SnpBinFile>, usize)> {
    let mut chunks = Vec::with_capacity(max_num_chunks);
    let mut writers = Vec::with_capacity(max_num_chunks);
    let mut max_snp_idx = 0;
    for (idx, mini_chunk) in records.chunks(min_chunk_size).into_iter().enumerate() {
        let file_idx = idx % max_num_chunks;
        if file_idx == chunks.len() {
            // We need to create a chunk
            let chunk_file: SnpBinFile = rover.make_path(format!("chunk_{file_idx}.bin"));
            let writer = chunk_file.lazy_writer()?;
            chunks.push(chunk_file);
            writers.push(writer);
        }
        let writer = &mut writers[file_idx];
        for record in mini_chunk {
            writer.write_item(&record)?;
            max_snp_idx = max_snp_idx.max(record.idx);
        }
    }
    // Close all writers
    for writer in writers {
        writer.finish()?;
    }
    Ok((chunks, max_snp_idx))
}
