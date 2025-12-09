//! Martian stage CALL_TAGS_GENETIC
//! Assign cells to samples using genetic demultiplexing (Souporcell).
use anyhow::Result;
use cr_types::filtered_barcodes::FilteredBarcodesCsv;
use cr_types::{GeneticDemuxParams, Mtx};
use itertools::Itertools;
use martian::prelude::{MartianRover, MartianStage};
use martian::{MartianVoid, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeWrite, LazyFileTypeIO};
use metric::TxHashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use souporcell::{
    SouporcellChunkResult, SouporcellParams, chunk_souporcell, collapse_probability_matrix,
    get_best_chunk_result, run_chunk,
};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct CallTagsGeneticStageInputs {
    /// Filtered barcodes CSV file.
    pub filtered_barcodes: FilteredBarcodesCsv,
    /// Path to the reference matrix.
    pub ref_mtx: Option<Mtx>,
    /// Path to the alternate matrix.
    pub alt_mtx: Option<Mtx>,
    /// Genetic demultiplexing parameters.
    pub genetic_demux_params: Option<GeneticDemuxParams>,
}

/// The Martian stage outputs.
#[derive(Serialize, MartianStruct)]
pub struct CallTagsGeneticStageOutputs {
    /// Path to the cell tags TSV file.
    // TODO: Use `BarcodeId` instead of `String` for the barcodes.
    pub barcodes_per_tag: Option<JsonFile<TxHashMap<usize, Vec<String>>>>,
}

/// Martian stage CALL_TAGS_GENETIC
/// Assign cells to samples using genetic demultiplexing (Souporcell).
pub struct CallTagsGenetic;

#[make_mro(volatile = strict)]
impl MartianStage for CallTagsGenetic {
    type StageInputs = CallTagsGeneticStageInputs;
    type StageOutputs = CallTagsGeneticStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        if args.genetic_demux_params.is_none() {
            return Ok(StageDef::new());
        }
        Ok(StageDef::with_join_resource(
            Resource::with_mem_gb(32).threads(16),
        ))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    /// Run the Martian stage CALL_TAGS_GENETIC.
    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let (Some(genetic_demux_params), Some(ref_mtx), Some(alt_mtx)) =
            (&args.genetic_demux_params, &args.ref_mtx, &args.alt_mtx)
        else {
            return Ok(Self::StageOutputs {
                barcodes_per_tag: None,
            });
        };
        let params = SouporcellParams::new(genetic_demux_params, ref_mtx, alt_mtx);
        // TODO: Use martian chunking
        let chunks = chunk_souporcell(
            rover.get_threads(),
            genetic_demux_params.restarts.unwrap_or(100),
            0,
        );
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(rover.get_threads())
            .build()
            .unwrap();
        // TODO: Can we iterate over the results instead of collecting them all?
        let chunks_results: Result<Vec<_>, _> = pool.install(|| {
            chunks
                .par_iter()
                .map(|chunk| run_chunk(chunk, &params))
                .collect()
        });
        let chunks_results = chunks_results?;
        let SouporcellChunkResult {
            log_loss,
            log_probs,
        } = get_best_chunk_result(&chunks_results);

        println!("Best probability: {log_loss}");
        let barcode_to_cluster_assignment = collapse_probability_matrix(&log_probs);
        let sorted_filtered_cbs = args
            .filtered_barcodes
            .lazy_reader()?
            .process_results(|iter| iter.map(|x| x.barcode).sorted())?
            .dedup();

        let barcodes_per_tag = sorted_filtered_cbs.zip(barcode_to_cluster_assignment).fold(
            TxHashMap::default(),
            |mut acc: TxHashMap<usize, Vec<String>>, (b, c_idx)| {
                acc.entry(c_idx).or_default().push(b.to_string());
                acc
            },
        );
        let barcodes_per_tag_file: JsonFile<_> = rover.make_path("barcodes_per_tag");
        barcodes_per_tag_file.write(&barcodes_per_tag)?;
        Ok(Self::StageOutputs {
            barcodes_per_tag: Some(barcodes_per_tag_file),
        })
    }
}
