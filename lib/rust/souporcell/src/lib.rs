//! Crate for Souporcell.
//!
//! Genetic demultiplexing using cell barcode by ref/alt count matrices.
#![deny(missing_docs)]

mod souporcell;

use anyhow::Result;
use cr_types::{GeneticDemuxParams, Mtx};
use itertools::Itertools;
use metric::TxHashMap;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use souporcell::{
    expectation_maximization, init_cluster_centers, load_cell_data, CellData, Params,
};

type ClusterIdx = usize;
type ProbMatrix = Vec<Vec<f32>>;
/// Struct for Souporcell chunk result
#[derive(Clone)]
pub struct SouporcellChunkResult {
    /// Log loss for the chunk
    pub log_loss: f32,
    /// Log probability matrix (log_probs[i][j]) for cell barcode (i) by cluster (j)
    pub log_probs: ProbMatrix,
}

/// Souporcell parameters
pub struct SouporcellParams {
    params: Params,
    loci_used: usize,
    cell_data: Vec<CellData>,
    locus_to_index: TxHashMap<usize, usize>,
}

/// Hold data for a Souporcell chunk
pub struct SouporcellChunk {
    seed: u64,
    epoch_count: usize,
    chunk_idx: ClusterIdx,
}

impl SouporcellParams {
    /// Create a new SouporcellParams instance.
    pub fn new(genetic_demux_params: &GeneticDemuxParams, ref_mtx: &Mtx, alt_mtx: &Mtx) -> Self {
        let params = {
            let mut params = Params {
                ref_mtx: ref_mtx.to_str().unwrap().to_string(),
                alt_mtx: alt_mtx.to_str().unwrap().to_string(),
                num_clusters: genetic_demux_params.num_clusters,
                ..Default::default()
            };
            if let Some(min_alt) = genetic_demux_params.min_alt_cells {
                params.min_alt = min_alt as u32;
            }
            if let Some(min_ref) = genetic_demux_params.min_ref_cells {
                params.min_ref = min_ref as u32;
            }
            if let Some(min_alt_umis) = genetic_demux_params.min_alt_umi {
                params.min_alt_umis = min_alt_umis as u32;
            }
            if let Some(min_ref_umis) = genetic_demux_params.min_ref_umi {
                params.min_ref_umis = min_ref_umis as u32;
            }
            if let Some(known_genotypes) = &genetic_demux_params.known_genotypes {
                params.known_genotypes = Some(known_genotypes.vcf.to_str().unwrap().to_string());
                params.known_genotypes_sample_names = known_genotypes.sample_names.clone();
            }
            params
        };
        let (loci_used, _total_cells, cell_data, _index_to_locuss, locus_to_index) =
            load_cell_data(&params);
        SouporcellParams {
            params,
            loci_used,
            cell_data,
            locus_to_index,
        }
    }
}

/// Create a vector of `SouporcellChunk`
pub fn chunk_souporcell(chunk_count: usize, restarts: usize, seed: u64) -> Vec<SouporcellChunk> {
    let epoch_count = restarts.div_ceil(chunk_count);
    let mut rng: SmallRng = SeedableRng::seed_from_u64(seed);
    (0..chunk_count)
        .map(|chunk_idx| SouporcellChunk {
            seed: rng.random::<u64>(),
            epoch_count,
            chunk_idx,
        })
        .collect()
}

/// Initialize cluster centers and run the expectation-maximization algorithm
fn run_epoch(
    rng: &mut SmallRng,
    chunk_idx: usize,
    epoch_idx: usize,
    params: &SouporcellParams,
) -> Result<SouporcellChunkResult> {
    let cluster_centers = init_cluster_centers(
        params.loci_used,
        &params.cell_data,
        &params.params,
        rng,
        &params.locus_to_index,
    )?;
    Ok(expectation_maximization(
        params.loci_used,
        cluster_centers,
        &params.cell_data,
        &params.params,
        epoch_idx,
        chunk_idx,
    ))
}

/// Run a chunk of Souporcell epochs.
/// Return the best log probabilities and the best total log probability among all epochs.
pub fn run_chunk(
    chunk: &SouporcellChunk,
    params: &SouporcellParams,
) -> Result<SouporcellChunkResult> {
    let mut rng = SmallRng::seed_from_u64(chunk.seed);
    (0..chunk.epoch_count)
        .map(|epoch_idx| run_epoch(&mut rng, chunk.chunk_idx, epoch_idx, params))
        .process_results(|iter| iter.max_by(cmp_results).unwrap())
}

/// Get best chunk result
pub fn get_best_chunk_result(chunk_results: &[SouporcellChunkResult]) -> SouporcellChunkResult {
    chunk_results
        .iter()
        .max_by(|a, b| cmp_results(a, b))
        .unwrap()
        .clone()
}

fn cmp_results(a: &SouporcellChunkResult, b: &SouporcellChunkResult) -> std::cmp::Ordering {
    a.log_loss
        .partial_cmp(&b.log_loss)
        .unwrap_or(std::cmp::Ordering::Equal)
}

/// Get the best cluster idx for each cell barcode
pub fn collapse_probability_matrix(best_log_probs: &[Vec<f32>]) -> Vec<ClusterIdx> {
    best_log_probs
        .iter()
        .map(|probs| {
            probs
                .iter()
                .position_max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .collect()
}
