//! high-level louvain implementation

use leiden::clustering::{Clustering, SimpleClustering};
use leiden::louvain::Louvain;
use leiden::louvain_parallel::ParallelLouvain;
use leiden::objective::{cpm, par_cpm};
use leiden::Network;
use log::info;
use ndarray::Array2;
use scan_rs::merge_clusters::relabel_by_size;
use std::time::Instant;

const TOLERANCE: f64 = 1e-6;

fn build_network(neighbors: &Array2<u32>) -> (Network, SimpleClustering) {
    let n_nodes = neighbors.shape()[0];
    let n_edges = n_nodes * neighbors.shape()[1];
    let adj = |i: usize, j: usize| -> Option<(u32, u32)> {
        if neighbors[[i, j]] == u32::MAX {
            None
        } else {
            Some((i as u32, neighbors[[i, j]]))
        }
    };
    let adjacency = (0..neighbors.shape()[0])
        .flat_map(|i| (0..neighbors.shape()[1]).filter_map(move |j| adj(i, j)));
    info!(
        "building network for neighbor matrix of shape: ({},{})",
        neighbors.shape()[0],
        neighbors.shape()[1]
    );
    let network = Louvain::build_network(n_nodes, n_edges, adjacency);
    let clustering: SimpleClustering = Clustering::init_different_clusters(n_nodes);
    (network, clustering)
}

fn do_relabeling(clustering: &SimpleClustering, n_nodes: usize) -> Vec<i16> {
    //! Relabel clusters by size
    let mut labels: Vec<i16> = Vec::with_capacity(n_nodes);
    let mut max_label = 0;
    for i in 0..n_nodes {
        let label = clustering.get(i) as i16;
        if max_label < label {
            max_label = label;
        }
        labels.push(label);
    }
    relabel_by_size(labels)
}

pub(crate) fn run_louvain(
    neighbors: &Array2<u32>,
    resolution: f64,
    seed: Option<usize>,
) -> Vec<i16> {
    //! Run louvain clustering
    let n_nodes = neighbors.shape()[0];
    let (network, mut clustering) = build_network(neighbors);

    let mut louvain = Louvain::new(resolution, seed);

    let mut score = cpm(resolution, &network, &clustering);
    info!("louvain starting cpm score: {:.6}", score);
    let now = Instant::now();
    for i in 0.. {
        let updated = louvain.iterate_one_level(&network, &mut clustering);
        info!("  iteration {}", i);
        let new_score = cpm(resolution, &network, &clustering);
        info!("    cpm score: {:.6}", new_score);
        if !updated || (new_score - score) <= TOLERANCE {
            info!("louvain optimized in {} iterations", i + 1);
            info!("louvain final cpm score: {:.6}", new_score);
            break;
        }
        score = new_score;
    }

    let elapsed = now.elapsed();
    println!("optimization time: {elapsed:.2?}");

    do_relabeling(&clustering, n_nodes)
}

pub(crate) fn run_louvain_parallel(neighbors: &Array2<u32>, resolution: f64) -> Vec<i16> {
    //! Run louvain clustering, parallel variant
    let n_nodes = neighbors.shape()[0];
    let (network, mut clustering) = build_network(neighbors);

    let mut louvain = ParallelLouvain::new(resolution);

    let mut score = par_cpm(resolution, &network, &clustering);
    info!("louvain starting cpm score {:.6}", score);
    let now = Instant::now();
    for i in 0.. {
        let updated = louvain.iterate_one_level(&network, &mut clustering);
        info!("  iteration {}", i);
        let new_score = par_cpm(resolution, &network, &clustering);
        info!("    cpm score: {:.6}", new_score);
        if !updated || ((new_score - score) < TOLERANCE) {
            info!("louvain optimized in {} iterations", i + 1);
            info!("louvain final cpm score: {:.6}", new_score);
            break;
        }
        score = new_score;
    }

    let elapsed = now.elapsed();
    println!("optimization time: {elapsed:.2?}");

    do_relabeling(&clustering, n_nodes)
}
