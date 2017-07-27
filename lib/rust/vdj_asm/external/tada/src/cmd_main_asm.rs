//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Final assembly command
//!
//! Assemble short contigs into larger sequences, terminated by branching extensions or hard stops
//! Compute the set of barcodes that touched each sequence in the graph.

use std::path::{Path, PathBuf};
use utils;
use debruijn;
use std::iter;
use debruijn::TempGraph;
use utils::MultiVec;
use kmer::{Lmer, Exts};
use martian::{JsonDict, MartianStage};

pub fn main_main_asm(shard_asms: Vec<PathBuf>, shard_bcs: Vec<PathBuf>, out_graph: &Path, out_edge_bcs: &Path)
{
    // Generate a graph from the shard assemblies
    fn build_graph(shard_asms : Vec<PathBuf>) -> (TempGraph, usize)
    {
        let mut sedges : Vec<(Lmer, Exts)> = Vec::new();

        info!("Loading sedge data");
        for p in shard_asms {
            let mut shard_sedges = utils::read_obj(&p).expect("loading shard");
            sedges.append(&mut shard_sedges);
        }

        let n_sedges = sedges.len();
        info!("Loaded {} sedges", n_sedges);

        info!("Building full graph");
        let graph = debruijn::build_edges(sedges);
        info!("Built graph. Got {} edges", graph.len());

        (graph, n_sedges)
    }

    // Build the graph and discard the sedge data
    let (mut graph, n_sedges) = build_graph(shard_asms);

    // Make a map from sedge index to graph edge index
    let mut sedge_to_edge : Vec<isize> = iter::repeat(-1).take(n_sedges).collect();
    for (edge_idx, edge) in Iterator::enumerate(graph.iter()) {
        for sedge_idx in edge.sedges.unwrap()
        {
            sedge_to_edge[*sedge_idx as usize] = edge_idx as isize;
        }
    }

    // Dont save the sedge list of the graph
    graph.edge_sedges = None;
    info!("Writing graph");
    utils::write_obj(&graph, out_graph).expect("write failed");


    info!("Collating BCs for edges");
    // Keep the BC list for each edge in an array of HashSets
    let mut bcs_sets : Vec<Vec<u32>> = Vec::with_capacity(graph.len());
    for _ in 0..graph.len() {
        bcs_sets.push(Vec::new());
    }

    // Iterate over shards, add the BCS for each sedge to the set for each edge
    let mut shard_idx_base = 0;
    for (i, shard_file) in shard_bcs.iter().enumerate() {
        if i % 100 == 0 { info!("Coalating BC for shard: {}", i); }

        let shard_bcs : MultiVec<u32> = utils::read_obj(&shard_file).expect("read");

        for sedge_idx in 0..shard_bcs.len() {
            let global_sedge_idx = sedge_idx + shard_idx_base;
            let edge_idx = sedge_to_edge[global_sedge_idx];

            let mut bcs_set = bcs_sets.get_mut(edge_idx as usize).unwrap();

            if bcs_set.len() < 1000 {
                bcs_set.extend(shard_bcs.get_slice(sedge_idx));
            }
        }

        // Occasionally compact-ify the bcs_sets
        if i % 64 == 0
        {
            info!("Compacting BC sets");
            // Occasionally compact-ify the bcs_sets
            for bcs_set in bcs_sets.iter_mut()
            {
                bcs_set.sort();
                bcs_set.dedup();
            }
        }

        shard_idx_base += shard_bcs.len();
    }


    for bcs_set in bcs_sets.iter_mut()
    {
        bcs_set.sort();
        bcs_set.dedup();
    }

    // Pull out per-edge BC lists into a MultiVec
    let mut edge_bcs = MultiVec::new();
    let mut temp_vec = Vec::new();

    for bc_set in bcs_sets {
        temp_vec.clear();
        temp_vec.extend(bc_set);
        temp_vec.sort();
        edge_bcs.add_slice(temp_vec.as_slice());
    }

    utils::write_obj(&edge_bcs, out_edge_bcs).expect("write failed");
}


pub struct MainAsmMartian;

impl MartianStage for MainAsmMartian {
    fn split(&self, _: JsonDict) -> JsonDict {
        panic!("non-splitting stage");
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {
        let sedge_asm = args["sedge_asm"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();
        let sedge_bcs = args["sedge_bcs"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();


        let out_asm = Path::new({outs.get("asm_graph").unwrap().as_string().unwrap().clone()});
        let out_bcs = Path::new({outs.get("asm_bcs").unwrap().as_string().unwrap().clone()});

        main_main_asm(sedge_asm, sedge_bcs, &out_asm, &out_bcs);
        outs.clone()
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, _: Vec<JsonDict>) -> JsonDict {
        panic!("Non-splitting stage");
    }
}
