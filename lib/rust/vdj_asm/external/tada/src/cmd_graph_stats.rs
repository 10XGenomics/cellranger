//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use utils;
use std::path::Path;
use utils::MultiVec;
use debruijn::TempGraph;
use csv;


#[derive(RustcEncodable)]
struct GraphNode {
    len: usize,
    num_bcs: usize,
    exts_left: u8,
    exts_right: u8,
}


pub fn main_graph_stats(graph: &Path, bcs: &Path, stats: &Path)
{
    let graph : TempGraph = utils::read_obj(graph).expect("read graph");
    let bcs : MultiVec<u32> = utils::read_obj(bcs).expect("read bcs");

    let mut wtr = csv::Writer::from_file(stats).expect("open csv").delimiter(b'\t');

    for e in graph.iter()
    {
        let record = GraphNode {
            len: e.sequence.length,
            num_bcs: bcs.vec_len[e.id] as usize,
            exts_left: e.exts.num_exts_l(),
            exts_right: e.exts.num_exts_r(),
        };

        wtr.encode(record).expect("csv write error");
    }
}

/*
pub fn main_show_permutation(perm: &Path)
{
    let permutation : Vec<usize> = utils::read_obj(perm).expect("couldn't read permutation");

    let graph : TempGraph = utils::read_obj(graph).expect("read graph");
    let bcs : MultiVec<u32> = utils::read_obj(bcs).expect("read bcs");

    let mut wtr = csv::Writer::from_file(stats).expect("open csv").delimiter(b'\t');

    for e in graph.iter()
    {
        let record = GraphNode {
            len: e.sequence.length,
            num_bcs: bcs.vec_len[e.id] as usize,
            exts_left: e.exts.num_exts_l(),
            exts_right: e.exts.num_exts_r(),
        };

        wtr.encode(record).expect("csv write error");
    }
}
*/
