//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use kmer;
use kmer::{K, Kmer, Dir};
use std::collections::{VecDeque, HashSet};
use bitenc;
use utils::DMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::iter::FromIterator;

use debruijn::{TempGraph, Edge, EdgeDb};

#[derive(Debug)]
pub struct PhasedEdge {
    pub id: u32,
    pub sequence: bitenc::BitEnc,
    pub l_exts: Vec<(usize, kmer::Dir, bool)>,
    pub r_exts: Vec<(usize, kmer::Dir, bool)>,
    pub covs: [usize; 3],
}

impl PhasedEdge {
    pub fn get_exts(&self, dir: kmer::Dir) -> &Vec<(usize, kmer::Dir, bool)> {
        match dir {
            Dir::Left => &self.l_exts,
            Dir::Right => &self.r_exts,
        }
    }
}

/// Special case for phased local assembly
pub fn edge_hap_summary<'a>(e: &Edge<'a>, kmer_barcodes: &DMap<Kmer, Vec<u32>>) -> [usize; 3]
{
    let mut read_counts: [HashSet<u32>; 3] = [HashSet::new(), HashSet::new(), HashSet::new()];

    for k in e.sequence.kmers() {
        let min_kmer = k.min_rc().0;
        match kmer_barcodes.get(&min_kmer) {
            Some(bc_list) => {
                for bc in bc_list {
                    let hap = bc >> 16;
                    read_counts.get_mut(hap as usize).expect("asdf").insert(*bc);
                }
            },
            None => panic!("can't find BC"),
        }
    }

    [read_counts[0].len(), read_counts[1].len(), read_counts[2].len()]
}


/// Special case for phased local assembly
pub fn path_hap_summary(path: &Vec<(usize, Dir)>, edges: &Vec<PhasedEdge>, kmer_barcodes: &DMap<Kmer, Vec<u32>>) -> [usize; 3]
{
    let mut read_counts: [HashSet<u32>; 3] = [HashSet::new(), HashSet::new(), HashSet::new()];

    for &(edge_id, _) in path.iter() {
        let edge = edges.get(edge_id).expect("oob");
        for k in edge.sequence.kmers() {
            let min_kmer = k.min_rc().0;
            match kmer_barcodes.get(&min_kmer) {
                Some(bc_list) => {
                    for bc in bc_list {
                        let hap = bc >> 16;
                        read_counts.get_mut(hap as usize).expect("asdf").insert(*bc);
                    }
                },
                None => panic!("can't find BC"),
            }
        }
    }

    [read_counts[0].len(), read_counts[1].len(), read_counts[2].len()]
}


pub fn get_phased_edges(graph: &TempGraph, valid_kmers: &DMap<Kmer, Vec<u32>>) -> Vec<PhasedEdge> {
    // Generate edge DB
    let db = EdgeDb::new(graph);
    let mut edges = Vec::new();

    for edge in graph.iter() {
        let exts = edge.exts;
        let l_kmer = edge.sequence.first_kmer();
        let r_kmer = edge.sequence.last_kmer();
        let mut l_exts = Vec::new();
        let mut r_exts = Vec::new();

        let covs = edge_hap_summary(&edge, valid_kmers);

        for i in 0..4
        {
            if exts.has_ext(Dir::Left, i) {
                let link = db.find_link(l_kmer.extend_left(i), Dir::Left).expect("missing link");
                l_exts.push(link);
            }

            if exts.has_ext(Dir::Right, i) {
                let link = db.find_link(r_kmer.extend_right(i), Dir::Right).expect("missing link");
                r_exts.push(link);
            }
        }

        let e =
            PhasedEdge {
                id: edge.id as u32,
                sequence: edge.sequence.to_bitenc(),
                l_exts: l_exts,
                r_exts: r_exts,
                covs: covs,
            };
        edges.push(e);
    }

    edges
}


pub fn max_path<F, F2>(edges: &Vec<PhasedEdge>, score: F, solid_path: F2) -> Vec<(usize, kmer::Dir)>
    where F: Fn(&PhasedEdge) -> f32, F2: Fn(&PhasedEdge) -> bool {

    let mut best = 0;
    for e in 0..edges.len() {
        if score(&edges[e]) > score(&edges[best]) {
            best = e;
        }
    }

    let oscore = |state| {
        match state  {
            None => 0.0,
            Some((id, _)) => score(&edges[id]),
        }
    };

    let osolid_path = |state| {
        match state  {
            None => false,
            Some((id, _)) => solid_path(&edges[id]),
        }
    };


    // Now expand in each direction, greedily taking the best path. Stop if we hit a node we've
    // already put into the path
    let mut used_nodes = HashSet::new();
    let mut path = VecDeque::new();


    // Start w/ initial state
    used_nodes.insert(best);
    path.push_front((best, Dir::Left));

    for init in [(best, Dir::Left, false), (best, Dir::Right, true)].iter() {

        let &(start_node, dir, do_flip) = init;
        let mut current = (start_node, dir);

        loop {
            let mut next = None;
            let (cur_id, incoming_dir) = current;

            let mut solid_paths = 0;
            for &(id, dir, _) in edges[cur_id].get_exts(incoming_dir.flip()) {
                let cand = Some((id, dir));
                if osolid_path(cand) {
                    solid_paths += 1;
                }

                if oscore(cand) > oscore(next) {
                    next = cand;
                }
            }

            if solid_paths > 1 {
                break;
            }

            match next {
                Some((next_id, next_incoming)) if !used_nodes.contains(&next_id) => {

                    if do_flip {
                        path.push_front((next_id, next_incoming.flip()));
                    } else {
                        path.push_back((next_id, next_incoming));
                    }

                    used_nodes.insert(next_id);
                    current = (next_id, next_incoming);
                },
                _ => break,
            }
        }
    }

    for p in path.iter() {
        println!("{:?}", p);
    }

    Vec::from_iter(path)
}

pub fn sequence_of_path(edges: &Vec<PhasedEdge>, path: &Vec<(usize, kmer::Dir)>) -> bitenc::BitEnc {
    let mut seq = bitenc::BitEnc::new();

    for (idx, &(edge_id, dir)) in path.iter().enumerate() {
        let start = if idx == 0 { 0 } else { K - 1 };

        let edge_seq = match dir {
            Dir::Left => edges[edge_id].sequence.clone(),
            Dir::Right => edges[edge_id].sequence.rc(),
        };

        for p in start..edge_seq.len() {
            seq.push(edge_seq.get(p).expect("oob"))
        }
    }

    seq
}




pub fn write_edge(edge: &PhasedEdge, f: &mut Write) {

    let covs = edge.covs;
    let color =
        if covs[0] >= 3 && covs[1] <= 1 {
            ",color=red"
        } else if covs[1] >= 3 && covs[0] <= 1 {
            ",color=green"
        } else {
            ""
        };

    writeln!(f, "n{} [label=\"{}  -- {} -- {:?}\"{},style=filled]", edge.id, edge.id, edge.sequence.len(), edge.covs, color).unwrap();


    for &(id, incoming_dir, _) in edge.l_exts.iter() {
        let color = match incoming_dir { Dir::Left => "blue", Dir::Right => "red"};
        writeln!(f, "n{} -> n{} [color={}]", id, edge.id, color).unwrap();
    }

    for &(id, incoming_dir, _) in edge.r_exts.iter() {
        let color = match incoming_dir { Dir::Left => "blue", Dir::Right => "red"};
        writeln!(f, "n{} -> n{} [color={}]", edge.id, id, color).unwrap();
    }
}


/// Write the POA to a dot file.
pub fn to_dot(edges: &Vec<PhasedEdge>, path: PathBuf) {

    let mut f = File::create(path).expect("couldn't open file");

    writeln!(&mut f, "digraph {{").unwrap();
    for e in edges.iter() {
        write_edge(e, &mut f);
    }
    writeln!(&mut f, "}}").unwrap();
}
