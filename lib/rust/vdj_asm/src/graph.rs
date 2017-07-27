//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

extern crate tada;
use tada::bitenc;
use tada::kmer;
use tada::debruijn;
use std::fs::File;
use std::io::Write;
use std::path::{PathBuf};
use std::collections::{HashSet, HashMap, VecDeque};
use constants::{P, K, UmiType, ReadType};
use sw;
use graph_read;
use std;
use tada::msp;
use itertools::{Itertools};
use std::iter::FromIterator;

/// BSPs from a vector of seqs
/// seqs is a vector of tuples (seq, qual, bc/UMI, read id)
pub fn build_bsps(reads: &Vec<&graph_read::Read>) -> Vec<tada::kmer::Bsp> {
    let mut bsps = Vec::new();
    let permutation = (0..1<<(P*2)).collect();

    for read in reads.iter() {
        // MSP read
        let seq = read.seq.iter().collect::<Vec<u8>>();
        let msp_parts = msp::simple_scan(K, P, seq.as_slice(), &permutation, false);
        for (_, _, start_pos, slice_length) in msp_parts {
            //assert!(slice_length >= k);
            // WARNING: THIS WON'T MAINTAIN THE READ ID because we convert to a
            // smaller type.
            let b = kmer::Bsp::new(seq.as_slice(), start_pos, slice_length, read.umi, read.id as u16);
            bsps.push(b);

            if (b.len() as usize) < K || (b.len() as usize) != b.sequence.len()
            {
                println!("bad bsp: {:?}", b);
                panic!("bad bsp: {:?}", b);
            }
        }
    }

    bsps
}


/// Read a shard and determine the valid kmers
/// Returns a tuple: First element is a map from kmers to all the partitions that
/// contain them. Only passing kmers included.
/// Second element of the tuple is a map from kmers to the number of occurrences.
/// All kmers are included.
#[inline(never)]
pub fn process_kmer_shard(bsps: &Vec<kmer::Bsp>, min_kmer_count: u32) ->
    (tada::utils::DMap<tada::kmer::Kmer, Vec<u32>>,
    tada::utils::DMap<kmer::Kmer, u32>) {

    let mut kmers = Vec::new();

    for b in bsps {
        for k in b.kmers() {
            kmers.push((k, b.partition));
        }
    }

    let mut all_kmer_count = tada::utils::new_dmap();
    let mut final_kmers: tada::utils::DMap<kmer::Kmer, Vec<u32>> = tada::utils::new_dmap();
    let mut observed_kmers = 0;

    kmers.sort();

    // kmers is a vector of tuples (kmer, partition). Group by kmer.
    // final_kmers will be a map from kmer to a list of partitions.
    for (kmer, group) in &kmers.into_iter().group_by_lazy(|elt| elt.0) {
        // set of barcodes
        let mut pris = HashSet::new();
        observed_kmers = observed_kmers + 1;
        for (_, bc) in group {
            pris.insert(bc);
        }

        all_kmer_count.insert(kmer, pris.len() as u32);
        if pris.len() > min_kmer_count as usize {
            final_kmers.insert(kmer, pris.into_iter().collect());
        }
    }

    println!("Kmers observed: {}. Kmers accepted: {}", observed_kmers, final_kmers.len());

    (final_kmers, all_kmer_count)
}

/// Build a graph from a vector of graph_read::Read objects.
/// min_kmers: min number of kmer occurrences to consider a kmer
/// max_kmers: maximum number of kmers to consider. If the total number of observed 
/// kmers (that pass the min_kmers cutoff) is greater than max_kmers, then we'll
/// get the top max_kmers kmers by frequency.
pub fn build_graph(reads: &Vec<&graph_read::Read>,
                   min_kmers:usize, max_kmers: usize) -> EdgeGraph {

    let bsps = build_bsps(reads);

    // This will chose between the kmer and its RC
    let (mut valid_kmers, all_kmers) = process_kmer_shard(&bsps, min_kmers as u32);

    if valid_kmers.len() > max_kmers {
        println!("Sampling {} out of {} valid_kmers", max_kmers, valid_kmers.len());
        valid_kmers = sample_kmers(&valid_kmers, &all_kmers, max_kmers);
    }

    //Returns list of tuples (kmer, extensions)
    //This will compute valid and invalid extension for all valid kmers.
    let kmer_extensions = tada::debruijn::kmer_extensions(&valid_kmers, &all_kmers, &bsps, false);

    //Returns list of tuples (Lmer, extensions). This is similar to kmer_extensions
    //but using a larger K. This is done by extending kmers with unique extensions.
    //This won't use invalid kmers to build the sedges but the extensions of the
    //sedges might correspond to invalid kmers.
    let sedges = tada::debruijn::build_sedges(&kmer_extensions, false);

    let graph = tada::debruijn::build_edges(sedges);
    println!("# Edges: {}", graph.start_pos.len());

    let barcoded_graph = EdgeGraph::new(&graph, &valid_kmers);

    barcoded_graph
}

fn sample_kmers(valid_kmers: &tada::utils::DMap<tada::kmer::Kmer, Vec<u32>>,
                all_kmers: &tada::utils::DMap<kmer::Kmer, u32>,
                max_kmers: usize) -> tada::utils::DMap<tada::kmer::Kmer, Vec<u32>>{
    // Sort by decreasing number of occurrences
    let mut sorted_kmers : Vec<tada::kmer::Kmer> = valid_kmers.keys().map(|x| x.clone()).collect();
    sorted_kmers.sort_by_key(|x| all_kmers.get(x).unwrap());
    sorted_kmers.reverse();
    // Take the first max_kmers
    let sel_kmers : HashSet<tada::kmer::Kmer> = sorted_kmers.iter().take(max_kmers).map(|x| *x).collect();
    HashMap::from_iter(valid_kmers.iter().filter(|&(x, _)| sel_kmers.contains(x)).map(|(x, y)| (*x, y.clone())))
}


pub struct ReadDb {
    // For each edge: read_id -> (edge_pos, read_pos, align_len)
    pub mappings: Vec<tada::utils::DMap<ReadType, (usize, usize, usize)>>,
    // Read_id -> Read
    pub reads: tada::utils::DMap<ReadType, graph_read::Read>
}

impl ReadDb {
    pub fn with_capacity(size: usize) -> ReadDb {
        let mut mappings = Vec::with_capacity(size);
        for _ in 0..size {
            let edge_mappings = tada::utils::new_dmap();
            mappings.push(edge_mappings);
        }
        ReadDb {
            reads: tada::utils::new_dmap(),
            mappings: mappings
        }
    }


    pub fn is_extension(edge_pos1: usize,
                        read_pos1: usize,
                        len1: usize,
                        edge_pos2: usize,
                        read_pos2: usize,
                        len2: usize,
                        dir: kmer::Dir)
                        -> bool {
        match dir {
            kmer::Dir::Right => edge_pos2 == 0 && read_pos2 == read_pos1 + len1 - K + 1,
            kmer::Dir::Left => edge_pos1 == 0 && read_pos1 == read_pos2 + len2 - K + 1,
        }
    }

    pub fn get_read(&self, read_id: ReadType) -> &graph_read::Read {
        self.reads.get(&read_id).unwrap()
    }
}

pub struct BarcodedEdge {
    pub id: u32,
    pub sequence: bitenc::BitEnc,
    /// Edges to the left of this edge. Vector of tuples
    /// (index_of_neighbor, direction of extension, reverse_complemented?)
    pub l_exts: Vec<(usize, kmer::Dir, bool)>,
    pub r_exts: Vec<(usize, kmer::Dir, bool)>,
    pub bcs: HashSet<u32>,
}

impl BarcodedEdge {
    pub fn new(vedge: &debruijn::VEdge, support: &HashSet<u32>) -> BarcodedEdge {
        BarcodedEdge {
            id: vedge.id,
            sequence: vedge.sequence.clone(),
            l_exts: vedge.l_exts.clone(),
            r_exts: vedge.r_exts.clone(),
            bcs: support.clone(),
        }
    }

    pub fn get_id(&self) -> u32 {
        return self.id;
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn support(&self) -> u32 {
        self.bcs.len() as u32
    }

    pub fn first_kmer(&self) -> kmer::Kmer {
        self.sequence.first_kmer()
    }

    pub fn last_kmer(&self) -> kmer::Kmer {
        self.sequence.last_kmer()
    }

    pub fn get_kmer(&self, pos: usize) -> kmer::Kmer {
        self.sequence.get_kmer(pos)
    }

    pub fn get(&self, pos: usize) -> Option<u8> {
        self.sequence.get(pos)
    }

    pub fn write_edge(&self, f: &mut Write, print_left: bool) {

        let color = ",color=green";

        writeln!(f,
                 "n{} [label=\"{} - {} - {}\"{},style=filled]",
                 self.id,
                 self.id,
                 self.sequence.len(),
                 self.bcs.len(),
                 color)
            .unwrap();

        for &(id, incoming_dir, _) in self.r_exts.iter() {
            let color = "red";
            let style = match incoming_dir {
                kmer::Dir::Left => "solid",
                kmer::Dir::Right => "dashed",
            };
            writeln!(f,
                     "n{} -> n{} [color={},style={}]",
                     self.id,
                     id,
                     color,
                     style)
                .unwrap();
        }

        if print_left {
            for &(id, incoming_dir, _) in self.l_exts.iter() {
                let color = "blue";
                let style = match incoming_dir {
                    kmer::Dir::Right => "solid",
                    kmer::Dir::Left => "dashed",
                };
                writeln!(f,
                         "n{} -> n{} [color={},style={}]",
                         self.id,
                         id,
                         color,
                         style)
                    .unwrap();
            }
        }
    }
}

pub struct EdgeGraph {
    pub edges: Vec<BarcodedEdge>,
    // From kmer to edge index and position within the edge
    pub kmer_map: tada::utils::DMap<kmer::Kmer, (usize, usize)>,
}

impl EdgeGraph {

    pub fn new(graph: &debruijn::TempGraph,
               valid_kmers: &tada::utils::DMap<tada::kmer::Kmer, Vec<u32>>)
               -> EdgeGraph {
        let edgedb = tada::debruijn::EdgeDb::new(graph);
        let mut edges = Vec::new();

        let mut kmer_map = tada::utils::new_dmap();

        // e is an Edge object: not very useful, we want to convert to VEdge which has pointers
        // to left and right edges.
        for (idx, _) in graph.iter().enumerate() {
            // Explore the graph at a depth of 1
            let mut ex = tada::debruijn::explore_graph(graph, &edgedb, idx, 1);

            for (_, vedge) in ex.drain() {
                let mut support = HashSet::new();

                for (kmer_idx, kmer) in vedge.sequence.kmers().iter().enumerate() {
                    match valid_kmers.get(&kmer) {
                        Some(id_list) => {
                            for id in id_list {
                                support.insert(*id);
                            }
                        }
                        None => panic!("kmer missing from valid_kmers"),
                    }
                    kmer_map.insert(*kmer, (idx, kmer_idx));
                }
                edges.push(BarcodedEdge::new(&vedge, &support));
            }
        }
        EdgeGraph {
            edges: edges,
            kmer_map: kmer_map,
        }
    }

    /// Get the list of edges that spell out the given sequence.
    ///
    /// * Return value:
    ///
    /// A tuple (path, match_start, match_len). Path is a vector
    /// (edge idx, start position in edge, start position in seq, alignment length)
    /// match_start is the position in the seq where matching started
    /// (this is the same as path[0].2  )
    /// match_len is the number of letters of sequence that were matched.
    pub fn trace_sequence(&self, seq: &bitenc::BitEnc)
                          -> (Vec<(usize, usize, usize, usize)>, usize, usize) {
        let mut path = Vec::new();
        let mut edge: usize = 0;
        let mut edge_pos: usize = 0;
        let mut seq_pos: usize = 0;
        let mut align_len: usize = 0;
        let mut total_match_len: usize = 0;

        let mut match_start = 0;
        let mut start_found = false;

        while !start_found && match_start < seq.len() - K {
            match self.kmer_map.get(&seq.get_kmer(match_start)) {
                Some(&(edge_idx, pos)) => {
                    edge = edge_idx;
                    edge_pos = pos;
                    align_len = K;
                    seq_pos = match_start;
                    total_match_len = K;
                    start_found = true;
                }
                None => {
                    match_start += 1;
                }
            };
        }

        if !start_found {
            return (path, seq.len(), 0);
        }

        while total_match_len + match_start < seq.len() {
            // We've matched up this whole edge. Test edge extensions.
            if edge_pos + align_len == self.edges[edge].len() {
                path.push((edge, edge_pos, seq_pos, align_len));
                let mut found = false;
                for &(ext, _, _) in self.edges[edge as usize].r_exts.iter() {
                    // Compare the first kmer of the extension to the last non-tested
                    // kmer of the sequence.
                    if self.edges[ext].first_kmer() == seq.get_kmer(seq_pos + align_len + 1 - K) {
                        edge = ext;
                        edge_pos = 0;
                        assert!(seq_pos + align_len + 1 >= K);
                        seq_pos = seq_pos + align_len + 1 - K;
                        align_len = K;
                        found = true;
                        // Match only grew by one.
                        total_match_len += 1;
                        break;
                    }
                }
                if !found {
                    break;
                }
            } else {
                let next_letter = self.edges[edge].get(edge_pos + align_len).unwrap();
                if seq.get(seq_pos + align_len).unwrap() == next_letter {
                    align_len += 1;
                    total_match_len += 1;
                } else {
                    // Mismatch found within the current edge.
                    path.push((edge, edge_pos, seq_pos, align_len));
                    break;
                }
            }
        }
        if total_match_len + match_start == seq.len() {
            path.push((edge, edge_pos, seq_pos, align_len));
        }
        if path.len() > 0 {
            assert!(match_start == path[0].2);
        }
        (path, match_start, total_match_len)
    }

    pub fn build_read_db(&self, reads: &Vec<&graph_read::Read>) -> ReadDb  {
        let mut db = ReadDb::with_capacity(self.edges.len());

        let mut nseqs = 0;
        for ref mut read in reads.iter() {
            let (new_path, _, match_len) = self.trace_sequence(&read.seq);
            if match_len < K {
                continue;
            }
            nseqs += 1;
            for &(edge, edge_pos, read_pos, align_len) in new_path.iter() {
                assert!(read_pos + align_len <= read.len());
                db.mappings[edge].insert(read.id, (edge_pos, read_pos, align_len));
            }

            assert!(!db.reads.contains_key(&(read.id)));
            db.reads.insert(read.id, (*(*read)).clone());
        }
        println!("Created ReadDb from {} reads", nseqs);
        db
    }

    pub fn len(&self) -> usize {
        self.edges.len()
    }

    fn get_edge(&self, edge_idx: u32) -> &BarcodedEdge {
        &self.edges[edge_idx as usize]
    }

    fn is_edge_source(&self, edge_idx: u32) -> bool {
        self.edges[edge_idx as usize].l_exts.is_empty()
    }

    pub fn extend_paths(&self,
                        starting_paths: &Vec<Vec<(u32, usize)>>,
                        dir: kmer::Dir,
                        max_depth: usize)
                        -> Vec<Vec<(u32, usize)>> {

        // Partial paths, not fully explored yet.
        let mut paths = Vec::new();
        for path in starting_paths {
            let mut path_nodes = HashSet::new();
            for &(node, _) in path {
                path_nodes.insert(node);
            }
            paths.push((path.clone(), path_nodes));
        }

        // Will collect all possible paths here.
        let mut final_paths = Vec::new();

        while !paths.is_empty() {
            let (path_to_explore, path_nodes) = paths.pop().unwrap();
            let (last_node, depth) = path_to_explore[path_to_explore.len() - 1];

            let mut path_extended = false;
            if depth < max_depth {
                let node_iter = match dir {
                    kmer::Dir::Left => self.edges[last_node as usize].l_exts.iter(),
                    kmer::Dir::Right => self.edges[last_node as usize].r_exts.iter(),
                };
                for &(ext, _, _) in node_iter {
                    let ext2 = ext as u32;
                    if !path_nodes.contains(&ext2) {
                        let mut new_path = path_to_explore.clone();
                        // The new node is always put at the end of the path (irrespective of dir).
                        // So the end of the path will always have the node most distant from the start.
                        new_path.push((ext as u32, depth + 1));
                        let mut new_path_nodes = path_nodes.clone();
                        new_path_nodes.insert(ext as u32);
                        paths.push((new_path, new_path_nodes));
                        path_extended = true;
                    }
                }
            }
            if !path_extended {
                final_paths.push(path_to_explore);
            }
        }
        final_paths
    }

    /// Gets all distinct paths from a given node in the given direction.
    ///
    /// # Return value
    ///
    /// A vector of paths. Each path is a vector of tuples (edge_idx, depth), where depth is the
    /// distance from the starting node. The first element of each path is always (start, 0).
    pub fn paths_from_node(&self,
                           start: u32,
                           dir: kmer::Dir,
                           max_depth: usize)
                           -> Vec<Vec<(u32, usize)>> {

        // Path initially consists of starting node.
        let new_path = vec![(start, 0 as usize)];
        let new_paths = vec![new_path];
        self.extend_paths(&new_paths, dir, max_depth)
    }
    
    
    fn get_read_support(&self,
                        new_path: &Vec<u32>,
                        read_db: &ReadDb,
                        dir: kmer::Dir,
                        umis: &HashSet<UmiType>,
                        align_params: &sw::AlignParams)
                        -> (tada::utils::DMap<ReadType, (usize, f64)>, f64) {

        let min_score = -1e10;
        let max_score = 1e10;

        let mut support_reads: tada::utils::DMap<ReadType, (usize, f64)> = tada::utils::new_dmap();

        let mut last_edge_idx = 0;
        let last_edge = new_path[last_edge_idx] as usize;
        for (&read_id, &(_, read_pos1, len1)) in read_db.mappings[last_edge].iter() {
            let read_umi = read_db.reads.get(&read_id).unwrap().umi;
            if umis.is_empty() || umis.contains(&read_umi) {
                let penalty_other_side = (align_params.clip as f64) * match dir {
                    kmer::Dir::Right => read_pos1 as f64,
                    kmer::Dir::Left => (read_db.get_read(read_id).len() - read_pos1 - len1) as f64,
                };
                support_reads.insert(read_id, (last_edge_idx, (align_params.match_score as f64) * (len1 as f64)  - penalty_other_side));
            }
        }

        while last_edge_idx < new_path.len() - 1 {

            let last_edge = new_path[last_edge_idx] as usize;
            let next_edge = new_path[last_edge_idx + 1] as usize;

            for (&read_id, &(edge_pos1, read_pos1, len1)) in read_db.mappings[last_edge].iter() {

                assert!(read_pos1 <= read_db.get_read(read_id).len());
                assert!(read_db.get_read(read_id).len() >= read_pos1 + len1);

                let read_umi = read_db.reads.get(&read_id).unwrap().umi;
                if !umis.is_empty() && !umis.contains(&read_umi) {
                    continue;
                }
                // Potential penalty for having to trim read.
                let penalty = (align_params.clip as f64) * match dir {
                    kmer::Dir::Left => read_pos1 as f64,
                    kmer::Dir::Right => (read_db.get_read(read_id).len() - read_pos1 - len1) as f64,
                };
                assert!(penalty >= 0.0 && penalty <= read_db.get_read(read_id).len() as f64);

                let penalty_other_side = (align_params.clip as f64) * match dir {
                    kmer::Dir::Right => read_pos1 as f64,
                    kmer::Dir::Left => (read_db.get_read(read_id).len() - read_pos1 - len1) as f64,
                };
                assert!(penalty_other_side >= 0.0 && penalty_other_side <= read_db.get_read(read_id).len() as f64);

                if read_db.mappings[next_edge].contains_key(&read_id) {
                    let &(edge_pos2, read_pos2, len2) =
                        read_db.mappings[next_edge].get(&read_id).unwrap();

                    if ReadDb::is_extension(edge_pos1,
                                            read_pos1,
                                            len1,
                                            edge_pos2,
                                            read_pos2,
                                            len2,
                                            dir) {
                        // if !support_reads.contains_key(&read) {
                        //    support_reads.insert(read, (last_edge_idx, len1 as f32 - penalty_other_side));
                        // }
                        // let curr_val = support_reads.get_mut(&read).unwrap();
                        // *curr_val = ((*curr_val).0, (*curr_val).1 + (len2 - K + 1) as f32);
                        if !support_reads.contains_key(&read_id) {
                            support_reads.insert(read_id, (last_edge_idx, (align_params.match_score as f64) * (len1 as f64) - penalty_other_side));
                        }
                        let curr_val = support_reads.get_mut(&read_id).unwrap();
                        if (*curr_val).1 > min_score && (*curr_val).1 < max_score {
                            assert!(len2 + 1 >= K);
                            *curr_val = ((*curr_val).0, (*curr_val).1 + (align_params.match_score as f64) * ((len2 + 1 - K) as f64));
                        }
                    } else {
                        // Assume that read was trimmed at the last position
                        if !support_reads.contains_key(&read_id) {
                            support_reads.insert(read_id, (last_edge_idx, (align_params.match_score as f64) * (len1 as f64) - penalty_other_side));
                        }
                        let curr_val = support_reads.get_mut(&read_id).unwrap();
                        if (*curr_val).1 > min_score && (*curr_val).1 < max_score {
                            *curr_val = ((*curr_val).0, (*curr_val).1 - penalty);
                        }
                    }
                } else {
                    // Assume that read was trimmed at the last position
                    if !support_reads.contains_key(&read_id) {
                        support_reads.insert(read_id, (last_edge_idx, (align_params.match_score as f64) * (len1 as f64) - penalty_other_side));
                    }
                    let curr_val = support_reads.get_mut(&read_id).unwrap();
                    if (*curr_val).1 > min_score && (*curr_val).1 < max_score {
                        *curr_val = ((*curr_val).0, (*curr_val).1 - penalty);
                    }
                }
            }
            last_edge_idx += 1;
        }
        let mut total_score = 0.0;
        for &align in support_reads.values() {
            if total_score > min_score && total_score < max_score {
                total_score = total_score + align.1;
            } else {
                break;
            }
        }

        (support_reads, total_score)
    }


    fn extension_qualities(&self, node: usize, read_db: &ReadDb, dir: kmer::Dir, rt_err: f64) -> Vec<(u32, u8)> {

        let node_iter = match dir {
            kmer::Dir::Left => self.edges[node].l_exts.iter(),
            kmer::Dir::Right => self.edges[node].r_exts.iter(),
        };

        let mut umi_reads : Vec<(UmiType, u8, u8)> = Vec::new();
        let mut ext_nodes = Vec::new();

        for (ext_idx, &(ext, _, _)) in node_iter.enumerate() {
            ext_nodes.push(ext as u32);

            assert!(ext_idx < 4);

            for (read_id, &(edge_pos1, read_pos1, len1)) in read_db.mappings[node].iter() {
                match read_db.mappings[ext].get(&read_id) {
                    Some(&(edge_pos2, read_pos2, len2)) => {
                        let read = read_db.reads.get(&read_id).unwrap();
                        let next_base_qual = read.quals[read_pos2];

                        if ReadDb::is_extension(edge_pos1, read_pos1, len1,
                                                edge_pos2, read_pos2, len2,
                                                dir){
                            umi_reads.push((read.umi, ext_idx as u8, next_base_qual));
                        }
                    },
                    None => {},
                }
            }
        }

        // Sort by UMI
        umi_reads.sort_by_key(|x| x.0);

        let quals = sw::pos_base_quals(&umi_reads, rt_err);
        let mut res = Vec::new();
        for idx in 0..ext_nodes.len() {
            res.push((ext_nodes[idx], quals[idx]));
        }
        res
    }

    /// Get paths from a given node based on qualities of branches.
    /// 
    /// Start from a node and extend in both directions. At each branch compute 
    /// the base qualities of all extensions. Follow all extensions that have 
    /// base qualities within qual_factor of the top extension and extend all
    /// paths constructed so far by all possible extensions. 
    ///
    /// Return value:
    /// A vector with tuples (path, support), one tuple for each computed path.
    /// - path is a vector of graph edges representing the path (which must contain "start").
    /// - support is a hashmap read_id -> (starting_node, score). Keys are read ids
    /// for reads on the corresponding path. starting_node is the node where the mapping
    /// of the read starts and score is (a rough approximation of the) SW score 
    /// of the read against the path.
    pub fn max_paths_by_qual(&self, start: u32, read_db: &ReadDb,
                            qual_factor: f32, rt_err: f64, min_qual: u8,
                            align_params: &sw::AlignParams)
                            -> Vec<(Vec<u32>, tada::utils::DMap<ReadType, (usize, f64)>)> {

        // Partial paths, not fully explored yet.
        let mut paths = Vec::new();
        let mut new_path = VecDeque::new();
        new_path.push_front(start);
        let mut new_path_nodes = HashSet::new();
        new_path_nodes.insert(start);

        // Path, set of path nodes, whether we can still extend to the left/right
        paths.push((new_path, new_path_nodes, true, true));

        let mut final_paths = Vec::new();

        while !paths.is_empty() {

            let (path_to_explore, path_nodes, can_extend_left, can_extend_right) = paths.pop().unwrap();

            // Get valid extensions in each direction
            let mut good_left_extensions = Vec::new();
            let mut good_right_extensions = Vec::new();

            if can_extend_left {
                let first_node = path_to_explore[0];

                let node_quals = self.extension_qualities(first_node as usize, read_db, kmer::Dir::Left, rt_err);
                let mut max_node_qual = std::f32::MIN;
                for (_, &(node, qual)) in node_quals.iter().enumerate() {
                    if !path_nodes.contains(&node) && qual as f32 > max_node_qual {
                        max_node_qual = qual as f32;
                    }
                }

                if max_node_qual >= (min_qual as f32) {
                    for (_, &(node, qual)) in node_quals.iter().enumerate() {
                        if !path_nodes.contains(&node) && qual as f32 >= qual_factor * (max_node_qual as f32) {
                            good_left_extensions.push(node);
                        }
                    }
                }
            }

            if can_extend_right {
                let last_node = path_to_explore[path_to_explore.len() - 1];
                let node_quals = self.extension_qualities(last_node as usize, read_db, kmer::Dir::Right, rt_err);
                let mut max_node_qual = std::f32::MIN;
                for (_, &(node, qual)) in node_quals.iter().enumerate() {
                    if !path_nodes.contains(&node) && qual as f32 > max_node_qual {
                        max_node_qual = qual as f32;
                    }
                }

                if max_node_qual >= (min_qual as f32) {
                    for (_, &(node, qual)) in node_quals.iter().enumerate() {
                        if !path_nodes.contains(&node) &&  qual as f32 >= qual_factor * (max_node_qual as f32) {
                            good_right_extensions.push(node);
                        }
                    }
                }
            }

            // Paths extended in one direction
            let mut tmp_paths = Vec::new();

            // If we determined that there are no good left extensions, then there's no point
            // trying to extend the starting node to the left next time we see it.
            let can_still_extend_left = !good_left_extensions.is_empty();
            let can_still_extend_right = !good_right_extensions.is_empty();

            if good_left_extensions.is_empty() {
                // Last element of the tuple is whether this is an extended path (true) or just a copy of
                // the initial path (false).
                tmp_paths.push((path_to_explore.clone(), path_nodes.clone(), can_still_extend_left, can_still_extend_right, false));
            } else {
                for &left_node in good_left_extensions.iter() {
                    let mut new_path = path_to_explore.clone();
                    new_path.push_front(left_node);
                    let mut new_path_nodes = path_nodes.clone();
                    new_path_nodes.insert(left_node);
                    tmp_paths.push((new_path, new_path_nodes, can_still_extend_left, can_still_extend_right, true));
                }
            }

            for &(ref tmp_path, ref tmp_nodes, l, r, extended) in tmp_paths.iter() {
                if good_right_extensions.is_empty() {
                    if extended {
                        paths.push((tmp_path.clone(), tmp_nodes.clone(), l, r));
                    } else {
                        // This will happen if there were no left extensions at all or no left
                        // extensions that didn't hit a cycle.
                        final_paths.push(tmp_path.iter().map(|x| *x).collect());
                    }
                } else {
                    for &right_node in good_right_extensions.iter() {
                        let mut new_path = tmp_path.clone();
                        new_path.push_back(right_node);
                        let mut new_path_nodes = tmp_nodes.clone();
                        new_path_nodes.insert(right_node);
                        paths.push((new_path, new_path_nodes, l, r));
                    }
                }
            }

        }

        let mut paths_with_scores = Vec::new();

        for ref path in final_paths.iter() {
            let empty = HashSet::new();
            let final_scores = self.get_read_support(&path,
                              &read_db,
                              kmer::Dir::Right, &empty, align_params).0;
            paths_with_scores.push(((*path).clone(), final_scores));
        }
        paths_with_scores
    }

    pub fn dfs_from_node(&self, start: u32, dir: kmer::Dir, max_depth: usize) -> Vec<u32> {
        let mut nodes = vec![start];
        let mut visited = HashSet::new();

        let mut nodes_to_visit = vec![(start, 0)];

        while !nodes_to_visit.is_empty() {
            let node_tuple = nodes_to_visit.pop().unwrap();
            let node = node_tuple.0;
            let depth = node_tuple.1;

            if !visited.contains(&node) {
                nodes.push(node);
                visited.insert(node);
                if depth < max_depth {
                    let edge_iter = match dir {
                        kmer::Dir::Right => self.edges[node as usize].r_exts.iter(),
                        kmer::Dir::Left => self.edges[node as usize].l_exts.iter(),
                    };
                    for &(ext, _, _) in edge_iter {
                        let ext2 = ext as u32;
                        if !visited.contains(&ext2) {
                            nodes_to_visit.push((ext2, depth + 1));
                        }
                    }
                }
            }
        }
        nodes
    }

    /// Get all nodes reachable from the given node, in either direction.
    /// If there are no cycles, the nodes will be returned in a topologically
    /// sorted order.
    pub fn search_from_node(&self, start: u32) -> Vec<u32> {
        // Each node is "seen" twice. The first time, we make sure its parents
        // are added. The second time we see the node, we know that its parents
        // have already been processed, so we can output it and continue with
        // the nodes children.
        let mut nodes = Vec::new();
        let mut seen = HashSet::new();
        let mut visited = HashSet::new();

        let mut nodes_to_visit = VecDeque::new();
        nodes_to_visit.push_front(start);

        while !nodes_to_visit.is_empty() {
            // peak at the node at the top of the deque.
            // If the node has already been seen, this means we've already
            // considered its parents. So remove and add to the output.
            let top = nodes_to_visit.pop_front().unwrap();
            if seen.contains(&top) {
                if !visited.contains(&top) {
                    nodes.push(top);
                }
                for &(ext, _, _) in self.edges[top as usize].r_exts.iter() {
                    let ext2 = ext as u32;
                    if !seen.contains(&ext2){
                        nodes_to_visit.push_back(ext2);
                    }
                }
                visited.insert(top);

            } else {
                // Do not remove the node. Instead add all its parents in
                // front of it and mark it as seen.
                nodes_to_visit.push_front(top);
                for &(ext, _, _) in self.edges[top as usize].l_exts.iter() {
                    let ext2 = ext as u32;
                    if !seen.contains(&ext2) {
                        nodes_to_visit.push_front(ext2);
                    }
                }
                seen.insert(top);
            }
        }
        nodes
    }


    pub fn get_all_paths(&self) -> Vec<Vec<u32>> {
        let mut paths = Vec::new();

        for (edge_idx, _) in self.edges.iter().enumerate() {
            if self.is_edge_source(edge_idx as u32) {
                let new_paths =
                    self.paths_from_node(edge_idx as u32, kmer::Dir::Right, self.edges.len());
                paths.append(&mut new_paths.iter()
                    .map(|x| x.iter().map(|y| y.0).collect())
                    .collect());
            }
        }
        paths
    }

    pub fn get_path_sequence(&self, path: &Vec<u32>) -> String {
        let mut seq: String = "".to_string();
        for (edge_idx, &edge) in path.iter().enumerate() {
            let mut new_seq = self.edges[edge as usize].sequence.to_dna_string();
            if edge_idx > 0 {
                new_seq = new_seq.chars().skip(K - 1).collect();
            }
            seq = seq + &new_seq;
        }
        seq
    }

    /// Gets connected components in the graph (ignoring direction of edges).
    pub fn connected_components(&self) -> Vec<Vec<u32>> {
        let mut components = Vec::new();
        let mut visited: HashSet<u32> = HashSet::new();

        for (_, edge) in self.edges.iter().enumerate() {
            if !visited.contains(&edge.id) {
                let new_edges = self.search_from_node(edge.id);
                for &new_edge in new_edges.iter() {
                    visited.insert(new_edge as u32);
                }
                components.push(new_edges);
            }
        }
        components
    }

    pub fn to_dot(&self, path: PathBuf, print_left: bool) {

        let mut f = File::create(path).expect("couldn't open file");

        writeln!(&mut f, "digraph {{").unwrap();
        for edge in self.edges.iter() {
            edge.write_edge(&mut f, print_left);
        }
        writeln!(&mut f, "}}").unwrap();
    }
}
