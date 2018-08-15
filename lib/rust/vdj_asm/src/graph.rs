//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//
use fxhash::FxHashMap;

use debruijn::{filter, kmer};
use debruijn::dna_string::DnaString;
use debruijn::{Exts, Dir, Kmer, Mer};
use debruijn::graph::{DebruijnGraph, Node};
use debruijn::compression::{SimpleCompress, compress_kmers};

use debruijn::Vmer;

use std::fmt::Debug;
use std::collections::{HashSet, VecDeque};
use constants::{UmiType, ReadType};
use sw;
use graph_read;
use std;

pub type Kmer1 = kmer::Kmer32;

/// Build a graph from a vector of graph_read::Read objects.
/// min_kmers: min number of kmer occurrences to consider a kmer
/// max_kmers: maximum number of kmers to consider. If the total number of observed
/// kmers (that pass the min_kmers cutoff) is greater than max_kmers, then we'll
/// get the top max_kmers kmers by frequency.
pub fn build_graph(reads: &Vec<&graph_read::Read>,
                   min_kmers:usize, max_kmers: usize) -> IndexedGraph<Kmer1, Vec<u32>> {

    let mut seqs = Vec::new();
    for r in reads {
        seqs.push((r.seq.clone(), Exts::empty(), r.umi));
    }

    let summarizer = filter::CountFilterSet::new(min_kmers);
    let mut valid_kmers : Vec<(Kmer1, (Exts, Vec<u32>))> = {
        let (kmer_hashmap, _) = filter::filter_kmers::<Kmer1, _, _, _, _>(&seqs, &Box::new(summarizer), true, false, 4);
        // convert hashmap back to vec
        kmer_hashmap.iter().map(|(k,e,d)| (*k, (*e, d.clone()))).collect()
    };

    println!("Kmers accepted: {}", valid_kmers.len());

    // if we have a ton of valid kmers, drop the ones with least coverage
    if valid_kmers.len() > max_kmers {
        println!("Observed > {} kmers. Truncating", max_kmers);
    }

    valid_kmers.sort_by_key(|x| -((x.1).1.len() as isize));
    valid_kmers.truncate(max_kmers);
    valid_kmers.sort();
    filter::remove_censored_exts(true, &mut valid_kmers);

    let cmp = SimpleCompress::new(|mut a: Vec<u32>, b: &Vec<u32>| { a.extend(b); a.sort(); a.dedup(); a });
    let graph = compress_kmers(true, cmp, &valid_kmers).finish();

    println!("# Edges: {}", graph.len());
    IndexedGraph::new(graph)
}

/*
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
*/

pub struct ReadDb {
    // For each edge: read_id -> (edge_pos, read_pos, align_len)
    pub mappings: Vec<FxHashMap<ReadType, (usize, usize, usize)>>,
    // Read_id -> Read
    pub reads: FxHashMap<ReadType, graph_read::Read>
}

impl ReadDb {
    pub fn with_capacity(size: usize) -> ReadDb {
        let mut mappings = Vec::with_capacity(size);
        for _ in 0..size {
            let edge_mappings = FxHashMap::default();
            mappings.push(edge_mappings);
        }
        ReadDb {
            reads: FxHashMap::default(),
            mappings: mappings
        }
    }


    pub fn is_extension(k: usize,
                        edge_pos1: usize,
                        read_pos1: usize,
                        len1: usize,
                        edge_pos2: usize,
                        read_pos2: usize,
                        len2: usize,
                        dir: Dir)
                        -> bool {
        match dir {
            Dir::Right => edge_pos2 == 0 && read_pos2 == read_pos1 + len1 - k + 1,
            Dir::Left => edge_pos1 == 0 && read_pos1 == read_pos2 + len2 - k + 1,
        }
    }

    pub fn get_read(&self, read_id: ReadType) -> &graph_read::Read {
        self.reads.get(&read_id).unwrap()
    }
}


pub struct IndexedGraph<K: Kmer, D> {
    pub graph: DebruijnGraph<K, D>,
    // From kmer to edge index and position within the edge
    pub kmer_map: FxHashMap<K, (usize, usize)>,
}

impl<K: Kmer, D: Debug> IndexedGraph<K, D> {

    pub fn new(graph: DebruijnGraph<K,D>) -> Self {

        let mut kmer_map: FxHashMap<K, (usize, usize)> = FxHashMap::default();

        // e is an Edge object: not very useful, we want to convert to VEdge which has pointers
        // to left and right edges.
        for idx in 0 .. graph.len() {
            // Explore the graph at a depth of 1
            let node = graph.get_node(idx);
            
            for (kmer_idx, kmer) in node.sequence().iter_kmers().enumerate() {
                kmer_map.insert(kmer, (idx, kmer_idx));
            }
        }

        IndexedGraph {
            graph: graph,
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
    pub fn trace_sequence(&self, seq: &DnaString)
                          -> (Vec<(usize, usize, usize, usize)>, usize, usize) {
        let mut path = Vec::new();
        let mut node_id: usize = 0;
        let mut node_pos: usize = 0;
        let mut seq_pos: usize = 0;
        let mut align_len: usize = 0;
        let mut total_match_len: usize = 0;

        let mut match_start = 0;
        let mut start_found = false;

        while !start_found && match_start < seq.len() - K::k() {
            match self.kmer_map.get(&seq.get_kmer(match_start)) {
                Some(&(node_idx, pos)) => {
                    node_id = node_idx;
                    node_pos = pos;
                    align_len = K::k();
                    seq_pos = match_start;
                    total_match_len = K::k();
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
            let node = self.graph.get_node(node_id);
            // We've matched up this whole edge. Test edge extensions.
            if node_pos + align_len == node.len() {
                path.push((node_id, node_pos, seq_pos, align_len));
                let mut found = false;
                for (ext_node_id, _, _) in node.r_edges() {
                    let ext_node = self.graph.get_node(ext_node_id);
                    // Compare the first kmer of the extension to the last non-tested
                    // kmer of the sequence.
                    let ext_kmer: K = ext_node.sequence().first_kmer();
                    let seq_kmer: K = seq.get_kmer(seq_pos + align_len + 1 - K::k());
                    
                    if ext_kmer == seq_kmer {
                        node_id = ext_node_id;

                        node_pos = 0;
                        assert!(seq_pos + align_len + 1 >= K::k());
                        seq_pos = seq_pos + align_len + 1 - K::k();
                        align_len = K::k();
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
                let next_letter = node.sequence().get(node_pos + align_len);
                if seq.get(seq_pos + align_len) == next_letter {
                    align_len += 1;
                    total_match_len += 1;
                } else {
                    // Mismatch found within the current edge.
                    path.push((node_id, node_pos, seq_pos, align_len));
                    break;
                }
            }
        }
        if total_match_len + match_start == seq.len() {
            path.push((node_id, node_pos, seq_pos, align_len));
        }
        if path.len() > 0 {
            assert!(match_start == path[0].2);
        }
        (path, match_start, total_match_len)
    }

    pub fn build_read_db(&self, reads: &Vec<&graph_read::Read>) -> ReadDb  {
        let mut db = ReadDb::with_capacity(self.graph.len());

        let mut nseqs = 0;
        for ref mut read in reads.iter() {
            let (new_path, _, match_len) = self.trace_sequence(&read.seq);
            if match_len < K::k() {
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
        self.graph.len()
    }

    pub fn get_node(&self, node_idx: usize) -> Node<K, D> {
        self.graph.get_node(node_idx)
    }

    fn is_edge_source(&self, node_idx: usize) -> bool {
        self.graph.get_node(node_idx).l_edges().len() == 0
    }

    pub fn extend_paths(&self,
                        starting_paths: &Vec<Vec<(usize, usize)>>,
                        dir: Dir,
                        max_depth: usize)
                        -> Vec<Vec<(usize, usize)>> {

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
                    Dir::Left => self.graph.get_node(last_node).l_edges(),
                    Dir::Right => self.graph.get_node(last_node).r_edges(),
                };
                for (ext, _, _) in node_iter {
                    if !path_nodes.contains(&ext) {
                        let mut new_path = path_to_explore.clone();
                        // The new node is always put at the end of the path (irrespective of dir).
                        // So the end of the path will always have the node most distant from the start.
                        new_path.push((ext, depth + 1));
                        let mut new_path_nodes = path_nodes.clone();
                        new_path_nodes.insert(ext);
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
                           start: usize,
                           dir: Dir,
                           max_depth: usize)
                           -> Vec<Vec<(usize, usize)>> {

        // Path initially consists of starting node.
        let new_path = vec![(start, 0)];
        let new_paths = vec![new_path];
        self.extend_paths(&new_paths, dir, max_depth)
    }


    fn get_read_support(&self,
                        new_path: &Vec<usize>,
                        read_db: &ReadDb,
                        dir: Dir,
                        umis: &HashSet<UmiType>,
                        scoring: &sw::Scoring)
                        -> (FxHashMap<ReadType, (usize, f64)>, f64) {

        let min_score = -1e10;
        let max_score = 1e10;

        let match_score = scoring.match_fn.match_score as f64;
        let clip_penalty = -scoring.xclip_prefix as f64;

        let mut support_reads: FxHashMap<ReadType, (usize, f64)> = FxHashMap::default();

        let mut last_edge_idx = 0;
        let last_edge = new_path[last_edge_idx];
        for (&read_id, &(_, _, len1)) in read_db.mappings[last_edge].iter() {
            let read_umi = read_db.reads.get(&read_id).unwrap().umi;
            if umis.is_empty() || umis.contains(&read_umi) {
                support_reads.insert(read_id, (last_edge_idx, match_score * (len1 as f64)  - clip_penalty));
            }
        }

        while last_edge_idx < new_path.len() - 1 {

            let last_edge = new_path[last_edge_idx];
            let next_edge = new_path[last_edge_idx + 1];

            for (&read_id, &(edge_pos1, read_pos1, len1)) in read_db.mappings[last_edge].iter() {

                assert!(read_pos1 <= read_db.get_read(read_id).len());
                assert!(read_db.get_read(read_id).len() >= read_pos1 + len1);

                let read_umi = read_db.reads.get(&read_id).unwrap().umi;
                if !umis.is_empty() && !umis.contains(&read_umi) {
                    continue;
                }

                if read_db.mappings[next_edge].contains_key(&read_id) {
                    let &(edge_pos2, read_pos2, len2) =
                        read_db.mappings[next_edge].get(&read_id).unwrap();

                    if ReadDb::is_extension(K::k(),
                                            edge_pos1,
                                            read_pos1,
                                            len1,
                                            edge_pos2,
                                            read_pos2,
                                            len2,
                                            dir) {

                        if !support_reads.contains_key(&read_id) {
                            support_reads.insert(read_id, (last_edge_idx, match_score * (len1 as f64) - clip_penalty));
                        }
                        let curr_val = support_reads.get_mut(&read_id).unwrap();
                        if (*curr_val).1 > min_score && (*curr_val).1 < max_score {
                            assert!(len2 + 1 >= K::k());
                            *curr_val = ((*curr_val).0, (*curr_val).1 + match_score * ((len2 + 1 - K::k()) as f64));
                        }
                    } else {
                        // Assume that read was trimmed at the last position
                        if !support_reads.contains_key(&read_id) {
                            support_reads.insert(read_id, (last_edge_idx, match_score * (len1 as f64) - clip_penalty));
                        }
                        let curr_val = support_reads.get_mut(&read_id).unwrap();
                        if (*curr_val).1 > min_score && (*curr_val).1 < max_score {
                            *curr_val = ((*curr_val).0, (*curr_val).1 - clip_penalty);
                        }
                    }
                } else {
                    // Assume that read was trimmed at the last position
                    if !support_reads.contains_key(&read_id) {
                        support_reads.insert(read_id, (last_edge_idx, match_score * (len1 as f64) - clip_penalty));
                    }
                    let curr_val = support_reads.get_mut(&read_id).unwrap();
                    if (*curr_val).1 > min_score && (*curr_val).1 < max_score {
                        *curr_val = ((*curr_val).0, (*curr_val).1 - clip_penalty);
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


    fn extension_qualities(&self, node: usize, read_db: &ReadDb, dir: Dir, rt_err: f64,
            ext_qual_lookup: &mut FxHashMap<(usize, i8), Vec<(usize, u8)> >) -> Vec<(usize, u8)> {

        let key = match dir {
            Dir::Left => (node, -1),
            Dir::Right => (node, 1),
        };

        if ext_qual_lookup.contains_key(&key) {
            return ext_qual_lookup.get(&key).unwrap().clone();
        }

        let node_iter = match dir {
            Dir::Left => self.graph.get_node(node).l_edges(),
            Dir::Right => self.graph.get_node(node).r_edges(),
        };

        let mut umi_reads : Vec<(UmiType, u8, u8)> = Vec::new();
        let mut ext_nodes = Vec::new();

        for (ext_idx, (ext, _, _)) in node_iter.into_iter().enumerate() {
            ext_nodes.push(ext);

            assert!(ext_idx < 4);

            for (read_id, &(edge_pos1, read_pos1, len1)) in read_db.mappings[node].iter() {
                match read_db.mappings[ext].get(&read_id) {
                    Some(&(edge_pos2, read_pos2, len2)) => {
                        let read = read_db.reads.get(&read_id).unwrap();
                        let next_base_qual = read.quals[read_pos2];

                        if ReadDb::is_extension(K::k(), edge_pos1, read_pos1, len1,
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

        let quals = sw::pos_base_quals(&umi_reads, rt_err, true);
        let mut res = Vec::new();
        for idx in 0..ext_nodes.len() {
            res.push((ext_nodes[idx], quals[idx]));
        }
        ext_qual_lookup.insert(key, res.clone());
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
    pub fn max_paths_by_qual(&self, start: usize, read_db: &ReadDb,
                            qual_factor: f32, rt_err: f64, min_qual: u8,
                            scoring: &sw::Scoring, mut ext_qual_lookup: &mut FxHashMap<(usize, i8), Vec<(usize, u8)> >)
                            -> Vec<(Vec<usize>, FxHashMap<ReadType, (usize, f64)>)> {

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

                let node_quals = self.extension_qualities(first_node, read_db, Dir::Left, rt_err, &mut ext_qual_lookup);
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
                let node_quals = self.extension_qualities(last_node, read_db, Dir::Right, rt_err, &mut ext_qual_lookup);
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
                              Dir::Right, &empty, scoring).0;
            paths_with_scores.push(((*path).clone(), final_scores));
        }
        paths_with_scores
    }

    pub fn dfs_from_node(&self, start: usize, dir: Dir, max_depth: usize) -> Vec<usize> {
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
                        Dir::Left => self.graph.get_node(node).l_edges(),
                        Dir::Right => self.graph.get_node(node).r_edges(),
                    };
                    for (ext, _, _) in edge_iter {
                        if !visited.contains(&ext) {
                            nodes_to_visit.push((ext, depth + 1));
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
    pub fn search_from_node(&self, start: usize) -> Vec<usize> {
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
                for (ext, _, _) in self.graph.get_node(top).r_edges() {
                    if !seen.contains(&ext){
                        nodes_to_visit.push_back(ext);
                    }
                }
                visited.insert(top);

            } else {
                // Do not remove the node. Instead add all its parents in
                // front of it and mark it as seen.
                nodes_to_visit.push_front(top);
                for (ext, _, _) in self.graph.get_node(top).l_edges() {
                    if !seen.contains(&ext) {
                        nodes_to_visit.push_front(ext);
                    }
                }
                seen.insert(top);
            }
        }
        nodes
    }

    pub fn get_all_paths(&self) -> Vec<Vec<usize>> {
        let mut paths = Vec::new();

        for node_idx in 0 .. self.graph.len() {
            if self.is_edge_source(node_idx) {
                let new_paths = self.paths_from_node(node_idx, Dir::Right, self.graph.len());
                paths.append(&mut new_paths.iter()
                    .map(|x| x.iter().map(|y| y.0).collect())
                    .collect());
            }
        }
        paths
    }

    pub fn get_path_sequence(&self, path: &Vec<usize>) -> String {
        let mut seq: String = "".to_string();
        for (node_idx, &node_id) in path.iter().enumerate() {
            let mut new_seq = self.get_node(node_id).sequence().to_dna_string();
            if node_idx > 0 {
                new_seq = new_seq.chars().skip(K::k() - 1).collect();
            }
            seq = seq + &new_seq;
        }
        seq
    }

    /// Gets connected components in the graph (ignoring direction of edges).
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let mut components = Vec::new();
        let mut visited: HashSet<usize> = HashSet::new();

        for node_id in 0 .. self.graph.len() {

            if !visited.contains(&node_id) {
                let new_edges = self.search_from_node(node_id);
                for &new_edge in new_edges.iter() {
                    visited.insert(new_edge);
                }
                components.push(new_edges);
            }
        }
        components
    }

    /*
    pub fn to_dot(&self, path: PathBuf, print_left: bool) {

        let mut f = File::create(path).expect("couldn't open file");

        writeln!(&mut f, "digraph {{").unwrap();
        for edge in self.edges.iter() {
            edge.write_edge(&mut f, print_left);
        }
        writeln!(&mut f, "}}").unwrap();
    }
    */
}
