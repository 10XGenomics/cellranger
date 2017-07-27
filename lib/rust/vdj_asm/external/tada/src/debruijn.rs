//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Utilities for building DeBruijn graphs efficiently
use utils;
use kmer::{self, K, Kmer, Dir, Exts, Bsp, Lmer, complement};
use std::collections::{VecDeque, HashSet};
use bitenc;
use bwt;
use utils::{MultiVec, DMap, new_dmap, LMap, new_lmap};

//use linked_hash_map::LinkedHashMap;
use std::hash::BuildHasherDefault;
use std::default::Default;
use std::hash::{Hash, Hasher};
use utils::FnvHasher;

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::iter::FromIterator;
use rayon::prelude::*;



#[inline]
fn flip_ext(flip: bool, base: u8, side: Dir) -> (u8, Dir) {
    if flip {
        (3 - base, side.flip())
    } else {
        (base, side)
    }
}

#[inline(never)]
fn extend(_kmer_exts: &mut DMap<Kmer, Exts>, kmer: Kmer, _base: u8, _dir: Dir, allow_rc: bool) {
    let mut min_kmer = kmer;
    let mut flip = false;
    if allow_rc {
        flip = kmer.min_rc().1;
        min_kmer = kmer.min_rc().0;
    }
    let (base, dir) = flip_ext(flip, _base, _dir);

    match _kmer_exts.get_mut(&min_kmer) {
        Some(exts) => *exts = exts.set(dir, base),
        None => {}
    }
}

#[inline(never)]
fn add_extensions(kmer_exts: &mut DMap<Kmer, Exts>, censor: &HashSet<Kmer>, k1: Kmer, k2: Kmer, allow_rc: bool) {

    if censor.contains(&k1) || censor.contains(&k2) {
        return
    }

    let ext_right = k2.get(kmer::K - 1);
    extend(kmer_exts, k1, ext_right, Dir::Right, allow_rc);

    let ext_left = k1.get(0);
    extend(kmer_exts, k2, ext_left, Dir::Left, allow_rc);
}


/// Find all possible extensions of each Kmer
///
/// Annotate each kmer in this partition with the set of edges
/// from it to other nodes in the full DeBruijn graph
/// If the following node is in the shard, we will encounter it locally,
/// if it's in a different shard, it will be encoded in the
/// Exts of the Bsp's that generated the kmer
#[inline(never)]
pub fn kmer_extensions<T, S>(valid_kmers: &DMap<Kmer, T>, all_kmers: &DMap<Kmer, S>, bsps: &Vec<Bsp>, allow_rc: bool) -> DMap<Kmer, Exts> {
    let mut kmer_exts: DMap<Kmer, Exts> = new_dmap();
    kmer_exts.reserve(valid_kmers.len() + 1000);

    // Censored kmers exist in this partition but were filtered
    let censored = {
        let mut cs = HashSet::new();
        for k in all_kmers.keys() {
            if !valid_kmers.contains_key(k) {
                cs.insert(*k);
            }
        }
        cs
    };

    for (&kmer, _) in valid_kmers.iter() {
        kmer_exts.insert(kmer, Exts::empty());
    }

    debug!("Processing BSPs: {}", bsps.len());
    let mut nbsps = 0;

    for b in bsps {
        let kmers = b.kmers();

        match b.extend_kmer(Dir::Left) {
            Some(ext_left) => add_extensions(&mut kmer_exts, &censored, ext_left, kmers[0], allow_rc),
            None => {}
        }

        match b.extend_kmer(Dir::Right) {
            Some(ext_right) => {
                add_extensions(&mut kmer_exts, &censored, (kmers[(kmers.len() - 1)]), ext_right, allow_rc)
            }
            None => {}
        }

        for i in 0..(kmers.len() - 1) {
            add_extensions(&mut kmer_exts, &censored, kmers[i], kmers[i + 1], allow_rc)
        }

        if nbsps % 1000000 == 0 {
            debug!("Processing BSPs: {} done", nbsps);
        }
        nbsps += 1;
    }

    kmer_exts
}

#[inline(never)]
pub fn barcodes_for_sedge(bc_vec: &mut Vec<u32>,
                          sedge: Lmer,
                          kmer_barcodes: &DMap<Kmer, Vec<u32>>) {

    bc_vec.clear();
    for kmer in sedge.iter_kmers() {
        match kmer_barcodes.get(&kmer.min_rc().0) {
            Some(v) => {
                for &bc in v {
                    bc_vec.push(bc);
                }
            }
            None => (),
        }
    }
}


#[derive(Copy, Clone)]
pub enum ExtMode {
    Unique(Kmer, Dir, Exts),
    Terminal(Exts),
}

/// Explore the graph from a Kmer
///
/// Attempt to extend kmer v in direction dir. Return:
///  - Unique(nextKmer, nextDir) if a single unique extension
///    is possible.  nextDir indicates the direction to extend nextDir
///    to preserve the direction of the extension.
/// - Term(ext) indicating the extensions at this end of the line
fn try_extend_kmer(available_kmers: &LMap<Kmer, ()>,
                   kmer_exts: &DMap<Kmer, Exts>,
                   v: Kmer,
                   dir: Dir, allow_rc: bool)
                   -> ExtMode {
    let exts = kmer_exts.get(&v).expect("didn't have kmer");
    if exts.num_ext_dir(dir) != 1 {
        ExtMode::Terminal(exts.single_dir(dir))
    } else {
        // Get the next kmer
        let ext_base = exts.get_unique_extension(dir).expect("should be unique");
        let mut next_kmer = v.extend(ext_base, dir);
        let mut do_flip = false;

        if allow_rc {
            do_flip = next_kmer.min_rc().1;
            next_kmer = next_kmer.min_rc().0;
        }

        let next_dir = dir.cond_flip(do_flip);

        // We can include this kmer in the line if:
        // a) it exists in the partition, and is still unused
        // b) the kmer we go to has a unique extension back in our direction

        if !available_kmers.contains_key(&next_kmer) {
            // This kmer isn't in this partition, or we've already used it
            return ExtMode::Terminal(exts.single_dir(dir));
        }

        // Direction we're approaching the new kmer from
        let new_incoming_dir = dir.flip().cond_flip(do_flip);
        let next_kmer_exts = kmer_exts.get(&next_kmer).expect("must have kmer");
        let incoming_count = next_kmer_exts.num_ext_dir(new_incoming_dir);
        let outgoing_exts = next_kmer_exts.single_dir(new_incoming_dir.flip());

        if incoming_count == 0 {
            panic!("unreachable");
        } else if incoming_count == 1 {
            // We have a unique path to next_kmer -- include it
            ExtMode::Unique(next_kmer, next_dir, outgoing_exts)
        } else {
            // there's more than one path
            // into the target kmer - don't include it
            ExtMode::Terminal(exts.single_dir(dir))
        }
    }
}



/// Build the maximal line starting at kmer in direction dir, at most max_dist long.
///
/// Also return the extensions at the end of this line.
/// Sub-lines break if their extensions are not available in this shard
#[inline(never)]
fn extend_sedge(kmer_exts: &DMap<Kmer, Exts>,
                available_kmers: &mut LMap<Kmer, ()>,
                kmer: Kmer,
                start_dir: Dir,
                max_dist: usize, allow_rc: bool)
                -> (Vec<(Kmer, Dir)>, Exts) {
    let mut current_dir = start_dir;
    let mut current_kmer = kmer;
    let mut path = Vec::new();
    let mut final_exts: Exts; // must get set below

    if max_dist == 0 {
        let first_exts = kmer_exts.get(&current_kmer).expect("didn't have kmer");
        return (path, first_exts.single_dir(start_dir));
    }

    let _ = available_kmers.remove(&kmer);

    loop {
        let ext_result = try_extend_kmer(available_kmers, kmer_exts, current_kmer, current_dir, allow_rc);

        match ext_result {
            ExtMode::Unique(next_kmer, next_dir, next_ext) => {
                path.push((next_kmer, next_dir));
                available_kmers.remove(&next_kmer);
                current_kmer = next_kmer;
                current_dir = next_dir;
                final_exts = next_ext
            }
            ExtMode::Terminal(ext) => {
                final_exts = ext;
                break;
            }
        }

        if path.len() >= max_dist {
            break;
        }
    }

    (path, final_exts)
}

/// Build the edge surrounding a kmer
#[inline(never)]
fn build_sedge(kmer_exts: &DMap<Kmer, Exts>,
               available_kmers: &mut LMap<Kmer, ()>,
               seed: Kmer, allow_rc: bool)
               -> (Lmer, Exts) {
    let (l_path, l_ext) = extend_sedge(kmer_exts, available_kmers, seed, Dir::Left, 92 - K, allow_rc);
    let (r_path, r_ext) = extend_sedge(kmer_exts,
                                       available_kmers,
                                       seed,
                                       Dir::Right,
                                       92 - K - l_path.len(), allow_rc);

    let mut edge_seq = VecDeque::new();
    for i in 0..K {
        edge_seq.push_back(seed.get(i));
    }

    // Add on the left path
    for &(next_kmer, dir) in l_path.iter() {
        let kmer = match dir {
            Dir::Left => next_kmer,
            Dir::Right => next_kmer.rc(),
        };

        edge_seq.push_front(kmer.get(0))
    }

    // Add on the right path
    for &(next_kmer, dir) in r_path.iter() {
        let kmer = match dir {
            Dir::Left => next_kmer.rc(),
            Dir::Right => next_kmer,
        };

        edge_seq.push_back(kmer.get(K - 1))
    }

    let left_extend = match l_path.last() {
        None => l_ext,
        Some(&(_, Dir::Left)) => l_ext,
        Some(&(_, Dir::Right)) => l_ext.complement(),
    };

    let right_extend = match r_path.last() {
        None => r_ext,
        Some(&(_, Dir::Left)) => r_ext.complement(),
        Some(&(_, Dir::Right)) => r_ext,
    };

    let edge = Lmer::new(edge_seq.iter());
    (edge, Exts::from_single_dirs(left_extend, right_extend))
}

/// Build all sedges until all kmers are exhausted
#[inline(never)]
pub fn build_sedges(kmer_exts: &DMap<Kmer, Exts>, allow_rc: bool) -> Vec<(Lmer, Exts)> {
    let mut edges: Vec<(Lmer, Exts)> = Vec::with_capacity(kmer_exts.len() / 16);

    let h = BuildHasherDefault::<FnvHasher>::default();
    let mut available_kmers: LMap<Kmer, ()> = LMap::with_capacity_and_hasher(kmer_exts.len(), h);
    available_kmers.extend(kmer_exts.keys().map(|k| (k.clone(), ())));

    loop {

        if available_kmers.len() % 100000 == 0
        {
            debug!("Building sedges: kmers left: {}", available_kmers.len());
        }

        match available_kmers.front() {
            Some((&k, _)) => {
                let sedge: (Lmer, Exts) = build_sedge(kmer_exts, &mut available_kmers, k, allow_rc);
                edges.push(sedge)
            }
            None => break,
        }
    }

    edges
}

/// Useful structure searching kmers in a graph.
pub struct SedgeDb {
    sedges: Vec<(Lmer, Exts)>,
    /// sedge indices sorted by the first kmer
    left_order: Vec<u32>,
    right_order: Vec<u32>,
}

impl SedgeDb {

    #[inline(never)]
    pub fn new(sedges: Vec<(Lmer, Exts)>) -> SedgeDb {

        let mut left_sort: Vec<u32> = Vec::with_capacity(sedges.len());
        let mut right_sort: Vec<u32> = Vec::with_capacity(sedges.len());
        for idx in 0..sedges.len()
        {
            left_sort.push(idx as u32);
            right_sort.push(idx as u32);
        }

        left_sort.sort_by_key(|idx| sedges[*idx as usize].0.first_kmer());
        right_sort.sort_by_key(|idx| sedges[*idx as usize].0.last_kmer());

        SedgeDb {
            sedges: sedges,
            left_order: left_sort,
            right_order: right_sort,
        }
    }

    fn get(&self, idx: usize) -> &(Lmer, Exts) {
        self.sedges.get(idx).expect("must have sedge")
    }

    pub fn search_kmer(&self, kmer: Kmer, side: Dir) -> Option<usize> {
        match side {
            Dir::Left => {
                let pos = self.left_order.binary_search_by_key(&kmer, |idx| { self.sedges[*idx as usize].0.first_kmer() });
                match pos {
                    Ok(idx) => Some(self.left_order[idx] as usize),
                    _ => None,
                }
            },
            Dir::Right => {
                let pos = self.right_order.binary_search_by_key(&kmer, |idx| { self.sedges[*idx as usize].0.last_kmer() });
                match pos {
                    Ok(idx) => Some(self.right_order[idx] as usize),
                    _ => None,
                }
            }
        }
    }

    pub fn find_link(&self, kmer: Kmer, dir: Dir) -> Option<(usize, Dir, bool)> {
        let rc = kmer.rc();

        // Only test self-consistent paths through
        // the edges
        // Avoids problems due to single kmer edges
        // (dir, next_side_incoming, rc)
        // (Dir::Left, Dir::Right, false) => true,
        // (Dir::Left, Dir::Left,  true) => true,
        // (Dir::Right, Dir::Left, false) => true,
        // (Dir::Right, Dir::Right, true) => true,

        match dir {
            Dir::Left => {
                match self.search_kmer(kmer, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, false)),
                    _ => (),
                }

                match self.search_kmer(rc, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, true)),
                    _ => (),
                }
            }

            Dir::Right => {
                match self.search_kmer(kmer, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, false)),
                    _ => (),
                }

                match self.search_kmer(rc, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, true)),
                    _ => (),
                }
            }
        }

        return None;
    }

    pub fn print(&self) {
        for l in 0..self.left_order.len() {
            let (bases, exts) = self.sedges[l];
            let l = bases.term_kmer(Dir::Left);
            let r = bases.term_kmer(Dir::Right);
            println!("L: {:?}, R {:?}, exts: {:?}", l, r, exts);
        }

        println!("left sort: {:?}", self.left_order);
        println!("rght sort: {:?}", self.right_order);
    }
}


// FIXME -- make generic w/ SedgeDb
pub struct EdgeDb<'a> {
    left_order: Vec<u32>,
    right_order: Vec<u32>,
    graph: &'a TempGraph
}

impl<'a> EdgeDb<'a> {

    #[inline(never)]
    pub fn new(graph: &TempGraph) -> EdgeDb {
        let mut left_sort: Vec<u32> = Vec::with_capacity(graph.start_pos.len());
        let mut right_sort: Vec<u32> = Vec::with_capacity(graph.start_pos.len());
        for i in 0 .. graph.start_pos.len()
        {
            left_sort.push(i as u32);
            right_sort.push(i as u32);
        }

        left_sort.sort_by_key( |idx| graph.sequence.get_kmer(graph.start_pos[*idx as usize]));
        right_sort.sort_by_key(|idx| graph.sequence.get_kmer(graph.start_pos[*idx as usize] + (graph.edge_len[*idx as usize] as usize) - K));

        EdgeDb {
            left_order: left_sort,
            right_order: right_sort,
            graph: graph,
        }
    }

    pub fn search_kmer(&self, kmer: Kmer, side: Dir) -> Option<usize> {
        match side {
            Dir::Left => {
                let pos = self.left_order.binary_search_by_key(&kmer,
                    |idx| self.graph.sequence.get_kmer(self.graph.start_pos[*idx as usize]));

                match pos {
                    Ok(idx) => Some(self.left_order[idx] as usize),
                    _ => None,
                }
            },
            Dir::Right => {
                let pos = self.right_order.binary_search_by_key(&kmer,
                    |idx| self.graph.sequence.get_kmer(self.graph.start_pos[*idx as usize] + (self.graph.edge_len[*idx as usize] as usize) - K));
                match pos {
                    Ok(idx) => Some(self.right_order[idx] as usize),
                    _ => None,
                }
            }
        }
    }

    pub fn find_link(&self, kmer: Kmer, dir: Dir) -> Option<(usize, Dir, bool)> {
        let rc = kmer.rc();

        // Only test self-consistent paths through
        // the edges
        // Avoids problems due to single kmer edges
        // (dir, next_side_incoming, rc)
        // (Dir::Left, Dir::Right, false) => true,
        // (Dir::Left, Dir::Left,  true) => true,
        // (Dir::Right, Dir::Left, false) => true,
        // (Dir::Right, Dir::Right, true) => true,
        //

        match dir {
            Dir::Left => {
                match self.search_kmer(kmer, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, false)),
                    _ => (),
                }

                match self.search_kmer(rc, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, true)),
                    _ => (),
                }
            }

            Dir::Right => {
                match self.search_kmer(kmer, Dir::Left) {
                    Some(idx) => return Some((idx, Dir::Left, false)),
                    _ => (),
                }

                match self.search_kmer(rc, Dir::Right) {
                    Some(idx) => return Some((idx, Dir::Right, true)),
                    _ => (),
                }
            }
        }

        return None;
    }
}



#[derive(Copy, Clone)]
enum ExtModeEdge {
    Unique(usize, Dir, Exts),
    Terminal(Exts),
}

/// Try to extend edge 'edge' in direction 'dir', returning:
/// - Unique(usize, dir, exts) if there a unique extension into a another sedge
/// - Terminal(exts) if there isn't -- exts denotes the neighboring kmers
#[inline(never)]
fn try_extend_edge(sedges: &SedgeDb,
                   available_edges: &LMap<usize, ()>,
                   edge: usize,
                   dir: Dir)
                   -> ExtModeEdge {
    let &(bases, exts) = sedges.sedges.get(edge).expect("Edge doesn't exist");

    if exts.num_ext_dir(dir) != 1 {
        ExtModeEdge::Terminal(exts.single_dir(dir))
    } else {
        // Get the next kmer
        let ext_base = exts.get_unique_extension(dir).expect("should be unique");
        let end_kmer = bases.term_kmer(dir);

        let next_kmer = end_kmer.extend(ext_base, dir);
        let (next_edge, next_side_incoming, rc) = match sedges.find_link(next_kmer, dir) {
            Some(e) => e,
            None => {
                println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                println!("end kmer: {:?}", end_kmer);
                println!("No kmer: {:?}", next_kmer);
                println!("rc: {:?}", next_kmer.min_rc().0);
                panic!(format!("No kmer: {:?}", next_kmer))
            }
        };

        let &(next_bases, next_exts) = sedges.sedges.get(next_edge).expect("must have edge");

        let consistent = (next_bases.len() == K) ||
                         match (dir, next_side_incoming, rc) {
            (Dir::Left, Dir::Right, false) => true,
            (Dir::Left, Dir::Left, true) => true,
            (Dir::Right, Dir::Left, false) => true,
            (Dir::Right, Dir::Right, true) => true,
            _ => {
                println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                println!("end kmer: {:?}", end_kmer);
                println!("next kmer: {:?}", next_kmer);
                println!("rc: {:?}", next_kmer.min_rc().0);
                println!("next bases: {:?}, next_side_incoming: {:?}, rc: {:?}",
                         next_bases,
                         next_side_incoming,
                         rc);
                false
            }
        };
        assert!(consistent);

        // We can include this kmer in the line if:
        // a) it exists in the partition, and is still unused
        // b) the kmer we go to has a unique extension back in our direction

        if !available_edges.contains_key(&next_edge) {
            // This kmer isn't in this partition, or we've already used it
            return ExtModeEdge::Terminal(exts.single_dir(dir));
        }

        // orientation of next edge
        let next_side_outgoing = next_side_incoming.flip();


        let incoming_count = next_exts.num_ext_dir(next_side_incoming);
        let outgoing_exts = next_exts.single_dir(next_side_outgoing);

        if incoming_count == 0 {
            println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
            println!("end kmer: {:?}", end_kmer);
            //panic!("unreachable");
            ExtModeEdge::Terminal(exts.single_dir(dir))
        } else if incoming_count == 1 {
            // We have a unique path to next_kmer -- include it
            ExtModeEdge::Unique(next_edge, next_side_outgoing, outgoing_exts)
        } else {
            // there's more than one path
            // into the target kmer - don't include it
            ExtModeEdge::Terminal(exts.single_dir(dir))
        }
    }
}

/// Generate complete unbranched edges
#[inline(never)]
fn extend_edge(edges: &SedgeDb,
               available_edges: &mut LMap<usize, ()>,
               edge: usize,
               start_dir: Dir)
               -> (Vec<(usize, Dir)>, Exts) {
    let mut current_dir = start_dir;
    let mut current_edge = edge;
    let mut path = Vec::new();
    let final_exts: Exts; // must get set below

    available_edges.remove(&edge);

    loop {
        let ext_result = try_extend_edge(edges, available_edges, current_edge, current_dir);

        match ext_result {
            ExtModeEdge::Unique(next_edge, next_dir_outgoing, _) => {
                let next_dir_incoming = next_dir_outgoing.flip();
                path.push((next_edge, next_dir_incoming));
                available_edges.remove(&next_edge);
                current_edge = next_edge;
                current_dir = next_dir_outgoing;
            }
            ExtModeEdge::Terminal(ext) => {
                final_exts = ext;
                break;
            }
        }
    }

    (path, final_exts)
}

// Determine the sequence and extensions of the maximal unbranched
// edge, centered around the given edge number
#[inline(never)]
fn build_edge(sedge_db: &SedgeDb,
              available_edges: &mut LMap<usize, ()>,
              seed_sedge: usize)
              -> (VecDeque<u8>, Exts, VecDeque<u32>) {
    let (l_path, l_ext) = extend_edge(sedge_db, available_edges, seed_sedge, Dir::Left);
    let (r_path, r_ext) = extend_edge(sedge_db, available_edges, seed_sedge, Dir::Right);

    // println!("lpath len: {}, rpl: {}", l_path.len(), r_path.len());

    // Stick together edge chunks to get full edge sequence
    let mut line_seq = VecDeque::new();
    let mut sedge_set = VecDeque::new();

    sedge_set.push_back(seed_sedge as u32);
    let &(start_edge_seq, _) = sedge_db.get(seed_sedge);
    for i in 0..start_edge_seq.len() {
        line_seq.push_back(start_edge_seq.get(i));
    }

    // Add on the left path
    for &(next_sedge, incoming_dir) in l_path.iter() {
        sedge_set.push_front(next_sedge as u32);
        let &(seq, _) = sedge_db.get(next_sedge);

        let ext_seq = match incoming_dir {
            Dir::Left => seq.rc(),
            Dir::Right => seq,
        };

        // println!("id: {}, in_dir: {:?}, add: {:?}", next_sedge, incoming_dir, ext_seq);
        // println!("raw: {:?}", seq);

        for i in (0..(ext_seq.len() - K + 1)).rev() {

            line_seq.push_front(ext_seq.get(i))
        }
    }

    // Add on the left path
    for &(next_sedge, incoming_dir) in r_path.iter() {
        sedge_set.push_back(next_sedge as u32);
        let &(seq, _) = sedge_db.get(next_sedge);

        let ext_seq = match incoming_dir {
            Dir::Left => seq,
            Dir::Right => seq.rc(),
        };

        // println!("id: {}, in_dir: {:?}, add: {:?}", next_sedge, incoming_dir, ext_seq);
        // println!("raw: {:?}", seq);

        for i in (K - 1)..(ext_seq.len()) {
            line_seq.push_back(ext_seq.get(i))
        }
    }

    let left_extend = match l_path.last() {
        None => l_ext,
        Some(&(_, Dir::Left)) => l_ext.complement(),
        Some(&(_, Dir::Right)) => l_ext,
    };

    let right_extend = match r_path.last() {
        None => r_ext,
        Some(&(_, Dir::Left)) => r_ext,
        Some(&(_, Dir::Right)) => r_ext.complement(),
    };

    // return sequence and extensions
    (line_seq, Exts::from_single_dirs(left_extend, right_extend), sedge_set)
}


/// Link together a DeBruijn graph from the given edges of the form:
/// (start_kmer, end_kmer, extension)
#[inline(never)]
pub fn build_edges(sedges: Vec<(Lmer, Exts)>) -> TempGraph {
    let mut graph = TempGraph::new(true);

    info!("Build LinkedHashMap...");
    let mut available_sedges: LMap<usize, ()> = new_lmap();
    available_sedges.reserve(sedges.len());
    available_sedges.extend((0..sedges.len()).map(|x| (x, ())));

    info!("Making SedgeDb");
    let mut sedge_db = SedgeDb::new(sedges);

    // Important!! Prune extensions to kmers that didn't pass the kmer-filtering
    // FIXME -- add a unit test that fails if this step doesn't happen
    info!("Cleaning dangling edges");
    fix_sedge_exts(&mut sedge_db);

    info!("Sedges to process: {}", available_sedges.len());
    loop {

        if available_sedges.len() % 1000000 == 0
        {
            debug!("Building graph: sedges left: {}", available_sedges.len());
        }

        match available_sedges.front() {
            Some((&idx, _)) => {
                let (seq, exts, sedges) = build_edge(&sedge_db, &mut available_sedges, idx);
                graph.add(seq, exts, sedges);
            }
            None => break,
        }
    }

    graph
}


// Determine the sequence and extensions of the maximal unbranched
// edge, centered around the given edge number
#[inline(never)]
fn build_edge2(edge_db: &EdgeDb,
              available_edges: &mut LMap<usize, ()>,
              seed_edge: usize)
              -> (VecDeque<u8>, Exts, VecDeque<u32>) {
    let (l_path, l_ext) = extend_edge2(edge_db, available_edges, seed_edge, Dir::Left);
    let (r_path, r_ext) = extend_edge2(edge_db, available_edges, seed_edge, Dir::Right);

    // println!("lpath len: {}, rpl: {}", l_path.len(), r_path.len());

    // Stick together edge chunks to get full edge sequence
    let mut line_seq = VecDeque::new();
    let mut sedge_set = VecDeque::new();

    let start_edge = edge_db.graph.get(seed_edge);
    sedge_set.extend(start_edge.sedges.unwrap());
    let seq = start_edge.sequence;
    for i in 0..seq.len() {
        line_seq.push_back(seq.get(i));
    }

    // Add on the left path
    for &(next_edge_id, incoming_dir) in l_path.iter() {
        // FIXME
        //sedge_set.push_front(next_sedge as u32);

        let next_edge = edge_db.graph.get(next_edge_id);
        let ext_seq = next_edge.sequence;

        let ext_seq = match incoming_dir {
            Dir::Left => {
                // RC
                for i in (0..(ext_seq.len() - K + 1)).rev() {
                    line_seq.push_front(complement(ext_seq.get(ext_seq.len() - 1 - i)))
                }
            },
            Dir::Right => {
                // Normal
                for i in (0..(ext_seq.len() - K + 1)).rev() {
                    line_seq.push_front(ext_seq.get(i))
                }
            },
        };


    }

    // Add on the left path
    for &(next_edge_id, incoming_dir) in r_path.iter() {
        //sedge_set.push_back(next_sedge as u32);

        let next_edge = edge_db.graph.get(next_edge_id);
        let ext_seq = next_edge.sequence;

        let ext_seq = match incoming_dir {
            Dir::Left => {
                // Normal
                for i in (K - 1)..(ext_seq.len()) {
                    line_seq.push_back(ext_seq.get(i));
                }
            },
            Dir::Right => {
                // RC
                for i in (K - 1)..(ext_seq.len()) {
                    line_seq.push_back(complement(ext_seq.get(ext_seq.len() - 1 - i)));
                }
            },
        };
    }

    let left_extend = match l_path.last() {
        None => l_ext,
        Some(&(_, Dir::Left)) => l_ext.complement(),
        Some(&(_, Dir::Right)) => l_ext,
    };

    let right_extend = match r_path.last() {
        None => r_ext,
        Some(&(_, Dir::Left)) => r_ext,
        Some(&(_, Dir::Right)) => r_ext.complement(),
    };

    // return sequence and extensions
    (line_seq, Exts::from_single_dirs(left_extend, right_extend), sedge_set)
}


/// Generate complete unbranched edges
#[inline(never)]
fn extend_edge2(edges: &EdgeDb,
               available_edges: &mut LMap<usize, ()>,
               edge: usize,
               start_dir: Dir)
               -> (Vec<(usize, Dir)>, Exts) {
    let mut current_dir = start_dir;
    let mut current_edge = edge;
    let mut path = Vec::new();
    let final_exts: Exts; // must get set below

    available_edges.remove(&edge);

    loop {
        let ext_result = try_extend_edge2(edges, available_edges, current_edge, current_dir);

        match ext_result {
            ExtModeEdge::Unique(next_edge, next_dir_outgoing, _) => {
                let next_dir_incoming = next_dir_outgoing.flip();
                path.push((next_edge, next_dir_incoming));
                available_edges.remove(&next_edge);
                current_edge = next_edge;
                current_dir = next_dir_outgoing;
            }
            ExtModeEdge::Terminal(ext) => {
                final_exts = ext;
                break;
            }
        }
    }

    (path, final_exts)
}



/// Try to extend edge 'edge' in direction 'dir', returning:
/// - Unique(usize, dir, exts) if there a unique extension into a another sedge
/// - Terminal(exts) if there isn't -- exts denotes the neighboring kmers
#[inline(never)]
fn try_extend_edge2(edges: &EdgeDb,
                   available_edges: &LMap<usize, ()>,
                   edge_id: usize,
                   dir: Dir)
                   -> ExtModeEdge {
    let edge = edges.graph.get(edge_id);
    let exts = edge.exts;
    let bases = edge.sequence;

    if exts.num_ext_dir(dir) != 1 {
        ExtModeEdge::Terminal(exts.single_dir(dir))
    } else {
        // Get the next kmer
        let ext_base = exts.get_unique_extension(dir).expect("should be unique");
        let end_kmer = bases.term_kmer(dir);

        let next_kmer = end_kmer.extend(ext_base, dir);
        let (next_edge_id, next_side_incoming, rc) = match edges.find_link(next_kmer, dir) {
            Some(e) => e,
            None => {
                println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                println!("end kmer: {:?}", end_kmer);
                println!("No kmer: {:?}", next_kmer);
                println!("rc: {:?}", next_kmer.min_rc().0);
                panic!(format!("No kmer: {:?}", next_kmer))
            }
        };

        let next_edge = edges.graph.get(next_edge_id);
        let next_bases = next_edge.sequence;
        let next_exts = next_edge.exts;

        let consistent = (next_bases.len() == K) ||
                         match (dir, next_side_incoming, rc) {
            (Dir::Left, Dir::Right, false) => true,
            (Dir::Left, Dir::Left, true) => true,
            (Dir::Right, Dir::Left, false) => true,
            (Dir::Right, Dir::Right, true) => true,
            _ => {
                println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                println!("end kmer: {:?}", end_kmer);
                println!("next kmer: {:?}", next_kmer);
                println!("rc: {:?}", next_kmer.min_rc().0);
                println!("next bases: {:?}, next_side_incoming: {:?}, rc: {:?}",
                         next_bases,
                         next_side_incoming,
                         rc);
                false
            }
        };
        assert!(consistent);

        // We can include this kmer in the line if:
        // a) it exists in the partition, and is still unused
        // b) the kmer we go to has a unique extension back in our direction

        if !available_edges.contains_key(&next_edge_id) {
            // This kmer isn't in this partition, or we've already used it
            return ExtModeEdge::Terminal(exts.single_dir(dir));
        }

        // orientation of next edge
        let next_side_outgoing = next_side_incoming.flip();


        let incoming_count = next_exts.num_ext_dir(next_side_incoming);
        let outgoing_exts = next_exts.single_dir(next_side_outgoing);

        if incoming_count == 0 {
            println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
            println!("end kmer: {:?}", end_kmer);
            panic!("unreachable");
        } else if incoming_count == 1 {
            // We have a unique path to next_kmer -- include it
            ExtModeEdge::Unique(next_edge_id, next_side_outgoing, outgoing_exts)
        } else {
            // there's more than one path
            // into the target kmer - don't include it
            ExtModeEdge::Terminal(exts.single_dir(dir))
        }
    }
}




// Kmers are initially decorated with left and right extensions
// that were observed at least once in the datasets.
// Kmer filtering means some of the target kmers won't actually exist
// The sedges now contain all valid kmers so remove extensions to Kmers
// that don't exist
#[inline(never)]
pub fn fix_sedge_exts(sedge_db: &mut SedgeDb)
{
    let mut new_exts: Vec<(usize, Exts)> = Vec::new();

    (0 .. sedge_db.sedges.len()).into_par_iter().map(|i| {
        let (seq, exts) = sedge_db.sedges[i];
        (i, get_true_exts(seq, exts, sedge_db))
    }).collect_into(&mut new_exts);

    for &(idx, ne) in new_exts.iter() {
        sedge_db.sedges[idx].1 = ne;
    }
}


pub fn get_true_exts(seq: Lmer, exts: Exts, sedge_db: &SedgeDb) -> Exts
{
    let mut new_exts = Exts::empty();
    let l_kmer = seq.first_kmer();
    let r_kmer = seq.last_kmer();

    for i in 0..4
    {
        if exts.has_ext(Dir::Left, i) {
            if sedge_db.find_link(l_kmer.extend_left(i), Dir::Left).is_some() {
                new_exts = new_exts.set(Dir::Left, i)
            }
        }

        if exts.has_ext(Dir::Right, i) {
            if sedge_db.find_link(r_kmer.extend_right(i), Dir::Right).is_some() {
                new_exts = new_exts.set(Dir::Right, i)
            }
        }
    }

    new_exts
}




#[derive(Clone, RustcEncodable, RustcDecodable)]
pub struct TempGraph {
    /// Concatenated sequence of all edges
    pub sequence: bitenc::BitEnc,
    /// Starting position of each edge in the graph's sequence
    pub start_pos: Vec<usize>,
    pub edge_len: Vec<u32>,
    /// Edge extensions
    pub exts: Vec<Exts>,
    /// Optional: for each edge, it's "sedges". A sedge is a tuple (Lmer, extensions)
    pub edge_sedges: Option<MultiVec<u32>>,
}

impl TempGraph {
    pub fn new(track_sedges: bool) -> TempGraph {

        let mv = if track_sedges { Some(MultiVec::new()) } else { None };

        TempGraph {
            sequence: bitenc::BitEnc::new(),
            start_pos: Vec::new(),
            edge_len: Vec::new(),
            exts: Vec::new(),
            edge_sedges: mv,
        }
    }

    pub fn trim_hairy_edges(&self, edge_db: &EdgeDb, max_len: usize) -> Vec<(usize, usize, Dir, u8)> {

        let mut edges_to_trim: Vec<Option<(usize, usize, Dir, u8)>> = Vec::new();

        (0..self.len()).into_par_iter().map(|i| {
            self.find_hairy_edge(edge_db, max_len, i)
        }).collect_into(&mut edges_to_trim);

        //for i in (0..self.len())  {
        //    edges_to_trim.push(self.find_hairy_edge(edge_db, max_len, i));
        //}

        let he: Vec<(usize, usize, Dir, u8)> = edges_to_trim.into_iter().filter_map(|x| x).collect();
        he
    }

    pub fn find_hairy_edge(&self, edge_db: &EdgeDb, max_len: usize, idx: usize) -> Option<(usize, usize, Dir, u8)> {
        let exts = self.exts[idx];
        if exts.num_exts_r() > 0 && exts.num_exts_l() > 0 {
            return None;
        }

        if self.edge_len[idx] > (max_len as u32) {
            return None;
        }

        if exts.num_exts_l() == 0 && exts.num_exts_r() == 1 {
            let next = exts.get_unique_extension(Dir::Right).unwrap();
            let bound_kmer = self.get(idx).sequence.last_kmer();
            let next_kmer = bound_kmer.extend_right(next);

            edge_db.find_link(next_kmer, Dir::Right).map(|(src,dir,flip)| {
                let ext_base_to_hair = if !flip { next_kmer.get(0) } else { complement(next_kmer.get(0)) };
                (idx, src, dir, ext_base_to_hair)
            })

        } else if exts.num_exts_r() == 0 && exts.num_exts_l() == 1 {
            let next = exts.get_unique_extension(Dir::Left).unwrap();
            let bound_kmer = self.get(idx).sequence.first_kmer();
            let next_kmer = bound_kmer.extend_left(next);

            edge_db.find_link(next_kmer, Dir::Left).map(|(src,dir,flip)| {
                let ext_base_to_hair = if !flip { next_kmer.get(0) } else { complement(next_kmer.get(0)) };
                (idx, src, dir, ext_base_to_hair)
            })

        } else {
            None
        }
    }


    pub fn add(&mut self, sequence: VecDeque<u8>, exts: Exts, sedges: VecDeque<u32>) {
        let start_pos = self.sequence.len();
        self.start_pos.push(start_pos);

        let edge_len = sequence.len();
        self.edge_len.push(edge_len as u32);

        self.exts.push(exts);

        match self.edge_sedges.as_mut() {
            Some(e) => e.add(sedges.into_iter()),
            None => ()
        }

        for b in sequence {
            self.sequence.push(b);
        }

    }

    /// Iterating through the graph returns Edge objects.
    /// These are "self-contained" objects that have the Edge sequence but not much extra information.
    /// "Exploring" the graph (with explore_graph), returns VEdge objects, which have pointers to
    /// left and right edges and make, well, exploration easier.
    pub fn iter(&self) -> TempGraphIterator {
        TempGraphIterator {
            pos: 0,
            graph: self,
        }
    }

    pub fn len(&self) -> usize {
        self.start_pos.len()
    }

    /// Get the pos-th edge.
    pub fn get<'a>(&'a self, pos: usize) -> Edge<'a> {
        if pos < self.start_pos.len() {
            let seq = bitenc::BitEncSlice {
                bitenc: &self.sequence,
                start: self.start_pos[pos],
                length: self.edge_len[pos] as usize,
            };

            let sedges = self.edge_sedges.as_ref().map(|e| e.get_slice(pos));

            Edge {
                id: pos,
                sequence: seq,
                exts: self.exts[pos],
                sedges: sedges,
            }
        } else {
            panic!("oob");
        }
    }
}


// Unbranched edge in the DeBruijn graph
#[derive(Debug)]
pub struct Edge<'a> {
    pub id: usize,
    pub sequence: bitenc::BitEncSlice<'a>,
    pub exts: Exts,
    pub sedges: Option<&'a[u32]>,
}

pub struct TempGraphIterator<'a> {
    pos: usize,
    graph: &'a TempGraph,
}

impl<'a> Iterator for TempGraphIterator<'a> {
    type Item = Edge<'a>;
    fn next(&mut self) -> Option<Edge<'a>> {
        if self.pos < self.graph.start_pos.len() {
            let e = self.graph.get(self.pos);
            self.pos += 1;
            Some(e)
        } else {
            None
        }
    }
}

/// Graph edge with pointers to left and right edges.
#[derive(Debug)]
pub struct VEdge {
    pub id: u32,
    pub sequence: bitenc::BitEnc,
    /// Edges to the left of this edge. Vector of tuples
    /// (index_of_neighbor, direction of extension, reverse_complemented?)
    pub l_exts: Vec<(usize, kmer::Dir, bool)>,
    pub r_exts: Vec<(usize, kmer::Dir, bool)>,
}

impl VEdge {

    pub fn write_edge(&self, f: &mut Write) {

        let color = ",color=green";

        writeln!(f, "n{} [label=\"{} -- {}\"{},style=filled]", self.id, self.id, self.sequence.len(), color).unwrap();


        for &(id, incoming_dir, _) in self.l_exts.iter() {
            let color = "blue";
            let style =
                match incoming_dir { Dir::Right => "solid", Dir::Left => "dashed" };

            writeln!(f, "n{} -> n{} [color={},style={}]", self.id, id, color, style).unwrap();
        }

        for &(id, incoming_dir, _) in self.r_exts.iter() {
            let color = "red";
            let style =
                match incoming_dir { Dir::Left => "solid", Dir::Right => "dashed" };
            writeln!(f, "n{} -> n{} [color={},style={}]", self.id, id, color, style).unwrap();
        }
    }
}

pub fn explore_graph(graph: &TempGraph, db: &EdgeDb, start: usize, distance: usize) -> DMap<usize, VEdge> {
    // BFS to some depth around start

    let mut visited_edges: DMap<usize, VEdge> = new_dmap();
    let mut to_visit = vec![start];

    for i in 0 .. distance {
        let mut next_to_visit = Vec::new();
        for edge_id in to_visit.iter() {

            if visited_edges.contains_key(&edge_id) {
                continue
            }

            let edge = graph.get(*edge_id);
            let exts = edge.exts;
            let l_kmer = edge.sequence.first_kmer();
            let r_kmer = edge.sequence.last_kmer();
            let mut l_exts = Vec::new();
            let mut r_exts = Vec::new();

            for i in 0..4
            {
                if exts.has_ext(Dir::Left, i) {
                    let link = db.find_link(l_kmer.extend_left(i), Dir::Left).expect("missing link");
                    l_exts.push(link);
                    if !visited_edges.contains_key(&link.0) {
                        next_to_visit.push(link.0)
                    }
                }

                if exts.has_ext(Dir::Right, i) {
                    let link = db.find_link(r_kmer.extend_right(i), Dir::Right).expect("missing link");
                    r_exts.push(link);
                    if !visited_edges.contains_key(&link.0) {
                        next_to_visit.push(link.0)
                    }
                }
            }

            let ve =
                VEdge {
                    id: *edge_id as u32,
                    sequence: edge.sequence.to_bitenc(),
                    l_exts: l_exts,
                    r_exts: r_exts,
                };
            visited_edges.insert(*edge_id, ve);
        }

        to_visit = next_to_visit;
    }

    visited_edges
}



pub fn do_explore(graph: String, edge_db_file: Option<String>) {
    // Load graph
    info!("Loading graph...");
    let g = utils::read_obj::<TempGraph>(Path::new(&graph)).unwrap();
    info!("Loaded.");
    info!("Seq length: {}", g.sequence.len());
    info!("NEdges: {}", g.start_pos.len());


    let edge_db = match edge_db_file {
        Some(infile) => {
            info!("Loading EdgeDb");
            let (left, right) = utils::read_obj::<(Vec<u32>, Vec<u32>)>(Path::new(&infile)).unwrap();
            EdgeDb {
                left_order: left,
                right_order: right,
                graph: &g
            }
        },
        None => {
            info!("Making EdgeDb");
            let edge_db = EdgeDb::new(&g);
            utils::write_obj(&(edge_db.left_order.clone(), edge_db.right_order.clone()), Path::new("edge_db.bin"));
            edge_db
        },
    };

    let he = g.trim_hairy_edges(&edge_db, 80);
    for h in he.iter() {
        println!("{:?}", h);
    }


    // Find some long-ish edges and explore the region
    for (idx, _, _, _) in he.into_iter().take(20){

        let mut ex = explore_graph(&g, &edge_db, idx, 4);

        let mut edges = Vec::new();
        for (k,v) in ex.drain() {
            edges.push(v);
        }

        let f = format!("he{}.dot", idx);
        to_dot(&edges, PathBuf::from(f));

    }


    let mut n = 0;

    // Find some long-ish edges and explore the region
    for (idx, e) in g.iter().enumerate() {
        if n > 20 {
            break;
        }

        if e.sequence.len() > 100 {
            n += 1;
            let mut ex = explore_graph(&g, &edge_db, idx, 4);

            let mut edges = Vec::new();
            for (k,v) in ex.drain() {
                edges.push(v);
            }

            let f = format!("g{}.dot", n);
            to_dot(&edges, PathBuf::from(f));
        }
    }

    /*
    let sa_val = 16;
    info!("Computing BWT...");
    let bwt1 = bwt::compute_bwt(&g.sequence, 3, sa_val);

    info!("Computing FMIndex");
    let fm = bwt::FMIndex::new(bwt1.bwt, bwt1.sa, sa_val, 3);

    let pos = fm.backward_search(seq.into_bytes().iter());
    info!("pos: {:?}", pos);
    */
}

/// Write the POA to a dot file.
pub fn to_dot(edges: &Vec<VEdge>, path: PathBuf) {

    let mut f = File::create(path).expect("couldn't open file");

    writeln!(&mut f, "digraph {{").unwrap();
    for e in edges.iter() {
        e.write_edge(&mut f);
    }
    writeln!(&mut f, "}}").unwrap();
}



#[cfg(test)]
mod tests {
    use kmer::*;
    use super::*;
    use sim_tests::random_dna;
    use utils::new_dmap;

    #[test]
    fn kmerize() {
        // Take a random input string, kmerize it, and put it back together
        let seq = random_dna(10000);

        let p = 6;
        let permutation = (0..(1 << 2 * p)).collect();

        let shard_bsps = Bsp::msp_read(K, p, 0, 0, &seq[..], &permutation, true);
        let mut all_kmers = new_dmap();
        let mut bsps = Vec::new();

        for (_, bsp) in shard_bsps {
            for km in bsp.kmers() {
                let (min_kmer, _) = km.min_rc();
                all_kmers.insert(min_kmer, true);
            }

            bsps.push(bsp);
        }

        let _main = kmer_extensions(&all_kmers, &all_kmers, &bsps, true);

        // pub fn kmer_extensions<T>(valid_kmers: HashMap<Kmer, T>, bsps: &Vec<Bsp>) -> HashMap<Kmer, Exts>
        // fn build_edge_from_kmer(
        //    kmer_exts: &HashMap<Kmer, Exts>,
        //    available_kmers: &mut HashSet<Kmer>,
        //    kmer: Kmer) -> (Lmer, Exts)

    }
}
