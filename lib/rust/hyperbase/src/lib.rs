// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.

// This file defines a HyperBasevector structure, and a Hyper structure, which
// is a HyperBasevector plus an involution and read ids for each edge.

use debruijn::compression::{compress_kmers, SimpleCompress};
use debruijn::dna_string::DnaString;
use debruijn::graph::DebruijnGraph;
use debruijn::kmer::Kmer20;
use debruijn::{filter, kmer, Exts, Kmer, Mer, Vmer};
use equiv::EquivRel;
use graph_simple::GraphSimple;
use kmer_lookup::make_kmer_lookup_20_single;
use petgraph::prelude::*;
use std::cmp::max;
use vector_utils::{
    bin_position, bin_position1_3, next_diff, next_diff1_3, reverse_sort, unique_sort,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// CONVERT FROM DEBRUIJN GRAPH TO PET GRAPH
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn debruijn_to_petgraph_hyperbasevector<K: Kmer>(
    g_in: &DebruijnGraph<K, Vec<u32>>,
    g_out: &mut Graph<u32, DnaString, Directed, u32>,
    inv: &mut Vec<u32>,
) {
    // Find the edges in the transformed graph, and build inv.

    let num_nodes = g_in.len();
    inv.clear();
    inv.reserve(2 * num_nodes);
    let mut edges: Vec<(i32, bool, DnaString)> = Vec::with_capacity(2 * num_nodes);
    let mut edge_to_fw: Vec<i32> = Vec::with_capacity(num_nodes);
    let mut edge_to_rc: Vec<i32> = Vec::with_capacity(num_nodes);

    // Dump the deBruijn graph.

    /*
    for i in 0 .. num_nodes {
        let node = g_in.get_node(i);
        let seq = node.sequence();
        println!( "\n{} = {}", i, seq.to_string() );
        for j in 0..node.r_edges().len() {
            let rnode = node.r_edges()[j];
            let e2 = rnode.0;
            println!( "==> {}, {:?}, {}", e2, rnode.1, rnode.2 );
        }
    }
    */

    // Proceed.

    let mut pal = Vec::<bool>::new();
    for i in 0..num_nodes {
        let node = g_in.get_node(i);
        let seq = node.sequence();

        // Update involution.

        let palindrome = seq == seq.rc();
        if !palindrome {
            inv.push((edges.len() + 1) as u32);
        }
        inv.push(edges.len() as u32);

        // Save edge.

        edge_to_fw.push(edges.len() as i32);
        if palindrome {
            edge_to_rc.push(edges.len() as i32);
        }
        edges.push((i as i32, true, seq.to_owned()));
        pal.push(palindrome);
        if !palindrome {
            edge_to_rc.push(edges.len() as i32);
            edges.push((i as i32, false, seq.rc().to_owned()));
            pal.push(palindrome);
        }
    }

    // Find pairs of adjacent edges in the transformed graph.

    let mut adj = Vec::<(i32, i32)>::new();
    for x in 0..edges.len() {
        if edges[x].1 {
            let e1 = edges[x].0;
            let node = g_in.get_node(e1 as usize);
            for j in 0..node.r_edges().len() {
                let rnode = node.r_edges()[j];
                let e2 = rnode.0;
                if !rnode.2 {
                    adj.push((x as i32, edge_to_fw[e2]));
                } else {
                    adj.push((x as i32, edge_to_rc[e2]));
                }
            }
        }
        if !edges[x].1 || pal[x] {
            let e1 = edges[x].0;
            let node = g_in.get_node(e1 as usize);
            for j in 0..node.l_edges().len() {
                let lnode = node.l_edges()[j];
                let e2 = lnode.0;
                if !lnode.2 {
                    adj.push((x as i32, edge_to_rc[e2]));
                } else {
                    adj.push((x as i32, edge_to_fw[e2]));
                }
            }
        }
    }

    // Find nodes in the transformed graph.  They are sets of edge ends under
    // the natural equivalence relation.

    let mut eq: EquivRel = EquivRel::new(2 * edges.len() as u32);
    for (left, right) in adj {
        eq.join((2 * left + 1) as usize, (2 * right) as usize);
    }
    let reps = eq.set_reps();

    // Now actually create the transformed graph.

    g_out.clear();
    g_out.reserve_exact_nodes(reps.len());
    g_out.reserve_exact_edges(edges.len());
    for i in 0..reps.len() {
        g_out.add_node(i as u32);
    }
    for (e, edge) in edges.into_iter().enumerate() {
        let v = bin_position(&reps, &eq.set_id(2 * e));
        let w = bin_position(&reps, &eq.set_id(2 * e + 1));
        g_out.add_edge(
            NodeIndex::<u32>::new(v as usize),
            NodeIndex::<u32>::new(w as usize),
            edge.2,
        );
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// HYPERBASEVECTOR DEFINITION
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// A HyperBasevector is a constant k and a digraph whose edges represent DNA
// sequences of length >= k, such that if two edges e1 and e2 abut at a vertex v,
// like so, u --e1--> v --e2--> w, then the last k-1 bases of e1 equal the first
// k-1 bases of e2.  We implement a HyperBasevector here using the graph structure
// in the petgraph crate.

pub struct HyperBasevector {
    pub k: i32,
    pub g: Graph<u32, DnaString, Directed, u32>,
}

impl HyperBasevector {
    pub fn new() -> HyperBasevector {
        HyperBasevector {
            k: 0,
            g: Graph::new(),
        }
    }

    // Return the number of kmers in an edge.

    pub fn kmers(&self, e: u32) -> usize {
        self.g[EdgeIndex::<u32>::new(e as usize)].len() - self.k as usize + 1
    }

    // Return the number of bases in an edge.

    pub fn bases(&self, e: u32) -> usize {
        self.g[EdgeIndex::<u32>::new(e as usize)].len()
    }

    // =============================================================================
    // Test for uniqueness of kmers.  This tests to see if a kmer appears at most
    // once in the graph.  This is not a requirement for a HyperBasevector.
    // =============================================================================

    pub fn test_unique(&self) {
        let mut edges = Vec::<DnaString>::new();
        for e in 0..self.g.edge_count() {
            edges.push(self.g.edge_obj(e as u32).clone());
        }
        let mut kmers_plus = Vec::<(Kmer20, i32, i32)>::new();
        make_kmer_lookup_20_single(&edges, &mut kmers_plus);
        let mut i: i64 = 0;
        let mut dups = 0;
        while i < kmers_plus.len() as i64 {
            let j = next_diff1_3(&kmers_plus, i as usize) as i64;
            if j - i > 1 {
                if dups == 0 {
                    println!("\ntest_unique failed");
                    println!(
                        "the kmer {} appears in edges {} and {}",
                        kmers_plus[i as usize].0.to_string(),
                        kmers_plus[i as usize].1,
                        kmers_plus[i as usize + 1].1
                    );
                }
                dups += j - i - 1;
            }
            i = j;
        }
        if dups > 0 {
            println!("of {} kmers, {} are duplicated", kmers_plus.len(), dups);
            panic!("bailing because test failed");
        }
    }
}

impl Default for HyperBasevector {
    fn default() -> Self {
        Self::new()
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// HYPER DEFINITION
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// A "Hyper" is a HyperBasevector, together with an involution "inv", mapping edges
// of the graph to reverse complement edges, and an assignment of forward oriented
// reads to each edge, represented by their "ids".  One can also keep track of
// actual paths for each read, but we do not do that here because tracking ids is
// sufficient for our purposes.

pub struct Hyper {
    pub h: HyperBasevector,
    pub inv: Vec<u32>,
    pub ids: Vec<Vec<u32>>,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// MEMBER FUNCTIONS FOR HYPER
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

impl Hyper {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // FUNCTIONS THAT DEPEND ONLY ON THE GRAPH
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // =============================================================================
    // Return the sequence associated to a path.
    // =============================================================================

    pub fn cat(&self, p: &[i32]) -> DnaString {
        let mut s = DnaString::new();
        if p.is_empty() {
            return s;
        }
        s = self.h.g.edge_obj(p[0] as u32).clone();
        for j in 1..p.len() as i32 {
            let b = &self.h.g.edge_obj(p[j as usize] as u32);
            for i in self.h.k - 1..b.len() as i32 {
                s.push(b.get(i as usize));
            }
        }
        s
    }

    // =============================================================================
    // Print the graph, component by component.
    // =============================================================================

    pub fn print(&self) {
        let mut comp = Vec::<Vec<u32>>::new();
        self.h.g.components_e(&mut comp);
        for (j, comp_j) in comp.into_iter().enumerate() {
            println!("\nCOMPONENT {}", j + 1);
            for e in comp_j {
                let v = self.h.g.to_left(e);
                let w = self.h.g.to_right(e);
                let b: DnaString = self.h.g[EdgeIndex::<u32>::new(e as usize)].clone();
                println!(
                    "\n{} ==(e={},len={},supp={})==> {}",
                    v,
                    e,
                    b.len() - self.h.k as usize + 1,
                    self.supp(e as usize),
                    w
                );
                println!("{}", b.to_string());
            }
        }
    }

    // ============================================================================================
    // Print the graph, component by component, with edge annotations.
    //
    // require_ann: if true, only show components having annotations.
    // hide_seq: if true, and there is an annotation on an edge, don't print the sequence
    // ============================================================================================

    pub fn print_with_annotations(&self, ann: &[String], require_ann: bool, hide_seq: bool) {
        let mut comp = Vec::<Vec<u32>>::new();
        self.h.g.components_e(&mut comp);
        let mut n = 0;
        for comp_j in comp {
            let mut have_ann = false;
            for e in &comp_j {
                if !ann[*e as usize].is_empty() {
                    have_ann = true;
                }
            }
            if require_ann && !have_ann {
                continue;
            }
            n += 1;
            println!("\nCOMPONENT {n}");
            for e in comp_j {
                let e = e as usize;
                let v = self.h.g.to_left(e as u32);
                let w = self.h.g.to_right(e as u32);
                let b: DnaString = self.h.g[EdgeIndex::<u32>::new(e)].clone();
                println!(
                    "\n{} ==(e={},len={},supp={})==> {}\n",
                    v,
                    e,
                    b.len() - self.h.k as usize + 1,
                    self.supp(e),
                    w
                );
                if !ann[e].is_empty() {
                    print!("{}", ann[e]);
                    if !hide_seq {
                        println!();
                    }
                }
                if !hide_seq || ann[e].is_empty() {
                    println!("{}", b.to_string());
                }
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // UTILITY FUNCTIONS
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // =============================================================================
    // Create a new Hyper.
    // =============================================================================

    pub fn new() -> Hyper {
        Hyper {
            h: HyperBasevector::new(),
            inv: Vec::new(),
            ids: Vec::new(),
        }
    }

    // =============================================================================
    // Create a new Hyper from data.
    // =============================================================================

    pub fn build_from_reads(&mut self, k: i32, reads: &[DnaString]) {
        // Only works if k = 20.

        assert_eq!(k, 20);

        // Build deBruijn graph.

        pub type Kmer1 = kmer::Kmer20;
        let mut seqs = Vec::new();
        for r in reads {
            seqs.push((r.clone(), Exts::empty(), 0));
        }
        let summarizer = filter::CountFilterSet::new(1);
        let mut valid_kmers: Vec<(Kmer1, (Exts, Vec<u32>))> = {
            let (kmer_hashmap, _) = filter::filter_kmers::<Kmer1, _, _, _, _>(
                &seqs,
                &Box::new(summarizer),
                false,
                false,
                4,
            );
            drop(seqs);
            kmer_hashmap
                .iter()
                .map(|(k, e, d)| (*k, (*e, d.clone())))
                .collect()
        };
        // ◼ There are several instances of drop here and below.  Why is it that
        // ◼ rust doesn't drop objects at last usage?
        valid_kmers.sort();
        let cmp = SimpleCompress::new(|mut a: Vec<u32>, b: &Vec<u32>| {
            a.extend(b);
            a.sort_unstable();
            a.dedup();
            a
        });
        let graph = compress_kmers(false, &cmp, &valid_kmers).finish_serial();
        drop(valid_kmers);

        // Translate to graph in which edges are sequences.

        debruijn_to_petgraph_hyperbasevector(&graph, &mut self.h.g, &mut self.inv);
        drop(graph);
        self.h.k = k;

        // Build self.ids.

        let evec = Vec::<u32>::new();
        for _e in 0..self.h.g.edge_count() {
            self.ids.push(evec.clone());
        }
        let mut kmers_plus = Vec::<(Kmer20, i32, i32)>::new();
        let mut edges = Vec::<DnaString>::new();
        for e in 0..self.h.g.edge_count() {
            edges.push(self.h.g.edge_obj(e as u32).clone());
        }
        make_kmer_lookup_20_single(&edges, &mut kmers_plus);
        drop(edges);
        let mut maxread = 0;
        for read in reads {
            maxread = max(maxread, read.len());
        }
        if maxread < k as usize {
            return;
        }
        let mut next_rpos: Vec<i32> = vec![0; reads.len()];
        for pos in 0..maxread - (k as usize) + 1 {
            for (id, read) in reads.iter().enumerate() {
                if pos + k as usize > read.len() {
                    continue;
                }
                if pos < next_rpos[id] as usize {
                    continue;
                }
                let x: Kmer20 = read.get_kmer(pos);
                let p = bin_position1_3(&kmers_plus, &x);
                if p < 0 {
                    continue;
                }
                let mut q: Vec<i32> = vec![kmers_plus[p as usize].1];
                let mut rpos = pos + k as usize;
                let mut e = kmers_plus[p as usize].1;
                let mut epos = kmers_plus[p as usize].2 + k;
                self.ids[e as usize].push(id as u32);
                loop {
                    if rpos == read.len() {
                        break;
                    }
                    let mut next = false;
                    if epos == self.h.bases(e as u32) as i32 {
                        let v = self.h.g.to_right(e as u32);
                        for j in 0..self.h.g.n_from(v as usize) {
                            let f = self.h.g.e_from(v as usize, j);
                            if self.h.g.edge_obj(f as u32).get((k - 1) as usize) == read.get(rpos) {
                                e = f as i32;
                                self.ids[e as usize].push(id as u32);
                                epos = k - 1;
                                q.push(e);
                                next = true;
                                break;
                            }
                        }
                        if !next {
                            break;
                        }
                    }
                    if !next {
                        if read.get(rpos) != self.h.g.edge_obj(e as u32).get(epos as usize) {
                            break;
                        }
                        rpos += 1;
                        epos += 1;
                    }
                }
                next_rpos[id] = (rpos - k as usize + 1) as i32;
            }
        }
        for e in 0..self.h.g.edge_count() {
            unique_sort(&mut self.ids[e]);
        }
        // self.h.test_unique();   // turn on if you want sanity check
        // self.test_involution(); // turn on if you want sanity check
        // self.test_overlaps();   // turn on if you want sanity check
    }

    // =============================================================================
    // Return the read support for an edge.
    // =============================================================================

    pub fn supp(&self, e: usize) -> usize {
        self.ids[e].len() + self.ids[self.inv[e] as usize].len()
    }

    // =============================================================================
    // Kill a list 'dels' of edges in a Hyper.  This will kill both the edges and
    // involuted edges, and updates inv and ids appropriately.  The input list can
    // have duplicates and does not have to be in order or contain the involuted
    // edges.
    //
    // Normally you would want to call kill_edges instead.
    // =============================================================================

    pub fn kill_edges_raw(&mut self, dels: &[u32]) {
        // Symmetrize and unique sort dels.

        let mut dels2 = dels.to_vec();
        for e in dels {
            dels2.push(self.inv[*e as usize]);
        }
        unique_sort(&mut dels2);

        // Delete edges in dels, updating inv and ids accordingly.

        for j in (0..dels2.len() as i32).rev() {
            let e = dels2[j as usize] as usize;
            if e != self.h.g.edge_count() - 1 {
                if self.h.g.edge_count() - 1 != self.inv[self.h.g.edge_count() - 1] as usize {
                    let last_index_inv = self.inv[self.h.g.edge_count() - 1];
                    self.inv[e] = *self.inv.last().unwrap();
                    self.inv[last_index_inv as usize] = e as u32;
                } else {
                    self.inv[e] = e as u32;
                }
            }
            self.inv.pop();
            if e != self.h.g.edge_count() - 1 {
                let n = self.ids.len() - 1;
                self.ids.swap(e, n);
            }
            self.ids.pop();
            self.h.g.remove_edge(EdgeIndex::<u32>::new(e));
        }
    }

    // =============================================================================
    // Kill a list 'dels' of edges in a Hyper, and remove unneeded vertices.  This
    // will kill both the edges and involuted edges, and updates inv and ids
    // appropriately.  The input list can have duplicates and does not have to be
    // in order or contain the involuted edges.
    //
    // The code is essentially from RemoveUnneededVertices2 in
    // paths/long/large/GapToyTools.cc, in the supernova codebase, and that drives
    // from the same function in the same place in the BroadCRD codebase.  There is
    // documentation in those files that describes some of the logic.
    // =============================================================================

    pub fn kill_edges(&mut self, dels: &[u32]) {
        // Kill the edges.

        self.kill_edges_raw(dels);

        // Find max read id.

        let mut maxread: i32 = -1;
        for e in 0..self.h.g.edge_count() {
            for id in &self.ids[e] {
                maxread = max(maxread, *id as i32);
            }
        }

        // Find pairs of edges v --> x --> w, where there is exactly one edge
        // entering and exiting a vertex x, excluding the case v = w.

        let mut vertex_kill = vec![false; self.h.g.node_count()];
        let mut vertex_queue = Vec::<NodeIndex<u32>>::new();
        for v in self.h.g.node_indices() {
            if self.h.g.edges_directed(v, Outgoing).count() == 1
                && self.h.g.edges_directed(v, Incoming).count() == 1
                && self
                    .h
                    .g
                    .edges_directed(v, Outgoing)
                    .next()
                    .unwrap()
                    .target()
                    != self
                        .h
                        .g
                        .edges_directed(v, Incoming)
                        .next()
                        .unwrap()
                        .source()
            {
                vertex_kill[v.index()] = true;
                vertex_queue.push(v);
            }
        }

        // Merge these edge pairs, with some exceptions, because of the involution.
        // This is complicated.

        let mut bound = Vec::<(u32, u32)>::new();
        while let Some(v) = vertex_queue.pop() {
            if !vertex_kill[v.index()] {
                continue;
            }
            let mut vleft = v;
            let mut eleft: usize;
            loop {
                vertex_kill[vleft.index()] = false;
                eleft = self.h.g.first_edge(vleft, Incoming).unwrap().index();
                vleft = self
                    .h
                    .g
                    .edges_directed(vleft, Incoming)
                    .next()
                    .unwrap()
                    .source();
                if !vertex_kill[vleft.index()] {
                    break;
                }
            }
            let mut eright: usize;
            let mut vright = v;
            loop {
                vertex_kill[vright.index()] = false;
                eright = self.h.g.first_edge(vright, Outgoing).unwrap().index();
                vright = self
                    .h
                    .g
                    .edges_directed(vright, Outgoing)
                    .next()
                    .unwrap()
                    .target();
                if !vertex_kill[vright.index()] {
                    break;
                }
            }
            if eleft < self.inv[eright] as usize {
                bound.push((eleft as u32, eright as u32));
                bound.push((self.inv[eright], self.inv[eleft]));
            }
        }
        let mut new_edge_numbers = Vec::<usize>::new();
        let mut to_delete = Vec::<u32>::new();
        let mut have: Vec<u32>;
        let mut havex: Vec<bool> = vec![false; (maxread + 1) as usize];
        self.ids.reserve(bound.len());
        while let Some(bounds) = bound.pop() {
            let new_edge_no: usize = self.h.g.edge_count();
            let mut new_edge = self.h.g.edge_obj(bounds.0).clone();
            have = self.ids[bounds.0 as usize].clone();
            for j in 0..have.len() {
                havex[have[j] as usize] = true;
            }
            to_delete.push(bounds.0);
            let mut v = self
                .h
                .g
                .edge_endpoints(EdgeIndex::<u32>::new(bounds.0 as usize))
                .unwrap()
                .1;
            loop {
                if v == self
                    .h
                    .g
                    .edge_endpoints(EdgeIndex::<u32>::new(bounds.1 as usize))
                    .unwrap()
                    .1
                {
                    break;
                }
                let edge_id = self.h.g.first_edge(v, Outgoing).unwrap();
                let edge = self.h.g[edge_id].clone();
                for j in 0..self.ids[edge_id.index()].len() {
                    let id = self.ids[edge_id.index()][j];
                    if !havex[id as usize] {
                        have.push(id);
                        havex[id as usize] = true;
                    }
                }
                to_delete.push(edge_id.index() as u32);
                for p in (self.h.k - 1)..(edge.len() as i32) {
                    new_edge.push(edge.get(p as usize));
                }
                v = self
                    .h
                    .g
                    .edges_directed(v, Outgoing)
                    .next()
                    .unwrap()
                    .target();
            }
            let v = self
                .h
                .g
                .edge_endpoints(EdgeIndex::<u32>::new(bounds.0 as usize))
                .unwrap()
                .0;
            let w = self
                .h
                .g
                .edge_endpoints(EdgeIndex::<u32>::new(bounds.1 as usize))
                .unwrap()
                .1;
            self.h.g.add_edge(v, w, new_edge);
            have.sort_unstable();
            for j in 0..have.len() {
                havex[have[j] as usize] = false;
            }
            self.ids.push(have.clone());
            new_edge_numbers.push(new_edge_no);
        }
        for i in 0..new_edge_numbers.len() {
            if i % 2 == 1 {
                continue;
            }
            self.inv.push(new_edge_numbers[i + 1] as u32);
            self.inv.push(new_edge_numbers[i] as u32);
        }

        // Now again delete edges, as at the beginning.

        self.kill_edges_raw(&to_delete);

        // Remove edgeless vertices.

        for v in (0..self.h.g.node_count()).rev() {
            if self.h.g.n_to(v) == 0 && self.h.g.n_from(v) == 0 {
                self.h.g.remove_node(NodeIndex::<u32>::new(v));
            }
        }
    }

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
    // FUNCTIONS INTENDED MOSTLY FOR DEBUGGING
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // =============================================================================
    // Test the involution to see if it is valid.
    // =============================================================================

    #[allow(dead_code)]
    fn test_involution(&self) {
        assert_eq!(self.h.g.edge_count(), self.inv.len());
        for e in 0..self.h.g.edge_count() {
            if self.inv[e] >= self.h.g.edge_count() as u32 {
                println!(
                    "inv[{}] = {}, but graph has only {} edges.",
                    e,
                    self.inv[e],
                    self.h.g.edge_count()
                );
                panic!("Involution test failed.");
            }
            assert_eq!(self.inv[self.inv[e] as usize], e as u32);
            let t = self.h.g[EdgeIndex::<u32>::new(e)].rc().to_string();
            assert_eq!(
                self.h.g[EdgeIndex::<u32>::new(self.inv[e] as usize)].to_string(),
                t
            );
        }
        let mut homomorphism_fails = Vec::<(usize, usize)>::new();
        let mut oks = 0;
        for v in 0..self.h.g.node_count() {
            for j1 in 0..self.h.g.n_to(v) {
                let e1 = self.h.g.e_to(v, j1);
                let re1 = self.inv[e1];
                for j2 in 0..self.h.g.n_from(v) {
                    let e2 = self.h.g.e_from(v, j2);
                    let re2 = self.inv[e2];
                    if self.h.g.to_right(re2) != self.h.g.to_left(re1) {
                        homomorphism_fails.push((e1, e2));
                    } else {
                        oks += 1;
                    }
                }
            }
        }
        if !homomorphism_fails.is_empty() {
            println!("\ntest_involution failed at homomorphism condition");
            println!(
                "there were {} fails and {} oks",
                homomorphism_fails.len(),
                oks
            );
            let (e1, e2) = (homomorphism_fails[0].0, homomorphism_fails[0].1);
            let (re1, re2) = (self.inv[e1], self.inv[e2]);
            println!("first has e1 = {e1}, e2 = {e2}, re1 = {re1}, re2 = {re2}");
            panic!("bailing because of test_involution failure");
        }
    }

    // =============================================================================
    // Test to see if overlap condition is satisfied.
    // =============================================================================

    #[allow(dead_code)]
    fn test_overlaps(&self) {
        for v in 0..self.h.g.node_count() {
            for j1 in 0..self.h.g.n_to(v) {
                let b1 = self.h.g.edge_obj(self.h.g.e_to(v, j1) as u32);
                let n1 = b1.len();
                let k = self.h.k as usize;
                for j2 in 0..self.h.g.n_from(v) {
                    let b2 = self.h.g.edge_obj(self.h.g.e_from(v, j2) as u32);
                    assert_eq!(b1.slice(n1 - (k - 1), n1), b2.slice(0, k - 1));
                }
            }
        }
    }

    // =============================================================================
    // Compute a checksum.  This will produce the same answers for two isomorphic
    // graphs.  Optionally, this captures the info in the supporting reads.  This
    // is a totally arbitrary checksum but it is effective in finding changes and
    // can be easily implemented in another language.
    // =============================================================================

    pub fn checksum_main(&self, use_reads: bool) -> u64 {
        let mut s: u64 = 0;
        for v in 0..self.h.g.node_count() {
            let n = self.h.g.n_to(v);
            let mut r: u64 = 1;
            for j in 0..n {
                let e = self.h.g.e_to(v, j);
                let b = &self.h.g[EdgeIndex::<u32>::new(e)].to_string();
                let mut k: u64 = 0;
                for l in 0..b.len() {
                    let c = &b.chars().nth(l).unwrap();
                    let mut i = 0;
                    if *c == 'C' {
                        i = 1;
                    } else if *c == 'G' {
                        i = 2;
                    } else if *c == 'T' {
                        i = 3;
                    }
                    k = k.wrapping_mul(3) + i;
                }
                if use_reads {
                    for i in 0..self.ids[e].len() {
                        let id = self.ids[e][i];
                        k = k.wrapping_add(((i + 1) * id as usize * id as usize) as u64);
                    }
                }
                r = r.wrapping_mul(k);
            }
            s = s.wrapping_add(r);
            let n = self.h.g.n_from(v);
            let mut r: u64 = 1;
            for j in 0..n {
                let e = self.h.g.e_from(v, j);
                let b = &self.h.g[EdgeIndex::<u32>::new(e)].to_string();
                let mut k: u64 = 0;
                for l in 0..b.len() {
                    let c = &b.chars().nth(l).unwrap();
                    let mut i = 0;
                    if *c == 'C' {
                        i = 1;
                    } else if *c == 'G' {
                        i = 2;
                    } else if *c == 'T' {
                        i = 3;
                    }
                    k = k.wrapping_mul(5) + i;
                }
                if use_reads {
                    for i in 0..self.ids[e].len() {
                        let id = self.ids[e][i];
                        k = k.wrapping_add(((i + 1) * id as usize * id as usize) as u64);
                    }
                }
                r = r.wrapping_mul(k);
            }
            s = s.wrapping_add(r);
        }
        s
    }

    pub fn checksum(&self) -> u64 {
        self.checksum_main(true)
    }

    pub fn checksum_hbv_only(&self) -> u64 {
        self.checksum_main(false)
    }

    // =============================================================================
    // Print component sizes, in an abbreviated form, where the size of each
    // component is the number of edges in it.
    // =============================================================================

    #[allow(dead_code)]
    fn print_comp_sizes(&mut self) {
        let mut comp = Vec::<Vec<u32>>::new();
        self.h.g.components_e(&mut comp);
        let mut sizes = Vec::<usize>::new();
        for c in comp {
            sizes.push(c.len());
        }
        reverse_sort(&mut sizes);
        print!("component sizes = [");
        let mut j = 0;
        loop {
            if j == sizes.len() {
                break;
            }
            let k = next_diff(&sizes, j);
            if j > 0 {
                print!(", ");
            }
            print!("{}^{}", sizes[j], k - j);
            j = k;
        }
        println!("]");
    }

    // =============================================================================
    // Return total read support, intended as a sort of checksum.
    // =============================================================================

    #[allow(dead_code)]
    fn total_supp(&mut self) -> usize {
        let mut total = 0;
        for e in 0..self.h.g.edge_count() {
            total += self.ids[e].len();
        }
        total
    }

    // =============================================================================
    // Return number of zero-support edges.
    // =============================================================================

    #[allow(dead_code)]
    fn zero_supp_edges(&mut self) -> usize {
        let mut z = 0;
        for e in 0..self.h.g.edge_count() {
            if self.supp(e) == 0 {
                z += 1;
            }
        }
        z
    }
}

impl Default for Hyper {
    fn default() -> Self {
        Self::new()
    }
}
