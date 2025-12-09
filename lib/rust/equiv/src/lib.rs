// Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

//! This file contains an equivalence relation struct.  There is a literature on
//! this and multiple rust implementations of sophisticated algorithms, see:
//!
//! 1. https://en.wikipedia.org/wiki/Disjoint-set_data_structure
//! 2. https://crates.io/crates/union-find
//! 3. https://crates.io/crates/disjoint-sets
//! 4. https://crates.io/crates/disjoint-set [seems defunct]
//! 5. https://crates.io/crates/fera-unionfind
//!
//! The code here is an optimized and rustified version of the code in the 10X
//! supernova codebase, which was adopted from the BroadCRD codebase.  The code here
//! uses a very naive algorithm that should not be competitive with the sophisticated
//! algorithms, but for unknown reasons, it is.  There are some comparisons to the
//! disjoint-sets crate at the end of this file.  The implementations in other
//! crates were not tested.
//!
//! Computational performance of EquivRel:
//! - storage = 3N bytes, where N is the set size; storage is flat
//! - initialization time = O(N)
//! - time to make n joins = O( n * log(N) )
//! - time to find all set reps = O(N)
//! - time to find a set = O(size of set)
//! - time to find the size of a set = O(1)
//! - time to find the set id of an element = O(1).

/// The type used to represent an element in an EquivRel disjoint set.
/// 32 bits are plenty to enumerate the kinds of entities we deal with in
/// experimental data. If in the future we need more that 2^32 elements in
/// certain cases, we should parameterize the EquivRel type so that either
/// width can be used.
///
/// We generally paper over this fact in the API by accepting and returning
/// usize, since these values are essentially always used to index into
/// collections. Using u32 internally keeps the memeory footprint down, though
/// this may not be of practical importance.
type ElementId = u32;

/// The type used to enumerate the individual disjoint sets.
type SetId = u32;

/// Represent the partition of N items into M distinct sets where M <= N.
///
/// The elements of the overall collection are represented as numeric indices
/// spanning 0..N; the indices can be used to index into an external fixed
/// collection that actually contains the data for each element.
pub struct EquivRel {
    /// The entry v at position i in this vector is the identity of the next
    /// element in the set that currently contains the element i.
    next_element: Vec<ElementId>,
    /// For each element in the set i, the corresponding value in this vector
    /// is the index of the set that the element is currently in.
    set_id: Vec<SetId>,
    /// For each set s, the corresponding value in this vector is the current
    /// number of elements contained in that set.
    size: Vec<u32>,
}

impl EquivRel {
    /// Create disjoint sets of n elements.
    /// The elements are initially contained in n sets of one element each.
    pub fn new(n: u32) -> EquivRel {
        EquivRel {
            next_element: (0..n).collect(),
            set_id: (0..n as SetId).collect(),
            size: vec![1; n as usize],
        }
    }

    /// Join the sets which currently contain the two provided elements.
    pub fn join(&mut self, a: usize, b: usize) {
        // Always empty the smaller set.  This is critical as otherwise
        // complexity of join would be O( n * N ) and not O( n * log(N) ).
        let (a, b) = if self.size_of_set_containing(a) < self.size_of_set_containing(b) {
            (b, a)
        } else {
            (a, b)
        };
        // The ID of the set that we are merging the smaller set into.
        let joined_set_id = self.set_id(a);
        // The ID of the set that will be empty after this join.
        let emptied_set_id = self.set_id(b);
        if joined_set_id == emptied_set_id {
            return;
        }

        // Now do the move.

        let new_size = self.size_of_set_containing(a) + self.size_of_set_containing(b);
        self.next_element.swap(a, b);
        let mut next = self.next_element(a);
        loop {
            if self.set_id(next) == joined_set_id {
                break;
            }
            self.set_id[next] = joined_set_id as SetId;
            next = self.next_element(next);
        }

        // Update set sizes.
        self.size[joined_set_id] = new_size as u32;
        self.size[emptied_set_id] = 0;
    }

    /// Return an iterator over a single arbitrary member from each set.
    pub fn set_reps_iter(&self) -> impl Iterator<Item = usize> + '_ {
        // This implementation relies on several implementation details:
        // - we only join sets, never split them
        // - the ith set always starts out containing the element i
        // This implies that the ith set will always contain the element i, so
        // the ids of all the non-empty sets is equivalent to a collection of
        // one element from every set.
        self.size
            .iter()
            .enumerate()
            .filter_map(|(set_id, size)| (*size > 0).then_some(set_id))
    }

    /// Return a vector containing a single arbitrary member from each set.
    pub fn set_reps(&self) -> Vec<usize> {
        self.set_reps_iter().collect()
    }

    /// Return the number of non-empty sets.
    pub fn n_sets(&self) -> usize {
        self.size.iter().filter(|s| **s > 0).count()
    }

    /// Return the size of the set containing the element a.
    pub fn size_of_set_containing(&self, a: usize) -> usize {
        self.size[self.set_id(a)] as usize
    }

    /// Return an iterator over all members of the set that contains the element a.
    pub fn members_of_set_containing(&self, a: usize) -> SetIterator<'_> {
        SetIterator::new(self, a)
    }

    /// Return the membership of every set, as a nested iterator performing no
    /// allocations.
    pub fn all_sets(&self) -> impl Iterator<Item = SetIterator<'_>> {
        self.set_reps_iter()
            .map(move |r| self.members_of_set_containing(r))
    }

    /// Return the membership of every set, as an iterator over Vec<i32>.
    /// A lot of the existing codebase expects the elements to be i32.
    /// This is a convenience adapter.
    pub fn all_sets_iter_i32(&self) -> impl Iterator<Item = Vec<i32>> + '_ {
        self.all_sets()
            .map(|set| set.map(|e| e as i32))
            .map(Iterator::collect)
    }

    /// Return the ID of the set that currently contains element a.
    pub fn set_id(&self, a: usize) -> usize {
        self.set_id[a] as usize
    }

    fn next_element(&self, a: usize) -> usize {
        self.next_element[a] as usize
    }
}

/// An iterator over the elements of the disjoint set containing starting_element.
///
/// Since the iterator knows the starting element, it can return the number of
/// elements in the set. It can also produce a copy of the iterator that has
/// been rewound to the beginning of iteration.
pub struct SetIterator<'a> {
    eq: &'a EquivRel,
    starting_element: usize,
    last_element: Option<usize>,
    count: usize,
}

impl<'a> SetIterator<'a> {
    fn new(eq: &'a EquivRel, starting_element: usize) -> Self {
        Self {
            eq,
            starting_element,
            last_element: None,
            count: 0,
        }
    }

    /// Return the number of elements in the set being iterated over.
    pub fn size(&self) -> usize {
        self.eq.size_of_set_containing(self.starting_element)
    }

    /// Return a copy of this set iterator, rewound to the beginning.
    pub fn duplicate(&self) -> Self {
        Self::new(self.eq, self.starting_element)
    }
}

impl Iterator for SetIterator<'_> {
    type Item = usize;

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.eq.size_of_set_containing(self.starting_element) - self.count;
        (remaining, Some(remaining))
    }

    fn next(&mut self) -> Option<Self::Item> {
        match self.last_element {
            Some(last) => {
                let next = self.eq.next_element(last);
                if next == self.starting_element {
                    None
                } else {
                    self.last_element = Some(next);
                    self.count += 1;
                    Some(next)
                }
            }
            None => {
                self.last_element = Some(self.starting_element);
                self.count += 1;
                Some(self.starting_element)
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PERFORMANCE COMPARISONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Comparison to the disjoint-sets crate.  Note that with disjoint sets, it is not
// clear how to find just one set, or the size of one set.  Two comparisons
// were carried out.  The first comparison follows.  Briefly, it shows that
// disjoint-sets is faster for this use case, but uses more memory, and in short,
// the conclusion is a "toss-up".

/*

    // Setup.

    const N : usize = 10_000_000;
    const L : usize = 20_000_000;
    let mut x : Vec<i64> = vec![ 0; L ];
    make_random_vec( &mut x, L );
    for i in 0..L { x[i] = x[i].abs() % (N as i64); }
    let peak = peak_mem_usage_bytes();
    let t = Instant::now( );

    // EquivRel version (use this or the following)
    // there are 1618950 sets
    // 2.6 seconds, delta peak mem = 137 Mb

    let mut e = EquivRel::new(N as i32);
    for j in 0..L/2 { e.join( x[j] as i32, x[j+L/2] as i32 ); }
    let mut reps = Vec::<i32>::new();
    e.orbit_reps( &mut reps );
    println!( "there are {} sets", reps.len() );
    let mut o = Vec::<i32>::new();
    for i in 0..reps.len() { e.orbit( reps[i] as i32, &mut o ); }

    // UnionFind version (use this or the previous);
    // there are 1618950 orbits
    // 1.5 seconds, delta peak mem = 258 Mb
    // disjoint-sets = "0.4.2"

    use disjoint_sets::UnionFind;
    let mut uf = UnionFind::<u32>::new(N as usize);
    for j in 0..L/2 { uf.union( x[j] as u32, x[j+L/2] as u32 ); }
    let reps = uf.to_vec();
    let mut repsx = Vec::<(u32,u32)>::new();
    for i in 0..reps.len() { repsx.push( (reps[i],i as u32) ); }
    repsx.sort();
    let mut o = Vec::<u32>::new();
    let mut orbits = Vec::<Vec<u32>>::new();
    for i in 0..repsx.len() {
        if i > 0 && repsx[i].0 != repsx[i-1].0 {
            orbits.push( o.clone() );
            o.clear();
        }
        o.push( repsx[i].1 );
    }
    orbits.push( o.clone() );
    println!( "there are {} orbits", orbits.len() );

    // Summarize.

    let delta_peak = peak_mem_usage_bytes() - peak;
    println!(
        "{} seconds used, delta peak mem = {} bytes", t.elapsed().as_secs_f64(), delta_peak );

*/

// The second comparison involves a change to the code in hyper.rs.  Here is the
// relevant chunk of code, using EquivRel:

/*

    // Find nodes in the transformed graph.  They are orbits of edge ends under
    // the natural equivalence relation.

    let mut eq : EquivRel = EquivRel::new( 2 * edges.len() as i32 );
    for i in 0..adj.len() {
        let left = adj[i].0;
        let right = adj[i].1;
        eq.join( 2*left + 1, 2*right );
    }
    let mut reps = Vec::<i32>::new();
    eq.orbit_reps( &mut reps );

    // Now actually create the transformed graph.

    g_out.clear();
    g_out.reserve_exact_nodes( reps.len() );
    g_out.reserve_exact_edges( edges.len() );
    for i in 0..reps.len() { g_out.add_node(i as u32); }
    for e in 0..edges.len() {
        let v = bin_position( &reps, &eq.class_id((2*e) as i32) );
        let w = bin_position( &reps, &eq.class_id((2*e+1) as i32) );
        g_out.add_edge( NodeIndex::<u32>::new(v as usize),
            NodeIndex::<u32>::new(w as usize), edges[e].2.clone() );
    }
}

*/

// and here is the relevant chunk of code using disjoint-sets:

/*

    use disjoint_sets::UnionFind;

    // Find nodes in the transformed graph.  They are orbits of edge ends under
    // the natural equivalence relation.

    let mut eq = UnionFind::<u32>::new( 2 * edges.len() );
    for i in 0..adj.len() {
        let left = adj[i].0 as u32;
        let right = adj[i].1 as u32;
        eq.union( 2*left + 1, 2*right );
    }
    let mut reps = eq.to_vec(); // list of orbit representatives
    reps.sort();

    // Now actually create the transformed graph.

    g_out.clear();
    g_out.reserve_exact_nodes( reps.len() );
    g_out.reserve_exact_edges( edges.len() );
    for i in 0..reps.len() { g_out.add_node(i as u32); }
    for e in 0..edges.len() {
        let v = bin_position( &reps, &eq.find( (2*e) as u32) );
        let w = bin_position( &reps, &eq.find( (2*e+1) as u32) );
        g_out.add_edge( NodeIndex::<u32>::new(v as usize),
            NodeIndex::<u32>::new(w as usize), edges[e].2.clone() );

*/

// Performance comparison: the test consisted of running the assemblies for
// 14 VDJ samples.  The actual test is not described here, but the results are
// as follows:
//
// version         server seconds   peak mem GB
// EquivRel        1366.10           7.65
// disjoint-sets   1570.20          13.45
//
// Of course one ought to be able to define a reproducible test that exhibits this
// performance difference.
