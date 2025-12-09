//! graph_simple
// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Define generic digraph functions.
//
// Some of the functions here are unnecessarily quadratic in the vertex degree for
// petgraph as a result of calling v_from and related functions below.  The
// quadratic behavior might be avoided but might make the code a bit less readable.
//
// These functions seem unnecessarily specialized to u32.

use petgraph::EdgeType;
use petgraph::prelude::*;
use std::collections::HashSet;

/// A simple graph.
pub trait GraphSimple<T> {
    /// Return the object associated to an edge id.
    fn edge_obj(&self, e: u32) -> &T;

    /// Return the source of an edge.
    fn to_left(&self, e: u32) -> u32;

    /// Return the target of an edge.
    fn to_right(&self, e: u32) -> u32;

    /// Return the number of edges exiting a given vertex.
    fn n_from(&self, v: usize) -> usize;

    /// Return the number of edges entering a given vertex.
    fn n_to(&self, v: usize) -> usize;

    /// Return id of the nth vertex exiting a given vertex id.
    /// Note that this is O(n).
    fn v_from(&self, v: usize, n: usize) -> usize;

    /// Return id of the nth vertex entering a given vertex id.
    /// Note that this is O(n).
    fn v_to(&self, v: usize, n: usize) -> usize;

    /// Return id of the nth edge exiting a given vertex id.
    /// Note that this is O(n).
    fn e_from(&self, v: usize, n: usize) -> usize;

    /// Return id of the nth edge entering a given vertex id.
    /// Note that this is O(n).
    fn e_to(&self, v: usize, n: usize) -> usize;

    /// get_predecessors: find all vertices which have a directed path to a vertex
    /// in v. This includes the vertices in v by definition. Return a sorted list x.
    fn get_predecessors(&self, v: &[i32], x: &mut Vec<u32>);

    /// start from one vertex.
    fn get_predecessors1(&self, v: i32, x: &mut Vec<u32>);

    /// get_successors: go the other way.
    fn get_successors(&self, v: &[i32], x: &mut Vec<u32>);

    /// start from one vertex.
    fn get_successors1(&self, v: i32, x: &mut Vec<u32>);

    /// Find the connected components.  Each component is a sorted list of vertices.
    fn components(&self, comp: &mut Vec<Vec<u32>>);

    /// Find the connected components as lists of edges.  Each component is an
    /// UNSORTED list of edges.
    fn components_e(&self, comp: &mut Vec<Vec<u32>>);
}

impl<S, T, U, V> GraphSimple<T> for Graph<S, T, U, V>
where
    U: EdgeType,
    V: petgraph::csr::IndexType,
{
    fn edge_obj(&self, e: u32) -> &T {
        &self[EdgeIndex::<V>::new(e as usize)]
    }

    fn to_left(&self, e: u32) -> u32 {
        self.edge_endpoints(EdgeIndex::<V>::new(e as usize))
            .unwrap()
            .0
            .index() as u32
    }

    fn to_right(&self, e: u32) -> u32 {
        self.edge_endpoints(EdgeIndex::<V>::new(e as usize))
            .unwrap()
            .1
            .index() as u32
    }

    fn n_from(&self, v: usize) -> usize {
        self.neighbors(NodeIndex::<V>::new(v)).count()
    }

    fn n_to(&self, v: usize) -> usize {
        self.neighbors_directed(NodeIndex::<V>::new(v), Incoming)
            .count()
    }

    fn v_from(&self, v: usize, n: usize) -> usize {
        self.edges_directed(NodeIndex::<V>::new(v), Outgoing)
            .nth(n)
            .unwrap()
            .target()
            .index()
    }

    fn v_to(&self, v: usize, n: usize) -> usize {
        self.edges_directed(NodeIndex::<V>::new(v), Incoming)
            .nth(n)
            .unwrap()
            .source()
            .index()
    }

    fn e_from(&self, v: usize, n: usize) -> usize {
        let mut e: EdgeIndex<V> = self.first_edge(NodeIndex::<V>::new(v), Outgoing).unwrap();
        for _j in 0..n {
            let f = self.next_edge(e, Outgoing).unwrap();
            e = f;
        }
        e.index()
    }

    fn e_to(&self, v: usize, n: usize) -> usize {
        let mut e: EdgeIndex<V> = self.first_edge(NodeIndex::<V>::new(v), Incoming).unwrap();
        for _j in 0..n {
            let f = self.next_edge(e, Incoming).unwrap();
            e = f;
        }
        e.index()
    }

    fn get_predecessors(&self, v: &[i32], x: &mut Vec<u32>) {
        let mut check: Vec<u32> = Vec::new();
        let mut tov: HashSet<u32> = HashSet::new();
        for v_i in v {
            let s: u32 = *v_i as u32;
            check.push(s);
            tov.insert(s);
        }
        while let Some(x) = check.pop() {
            let n = self.n_to(x as usize);
            for i in 0..n {
                let y = self.v_to(x as usize, i);
                if tov.contains(&(y as u32)) {
                    continue;
                }
                check.push(y as u32);
                tov.insert(y as u32);
            }
        }
        x.clear();
        for v in tov {
            x.push(v);
        }
        x.sort_unstable();
    }

    fn get_predecessors1(&self, v: i32, x: &mut Vec<u32>) {
        let vs = vec![v];
        self.get_predecessors(&vs, x);
    }

    fn get_successors(&self, v: &[i32], x: &mut Vec<u32>) {
        let mut check: Vec<u32> = Vec::new();
        let mut fromv: HashSet<u32> = HashSet::new();
        for v_i in v {
            let s: u32 = *v_i as u32;
            check.push(s);
            fromv.insert(s);
        }
        while let Some(x) = check.pop() {
            let n = self.n_from(x as usize);
            for i in 0..n {
                let y = self.v_from(x as usize, i);
                if fromv.contains(&(y as u32)) {
                    continue;
                }
                check.push(y as u32);
                fromv.insert(y as u32);
            }
        }
        x.clear();
        for v in fromv {
            x.push(v);
        }
        x.sort_unstable();
    }

    fn get_successors1(&self, v: i32, x: &mut Vec<u32>) {
        let vs = vec![v];
        self.get_successors(&vs, x);
    }

    fn components(&self, comp: &mut Vec<Vec<u32>>) {
        comp.clear();
        let mut used: Vec<bool> = vec![false; self.node_count()];
        let mut c: Vec<u32> = Vec::new();
        let mut cnext: Vec<u32> = Vec::new();
        for v in 0..self.node_count() {
            if used[v] {
                continue;
            }
            c.clear();
            cnext.clear();
            cnext.push(v as u32);
            while let Some(w) = cnext.pop() {
                if used[w as usize] {
                    continue;
                }
                used[w as usize] = true;
                c.push(w);
                let n = self.n_from(w as usize);
                for j in 0..n {
                    cnext.push(self.v_from(w as usize, j) as u32);
                }
                let n = self.n_to(w as usize);
                for j in 0..n {
                    cnext.push(self.v_to(w as usize, j) as u32);
                }
            }
            c.sort_unstable();
            comp.push(c.clone());
        }
    }

    fn components_e(&self, comp: &mut Vec<Vec<u32>>) {
        self.components(comp);
        for comp_i in comp {
            let mut c = Vec::<u32>::new();
            for comp_i_j in comp_i.iter() {
                let v = *comp_i_j;
                let n = self.n_from(v as usize);
                for l in 0..n {
                    c.push(self.e_from(v as usize, l) as u32);
                }
            }
            *comp_i = c;
        }
    }
}
