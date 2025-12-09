// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// This file provides the single function graph_filter.

use super::BarcodeFilter;
use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::BarcodeContigs;
use graph_simple::GraphSimple;
use io_utils::fwriteln;
use petgraph::prelude::*;
use rayon::prelude::*;
use std::cmp::{max, min};
use std::io::Write;
use vector_utils::{bin_member, bin_position, erase_if, lower_bound, next_diff12_3, reverse_sort};

/// Create a digraph which has one vertex for each V..J that appears in a productive
/// pair, and for a given light chain and a given heavy chain vertex, a weighted edge
/// from the light to the heavy, so long as they co-occur in some cell; the weight is
/// a pair (numi, ncells), where numi is the sum over all cells for which the co-occur,
/// of the minimum of the number of UMIs that support the two chains.
/// seqs = { (V..J, is_igh) }
/// Also tracking CDR3_AAs for now for exploratory purposes at least.  This makes
/// things less efficient.
///
/// As part of this, identify weak links and kill them.
///
/// This code tiptoes around the fact that whereas these calculations should be at the clonotype
/// level, we compute here at the exact subclonotype level.  The use of V segments is part of this.
///
/// Hmm, seems like the edges go from heavy to light.
#[derive(Default)]
pub(crate) struct GraphFilter {
    pub(crate) log_graph: bool,
}

impl BarcodeFilter<BarcodeContigs> for GraphFilter {
    fn fate_keys(&self, item: &BarcodeContigs) -> impl Iterator<Item = (usize, String)> {
        std::iter::once((item[0].dataset_index, item[0].barcode.clone()))
    }

    fn filter(&self, tig_bc: &[BarcodeContigs]) -> Vec<Option<BarcodeFate>> {
        let mut seqs = Vec::<(&[u8], bool, &str, usize)>::new();
        for tigi in tig_bc {
            for x in tigi {
                seqs.push((x.seq(), x.left, x.cdr3_aa.as_str(), x.v_ref_id));
            }
        }
        seqs.par_sort();
        seqs.dedup();

        // If there are multiple seqs entries whose first three elements agree,
        // delete all but the first.

        let mut to_delete = vec![false; seqs.len()];
        let mut i = 0;
        while i < seqs.len() {
            let mut j = i + 1;
            while j < seqs.len() {
                if seqs[j].0 != seqs[i].0 || seqs[j].1 != seqs[i].1 || seqs[j].2 != seqs[i].2 {
                    break;
                }
                j += 1;
            }
            for d in &mut to_delete[i + 1..j] {
                *d = true;
            }
            i = j;
        }
        erase_if(&mut seqs, &to_delete);

        // Proceed.
        let results: Vec<_> = tig_bc
            .par_iter()
            .map(|tig_bc| {
                let mut edges = Vec::new();
                for x1 in tig_bc {
                    if x1.left {
                        let p1 =
                            lower_bound(&seqs, &(x1.seq(), false, x1.cdr3_aa.as_str(), 0)) as usize;
                        for x2 in tig_bc {
                            if !x2.left {
                                let p2 =
                                    lower_bound(&seqs, &(x2.seq(), false, x2.cdr3_aa.as_str(), 0))
                                        as usize;
                                edges.push((p1, p2, min(x1.umi_count, x2.umi_count)));
                            }
                        }
                    }
                }
                edges
            })
            .collect();
        let mut edges0: Vec<_> = results.into_iter().flatten().collect();
        edges0.sort_unstable();
        let mut edges1 = Vec::new();
        let mut i = 0;
        while i < edges0.len() {
            let j = next_diff12_3(&edges0, i);
            let mut weight = 0;
            for e in &edges0[i..j] {
                weight += e.2;
            }
            edges1.push((edges0[i].0, edges0[i].1, (weight, j - i)));
            i = j;
        }
        let mut g = Graph::<u32, (usize, usize), Directed, u32>::new();
        g.reserve_exact_nodes(seqs.len());
        g.reserve_exact_edges(edges1.len());
        for i in 0..seqs.len() {
            g.add_node(i as u32);
        }
        for (v, w, weight) in edges1 {
            g.add_edge(NodeIndex::<u32>::new(v), NodeIndex::<u32>::new(w), weight);
        }

        // Kill weak branches from light to heavy chains.  Also kill light chain onesies that
        // have too many heavy chain partners.
        //
        // ********************************************************************************************
        // THIS IS TURNED OFF.  Reason: if a light chain is ubiquitous, then this code would be
        // prejudiciously deleting pairs that use it.  The code kills a lot of real cells.
        // When we turned off this code, we got one additional false positive, but it seems like we
        // "should" have the false positive.  The code also resulted in the creation of a few more
        // 5-chain clonotypes, and lots more 4-chain clonotypes (both in the ~400 dataset run).
        // ********************************************************************************************

        /*
        let mut log = Vec::<u8>::new();
        fwriteln!(log, "\nBRANCHING FROM LIGHT CHAINS");
        const MIN_RATIO_KILL: usize = 8;
        const MAX_KILL: usize = 5;
        const MAX_KILL_CELLS: usize = 2;
        const MAX_PARTNERS: usize = 50;
        let mut kills = Vec::<(usize, usize)>::new();
        let mut badones = Vec::<usize>::new();
        let mut results = Vec::<(usize, Vec<(usize, usize)>, Vec<usize>, Vec<u8>)>::new();
        for v in 0..g.node_count() {
            results.push((v, Vec::new(), Vec::new(), Vec::new()));
        }
        results.par_iter_mut().for_each(|res| {
            let v = res.0;
            let log = &mut res.3;
            if g.n_to(v) > 1 {
                let mut stats = Vec::<(usize, usize, usize)>::new();
                fwriteln!(log, "\nlight chain {} = {}", v, seqs[v].2);
                for i in 0..g.n_to(v) {
                    let (w, e) = (g.v_to(v, i), g.e_to(v, i));
                    let numi = g.edge_obj(e as u32).0;
                    let ncells = g.edge_obj(e as u32).1;
                    fwriteln!(
                        log,
                        "• heavy chain {} = {}, weight = {}/{}",
                        w,
                        seqs[w].2,
                        numi,
                        ncells
                    );
                    stats.push((numi, ncells, w));
                }
                reverse_sort(&mut stats);
                for i in 1..stats.len() {
                    // Below, the part seqs[stats[i].2].3 != seqs[stats[0].2].3
                    // is requiring that the partner V segments are different.  See discussion at
                    // the top around the roles of clonotypes versus exact subclonotypes.

                    if seqs[stats[i].2].3 != seqs[stats[0].2].3 {
                        let numi = stats[i].0;
                        let ncells = stats[i].1;
                        let numi_best = stats[0].0;
                        let ncells_best = stats[0].1;
                        if numi_best >= MIN_RATIO_KILL * max(1, numi) && numi <= MAX_KILL {
                            res.1.push((v, stats[i].2));
                        } else if numi_best >= numi && ncells_best >= MIN_RATIO_KILL * max(1, ncells) {
                            if ncells <= MAX_KILL_CELLS {
                                let w = stats[i].2;
                                if graph {
                                    println!(
                                        "\nkill type 1, from {} to {}, ncells = {}, numi = {}",
                                        seqs[v].2, seqs[w].2, ncells, numi
                                    );
                                    println!(
                                        "killed by {} to {}, ncells = {}, numi = {}",
                                        seqs[v].2, seqs[stats[0].2].2, ncells_best, numi_best
                                    );
                                }
                                res.1.push((v, w));
                            } else {
                                let w = stats[i].2;
                                for j in 0..g.n_from(w) {
                                    let (vx, ex) = (g.v_from(w, j), g.e_from(w, j));
                                    let numix = g.edge_obj(ex as u32).0;
                                    let ncellsx = g.edge_obj(ex as u32).1;
                                    if vx != v
                                        && ncellsx >= MIN_RATIO_KILL * ncells
                                        && numix >= MIN_RATIO_KILL * numi
                                    {
                                        if graph {
                                            println!(
                                                "\nkill type 2, from {} to {}, ncells = {}, numi = {}",
                                                seqs[v].2, seqs[w].2, ncells, numi
                                            );
                                            println!(
                                                "killed by {} to {}, ncells = {}, numi = {}",
                                                seqs[v].2, seqs[stats[0].2].2, ncells_best, numi_best
                                            );
                                        }
                                        res.1.push((v, w));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if g.n_to(v) > MAX_PARTNERS {
                res.2.push(v);
            }
        });
        for i in 0..results.len() {
            kills.append(&mut results[i].1.clone());
            badones.append(&mut results[i].2.clone());
            log.append(&mut results[i].3.clone());
        }
        kills.sort();
        // presumably badly inefficient
        let mut to_delete = vec![false; tig_bc.len()];
        let mut results = Vec::<(usize, bool)>::new();
        for i in 0..tig_bc.len() {
            results.push((i, false));
        }
        results.par_iter_mut().for_each(|res| {
            let i = res.0;
            for j1 in 0..tig_bc[i].len() {
                if tig_bc[i][j1].left {
                    continue;
                }
                let x1 = &tig_bc[i][j1];
                let m1 = (x1.seq.clone(), x1.left, x1.cdr3_aa.clone(), x1.v_ref_id);
                let p1 = bin_position(&seqs, &m1) as usize;
                for j2 in 0..tig_bc[i].len() {
                    if !tig_bc[i][j2].left {
                        continue;
                    }
                    let x2 = &tig_bc[i][j2];
                    let m2 = (x2.seq.clone(), x2.left, x2.cdr3_aa.clone(), x2.v_ref_id);
                    let p2 = bin_position(&seqs, &m2) as usize;
                    if bin_member(&kills, &(p1, p2)) {
                        res.1 = true;
                    }
                }
            }
            if tig_bc[i].len() == 1 {
                let x0 = &tig_bc[i][0];
                let m0 = (x0.seq.clone(), x0.left, x0.cdr3_aa.clone(), x0.v_ref_id);
                let p = bin_position(&seqs, &m0) as usize;
                if bin_member(&badones, &p) {
                    res.1 = true;
                }
            }
        });
        for i in 0..tig_bc.len() {
            to_delete[i] = results[i].1;
            if to_delete[i] {
                ndels += 1;
            }
        }
        for i in 0..tig_bc.len() {
            if to_delete[i] {
                fate[tig_bc[i][0].dataset_index].insert(
                    tig_bc[i][0].barcode.clone(),
                    "failed GRAPH_FILTER filter".to_string(),
                );
            }
        }
        if !ctl.gen_opt.ngraph_filter {
            erase_if(&mut tig_bc, &to_delete);
        }
        if graph {
            fwriteln!(log, "");
            print!("{}", strme(&log));
        }
        */

        // Kill weak branches from heavy to light chains.

        let mut log = Vec::<u8>::new();
        fwriteln!(log, "\nBRANCHING FROM HEAVY CHAINS");
        const MIN_RATIO_KILL_HEAVY: usize = 8;
        const MAX_KILL_HEAVY: usize = 6;
        const MAX_KILL_HEAVY_CELLS: usize = 1;
        let mut kills = Vec::<(usize, usize)>::new();
        let mut results = Vec::<(usize, Vec<(usize, usize)>, Vec<u8>)>::new();
        for v in 0..g.node_count() {
            results.push((v, Vec::new(), Vec::new()));
        }
        results.par_iter_mut().for_each(|res| {
            let v = res.0;
            let log = &mut res.2;
            if g.n_from(v) > 1 {
                let mut stats = Vec::<((usize, usize), usize)>::new();
                fwriteln!(log, "\nheavy chain {} = {}", v, seqs[v].2);
                for i in 0..g.n_from(v) {
                    let w = g.v_from(v, i);
                    let e = g.e_from(v, i);
                    let weight = g.edge_obj(e as u32);
                    fwriteln!(
                        log,
                        "• light chain {} = {}, weight = {}/{}",
                        w,
                        seqs[w].2,
                        weight.1,
                        weight.0
                    );
                    stats.push((*weight, w));
                }
                reverse_sort(&mut stats);
                for i in 1..stats.len() {
                    if (stats[0].0).0 >= MIN_RATIO_KILL_HEAVY * max(1, (stats[i].0).0)
                        && (stats[i].0).0 <= MAX_KILL_HEAVY
                        && (stats[i].0).1 <= MAX_KILL_HEAVY_CELLS
                    {
                        if self.log_graph {
                            let w = stats[i].1;
                            println!(
                                "\nkill type 3, from {} to {}\nkilled by {} to {}",
                                seqs[v].2, seqs[w].2, seqs[v].2, seqs[stats[0].1].2
                            );
                        }
                        res.1.push((v, stats[i].1));
                    }
                }
            }
        });
        for (_, mut r1, mut r2) in results {
            kills.append(&mut r1);
            log.append(&mut r2);
        }
        kills.sort_unstable();
        // presumably badly inefficient
        let to_delete: Vec<_> = tig_bc
            .par_iter()
            .map(|tig_bc| {
                for tbc_0 in tig_bc {
                    if !tbc_0.left {
                        continue;
                    }
                    let x1 = &tbc_0;
                    let m1 = (x1.seq(), x1.left, x1.cdr3_aa.as_str(), x1.v_ref_id);
                    let p1 = bin_position(&seqs, &m1) as usize;
                    for tbc_1 in tig_bc {
                        if tbc_1.left {
                            continue;
                        }
                        let x2 = &tbc_1;
                        let m2 = (x2.seq(), x2.left, x2.cdr3_aa.as_str(), x2.v_ref_id);
                        let p2 = bin_position(&seqs, &m2) as usize;
                        if bin_member(&kills, &(p1, p2)) {
                            return Some(BarcodeFate::GraphFilter);
                        }
                    }
                }
                None
            })
            .collect();

        let ndels: usize = to_delete.iter().flatten().count();

        if self.log_graph {
            fwriteln!(log, "");
            print!("{}", std::str::from_utf8(&log).unwrap());
            println!("total graph filter deletions = {ndels}");
        }

        to_delete
    }
}
