// Print a hyper with annotations.

#![allow(clippy::many_single_char_names)]

use debruijn::dna_string::DnaString;
use debruijn::kmer::Kmer20;
use debruijn::{Mer, Vmer};
use graph_simple::GraphSimple;
use hyperbase::Hyper;
use kmer_lookup::make_kmer_lookup_20_single;
use petgraph::prelude::*;
use std::io::prelude::*;
use vdj_ann::annotate::{annotate_seq, print_annotations};
use vdj_ann::refx::RefData;
use vector_utils::{lower_bound1_3, unique_sort, upper_bound1_3};

// =============================================================================
// Print the graph, component by component, with reference annotations.
// Elide printing of rc components: if a component is not annotated and its
// rc is, don't print it; if neither a component nor its rc is annotated, only
// print one of them.  This uses components_e_pos_sorted, which may turn out
// to be too slow.  Optionally show read support for each edge (which produces
// voluminous output and is highly inefficient).
// =============================================================================

#[allow(clippy::too_many_arguments)]
pub fn print_with_annotations(
    free: bool,
    this: &Hyper,
    reads: &[DnaString],
    umi_id: &[i32],
    refdata: &RefData,
    show_supp: bool,
    print_seq_edges: bool,
    log: &mut Vec<u8>,
) {
    // Find components and index them.

    let mut comp = Vec::<Vec<u32>>::new();
    this.h.g.components_e_pos_sorted(&mut comp);
    let mut to_comp: Vec<i32> = vec![-1; this.h.g.edge_count()];
    for c in 0..comp.len() {
        for j in 0..comp[c].len() {
            to_comp[comp[c][j] as usize] = c as i32;
        }
    }

    // Determine which components have annotations.

    let mut is_ann: Vec<bool> = vec![false; comp.len()];
    if !free {
        for j in 0..comp.len() {
            for i in 0..comp[j].len() {
                let e = comp[j][i];
                let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
                annotate_seq(this.h.g.edge_obj(e), refdata, &mut ann, true, false, true);
                if !ann.is_empty() {
                    is_ann[j] = true;
                }
            }
        }
    }

    // Make lookup table for edges if showing support.

    let mut kmers_plus = Vec::<(Kmer20, i32, i32)>::new();
    if show_supp {
        let mut edges = Vec::<DnaString>::new();
        for e in 0..this.h.g.edge_count() {
            edges.push(this.h.g.edge_obj(e as u32).clone());
        }
        make_kmer_lookup_20_single(&edges, &mut kmers_plus);
    }

    // Print.

    let mut ccount = 0;
    for j in 0..comp.len() {
        // Elide rc components.

        let j_rc = to_comp[this.inv[comp[j][0] as usize] as usize] as usize;
        if !is_ann[j] && is_ann[j_rc] {
            continue;
        }
        if !is_ann[j] && !is_ann[j_rc] && j_rc < j {
            continue;
        }

        // Now actually print.

        ccount += 1;
        fwriteln!(log, "\nCOMPONENT {}", ccount);
        for i in 0..comp[j].len() {
            let e = comp[j][i];
            let v = this.h.g.to_left(e);
            let w = this.h.g.to_right(e);
            let b: DnaString = this.h.g[EdgeIndex::<u32>::new(e as usize)].clone();
            let re = this.inv[e as usize];
            let mut u = Vec::<i32>::new();
            for pass in 0..2 {
                let mut f = e;
                if pass == 1 {
                    f = re;
                }
                if !this.ids[f as usize].is_empty() {
                    u.push(umi_id[this.ids[f as usize][0] as usize]);
                }
                for j in 1..this.ids[f as usize].len() {
                    if umi_id[this.ids[f as usize][j] as usize]
                        != umi_id[this.ids[f as usize][j - 1] as usize]
                    {
                        u.push(umi_id[this.ids[f as usize][j] as usize]);
                    }
                }
            }
            unique_sort(&mut u);
            fwriteln!(
                log,
                "\n{} ==(e={},len={},s={},u={})==> {}",
                v,
                e,
                b.len() - this.h.k as usize + 1,
                this.supp(e as usize),
                u.len(),
                w
            );
            // writeln!( log, "" ); // should only do if there are annotations
            if print_seq_edges {
                fwriteln!(log, "{}", this.h.g.edge_obj(e).to_string());
            }
            if !free {
                let verbose = false;
                print_annotations(this.h.g.edge_obj(e), refdata, log, false, true, verbose);
            }

            // Print edge support.  Only shows fw support, and only shows
            // perfect match intervals, of length >= 20.
            // ◼ There are multiple places in the code where we find maximal
            // ◼ perfect match intervals.  Should make into a function.

            if show_supp {
                let mut matches = Vec::<(usize, usize, u32, i32)>::new();
                let k = this.h.k as usize;
                let r = &this.h.g.edge_obj(e);
                for j in 0..this.ids[e as usize].len() {
                    let id = this.ids[e as usize][j];
                    let b = &reads[id as usize];
                    for l in 0..b.len() - k + 1 {
                        let x: Kmer20 = b.get_kmer(l);
                        let low = lower_bound1_3(&kmers_plus, &x);
                        let high = upper_bound1_3(&kmers_plus, &x);
                        for m in low..high {
                            let t = kmers_plus[m as usize].1 as usize;
                            if t != e as usize {
                                continue;
                            }
                            let p = kmers_plus[m as usize].2 as usize;
                            if l > 0 && p > 0 && b.get(l - 1) == r.get(p - 1) {
                                continue;
                            }
                            let mut len = k;
                            while l + len < b.len() && p + len < r.len() {
                                if b.get(l + len) != r.get(p + len) {
                                    break;
                                }
                                len += 1;
                            }
                            // off is inferred read start on edge
                            let off = p as i32 - l as i32;
                            matches.push((p, len, id, off));
                        }
                    }
                }
                matches.sort_unstable();
                if !matches.is_empty() {
                    fwriteln!(log, "");
                }
                let mut i = 0;
                while i < matches.len() {
                    let mut j = i + 1;
                    while j < matches.len() {
                        if matches[j].0 != matches[i].0 {
                            break;
                        }
                        if matches[j].1 != matches[i].1 {
                            break;
                        }
                        if matches[j].3 != matches[i].3 {
                            break;
                        }
                        if umi_id[matches[j].2 as usize] != umi_id[matches[i].2 as usize] {
                            break;
                        }
                        j += 1;
                    }
                    let p = matches[i].0;
                    let len = matches[i].1;
                    let id = matches[i].2;
                    let u = umi_id[id as usize];
                    let off = matches[i].3;
                    if j - i == 1 {
                        fwriteln!(
                            log,
                            "{}-{}, off = {}, u = {}, id = {}",
                            p,
                            p + len,
                            off,
                            u,
                            id
                        );
                    } else {
                        let mut ids = Vec::<u32>::new();
                        for m in &matches[i..j] {
                            ids.push(m.2);
                        }
                        fwriteln!(
                            log,
                            "{}-{}, off = {}, u = {}, id = {:?}",
                            p,
                            p + len,
                            off,
                            u,
                            ids
                        );
                    }
                    i = j;
                }
            }
        }
    }
}
