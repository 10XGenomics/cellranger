// This file contains code to make contigs.

#![allow(clippy::many_single_char_names)]
// TODO: fix these.
#![allow(clippy::needless_range_loop)]

use self::equiv::EquivRel;
use self::graph_simple::GraphSimple;
use self::hyperbase::Hyper;
use self::kmer_lookup::make_kmer_lookup_20_single;
use crate::barcode_data::BarcodeData;
use crate::constants::UmiType;
use crate::heuristics::Heuristics;
use crate::log_opts::LogOpts;
use crate::ref_free::{alt_strong_paths, uber_strong_paths};
use crate::sw::pos_base_quals_helper;
use crate::utils;
use debruijn::dna_string::{ndiffs, DnaString};
use debruijn::kmer::Kmer20;
use debruijn::{Mer, Vmer};
use itertools::Itertools;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use string_utils::stringme;
use vdj_ann::annotate::{annotate_seq, get_cdr3, Annotation};
use vdj_ann::refx::RefData;
use vdj_ann::transcript::{is_productive_contig, junction_seq, junction_supp, junction_supp_core};
use vector_utils::{
    bin_member, bin_position, contains, erase_if, lower_bound1_3, next_diff, next_diff1_2,
    next_diff1_3, next_diff1_6, reverse_sort, sort_sync2, unique_sort, upper_bound1_3, VecUtils,
};
use {equiv, graph_simple, hyperbase, kmer_lookup};

// Given a vector x of length N > n, return n randomly selected elements, and sort the resulting
// vector.  This is deterministic and will not change when crates are updated.

pub fn random_subvec(x: &[String], n: usize) -> Vec<String> {
    assert!(x.len() > n);
    let mut s = HashSet::<usize>::new();
    let mut y = Vec::<String>::new();
    let mut r = 0_i64;
    for _ in 0..n {
        loop {
            r = 6_364_136_223_846_793_005i64
                .wrapping_mul(r)
                .wrapping_add(1_442_695_040_888_963_407);
            let p = (r as usize) % x.len();
            if !s.contains(&p) {
                s.insert(p);
                y.push(x[p].clone());
                break;
            }
        }
    }
    y.sort();
    y
}

// Make contigs and clean them up using the reference, then determine if there
// is a pairing.

#[allow(clippy::too_many_arguments)]
pub fn make_contigs(
    // INPUTS:
    single_end: bool,
    is_gd: Option<bool>,
    primers: &[Vec<u8>],
    reads: &[DnaString],
    quals: &[Vec<u8>],
    x: &mut Hyper,
    umi_id: &[i32],
    uu: &[String],
    refdata: &RefData,
    refdata_full: &RefData,
    rkmers_plus_full_20: &[(Kmer20, i32, i32)],
    heur: &Heuristics,
    min_contig_length: Option<usize>,

    // OUTPUTS:
    con: &mut Vec<DnaString>,                  // good contigs (or all?)
    jsupp: &mut Vec<(i32, i32)>,               // junction support for good contigs
    con2: &mut Vec<DnaString>,                 // reject contigs
    conxq: &mut Vec<Vec<u8>>,                  // qual score for all contigs
    cumi: &mut Vec<Vec<i32>>,                  // umis assigned to each contig
    cids: &mut Vec<Vec<i32>>,                  // reads assigned to each contig
    barcode_data: &mut BarcodeData,            // data on barcode including pairing
    validated_umis: &mut Vec<Vec<String>>,     // validated UMIs
    non_validated_umis: &mut Vec<Vec<String>>, // non-validated UMIs
    invalidated_umis: &mut Vec<Vec<String>>,   // invalidated UMIs

    // LOGGING:
    log: &mut Vec<u8>,
    log_opts: &LogOpts,
) {
    let gd_mode = is_gd.unwrap_or(false);
    // Unpack refdata.

    let refs_full = &refdata_full.refs;
    let rheaders_full = &refdata_full.rheaders;
    let rheaders = &refdata.rheaders;

    // Make rc primers.

    let mut rc_primers = Vec::<Vec<u8>>::new();
    for i in 0..primers.len() {
        let mut p = primers[i].clone();
        p.reverse();
        for j in 0..p.len() {
            if p[j] == b'A' {
                p[j] = b'T';
            } else if p[j] == b'C' {
                p[j] = b'G';
            } else if p[j] == b'G' {
                p[j] = b'C';
            } else if p[j] == b'T' {
                p[j] = b'A';
            }
        }
        rc_primers.push(p);
    }

    // Find simple paths, mimicking trace_umis.

    let mut simple_paths = Vec::<Vec<i32>>::new();
    let mut z = umi_id.to_vec();
    unique_sort(&mut z);
    let mut ue = vec![Vec::<usize>::new(); z.len()];
    for e in 0..x.h.g.edge_count() {
        for j in 0..x.ids[e].len() {
            let u = umi_id[x.ids[e][j] as usize];
            ue[u as usize].push(e);
        }
    }
    for u in 0..z.len() {
        unique_sort(&mut ue[u]);
        let es = &ue[u];
        let mut used = vec![false; es.len()];
        for i in 0..es.len() {
            if used[i] {
                continue;
            }
            used[i] = true;
            let mut p = vec![es[i]];
            loop {
                let mut exts = Vec::<usize>::new();
                for j in 0..es.len() {
                    if p.contains(&es[j]) || used[j] {
                        continue;
                    }
                    if x.h.g.to_right(es[j] as u32) == x.h.g.to_left(p[0] as u32) {
                        exts.push(j);
                    }
                }
                if exts.len() != 1 {
                    break;
                }
                used[exts[0]] = true;
                p.insert(0, es[exts[0]]);
            }
            loop {
                let mut exts = Vec::<usize>::new();
                for j in 0..es.len() {
                    if p.contains(&es[j]) || used[j] {
                        continue;
                    }
                    if x.h.g.to_left(es[j] as u32) == x.h.g.to_right(p[p.len() - 1] as u32) {
                        exts.push(j);
                    }
                }
                if exts.len() != 1 {
                    break;
                }
                used[exts[0]] = true;
                p.push(es[exts[0]]);
            }
            let mut pp = Vec::<i32>::new();
            for i in 0..p.len() {
                pp.push(p[i] as i32);
                // println!( "{}", p.iter().format(",")  ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            }
            simple_paths.push(pp);
        }
    }
    // println!(""); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Find prohibited bubble traversals.  Here we look for the following,
    // for which the data support exactly two traversals:
    //
    //     e1       f       g1
    // v1 ====> v2 ---> v3 ====> v4
    //     e2               g2
    //
    // where two traversals are strongly supported by umis and the other two are not.

    const MIN_RATIO: usize = 5;
    const MAX_KILL: usize = 4;
    let g = &x.h.g;
    let mut prohibited = Vec::<Vec<i32>>::new();
    for v1 in 0..g.node_count() {
        if g.n_from(v1) != 2 {
            continue;
        }
        let (e1, e2) = (g.e_from(v1, 0), g.e_from(v1, 1));
        let v2 = g.v_from(v1, 0);
        if v2 != g.v_from(v1, 1) || g.n_to(v2) != 2 || g.n_from(v2) != 1 {
            continue;
        }
        let f = g.e_from(v2, 0);
        let v3 = g.v_from(v2, 0);
        if g.n_from(v3) != 2 {
            continue;
        }
        let (g1, g2) = (g.e_from(v3, 0), g.e_from(v3, 1));
        let v4 = g.v_from(v3, 0);
        if v4 != g.v_from(v3, 1) || g.n_to(v4) != 2 {
            continue;
        }
        let p11 = [e1 as i32, f as i32, g1 as i32].to_vec();
        let p12 = [e1 as i32, f as i32, g2 as i32].to_vec();
        let p21 = [e2 as i32, f as i32, g1 as i32].to_vec();
        let p22 = [e2 as i32, f as i32, g2 as i32].to_vec();
        let (mut n11, mut n12) = (0, 0);
        let (mut n21, mut n22) = (0, 0);
        for x in &simple_paths {
            if contains(x, &p11) {
                n11 += 1;
            }
            if contains(x, &p12) {
                n12 += 1;
            }
            if contains(x, &p21) {
                n21 += 1;
            }
            if contains(x, &p22) {
                n22 += 1;
            }
        }
        if min(n11, n22) >= MIN_RATIO * max(n12, n21) && n11 >= 1 {
            if n12 <= MAX_KILL {
                prohibited.push(p12.clone());
            }
            if n21 <= MAX_KILL {
                prohibited.push(p21.clone());
            }
        }
        if min(n12, n21) >= MIN_RATIO * max(n11, n22) && n12 >= 1 {
            if n11 <= MAX_KILL {
                prohibited.push(p11.clone());
            }
            if n22 <= MAX_KILL {
                prohibited.push(p22.clone());
            }
        }
    }

    // Define another class of prohibited paths.  Here we look for a vertex v having two
    // edges e1, e2 entering it, and one edge f exiting it, // in which f has a V annotation
    // followed by a J annotation.  Then if e1,f or e2,f has twice as many umis as the other,
    // prohibit the weak ei,f path.
    //
    // A possible problem with this is that e1 and e2 don't necessarily have the same length,
    // and so the umi accounting isn/t completely fair.  To avoid this problem, we require
    // that the stronger edge is also shorter.
    //
    // This condition (or perhaps also subjecting simple paths to prohibition) lowered sensitivity
    // by about 0.1%.  It also may increase run time by ~1% because it annotates every edge.

    /*
    const MIN_RATIO_VJ : usize = 2;
    let y = &x.h.g;
    for v in 0..g.node_count() {
        if g.n_to(v) != 2 || g.n_from(v) != 1 {
            continue;
        }
        let (e1, e2) = (g.e_to(v, 0), g.e_to(v, 1));
        let (mut n1, mut n2) = (0, 0);
        for pass in 0..2 {
            let mut u = Vec::<i32>::new();
            let mut e = e1 as usize;
            if pass == 1 {
                e = e2 as usize;
            }
            if !x.ids[e].is_empty() {
                u.push(umi_id[x.ids[e][0] as usize]);
            }
            for j in 1..x.ids[e].len() {
                if umi_id[x.ids[e][j] as usize] != umi_id[x.ids[e][j - 1] as usize] {
                    u.push(umi_id[x.ids[e][j] as usize]);
                }
            }
            unique_sort(&mut u);
            if pass == 0 {
                n1 = u.len();
            } else {
                n2 = u.len();
            }
        }
        let (mut ok1, mut ok2) = (false, false);
        if n1 >= MIN_RATIO_VJ * n2 && n1 > 0
            && y.edge_obj(e1 as u32).len() <= y.edge_obj(e2 as u32).len() {
            ok1 = true;
        } else if n2 >= MIN_RATIO_VJ * n1 && n2 > 0
            && y.edge_obj(e2 as u32).len() <= y.edge_obj(e1 as u32).len() {
            ok2 = true;
        }
        if !ok1 && !ok2 {
            continue;
        }
        let f = g.e_from(v, 0);
        let b = &y.edge_obj(f as u32);
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(&b, refdata, &mut ann, true, false, false);
        let (mut vstart, mut jstart) = (-1 as i32, -1 as i32);
        for i in 0..ann.len() {
            let t = ann[i].2 as usize;
            if refdata.segtype[t] == "V" {
                vstart = ann[i].0 as i32;
            }
            if refdata.segtype[t] == "J" {
                jstart = ann[i].0 as i32;
            }
        }
        if vstart >= 0 && jstart > vstart {
            if ok1 {
                println!( "prohibiting2 {},{}", e2, f );
                printme!( n1, n2 );
                prohibited.push( vec![ e2 as i32, f as i32 ] );
            } else if ok2 {
                println!( "prohibiting1 {},{}", e1, f );
                printme!( n1, n2 );
                prohibited.push( vec![ e1 as i32, f as i32 ] );
            }
        }
    }
    */

    // Make a path for each UMI, and kill duplicates.

    let min_contig_length: usize = min_contig_length.unwrap_or(300);
    let mut strong = Vec::<(i32, Vec<i32>)>::new();

    uber_strong_paths(x, umi_id, &mut strong);
    let mut alt_strong = Vec::<Vec<i32>>::new();
    alt_strong_paths(x, umi_id, &mut alt_strong);

    /*
    // Print:
    println!( "\nuber strong paths:" );
    for i in 0..strong.len() {
        println!( "{} = {} from umi {}", i, strong[i].1.iter().format(","), strong[i].0 );
    }
    println!( "\nalt strong paths:" );
    for i in 0..alt_strong.len() {
        println!( "{} = {}", i, alt_strong[i].iter().format(",") );
    }
    println!( "\ntrace paths:" );
    */

    let mut all = Vec::<Vec<i32>>::new();
    for i in 0..strong.len() {
        all.push(strong[i].1.clone());
    }
    all.append(&mut alt_strong);

    // Add in the simple paths.

    all.append(&mut simple_paths.clone());

    // Remove prohibited paths.

    let mut to_delete = vec![false; all.len()];
    for i in 0..all.len() {
        for j in 0..prohibited.len() {
            if contains(&all[i], &prohibited[j]) {
                to_delete[i] = true;
            }
        }
    }
    erase_if(&mut all, &to_delete);

    // Collect all the strong edges.

    unique_sort(&mut all);
    if log_opts.print_strong {
        for i in 0..strong.len() {
            fwriteln!(log, "strong: {}", strong[i].1.iter().format(","));
        }
    }

    con.clear();
    con2.clear();
    let mut conp = Vec::<Vec<i32>>::new();
    for j in 0..all.len() {
        let c = x.cat(&all[j]);

        // Require that contig contains a primer or a CDR3.

        let mut ok = false;
        let cb = c.to_ascii_vec();
        for i in 0..rc_primers.len() {
            if contains(&cb, &rc_primers[i]) {
                ok = true;
                break;
            }
        }
        if !ok {
            let cdr3 = get_cdr3(&c.slice(0, c.len()));
            if !cdr3.is_empty() {
                ok = true;
            }
        }
        if !ok {
            continue;
        }

        // Save contig if long enough.

        if c.len() >= min_contig_length {
            con.push(c);
            conp.push(all[j].clone());
        }
    }
    // unique_sort_sync( &con, &conp );
    let mut cc = Vec::<(DnaString, Vec<i32>)>::new();
    for i in 0..con.len() {
        cc.push((con[i].clone(), conp[i].clone()));
    }
    cc.sort();
    con.clear();
    conp.clear();
    let mut i = 0;
    while i < cc.len() {
        let j = next_diff1_2(&cc, i);
        con.push(cc[i].0.clone());
        conp.push(cc[i].1.clone());
        i = j;
    }

    // Trim stuff before a UTR.  Also trim stuff after constant regions, which
    // is intended to address very rare cases where a contig is the concatenation
    // of two transcripts in reverse orientation.

    if !heur.free {
        for j in 0..con.len() {
            let ann = annotate_seq(&con[j], refdata, true, false, false);
            if !ann.is_empty() {
                let (mut start, mut stop) = (0, con[j].len());
                let l = ann.len() - 1;
                let (t1, t2) = (ann[0].ref_id as usize, ann[l].ref_id as usize);
                if refdata.segtype[t1] == "U" {
                    start = ann[0].tig_start as usize;
                }
                if refdata.segtype[t2] == "C" {
                    stop = ann[l].tig_start as usize + ann[l].match_len as usize;
                }
                if start > 0 || stop < con[j].len() {
                    con[j] = con[j].slice(start, stop).to_owned();
                }
            }
        }
    }

    // Trim the consensus sequences to stop right after the primer sites.  It would
    // seem like this would not have an effect, but it does, albeit small.

    if heur.prim_trim {
        let mut rc_primers = Vec::<String>::new(); // IS THIS DUPLICATED ABOVE?
        for i in 0..primers.len() {
            let mut p = primers[i].clone();
            p.reverse();
            for j in 0..p.len() {
                if p[j] == b'A' {
                    p[j] = b'T';
                } else if p[j] == b'C' {
                    p[j] = b'G';
                } else if p[j] == b'G' {
                    p[j] = b'C';
                } else if p[j] == b'T' {
                    p[j] = b'A';
                }
            }
            rc_primers.push(stringme(&p));
        }
        for j in 0..con.len() {
            let s: String = con[j].to_string();
            for i in 0..rc_primers.len() {
                let k = rc_primers[i].len();
                if s.len() >= k {
                    for p in 0..s.len() - k + 1 {
                        if p + k <= con[j].len() && s[p..p + k] == rc_primers[i] {
                            con[j] = con[j].slice(0, p + k).to_owned();
                            break;
                        }
                    }
                }
            }
        }
    }

    // Delete contigs having only C annotations.  Also test for V annotation.
    if !heur.free {
        let mut to_delete = vec![false; con.len()];
        for i in 0..con.len() {
            let ann = annotate_seq(&con[i], refdata, true, false, false);
            let mut c_only = true;
            for a in ann {
                let t = a.ref_id as usize;
                if !refdata.is_c(t) {
                    c_only = false;
                }
            }
            if c_only {
                to_delete[i] = true;
            }
        }
        erase_if(con, &to_delete);
        erase_if(&mut conp, &to_delete);
    }

    // If a contig has a single-base indel relative to the reference, that is
    // supported by only UMI, or only one UMI plus one additional read, fix it.

    if !heur.free {
        let k = 20;
        let mut to_edge = HashMap::<Kmer20, usize>::new();
        for e in 0..x.h.g.edge_count() {
            let b = &x.h.g.edge_obj(e as u32);
            for l in 0..b.len() - k + 1 {
                to_edge.insert(b.get_kmer(l), e);
            }
        }
        for i in 0..con.len() {
            let mut edited = false;
            let mut last_match = false;
            let mut next = 0_u8;
            let b = con[i].clone();
            // let mut last_locs = Vec::<(usize,usize)>::new();
            for l in 0..b.len() - k + 1 {
                if edited {
                    break;
                }
                let w: Kmer20 = b.get_kmer(l);
                let low = lower_bound1_3(rkmers_plus_full_20, &w);
                let high = upper_bound1_3(rkmers_plus_full_20, &w);
                let mut nexts = Vec::<u8>::new();
                // let mut locs = Vec::<(usize,usize)>::new();
                for r in low..high {
                    let t = rkmers_plus_full_20[r as usize].1 as usize;
                    if rheaders_full[t].contains("segment") {
                        continue;
                    }
                    let p = rkmers_plus_full_20[r as usize].2 as usize;
                    if p + k < refs_full[t].len() {
                        // locs.push( (t,p) );
                        nexts.push(refs_full[t].get(p + k));
                    }
                }
                unique_sort(&mut nexts);
                if nexts.solo() {
                    next = nexts[0];
                    last_match = true;
                } else if !last_match {
                    continue;
                } else if nexts.is_empty() {
                    last_match = false;
                }
                if !nexts.is_empty() {
                    continue;
                }
                let pos = l;

                // Two passes.  First pass is delete the last base.
                // Second pass is insert the reference base.

                for pass in 0..2 {
                    let mut c = DnaString::new();
                    for m in 0..pos + k - 1 {
                        c.push(con[i].get(m));
                    }
                    if pass == 0 {
                        for m in pos + k..con[i].len() {
                            c.push(con[i].get(m));
                        }
                    } else {
                        c.push(next);
                        for m in pos + k - 1..con[i].len() {
                            c.push(con[i].get(m));
                        }
                    }
                    let req = 10; // required extension
                    if l + k + req > c.len() {
                        continue;
                    }
                    let mut ok = true;
                    for m in l..l + req {
                        let q: Kmer20 = c.get_kmer(m);
                        let low = lower_bound1_3(rkmers_plus_full_20, &q);
                        let high = upper_bound1_3(rkmers_plus_full_20, &q);
                        if low == high {
                            ok = false;
                            break;
                        }
                    }
                    if !ok {
                        continue;
                    }

                    // Test umi support.

                    let mut e = -1_i32;
                    if to_edge.contains_key(&w) {
                        e = to_edge[&w] as i32;
                    }
                    let mut umis = Vec::<i32>::new();
                    let mut ucount = Vec::<usize>::new();
                    if e >= 0 {
                        for l in 0..x.ids[e as usize].len() {
                            let id = x.ids[e as usize][l] as usize;
                            let u = umi_id[id];
                            let numi = umis.len();
                            if numi > 0 && u == umis[numi - 1] && ucount[numi - 1] >= 2 {
                                continue;
                            }
                            let r = &reads[id];
                            if r.len() >= k {
                                for z in 0..r.len() - k + 1 {
                                    let y: Kmer20 = r.get_kmer(z);
                                    if y == w {
                                        if numi > 0 && u == umis[numi - 1] {
                                            ucount[numi - 1] += 1;
                                        } else {
                                            umis.push(u);
                                            ucount.push(1);
                                        }
                                        break;
                                    }
                                }
                            }
                            if ucount.len() > 2 {
                                break;
                            }
                        }
                    }
                    if ucount.len() > 2 {
                        continue;
                    }
                    if ucount.duo() && ucount[0] > 1 && ucount[1] > 1 {
                        continue;
                    }

                    // Make edit.

                    con[i] = c;
                    edited = true;
                    break;
                }
                if edited {
                    break;
                }
            }
        }
    }

    // Filter out contigs that are now too short.

    let mut to_delete = vec![false; con.len()];
    for i in 0..con.len() {
        if con[i].len() < min_contig_length {
            to_delete[i] = true;
        }
    }
    erase_if(con, &to_delete);
    erase_if(&mut conp, &to_delete);

    // Move invalid contigs to the reject pile.

    let mut to_delete: Vec<bool> = vec![false; con.len()];
    let mut ann_all = Vec::<Vec<Annotation>>::new();
    for i in 0..con.len() {
        let mut good = true;
        if !heur.free {
            let ann = annotate_seq(&con[i], refdata, true, false, true);
            if !is_productive_contig(&con[i], refdata, &ann).0 {
                good = false;
            }
            ann_all.push(ann);
        } else {
            let cdr3 = get_cdr3(&con[i].slice(0, con[i].len()));
            if cdr3.is_empty() {
                good = false;
            }
        }
        if !good {
            to_delete[i] = true;
            con2.push(con[i].clone());
        }
    }
    erase_if(con, &to_delete);
    erase_if(&mut conp, &to_delete);
    erase_if(&mut ann_all, &to_delete);

    /*
    // Print:
    println!( "\nCONTIGS BEFORE UNIQIFICATION" );
    for u in 0..con.len() {
        println!( "\ncontig{}", u+1 );
        let mut log = Vec::<u8>::new();
        print_annotations( &con[u], &refdata, &mut log, true, false, false );
        println!( "{}", strme(&log) );
    }
    */

    // Compute the number of umis assigned to each contig.  To do this, for each
    // read, we compute the number of edges it shares with each contig path.  Each
    // pair (umi,contig) is assigned the sum of edge counts for the umi.  For a
    // given umi, all the contigs sharing the top score get that umi assigned to
    // them.
    //
    // ◼ There is another assignment of UMIs below, independent of this.

    let mut conps = conp.clone();
    for i in 0..conps.len() {
        unique_sort(&mut conps[i]);
    }
    let mut umi_mcount_contig = Vec::<(i32, i32, i32)>::new();
    for j in 0..con.len() {
        let mut umi = Vec::<i32>::new();
        for e in 0..x.h.g.edge_count() {
            for i in 0..x.ids[e].len() {
                if bin_member(&conps[j], &(e as i32)) {
                    umi.push(umi_id[x.ids[e][i] as usize]);
                }
            }
        }
        umi.sort_unstable();
        let mut i = 0;
        while i < umi.len() {
            let k = next_diff(&umi, i);
            umi_mcount_contig.push((umi[i], i as i32 - k as i32, j as i32));
            i = k;
        }
    }
    umi_mcount_contig.sort_unstable();
    let mut umis = vec![0; con.len()];
    let mut i = 0;
    while i < umi_mcount_contig.len() {
        let j = next_diff1_3(&umi_mcount_contig, i);
        for k in i..j {
            if umi_mcount_contig[k].1 > umi_mcount_contig[i].1 {
                break;
            }
            umis[umi_mcount_contig[k].2 as usize] += 1;
        }
        i = j;
    }

    /*
    // Print:
    for i in 0..con.len() {
        println!( ">tig{} = {}", i, conp[i].iter().format(",") );
        println!( "{}", con[i].to_string() );
    }
    */

    // Compute V segment mismatches.

    let mut vmis = vec![1_000_000; con.len()];
    if !heur.free {
        for u in 0..con.len() {
            for m in 0..ann_all[u].len() {
                let t = ann_all[u][m].ref_id as usize;
                const EXCLUDE: usize = 15;
                let r = &refdata.refs[t];
                let start = ann_all[u][m].tig_start as usize;
                if refdata.is_v(t) && ann_all[u][m].ref_start == 0 && r.len() >= EXCLUDE {
                    let mut mis = 0;
                    for z in 0..r.len() - EXCLUDE {
                        if start + z >= con[u].len() || con[u].get(start + z) != r.get(z) {
                            mis += 1;
                        }
                    }
                    vmis[u] = mis;
                    break;
                }
            }
        }
    }

    // Compute junction sequences, and gather other info.

    let mut jseqs = Vec::<(DnaString, usize, i32, usize)>::new();
    for j in 0..con.len() {
        let mut jseq = DnaString::new();
        if !heur.free {
            junction_seq(&con[j], refdata, &ann_all[j], &mut jseq);
        } else {
            let cdr3 = get_cdr3(&con[j].slice(0, con[j].len()));
            if cdr3.solo() {
                jseq = con[j]
                    .slice(
                        cdr3[0].start_position_on_contig,
                        cdr3[0].start_position_on_contig + 3 * cdr3[0].aa_seq.len(),
                    )
                    .to_owned();
            }
        }
        jseqs.push((jseq, vmis[j], -umis[j], j));
    }

    // Define equivalence relation on junction sequences: they are equivalent if they differ
    // by <= 2 mismatches.

    const MAX_JUNCTION_MISMATCHES: usize = 2;
    let mut eq: EquivRel = EquivRel::new(jseqs.len() as u32);
    for i1 in 0..jseqs.len() {
        for i2 in i1 + 1..jseqs.len() {
            if eq.set_id(i1) != eq.set_id(i2) && jseqs[i1].0.len() == jseqs[i2].0.len() {
                let diffs = ndiffs(&jseqs[i1].0, &jseqs[i2].0);
                if diffs <= MAX_JUNCTION_MISMATCHES {
                    eq.join(i1, i2);
                }
            }
        }
    }

    // Rework jseqs, replacing junction sequence by its set id.

    let mut jseqs2 = Vec::<(usize, bool, usize, i32, usize, usize)>::new();
    for i in 0..jseqs.len() {
        let mut have_c = false;
        if !heur.free {
            for j in 0..ann_all[i].len() {
                let t = ann_all[i][j].ref_id as usize;
                if refdata.is_c(t) {
                    have_c = true;
                }
            }
        }
        jseqs2.push((
            eq.set_id(i),
            !have_c,
            jseqs[i].1,
            jseqs[i].2,
            jseqs[i].3,
            jseqs[i].0.len(),
        ));
    }
    jseqs2.sort_unstable();

    // Uniqify by junction: where we have multiple sequences sharing the same
    // junction, except for up to two mismatches, arbitrarily pick one that has the most UMIs.
    // MODIFIED: first look at V mismatches.
    // MODIFIED: first look at presence of C annotation
    //
    // ◼ The error rate in V segments might be reduced by exploiting global
    // ◼ kmer frequency here.
    //
    // ◼ In the denovo case, it seems like we should be using the cdr3 calculation
    // ◼ that uses annotation.
    let mut conx = Vec::<DnaString>::new();
    let mut conpx = Vec::<Vec<i32>>::new();
    let mut ann_allx = Vec::<Vec<Annotation>>::new();
    let mut i = 0;
    while i < jseqs2.len() {
        let j = next_diff1_6(&jseqs2, i);
        if jseqs2[i].5 > 0 {
            conx.push(con[jseqs2[i].4].clone());
            conpx.push(conp[jseqs2[i].4].clone());
            if !heur.free {
                ann_allx.push(ann_all[jseqs2[i].4].clone());
            }
        } else {
            for k in i..j {
                conx.push(con[jseqs2[k].4].clone());
                conpx.push(conp[jseqs2[k].4].clone());
                if !heur.free {
                    ann_allx.push(ann_all[jseqs2[k].4].clone());
                }
            }
        }
        i = j;
    }
    *con = conx;
    conp = conpx;
    ann_all = ann_allx;

    // Competitively delete TRAs and TRBs based on junction support.
    // Note that this is funkily generalized to IG.

    let mut tra = Vec::<ContigJunction>::new();
    let mut trb = Vec::<ContigJunction>::new();
    let mut to_delete: Vec<bool> = vec![false; con.len()];
    if heur.free {
        let mut all = Vec::<(i32, i32, usize)>::new();
        for j in 0..con.len() {
            let mut jseq_long = DnaString::new();
            let cdr3 = get_cdr3(&con[j].slice(0, con[j].len()));
            if cdr3.solo() {
                let jseq = con[j]
                    .slice(
                        cdr3[0].start_position_on_contig,
                        cdr3[0].start_position_on_contig + 3 * cdr3[0].aa_seq.len(),
                    )
                    .to_owned();
                if jseq.len() > 100 {
                    jseq_long = jseq.clone();
                } else {
                    let ext = 100 - jseq.len();
                    let start = if ext <= cdr3[0].start_position_on_contig {
                        cdr3[0].start_position_on_contig - ext
                    } else {
                        0
                    };
                    jseq_long = con[j]
                        .slice(
                            start,
                            cdr3[0].start_position_on_contig + 3 * cdr3[0].aa_seq.len(),
                        )
                        .to_owned();
                }
            }
            let mut jsupp_this: (i32, i32) = (0, 0);
            if !jseq_long.is_empty() {
                junction_supp_core(reads, x, umi_id, &jseq_long, &mut jsupp_this);
            }
            jsupp.push(jsupp_this);
            all.push((jsupp_this.0, jsupp_this.1, j));
        }
        reverse_sort(&mut all);
        for j1 in 0..all.len() {
            for j2 in 2..all.len() {
                let t2 = all[j2].2;
                let (u1, u2) = (all[j1].0, all[j2].0);
                let (n1, n2) = (all[j1].1, all[j2].1);
                let mut del = false;
                // ◼ ugly because the same thresholds are repeated below
                if u1 >= u2 && n1 >= 10 * n2 {
                    del = true;
                }
                if u1 >= 2 * u2 && u2 == 1 && n1 >= 2 * n2 {
                    del = true;
                }
                if u1 >= 3 * u2 && u2 == 1 && n1 >= n2 {
                    del = true;
                }
                if u1 >= 5 * u2 && u2 <= 2 && n1 >= 2 * n2 {
                    del = true;
                }
                if del {
                    to_delete[t2] = true;
                }
            }
        }
    } else {
        for j in 0..con.len() {
            let mut jsupp_this: (i32, i32) = (0, 0);
            junction_supp(
                &con[j],
                reads,
                x,
                umi_id,
                refdata,
                &ann_all[j],
                &mut jsupp_this,
            );
            jsupp.push(jsupp_this);
            let mut jseq = DnaString::new();
            junction_seq(&con[j], refdata, &ann_all[j], &mut jseq);
            let ann = annotate_seq(&con[j], refdata, true, false, false);
            let mut this_is_tra = false;
            for a in ann {
                // This definition of is_tra makes absolutely no sense, since
                // IG heavy chains have a D segment.  But it may be harmless.
                // Heavy chains do seem to have less UMI support than light chains.
                if rheaders[a.ref_id as usize].contains("TRAJ")
                    || (rheaders[a.ref_id as usize].contains("TRGJ") && gd_mode)
                    || rheaders[a.ref_id as usize].contains("IGHJ")
                {
                    this_is_tra = true;
                }
            }
            if this_is_tra {
                tra.push(ContigJunction {
                    jsupp_umi: jsupp[j].0,
                    jsupp_reads: jsupp[j].1,
                    con_idx: j,
                    junction_seq: jseq,
                    contig_seq: con[j].clone(),
                });
            } else {
                trb.push(ContigJunction {
                    jsupp_umi: jsupp[j].0,
                    jsupp_reads: jsupp[j].1,
                    con_idx: j,
                    junction_seq: jseq,
                    contig_seq: con[j].clone(),
                });
            }
        }
        reverse_sort(&mut tra);
        reverse_sort(&mut trb);
        for pass in 0..2 {
            let mut tr = &mut tra;
            if pass == 1 {
                tr = &mut trb;
            };
            let mut to_deletet: Vec<bool> = vec![false; tr.len()];
            for j in 1..tr.len() {
                let (u1, u2) = (tr[0].jsupp_umi, tr[j].jsupp_umi);
                let (n1, n2) = (tr[0].jsupp_reads, tr[j].jsupp_reads);
                let t = tr[j].con_idx;
                let mut del = false;
                if u1 >= u2 && n1 >= 10 * n2 {
                    del = true;
                }
                if u1 >= 2 * u2 && u2 == 1 && n1 >= 2 * n2 {
                    del = true;
                }
                if u1 >= 3 * u2 && u2 == 1 && n1 >= n2 {
                    del = true;
                }
                if u1 >= 5 * u2 && u2 <= 2 && n1 >= 2 * n2 {
                    del = true;
                }
                if del {
                    to_delete[t] = true;
                    to_deletet[j] = true;
                }
            }
            erase_if(tr, &to_deletet);
        }
    }
    let mut tocon: Vec<i32> = vec![-1; con.len()];
    let mut pos = 0;
    for i in 0..con.len() {
        if to_delete[i] {
            continue;
        }
        tocon[i] = pos;
        pos += 1;
    }
    erase_if(con, &to_delete);
    erase_if(&mut conp, &to_delete);
    erase_if(jsupp, &to_delete);
    erase_if(&mut ann_all, &to_delete);

    // Print paths.

    if log_opts.paths {
        fwriteln!(log, "\nFINAL GOOD CONTIG PATHS");
        for j in 0..conp.len() {
            fwriteln!(log, "{} = {}", j + 1, conp[j].iter().format(","));
        }
    }

    // If nearly all the kmers in one reject contig are contained in another longer
    // reject contig, delete the shorter one.

    let mut zkmers_plus = Vec::<(Kmer20, i32, i32)>::new();
    make_kmer_lookup_20_single(con2, &mut zkmers_plus);
    let mut to_delete2 = vec![false; con2.len()];
    let mut share = Vec::<Vec<i32>>::new();
    for _i in 0..con2.len() {
        let x: Vec<i32> = vec![0; con2.len()];
        share.push(x);
    }
    let mut i = 0;
    while i < zkmers_plus.len() {
        let j = next_diff1_3(&zkmers_plus, i);
        for k1 in i..j {
            for k2 in k1 + 1..j {
                let (t1, t2) = (zkmers_plus[k1].1, zkmers_plus[k2].1);
                if t1 != t2 {
                    share[t1 as usize][t2 as usize] += 1;
                    share[t2 as usize][t1 as usize] += 1;
                }
            }
        }
        i = j;
    }
    const MIN_COV: f64 = 0.75;
    for t1 in 0..con2.len() {
        for t2 in 0..con2.len() {
            let nkmers2 = (con2[t2].len() - x.h.k as usize + 1) as f64;
            if share[t1][t2] as f64 >= MIN_COV * nkmers2
                && (con2[t1].len() > con2[t2].len()
                    || (con2[t1].len() == con2[t2].len() && t1 < t2))
            {
                to_delete2[t2] = true;
            }
        }
    }
    erase_if(con2, &to_delete2);

    // If nearly all the annotated kmers in one reject contig are contained in a
    // good contig, delete the reject contig.

    if !heur.free {
        let mut ann_kmers2 = vec![0; con2.len()];
        let mut ann2s = Vec::<Vec<bool>>::new();
        for i in 0..con2.len() {
            let ann = annotate_seq(&con2[i], refdata, false, false, false);
            let mut ann2 = vec![false; con2[i].len() - x.h.k as usize + 1];
            for j in 0..ann.len() {
                let (start, len) = (ann[j].tig_start as usize, ann[j].match_len as usize);
                if start + len >= x.h.k as usize {
                    let stop = start + len - x.h.k as usize + 1;
                    for l in start..stop {
                        ann2[l] = true;
                    }
                }
            }
            for j in 0..ann2.len() {
                if ann2[j] {
                    ann_kmers2[i] += 1;
                }
            }
            ann2s.push(ann2);
        }
        let mut wkmers_plus = Vec::<(Kmer20, i32, i32)>::new();
        let mut con12 = con.clone();
        for b in con2.as_slice() {
            con12.push(b.clone());
        }
        make_kmer_lookup_20_single(&con12, &mut wkmers_plus);
        let mut to_delete12 = vec![false; con2.len()];
        let mut share = vec![vec![0_i32; con12.len()]; con12.len()];
        let mut i = 0;
        while i < wkmers_plus.len() {
            let j = next_diff1_3(&wkmers_plus, i);
            for k1 in i..j {
                for k2 in k1 + 1..j {
                    let t1 = wkmers_plus[k1].1;
                    let t2 = wkmers_plus[k2].1 - con.len() as i32;
                    if t1 >= con.len() as i32 || t2 < 0 {
                        continue;
                    }
                    let p2 = wkmers_plus[k2].2 as usize;
                    if !ann2s[t2 as usize][p2] {
                        continue;
                    }
                    share[t1 as usize][t2 as usize] += 1;
                }
            }
            i = j;
        }
        for t1 in 0..con.len() {
            for t2 in 0..con2.len() {
                if share[t1][t2] as f64 >= MIN_COV * ann_kmers2[t2] as f64 {
                    to_delete12[t2] = true;
                }
            }
        }
        erase_if(con2, &to_delete12);
    }

    // Move contigs.  Report total number.
    // ◼ This is totally confusing.  We move bad contigs, but then in practice they
    // ◼ appear to be deleted later.  It's not obvious from the code why the are
    // ◼ guaranteed to be deleted later.

    for x in con2.as_slice() {
        con.push(x.clone());
        jsupp.push((0, 0));
    }
    barcode_data.ncontigs = con.len() as i32;

    // Find maximal perfect matches of the reads to the contigs.  Only reads
    // that are placed on the graph are used.

    let mut placed = vec![false; reads.len()];
    for e in 0..x.h.g.edge_count() {
        for id in &x.ids[e] {
            placed[*id as usize] = true;
        }
    }
    let mut perf = Vec::<PerfectMatch>::new();
    const K: usize = 20;
    let mut bkmers_plus = Vec::<(Kmer20, i32, i32)>::new();
    let mut conx = con.clone();
    for b in con2.as_slice() {
        conx.push(b.clone());
    }
    make_kmer_lookup_20_single(&conx, &mut bkmers_plus);
    for id in 0..reads.len() {
        if !placed[id] {
            continue;
        }
        let b = &reads[id];
        for l in 0..(b.len() - K + 1) {
            let x: Kmer20 = b.get_kmer(l);
            let low = lower_bound1_3(&bkmers_plus, &x) as usize;
            for r in low..bkmers_plus.len() {
                if bkmers_plus[r].0 != x {
                    break;
                }
                let t = bkmers_plus[r].1 as usize; // contig that the kmer matches to
                let p = bkmers_plus[r].2 as usize; // start position of match on the contig
                if l > 0 && p > 0 && b.get(l - 1) == conx[t].get(p - 1) {
                    continue;
                }
                let mut len = K;
                while l + len < b.len() && p + len < conx[t].len() {
                    if b.get(l + len) != conx[t].get(p + len) {
                        break;
                    }
                    len += 1;
                }
                perf.push(PerfectMatch {
                    contig_id: t as i32,
                    read_id: id as i32,
                    match_len: -(len as i32),
                    start_pos_read: l as i32,
                    start_pos_contig: p as i32,
                });
            }
        }
    }
    perf.sort_unstable();

    // Assign UMIs to contigs.  We assign each UMI to at most one contig.
    // We let cumi[t] be the list of UMIs assigned to contig t.

    cumi.clear();
    for _i in 0..conx.len() {
        cumi.push(Vec::<i32>::new());
    }
    let mut numi = 0;
    for u in umi_id {
        numi = max(numi, u + 1);
    }
    let mut go = Vec::<Vec<i32>>::new();
    for _u in 0..numi {
        go.push(vec![0; conx.len()]);
    }
    for x in &perf {
        go[umi_id[x.read_id as usize] as usize][x.contig_id as usize] += -x.match_len;
    }
    for u in 0..numi {
        let mut y = Vec::<(i32, i32)>::new();
        for j in 0..conx.len() {
            y.push((go[u as usize][j], j as i32));
        }
        reverse_sort(&mut y);
        if !conx.is_empty() && y[0].0 > 0 {
            cumi[y[0].1 as usize].push(u);
        }
    }
    for t in 0..conx.len() {
        unique_sort(&mut cumi[t]);
    }

    // Set up to track validated and non-validated and invalidated UMIs.

    *validated_umis = vec![Vec::<String>::new(); conx.len()];
    *non_validated_umis = vec![Vec::<String>::new(); conx.len()];
    *invalidated_umis = vec![Vec::<String>::new(); conx.len()];

    // Form pileup and then compute quality scores.  Also compute cids.  Note that the pileup
    // uses only bases that agree with the contig, which does not necessarily make sense.  It also
    // does not respect the assignment of reads to contigs, which does not necessarily make sense.
    //
    // pilex: This is a different pileup:
    // - We only pile a read on a contig if its UMI is assigned to that contig.
    // - The offset of a read is determined by its longest match to the contig.
    // - Every base on the read (that is not off the end) is assigned to the contig.

    cids.clear();
    for _i in 0..conx.len() {
        cids.push(Vec::<i32>::new());
    }
    let mut perfp = 0;
    let mut contig_count = 0;
    for t in 0..conx.len() {
        let mut pile = Vec::<Vec<(UmiType, u8, u8)>>::new(); // {(umi, read base, qual score)}
        for _ in 0..conx[t].len() {
            pile.push(Vec::<(UmiType, u8, u8)>::new());
        }
        let mut pilex = Vec::<Vec<(UmiType, u8, u8)>>::new(); // {(umi, read base, qual score)}
        for _ in 0..conx[t].len() {
            pilex.push(Vec::<(UmiType, u8, u8)>::new());
        }
        let mut assigned = HashSet::<i32>::new();
        while perfp < perf.len() {
            let x = &perf[perfp];
            let t2 = x.contig_id;
            if t2 > t as i32 {
                break;
            }
            let id = x.read_id;
            let len = -x.match_len;
            let l = x.start_pos_read;
            let p = x.start_pos_contig;
            let u = umi_id[id as usize];
            if bin_member(&cumi[t], &u) {
                cids[t].push(id);
                if !assigned.contains(&id) {
                    assigned.insert(id);
                    // position p on contig maps to position l on read
                    let r = &reads[id as usize];
                    for m in 0..r.len() {
                        let cpos = m as i32 - l + p;
                        if cpos >= 0 && cpos < conx[t].len() as i32 {
                            pilex[cpos as usize].push((
                                u as UmiType,
                                r.get(m),
                                quals[id as usize][m],
                            ));
                        }
                    }
                }

                // Also save the partner of the read.

                if !single_end {
                    if id % 2 == 0 {
                        cids[t].push(id + 1);
                    } else {
                        cids[t].push(id - 1);
                    }
                }
            }
            for j in 0..len {
                pile[(p + j) as usize].push((
                    u as UmiType,
                    reads[id as usize].get((l + j) as usize),
                    quals[id as usize][(l + j) as usize],
                ));
            }
            perfp += 1;
        }
        for t in 0..conx.len() {
            unique_sort(&mut cids[t]);
        }

        // Validate UMIs.  For each productive contig, we create a list of UMIs which appear
        // to have complete and correct coverage from the beginning of the V segment to the end
        // of the J segment.

        if !cids[t].is_empty() && t < con.len() {
            fwriteln!(
                log,
                "\nUMI info for barcode {} contig {} = {}...",
                barcode_data.barcode,
                contig_count + 1,
                conx[t].slice(0, 10).to_string(),
            );
            contig_count += 1;
            let ann = annotate_seq(&conx[t], refdata, true, false, true);
            let mut vstart = None;
            let mut jstop = None;
            for i in 0..ann.len() {
                let u = ann[i].ref_id;
                if refdata.is_v(u as usize) {
                    vstart = Some(ann[i].tig_start);
                    break;
                }
            }
            for i in (0..ann.len()).rev() {
                let u = ann[i].ref_id;
                if refdata.is_j(u as usize) {
                    jstop = Some(ann[i].tig_start + ann[i].match_len);
                    break;
                }
            }
            if vstart.is_none() {
                fwriteln!(log, "vstart = None"); // should not be possible, and will cause crash
            }
            if jstop.is_none() {
                fwriteln!(log, "jstop = None"); // should not be possible, and will cause crash
            }
            let mut umis = Vec::<UmiType>::new();
            for i in 0..cids[t].len() {
                umis.push(umi_id[cids[t][i] as usize] as u32);
            }
            unique_sort(&mut umis);
            let mut count = vec![0; umis.len()]; // read count for each UMI
            for i in 0..cids[t].len() {
                let u = umi_id[cids[t][i] as usize] as u32;
                let p = bin_position(&umis, &u) as usize;
                count[p] += 1;
            }
            if let Some(vstart) = vstart {
                if let Some(jstop) = jstop {
                    if vstart < jstop {
                        for i in 0..umis.len() {
                            fwrite!(log, "umi {} = {} reads:", uu[umis[i] as usize], count[i]);
                            let vstart = vstart as usize;
                            let vj = &conx[t].slice(vstart, jstop as usize);
                            let mut good = vec![false; vj.len()];
                            let mut bad = vec![false; vj.len()];
                            let mut verybad = vec![false; vj.len()];
                            let mut val = true;
                            let mut inval = false;
                            for p in 0..vj.len() {
                                let mut qsums = vec![0_usize; 4];
                                for j in 0..pilex[vstart + p].len() {
                                    if pilex[vstart + p][j].0 == umis[i] {
                                        qsums[pilex[vstart + p][j].1 as usize] +=
                                            pilex[vstart + p][j].2 as usize;
                                    }
                                }
                                let mut ids = vec![0, 1, 2, 3];
                                sort_sync2(&mut qsums, &mut ids);
                                qsums.reverse();
                                ids.reverse();
                                let q = qsums[0] - qsums[1];
                                let b = ids[0];
                                if b == conx[t].get(vstart + p) && q >= 20 {
                                    good[p] = true;
                                } else {
                                    val = false;
                                }
                                if b != conx[t].get(vstart + p) && qsums[0] > 0 {
                                    bad[p] = true;
                                    inval = true;
                                }
                                if b != conx[t].get(vstart + p) && q >= 40 {
                                    verybad[p] = true;
                                }
                            }
                            let mut j = 0;
                            while j < good.len() {
                                let k = next_diff(&good, j);
                                if good[j] {
                                    fwrite!(log, " +");
                                } else {
                                    fwrite!(log, " -");
                                }
                                fwrite!(log, "{}", k - j);
                                if verybad[j] {
                                    fwrite!(log, "XX");
                                } else {
                                    let mut b = false;
                                    for m in j..k {
                                        if bad[m] {
                                            b = true;
                                        }
                                    }
                                    if b {
                                        fwrite!(log, "X");
                                    }
                                }
                                j = k;
                            }
                            if val {
                                validated_umis[t].push(uu[umis[i] as usize].clone());
                                fwrite!(log, " validated");
                            } else if !inval {
                                non_validated_umis[t].push(uu[umis[i] as usize].clone());
                                fwrite!(log, " non-validated");
                            } else {
                                invalidated_umis[t].push(uu[umis[i] as usize].clone());
                                fwrite!(log, " invalidated");
                            }
                            fwriteln!(log, "");
                        }
                    }
                }
            }
            const UMI_CAP: usize = 20;
            if validated_umis[t].len() > UMI_CAP {
                validated_umis[t] = random_subvec(&validated_umis[t], UMI_CAP);
            }
            if non_validated_umis[t].len() > UMI_CAP {
                non_validated_umis[t] = random_subvec(&non_validated_umis[t], UMI_CAP);
            }
            if invalidated_umis[t].len() > UMI_CAP {
                invalidated_umis[t] = random_subvec(&invalidated_umis[t], UMI_CAP);
            }
        }

        // Continue with computing quality scores.

        let rt_err = 1e-4;
        let mut x = Vec::<u8>::new();
        // let mut hards = 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        for j in 0..conx[t].len() {
            // Compute the quality score using an expensive method.  However we
            // use some heuristic shortcuts first that yield approximately the same
            // values, and are much faster.
            // ◼ It would be good to put this on a firmer basis, and perhaps
            // also have more shortcuts.

            let mut count = [0; 4];
            for k in 0..pile[j].len() {
                count[pile[j][k].1 as usize] += 1;
            }
            let mut calls = 0;
            for k in 0..4 {
                if count[k] > 0 {
                    calls += 1;
                }
            }
            let mut nmax = 0;
            let mut numax = 0;
            let mut numi = 0;
            let total = pile[j].len();
            let mut maxq = 0;
            let mut k = 0;
            while k < pile[j].len() {
                let mut top = 0;
                let mut l = k;
                while l < pile[j].len() {
                    if pile[j][l].0 != pile[j][k].0 {
                        break;
                    }
                    let q = pile[j][l].2;
                    top = max(top, q);
                    if q >= 32 {
                        nmax += 1;
                    }
                    maxq = max(maxq, q);
                    l += 1;
                }
                numi += 1;
                if top >= 32 {
                    numax += 1;
                }
                k = l;
            }

            /*
            let q4 = pos_base_quals_helper( &pile[j], rt_err, utils::logaddexp_arr );
            let qval = q4[ conx[t].get(j) as usize ];
            */

            if calls == 1 && numax >= 2 {
                x.push(60);
            } else if calls == 1 && numi == 1 && nmax >= 2 {
                x.push(40);
            }
            // weird but this is what happens:
            else if total == 1 && maxq >= 22 {
                x.push(29);
            }
            // Do the expensive calculation.
            else {
                /*
                // if calls == 1 && nu37 >= 2 && qval != 60 {
                    println!( "" );
                    for k in 0..pile[j].len() {
                        print!( "{},{},{} ",
                            pile[j][k].0, pile[j][k].1, pile[j][k].2 );
                    }
                    println!( "\n\n{:?}", count );
                    println!( "calls = {}", calls );
                    println!( "numax = {}", numax );
                    println!( "numi = {}\n", numi );
                    println!( "{}", qval );
                // }
                */

                // hards += 1; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                let q4 = pos_base_quals_helper(&pile[j], rt_err, utils::logaddexp_arr);
                x.push(q4[conx[t].get(j) as usize]);
            }
        }
        // println!( "{} / {}", hards, conx[t].len() ); // XXXXXXXXXXXXXXXXXXXXXXXX
        conxq.push(x);
    }

    // Delete contigs that have no reads assigned to them.  This completes the
    // creation of contigs.
    // ◼ This is likely doing unnecessary damage.  The underlying problem may be
    // ◼ that we have not sufficiently reduced contig redundancy.

    let mut dels1 = Vec::<i32>::new();
    let mut noreads = vec![false; cids.len()];
    let mut noreads1 = vec![false; con.len()];
    let mut noreads2 = vec![false; con2.len()];
    for i in 0..cids.len() {
        if cids[i].is_empty() {
            noreads[i] = true;
            if i < con.len() {
                noreads1[i] = true;
                dels1.push(i as i32);
            } else {
                noreads2[i - con.len()] = true;
            }
        }
    }
    let mut count1 = 0;
    let mut to_new1: Vec<i32> = vec![-1; con.len()];
    for i in 0..con.len() {
        if noreads1[i] {
            continue;
        }
        to_new1[i] = count1;
        count1 += 1;
    }
    let mut tra_del = vec![false; tra.len()];
    for i in 0..tra.len() {
        if noreads1[tocon[tra[i].con_idx] as usize] {
            tra_del[i] = true;
        }
    }
    erase_if(&mut tra, &tra_del);
    let mut trb_del = vec![false; trb.len()];
    for i in 0..trb.len() {
        if noreads1[tocon[trb[i].con_idx] as usize] {
            trb_del[i] = true;
        }
    }
    erase_if(&mut trb, &trb_del);
    erase_if(con, &noreads1);
    erase_if(jsupp, &noreads1);
    erase_if(con2, &noreads2);
    erase_if(conxq, &noreads);
    erase_if(cumi, &noreads);
    erase_if(cids, &noreads);
    erase_if(validated_umis, &noreads);
    erase_if(non_validated_umis, &noreads);
    erase_if(invalidated_umis, &noreads);
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct PerfectMatch {
    contig_id: i32,
    read_id: i32,
    match_len: i32,
    start_pos_read: i32,
    start_pos_contig: i32,
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct ContigJunction {
    jsupp_umi: i32,
    jsupp_reads: i32,
    con_idx: usize,
    junction_seq: DnaString,
    contig_seq: DnaString,
}
