// This file contains code to process a barcode.

#![allow(clippy::many_single_char_names)]

use crate::barcode_data::{BarcodeData, BarcodeDataBrief};
use crate::constants::{CHAIN_TYPES, PRIMER_EXT_LEN};
use crate::contigs::make_contigs;
use crate::heuristics::Heuristics;
use crate::hops::FlowcellContam;
use crate::log_opts::LogOpts;
use crate::print_hyper::print_with_annotations;
use crate::ref_free::{simplify_without_ref, uber_strong_paths, umis2};
use align_tools::{affine_align, complexity};
use amino::{have_start, have_stop};
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer12, Kmer20};
use debruijn::{Kmer, Mer, Vmer};
use graph_simple::GraphSimple;
use hyperbase::Hyper;
use itertools::{izip, Itertools};
use kmer_lookup::{make_kmer_lookup_20_single, match_12};
use stats_utils::percent_ratio;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::io::prelude::*;
use std::iter::zip;
use std::mem::swap;
use std::time::Instant;
use string_utils::{abbrev_list, stringme, TextUtils};
use superslice::Ext;
use tenkit2::io::print_compressed;
use tenkit2::pack_dna::{pack_bases_16x, reverse_complement, unpack_bases_10};
use vdj_ann::annotate::{
    annotate_seq, annotate_seq_core, chain_type, get_cdr3, get_cdr3_using_ann, print_annotations,
    print_cdr3, print_cdr3_using_ann, print_some_annotations, print_start_codon_positions,
    ContigAnnotation, JunctionSupport,
};
use vdj_ann::refx::RefData;
use vdj_ann::transcript::is_valid;
use vector_utils::{
    bin_member, intersection, lower_bound1_3, meet_size, next_diff, position, unique_sort,
    upper_bound1_3,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// REFERENCE-ASSISTED GRAPH CLEANING
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Pop bubbles using reference.

fn pop_bubbles_using_reference(x: &mut Hyper, umi_id: &[i32], refdata: &RefData) {
    let refs = &refdata.refs;
    let rkmers_plus = &refdata.rkmers_plus;
    let mut dels = Vec::<u32>::new();
    for v in 0..x.h.g.node_count() {
        if x.h.g.n_to(v) != 1 || x.h.g.n_from(v) != 2 {
            continue;
        }
        let w = x.h.g.v_from(v, 0);
        if x.h.g.v_from(v, 1) != w || x.h.g.n_to(w) != 2 || x.h.g.n_from(w) != 1 {
            continue;
        }
        let mut e1 = x.h.g.e_from(v, 0);
        let mut e2 = x.h.g.e_from(v, 1);
        for pass in 0..2 {
            if pass == 1 {
                swap(&mut e1, &mut e2);
            }

            // Compute properties of bubble branches.

            let (mut u1, mut u2) = (Vec::<i32>::new(), Vec::<i32>::new());
            umis2(x, umi_id, e1 as i32, &mut u1);
            umis2(x, umi_id, e2 as i32, &mut u2);
            unique_sort(&mut u1);
            unique_sort(&mut u2);
            let (nu1, nu2) = (u1.len(), u2.len());
            let mut ann1 = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(
                x.h.g.edge_obj(e1 as u32),
                refdata,
                &mut ann1,
                false,
                false,
                false,
            );
            let mut ann2 = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(
                x.h.g.edge_obj(e2 as u32),
                refdata,
                &mut ann2,
                false,
                false,
                false,
            );
            let (re1, re2) = (x.inv[e1], x.inv[e2]);
            let mut ann1rc = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(
                x.h.g.edge_obj(re1),
                refdata,
                &mut ann1rc,
                false,
                false,
                false,
            );
            let mut ann2rc = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(
                x.h.g.edge_obj(re2),
                refdata,
                &mut ann2rc,
                false,
                false,
                false,
            );
            let a1 = ann1.len() + ann1rc.len() > 0;
            let a2 = ann2.len() + ann2rc.len() > 0;
            let k1 = x.h.kmers(e1 as u32) as i32;
            let k2 = x.h.kmers(e2 as u32) as i32;
            let (s1, s2) = (x.supp(e1), x.supp(e2));
            let k = x.h.k;

            // Carry out tests.

            let mut ok = false;
            if k1 == k2 {
                if k1 <= 2 * k && nu1 >= nu2 && a1 && !a2 {
                    ok = true;
                }
                if nu1 >= 3 && nu2 == 1 && a1 {
                    ok = true;
                }
                if k1 == k && nu1 >= nu2 && s1 > s2 && !a2 {
                    ok = true;
                }
            }

            // One more test.  If the bubble is a SNP, and one branch is in the
            // reference but the other is not, and there are at least as many UMIs
            // on the 'reference' branch, delete the other one.
            //
            // And a final test after this: at a SNP if one branch has more UMIs,
            // it wins.  *** THIS IS OFF ***

            if k1 == k2 && k1 == k && nu1 >= nu2 {
                let mut diffs = 0;
                for i in 0..(k as usize) {
                    if x.h.g.edge_obj(e1 as u32).get(i) != x.h.g.edge_obj(e2 as u32).get(i) {
                        diffs += 1;
                    }
                }
                if diffs == 1 {
                    let m1 = match_12(x.h.g.edge_obj(e1 as u32), refs, rkmers_plus);
                    let m2 = match_12(x.h.g.edge_obj(e2 as u32), refs, rkmers_plus);
                    if m1 && !m2 {
                        ok = true;
                    }
                    // if nu1 > nu2 { ok = true; }
                }
            }
            if ok {
                dels.push(e2 as u32);
                break;
            }
        }
    }
    x.kill_edges(&dels);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// FIND BEST ALIGNMENT OF A V SEGMENT TO A CONTIG
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Requirements for INTERESTING:
//
// 1. Amongst V segments of length >= 300 bases, find one with a gap-free
//    alignment to tig, starting with a start codon on tig, and having no stop
//    codon, for which the percent identity is highest.
//
// 2. Percent identity must be between 50 and 85%.
//
// 3. Maximum perfect match to any reference sequence must be < 60 bases.

#[allow(clippy::too_many_arguments)]
fn analyze_vscore(
    label: &str,
    tig: &DnaString,
    refs: &[DnaString],
    rheaders: &[String],
    rkmers_plus: &[(Kmer12, i32, i32)],
    is_tcr: bool,
    is_bcr: bool,
    log: &mut Vec<u8>,
) {
    let mut best = -1_f64;
    let mut best_loc = 0;
    let mut tbest = 0;
    // let mut comments = String::new();
    for t in 0..refs.len() {
        if refs[t].len() < 300 {
            continue;
        }
        if tig.len() < refs[t].len() {
            continue;
        }
        if rheaders[t].contains("UTR") {
            continue;
        }
        let mut ok = false;
        if is_tcr
            && (rheaders[t].contains("TRAV")
                || rheaders[t].contains("TRBV")
                || rheaders[t].contains("TRDV"))
        {
            ok = true;
        }
        if is_bcr
            && (rheaders[t].contains("IGHV")
                || rheaders[t].contains("IGKV")
                || rheaders[t].contains("IGLV"))
        {
            ok = true;
        }
        if !ok {
            continue;
        }
        for i in 0..tig.len() - refs[t].len() + 1 {
            if !have_start(tig, i) {
                continue;
            }
            let mut matches = 0;
            for j in 0..refs[t].len() {
                if refs[t].get(j) == tig.get(i + j) {
                    matches += 1;
                }
            }
            let score = percent_ratio(matches, refs[t].len());
            if score > best {
                let mut stop = false;
                for j in 0..(refs[t].len() - 3) / 3 {
                    if have_stop(tig, i + 3 * j) {
                        stop = true;
                    }
                }
                if !stop {
                    best = score;
                    best_loc = i;
                    tbest = t;
                }
            }
        }
    }
    if best < 0 as f64 {
        return;
    }
    let mut mpm = 0; // max perfect match for all V segments
    let b = &tig.slice(best_loc, best_loc + refs[tbest].len());
    const K: usize = 12;
    for l in 0..(b.len() - K + 1) {
        let x: Kmer12 = b.get_kmer(l);
        let low = lower_bound1_3(rkmers_plus, &x);
        let high = upper_bound1_3(rkmers_plus, &x);
        for r in low..high {
            let t = rkmers_plus[r as usize].1 as usize;
            let p = rkmers_plus[r as usize].2 as usize;
            if l > 0 && p > 0 && b.get(l - 1) == refs[t].get(p - 1) {
                continue;
            }
            let mut len = K;
            while l + len < b.len() && p + len < refs[t].len() {
                if b.get(l + len) != refs[t].get(p + len) {
                    break;
                }
                len += 1;
            }
            mpm = max(mpm, len);
        }
    }
    let mut interesting = false;
    if (50_f64..=85_f64).contains(&best) && mpm < 60 {
        // Double check.  Is there a V segment that starts at the same place and
        // has passing identity?

        let mut better = false;
        for t in 0..refs.len() {
            let mut ok = false;
            if is_tcr && (rheaders[t].contains("TRAV") || rheaders[t].contains("TRBV")) {
                ok = true;
            }
            if is_bcr
                && (rheaders[t].contains("IGHV")
                    || rheaders[t].contains("IGKV")
                    || rheaders[t].contains("IGLV"))
            {
                ok = true;
            }
            if !ok {
                continue;
            }
            let mut mat = 0;
            let mut total = 0;
            for j in 0..refs[t].len() {
                if best_loc + j == tig.len() {
                    break;
                }
                total += 1;
                if refs[t].get(j) == tig.get(best_loc + j) {
                    mat += 1;
                }
            }
            let alt = percent_ratio(mat, total);
            if alt > 85_f64 {
                better = true;
            }
        }
        if !better {
            // comments = ", INTERESTING".to_string();
            interesting = true;
        }
    }
    // fwriteln!( log, "\nVSCORE = {:2.1}%, max perf match = {}{}",
    //     best, mpm, comments );
    if interesting {
        fwriteln!(log, ">vscore_{}_{:2.1}%\n{}", label, best, b.to_string());
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PRINT CONTIG INFO
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[allow(clippy::too_many_arguments)]
fn print_contig_info(
    tigname: &str,
    label: &str,
    log: &mut Vec<u8>,
    good: bool,
    id: usize,
    tig: &DnaString,
    ann: &[(i32, i32, i32, i32, i32)],
    jsupp: &(i32, i32),
    validated_umis: &[String],
    non_validated_umis: &[String],
    invalidated_umis: &[String],
    tigq: &[u8],
    cumi: &[i32],
    cids: &[i32],
    barcode_data: &BarcodeData,
    is_tcr: bool,
    is_bcr: bool,
    is_gd: Option<bool>,
    refdata: &RefData,
    refdatax: &RefData,
    free: bool,
    log_opts: &LogOpts,
) {
    let gd_mode = is_gd.unwrap_or(false);
    let refs = &refdata.refs;
    let rheaders = &refdata.rheaders;
    let rkmers_plus = &refdata.rkmers_plus;
    fwriteln!(log, "\nTIG {}[bases={}]", id + 1, tig.len());
    if log_opts.print_seq {
        fwriteln!(log, "{}", tig.to_string());
    }
    if log_opts.print_qual {
        print_compressed(log, tigq);
    }
    let verbose = false;
    if good {
        print_some_annotations(refdata, ann, log, verbose);
    } else {
        print_annotations(tig, refdatax, log, false, true, verbose);
    }
    let junction_support = match !free && good {
        true => Some(JunctionSupport {
            reads: jsupp.1,
            umis: jsupp.0,
        }),
        false => None,
    };

    // Print CDR3 locations.

    let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
    if !free {
        print_cdr3_using_ann(tig, refdata, ann, log);
        get_cdr3_using_ann(tig, refdata, ann, &mut cdr3);
    } else {
        print_cdr3(&tig.slice(0, tig.len()), log);
        get_cdr3(&tig.slice(0, tig.len()), &mut cdr3);
    }
    if good && cdr3.is_empty() {
        fwriteln!(log, "warning: contig labeled good but has no cdr3");
    }

    // Report umi and read assignments.

    let cumi_short = cumi.iter().take(10).copied().collect::<Vec<_>>();
    fwrite!(log, "umis assigned: {:?}", cumi_short);
    if cumi.len() > cumi_short.len() {
        fwrite!(log, " and {} others", cumi.len() - cumi_short.len());
    }
    let mut z = Vec::<i32>::new();
    let mut a = cumi.to_vec();
    let mut b = barcode_data.xuids.clone();
    a.sort_unstable();
    b.sort_unstable();
    intersection(&a, &b, &mut z);
    fwriteln!(log, "\nof which {} are surviving nonsolos", z.len());
    fwriteln!(log, "reads assigned: {:?}", cids.len());

    // Gather some annotation info.

    let (mut vs, mut js) = (Vec::<usize>::new(), Vec::<usize>::new());
    for i in ann {
        let t = i.2 as usize;
        let h = &rheaders[t];
        if h.contains("UTR") {
            continue;
        }
        if h.contains("TRAV")
            || h.contains("TRBV")
            || (h.contains("TRGV") && gd_mode)
            || (h.contains("TRDV") && gd_mode)
            || h.contains("IGHV")
            || h.contains("IGKV")
            || h.contains("IGLV")
        {
            vs.push(t);
        }
        if h.contains("TRAJ")
            || h.contains("TRBJ")
            || (h.contains("TRGJ") && gd_mode)
            || (h.contains("TRDJ") && gd_mode)
            || h.contains("IGHJ")
            || h.contains("IGKJ")
            || h.contains("IGLJ")
        {
            js.push(t);
        }
    }

    // Print start codon positions.

    print_start_codon_positions(tig, log);

    // Report certain split V annotations.

    let mut split_v = false;
    for i in 0..ann.len() {
        if i + 1 < ann.len()
            && rheaders[ann[i].2 as usize].contains("V-REGION")
            && ann[i].2 == ann[i + 1].2
        {
            if ann[i].0 + ann[i].1 == ann[i + 1].0
                && ann[i].3 + ann[i].1 < ann[i + 1].3
                && ann[i + 1].3 - ann[i].1 - ann[i].3 <= 27
                && (ann[i + 1].3 - ann[i].1 - ann[i].3) % 3 == 0
            {
                split_v = true;
                fwriteln!(
                    log,
                    "see deletion of {} bases at pos {} on {}",
                    ann[i + 1].3 - ann[i].1 - ann[i].3,
                    ann[i].1 + ann[i].3,
                    rheaders[ann[i].2 as usize]
                );
            }
            if ann[i].3 + ann[i].1 == ann[i + 1].3
                && ann[i].1 + ann[i].0 < ann[i + 1].0
                && ann[i + 1].0 - ann[i].0 - ann[i].1 <= 27
                && (ann[i + 1].0 - ann[i].0 - ann[i].1) % 3 == 0
            {
                split_v = true;
                fwriteln!(
                    log,
                    "see insertion of {} at pos {} on {}",
                    tig.slice((ann[i].0 + ann[i].1) as usize, ann[i + 1].0 as usize)
                        .to_string(),
                    ann[i + 1].3,
                    rheaders[ann[i].2 as usize]
                );
            }
        }
    }

    // Check for nonstandard stuff.  Exempt cases where we have a double V that
    // exhibits a 'standard' indel.

    let mut split_u = false;
    for i in 0..ann.len() {
        if i + 1 < ann.len()
            && rheaders[ann[i].2 as usize].contains("UTR")
            && ann[i].2 == ann[i + 1].2
        {
            if ann[i].0 + ann[i].1 == ann[i + 1].0
                && ann[i].3 + ann[i].1 < ann[i + 1].3
                && ann[i + 1].3 - ann[i].1 - ann[i].3 == 1
            {
                split_u = true;
            }
            if ann[i].3 + ann[i].1 == ann[i + 1].3
                && ann[i].1 + ann[i].0 < ann[i + 1].0
                && ann[i + 1].0 - ann[i].0 - ann[i].1 == 1
            {
                split_u = true;
            }
        }
    }
    if good {
        let mut funny = false;
        let mut u = Vec::<String>::new();
        let mut v = Vec::<String>::new();
        let mut d = Vec::<String>::new();
        let mut j = Vec::<String>::new();
        let mut c = Vec::<String>::new();
        for i in ann {
            let t = i.2 as usize;
            let name = rheaders[t].after("|").between("|", "|");
            if rheaders[t].contains("UTR") {
                u.push(name.to_string());
            }
            if rheaders[t].contains("V-REGION") {
                v.push(name.to_string());
            }
            if rheaders[t].contains("D-REGION") {
                d.push(name.to_string());
            }
            if rheaders[t].contains("J-REGION") {
                j.push(name.to_string());
            }
            if rheaders[t].contains("C-REGION") {
                c.push(name.to_string());
            }
        }
        if u.len() > 1 || v.len() > 1 || d.len() > 1 || j.len() > 1 || c.len() > 1 {
            funny = true;
        }
        if u.len() == 1 && v.len() == 1 && u[0] != v[0] {
            funny = true;
        }
        if u.len() <= 1 && v.len() == 2 && d.len() <= 1 && j.len() <= 1 && c.len() <= 1 && split_v {
            funny = false;
        }
        if u.len() == 2 && v.len() <= 1 && d.len() <= 1 && j.len() <= 1 && c.len() <= 1 && split_u {
            funny = false;
        }
        if funny {
            fwriteln!(log, "funny annotation");
        }
    }

    // Do some optional analyses.
    //
    // The vscore analysis is a bit extravagant.  It does the following:
    // For each contig, find the gap-free alignment between a V segment and the
    // contig having the highest percent identity, and report this percent; this
    // can be used to find putative novel V segments: apply this to output
    // | grep -A1 vscore_seq | grep -v ">" | grep -v '\-' | sort |
    // uniq -c | sort -n -r | less

    if log_opts.vis && vs.len() == 1 && js.len() == 1 {
        let mut r = refs[vs[0]].clone();
        let r2 = &refs[js[0]];
        for i in 0..r2.len() {
            r.push(r2.get(i));
        }
        let a = affine_align(&r, tig);
        fwriteln!(log, "\nalign complexity = {}", complexity(&a));
        fwriteln!(log, ">VJ\n{}", r.to_string());
        fwriteln!(log, ">tig\n{}", tig.to_string());
    }
    analyze_vscore(label, tig, refs, rheaders, rkmers_plus, is_tcr, is_bcr, log);

    // Print json contig annotations.

    if log_opts.json {
        let high_confidence = true; // might change
        let is_cell = true; // might change
        let bar = "===============================================================\
                   =====================";
        fwriteln!(log, "\n{}\n\nJSON ANNOTATION\n", bar);
        let json = ContigAnnotation::from_annotate_seq(
            tig,
            tigq,
            tigname,
            refdata,
            ann,
            cids.len(),
            cumi.len(),
            high_confidence,
            Some(validated_umis.to_vec()),
            Some(non_validated_umis.to_vec()),
            Some(invalidated_umis.to_vec()),
            is_cell,
            good,
            junction_support,
        );
        json.print(log);
        fwriteln!(log, "\n\n{}", bar);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PROCESS A BARCODE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Note that this trims the reads.

#[allow(clippy::too_many_arguments)]
pub fn process_barcode(
    // INPUTS:
    label: &str,
    single_end: bool,
    is_tcr: bool,
    is_bcr: bool,
    is_gd: Option<bool>,
    inner_primers: &[Vec<u8>],
    outer_primers: &[Vec<u8>],
    inner_primer_exts: &[Vec<Vec<u8>>],
    outer_primer_exts: &[Vec<Vec<u8>>],
    reads: &mut [DnaString],
    quals: &mut [Vec<u8>],
    umi_id: &[i32],
    uu: &[String],
    n50_n50_rpu: i32,
    refdata: &RefData,
    refdata_full: &RefData,
    rkmers_plus_full_20: &[(Kmer20, i32, i32)],
    refdatax: &RefData,
    lena: usize,
    contam: &FlowcellContam,
    // OUTPUTS:
    barcode_data: &mut BarcodeData,
    barcode_data_brief: &mut BarcodeDataBrief,
    conx: &mut Vec<DnaString>,
    conxq: &mut Vec<Vec<u8>>,
    productive: &mut Vec<bool>,
    validated_umis: &mut Vec<Vec<String>>,
    non_validated_umis: &mut Vec<Vec<String>>,
    invalidated_umis: &mut Vec<Vec<String>>,
    cids: &mut Vec<Vec<i32>>,
    cumi: &mut Vec<Vec<i32>>,
    nedges: &mut usize,
    junction_support: &mut Vec<Option<JunctionSupport>>,
    log: &mut Vec<u8>,
    // INPUTS:
    log_opts: &LogOpts,
    heur: &Heuristics,
) {
    let gd_mode = is_gd.unwrap_or(false);
    // Unpack refdata.

    let t = Instant::now();
    log_opts.report_perf_stats(log, &t, "upon entering process_barcode");
    let rheaders = &refdata.rheaders;
    let rheaders_full = &refdata_full.rheaders;

    // Get cross-flowcell contamination.

    let mut bad_umis = Vec::<String>::new();
    if !contam.entry.is_empty() {
        const MIN_RATIO_CONTAM: u32 = 10;
        let mut lid = 60000_u16;
        for i in 0..contam.lena_index.len() {
            if bin_member(&contam.lena_index[i], &lena) {
                lid = i as u16;
                break;
            }
        }
        if lid == 60000_u16 {
            panic!("Failed to find lena id in contam.lena_index.");
        }
        let bc = barcode_data_brief.barcode.before("-").as_bytes();
        let mut bc_packed = [0_u8; 4];
        pack_bases_16x(bc, &mut bc_packed);
        let low = contam.entry.lower_bound_by_key(&bc_packed, |a| a.bc);
        let high = contam.entry.upper_bound_by_key(&bc_packed, |a| a.bc);
        let mut i = low;
        while i < high {
            let mut j = i + 1;
            while j < high {
                if contam.entry[j].umi != contam.entry[i].umi {
                    break;
                }
                j += 1;
            }
            let mut max_mult = 0;
            for k in i..j {
                max_mult = max(max_mult, contam.entry[k].mult);
            }
            for k in i..j {
                if MIN_RATIO_CONTAM * contam.entry[k].mult <= max_mult && contam.entry[k].lid == lid
                {
                    let umi_packed = &contam.entry[k].umi;
                    let mut umi = [0_u8; 10];
                    unpack_bases_10(umi_packed, &mut umi);
                    bad_umis.push(stringme(&umi));
                }
            }
            i = j;
        }
        bad_umis.sort();
    }

    // Start log.

    fwriteln!(
        log,
        "\n==========================================================\
         ==========================\n"
    );

    // Trim the reads to remove sequence after inner primer sites in the reverse
    // (WRONG) orientation.  This has a large positive effect.
    // ◼ Should test: what happens if you remove them in the right orientation?
    // ◼ Or the wrong and the right orientation?

    let mut rc_inner_primers = Vec::<String>::with_capacity(inner_primers.len());
    let rc_inner_primers_bytes = inner_primers
        .iter()
        .map(|p| {
            let mut p = p.clone();
            reverse_complement(&mut p);
            rc_inner_primers.push(stringme(&p));
            p
        })
        .collect::<Vec<_>>();
    let k = 20;
    if heur.prim_trim {
        for i in 0..reads.len() {
            // In double end case, tried ignoring even-numbered reads, but that
            // appeared to yield somewhat worse results.

            let mut s = reads[i].to_ascii_vec();
            reverse_complement(&mut s);
            let mut trimmed = false;
            'outer: for j in 0..inner_primers.len() {
                let n = inner_primers[j].len();
                if s.len() >= n {
                    for p in (0..s.len() - n + 1).rev() {
                        if s[p..p + n] == rc_inner_primers_bytes[j][0..n] {
                            s = s[0..p + n].to_owned();
                            trimmed = true;
                            break 'outer;
                        }
                    }
                }
            }
            if trimmed {
                let trim = reads[i].len() - s.len();
                let b = DnaString::from_acgt_bytes(&s);
                let s = b.rc().to_string();
                reads[i] = DnaString::from_dna_string(&s);
                quals[i] = quals[i][trim..quals[i].len()].to_owned();
            }
        }
    }

    // Compute primer stats.
    // ◼ The code duplication here for inner/outer is really ugly.

    let rc_outer_primers_bytes = outer_primers
        .iter()
        .map(|p| {
            let mut p = p.clone();
            reverse_complement(&mut p);
            p
        })
        .collect::<Vec<_>>();
    barcode_data.inner_hit = vec![0; inner_primers.len()];
    barcode_data.inner_hit_good = vec![0; inner_primers.len()];
    barcode_data.outer_hit = vec![0; outer_primers.len()];
    barcode_data.outer_hit_good = vec![0; outer_primers.len()];
    const MAX_MIS: usize = 2;
    for (i, read) in reads.iter().enumerate() {
        // Only use the second read in each original pair and only use a fraction
        // of the reads.

        if single_end && i % 16 != 0 {
            continue;
        }
        if !single_end && i % 16 != 1 {
            continue;
        }

        barcode_data.primer_sample_count += 1;
        let s = read.to_ascii_vec();
        let mut trimmed = false;
        'outerx: for j in 0..inner_primers.len() {
            let n = inner_primers[j].len();
            if s.len() >= n {
                let p = s.len() - n;
                if s[p..p + n] == rc_inner_primers_bytes[j][0..n] {
                    trimmed = true;
                    barcode_data.inner_hit[j] += 1;
                    if p + n >= PRIMER_EXT_LEN {
                        let mut min_mis = 1000000;
                        for l in 0..inner_primer_exts[j].len() {
                            let x = &inner_primer_exts[j][l];
                            let mut mis = 0;
                            for m in 0..PRIMER_EXT_LEN {
                                if s[p + n - m - 1] != x[PRIMER_EXT_LEN - m - 1] {
                                    mis += 1;
                                }
                            }
                            min_mis = min(min_mis, mis);
                        }
                        if min_mis <= MAX_MIS {
                            barcode_data.inner_hit_good[j] += 1;
                        }
                    }
                    break 'outerx;
                }
            }
        }
        if !trimmed {
            'outerx2: for j in 0..outer_primers.len() {
                let n = outer_primers[j].len();
                if s.len() >= n {
                    let p = s.len() - n;
                    if s[p..p + n] == rc_outer_primers_bytes[j][0..n] {
                        barcode_data.outer_hit[j] += 1;
                        if p + n >= PRIMER_EXT_LEN {
                            let mut min_mis = 1000000;
                            for l in 0..outer_primer_exts[j].len() {
                                let x = &outer_primer_exts[j][l];
                                let mut mis = 0;
                                for m in 0..PRIMER_EXT_LEN {
                                    if s[p + n - m - 1] != x[PRIMER_EXT_LEN - m - 1] {
                                        mis += 1;
                                    }
                                }
                                min_mis = min(min_mis, mis);
                            }
                            if min_mis <= MAX_MIS {
                                barcode_data.outer_hit_good[j] += 1;
                            }
                        }
                        break 'outerx2;
                    }
                }
            }
        }
    }

    // Kill contamination.

    if !bad_umis.is_empty() {
        for (read, qual, &umi) in izip!(reads.iter_mut(), quals.iter_mut(), umi_id) {
            if bin_member(&bad_umis, &uu[umi as usize]) {
                read.clear();
                qual.clear();
            }
        }
    }

    // Build the graph.

    log_opts.report_perf_stats(log, &t, "before building graph");
    let mut x: Hyper = Hyper::new();
    x.build_from_reads(20, reads);
    let edges_initial = x.h.g.edge_count();
    log_opts.report_perf_stats(log, &t, "after building graph");

    // Simplify the graph without using the reference.

    let _t = Instant::now();
    simplify_without_ref(&mut x, umi_id, log, log_opts);

    // Simplify the graph using the reference.  First we pop bubbles, then we
    // kill certain branches that arise from alternate splicing e.g. hopping from
    // from TRBV18 to TRBV19, when the event is supported by only one UMI.

    if !heur.free {
        let mut dels = Vec::<u32>::new();
        pop_bubbles_using_reference(&mut x, umi_id, refdata);
        let mut ann = vec![Vec::<(i32, i32, i32, i32, i32)>::new(); x.h.g.edge_count()];
        for (e, anne) in ann.iter_mut().enumerate() {
            annotate_seq(x.h.g.edge_obj(e as u32), refdata, anne, false, false, false);
        }
        for (e1, (anne1, xids1)) in zip(&ann, &x.ids).take(x.h.g.edge_count()).enumerate() {
            if anne1.is_empty() {
                continue;
            }
            let t1 = anne1[anne1.len() - 1].2 as usize;
            if !rheaders[t1].contains("V-REGION") {
                continue;
            }
            let gene1 = rheaders[t1].after("|").between("|", "|");
            if !gene1.contains('V') {
                continue;
            }
            let n1 = match gene1.after("V").parse::<i32>() {
                Ok(n) => n,
                Err(_) => continue,
            };
            let right = x.h.g.to_right(e1 as u32);
            for (e2, anne2) in ann.iter().take(x.h.g.edge_count()).enumerate() {
                // oops quadratic
                if right != x.h.g.to_left(e2 as u32) {
                    continue;
                }
                if anne2.is_empty() {
                    continue;
                }
                let t2 = anne2[0].2 as usize;
                if !rheaders[t2].contains("V-REGION") {
                    continue;
                }
                let gene2 = rheaders[t2].after("|").between("|", "|");
                if !gene2.contains('V') || gene2.after("V").parse::<i32>().is_err() {
                    continue;
                }
                let n2 = gene2.after("V").force_i32();
                if n1 >= n2 || n1 < n2 - 2 {
                    continue;
                }
                let mut umis = Vec::<i32>::new();
                for &i in xids1 {
                    let u = umi_id[i as usize];
                    if umis.is_empty() || umis[0] != u {
                        umis.push(u);
                    }
                    if umis.len() > 1 {
                        break;
                    }
                }
                if umis.len() <= 1 {
                    dels.push(e1 as u32);
                }
            }
        }
        x.kill_edges(&dels);
    }
    log_opts.report_perf_stats(log, &t, "after simplifying graph");
    fwriteln!(
        log,
        "graph has {} edges initially, {} edges after simplification",
        edges_initial,
        x.h.g.edge_count()
    );

    // Function to form and print paths from a list of edges.

    fn print_paths_from_edges(x: &Hyper, es: &[usize], log: &mut Vec<u8>) {
        let mut used = vec![false; es.len()];
        for i in 0..es.len() {
            if used[i] {
                continue;
            }
            used[i] = true;
            let mut p = vec![es[i]];
            loop {
                let mut exts = Vec::<usize>::new();
                for (j, &esj) in es.iter().enumerate() {
                    if p.contains(&esj) || used[j] {
                        continue;
                    }
                    if x.h.g.to_right(esj as u32) == x.h.g.to_left(p[0] as u32) {
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
                for (j, &esj) in es.iter().enumerate() {
                    if p.contains(&esj) || used[j] {
                        continue;
                    }
                    if x.h.g.to_left(esj as u32) == x.h.g.to_right(p[p.len() - 1] as u32) {
                        exts.push(j);
                    }
                }
                if exts.len() != 1 {
                    break;
                }
                used[exts[0]] = true;
                p.push(es[exts[0]]);
            }
            fwriteln!(log, "path: {}", p.iter().format(","));
        }
    }

    // Trace umis.

    if log_opts.trace_umis {
        fwriteln!(log, "\nTRACING UMIS");
        let mut ue = vec![Vec::<usize>::new(); uu.len()];
        for e in 0..x.h.g.edge_count() {
            for &j in &x.ids[e] {
                let u = umi_id[j as usize];
                ue[u as usize].push(e);
            }
        }
        for (u, mut ue) in ue.into_iter().take(uu.len()).enumerate() {
            let n = ue.len();
            unique_sort(&mut ue);
            fwriteln!(log, "\nu = {}, nreads = {}", u, n);
            print_paths_from_edges(&x, &ue, log);
        }
    }

    // Trace sequence.

    if !log_opts.trace_seq.is_empty() {
        fwriteln!(log, "\nsearching for TRACE_SEQ\n");
        fwriteln!(
            log,
            "Showing edges through the graph that overlap TRACE_SEQ"
        );
        fwriteln!(
            log,
            "by at least one kmer; these are displayed in paths but"
        );
        fwriteln!(
            log,
            "these paths are not necessarily the only paths these edges lie in.\n"
        );
        fwriteln!(log, "Also showing missing kmers.\n");
        let mut lookup = Vec::<(Kmer20, i32, i32)>::new();
        let b = DnaString::from_dna_string(&log_opts.trace_seq);
        let bv = vec![b.clone()];
        make_kmer_lookup_20_single(&bv, &mut lookup);
        let mut es = Vec::new();
        let mut found = vec![false; b.len() - k + 1];
        for e in 0..x.h.g.edge_count() {
            let m = &x.h.g.edge_obj(e as u32);
            for i in 0..m.len() - k + 1 {
                let y: Kmer20 = m.get_kmer(i);
                let low = lower_bound1_3(&lookup, &y) as usize;
                let high = upper_bound1_3(&lookup, &y) as usize;
                if low < high {
                    es.push(e);
                    for j in low..high {
                        found[lookup[j].2 as usize] = true;
                    }
                }
            }
        }
        unique_sort(&mut es);
        for (i, found) in found.into_iter().enumerate() {
            if !found {
                let z: Kmer20 = b.get_kmer(i);
                fwriteln!(log, "missing kmer {} = {}", i, z.to_string());
            }
        }
        print_paths_from_edges(&x, &es, log);
    }

    // Compute TRAV, TRBV, TRGV, and TRDV read counts. It is possible that it would be better
    // to compute this on the unsimplified graph.
    // To add note here on why we're doing this.

    // [ av, bv,
    //   gv, dv,
    //   aj, bj,
    //   gj, dj ]
    let mut tr = vec![Vec::<u32>::default(); 8];
    if !heur.free {
        for e in 0..x.h.g.edge_count() {
            let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(
                x.h.g.edge_obj(e as u32),
                refdata,
                &mut ann,
                false,
                false,
                false,
            );
            // [ av, bv,
            //   gv, dv,
            //   aj, bj,
            //   gj, dj ]
            let mut is_tr = vec![false; 8];
            for i in 0..ann.len() {
                let h = &rheaders[ann[i].2 as usize];
                if h.contains("TRAV") {
                    is_tr[0] = true;
                }
                if h.contains("TRBV") {
                    is_tr[1] = true;
                }
                if gd_mode && h.contains("TRGV") {
                    is_tr[2] = true;
                }
                if gd_mode && h.contains("TRDV") {
                    is_tr[3] = true;
                }
                if h.contains("TRAJ") {
                    is_tr[4] = true;
                }
                if h.contains("TRBJ") {
                    is_tr[5] = true;
                }
                if gd_mode && h.contains("TRGJ") {
                    is_tr[6] = true;
                }
                if gd_mode && h.contains("TRDJ") {
                    is_tr[7] = true;
                }
            }
            if is_tr.iter().any(|&x| x) {
                for i in 0..x.ids[e].len() {
                    for j in 0..is_tr.len() {
                        if is_tr[j] {
                            tr[j].push(x.ids[e][i]);
                        }
                    }
                }
            }
        }
        for v in &mut tr {
            v.sort_unstable();
            v.dedup();
        }

        barcode_data.trav_reads = tr[0].len() as i32;
        barcode_data.trbv_reads = tr[1].len() as i32;
        barcode_data.trgv_reads = tr[2].len() as i32;
        barcode_data.trdv_reads = tr[3].len() as i32;
        barcode_data.travj_reads = meet_size(&tr[0], &tr[4]) as i32;
        barcode_data.trbvj_reads = meet_size(&tr[1], &tr[5]) as i32;
        barcode_data.trgvj_reads = meet_size(&tr[2], &tr[6]) as i32;
        barcode_data.trdvj_reads = meet_size(&tr[3], &tr[7]) as i32;
    }

    // Print the reads.

    if log_opts.umi_seq {
        let mut id = 0;
        while id < reads.len() {
            let id2 = next_diff(umi_id, id);
            let null = reads[id..id2].iter().all(DnaString::is_empty);
            if !null {
                fwriteln!(log, "\numi {} = {}", umi_id[id], uu[umi_id[id] as usize]);
                for (idx, read) in reads[id..id2].iter().enumerate() {
                    fwriteln!(
                        log,
                        ">{}={}.{}\n{}",
                        idx,
                        (idx - id) / 2 + 1,
                        (idx - id) % 2 + 1,
                        read.to_string()
                    );
                }
            }
            id = id2;
        }
    }

    // Compute ucounts, and umi_info.

    let mut ucounts = Vec::<i32>::new();
    let mut total_ucounts = 0;
    let mut id = 0;
    let none_pos = CHAIN_TYPES
        .into_iter()
        .enumerate()
        .find(|v| v.1 == "None")
        .map_or(-1_i32, |v| v.0 as i32);
    while id < reads.len() {
        let id2 = next_diff(umi_id, id);
        let null = reads[id..id2].iter().all(DnaString::is_empty);
        let mut n = (id2 - id) as i32;
        if !single_end {
            n /= 2;
        }
        if !null {
            total_ucounts += 1;
        }
        if (n > 1 || n50_n50_rpu == 2) && !null {
            ucounts.push(n);
        }
        barcode_data
            .umi_info
            .push((n, none_pos as i8, uu[umi_id[id] as usize].clone()));
        id = id2;
    }
    barcode_data.total_ucounts = total_ucounts;
    ucounts.sort_unstable();
    if !log_opts.nucounts {
        fwriteln!(log, "\ntotal ucounts = {}", total_ucounts);
        fwriteln!(
            log,
            "nonsolo ucounts = {}[{}]",
            ucounts.len(),
            abbrev_list(&ucounts)
        );
    }
    if log_opts.print_umi {
        let mut x = Vec::<(i32, String)>::new();
        let mut id = 0;
        while id < reads.len() {
            let id2 = next_diff(umi_id, id);
            let mut n = (id2 - id) as i32;
            if !single_end {
                n /= 2;
            }
            if (n > 1 || n50_n50_rpu == 2) && !reads[id..id2].iter().all(DnaString::is_empty) {
                let u = umi_id[id];
                x.push((n, uu[u as usize].clone()));
            }
            id = id2;
        }
        x.sort();
        fwrite!(log, "\numis = ");
        for (i, xi) in x.into_iter().enumerate() {
            if i > 0 {
                fwrite!(log, ", ");
            }
            fwrite!(log, "{}[{}]", xi.1, xi.0);
        }
        fwriteln!(log, "");
    }
    barcode_data.ucounts = ucounts;

    // Find surviving umis.
    //
    // 1. Find strong paths p that have an annotation, or in the free case,
    //    match to a primer.
    //    ◼ Why not require correct orientation?
    //
    // 2. Find all the edges on a strong path.
    //
    // 3. Find all the reads on such a good edge.
    //
    // 4. Find the umis for these reads.
    //
    // 5. Restrict to those umis for which 50% of their kmers are contained in
    //    good edges.
    //
    // 6. In the non-free case, if no strong path had a V annotation, kill all the
    //    surviving umis.
    //    ◼ Why not require that there is a V annotation in the correct orientation?
    //
    // WARNING!  This code is not robust to datasets in which errors have been
    // artificially and randomly introduced.  At a sufficiently high error rate,
    // what you'll see is no surviving umis.  There is no obvious way to fix this,
    // but also no obvious need.  If we wanted to make the cost robust to naturally
    // occurring errors at a high rate, we might restrict attention to kmers having
    // only high quality scores.

    let mut strong = Vec::<(i32, Vec<i32>)>::new();
    // ◼ Below, needed x to be mutable, why?
    uber_strong_paths(&mut x, umi_id, &mut strong);
    let mut es = Vec::<i32>::new();
    let mut have_v = false;
    let mut strongx = strong.iter().map(|i| &i.1).collect::<Vec<_>>();
    unique_sort(&mut strongx);
    for p in strongx {
        let c = x.cat(p);
        let mut anns = 0;
        if !heur.free {
            let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq(&c, refdata, &mut ann, false, false, false);
            for i in &ann {
                let h = &refdata.rheaders[i.2 as usize];
                // ◼: Should just look for V region, and below too.
                if h.contains("TRAV")
                    || h.contains("TRBV")
                    || h.contains("IGHV")
                    || h.contains("IGKV")
                    || h.contains("IGLV")
                {
                    have_v = true;
                }
            }
            anns += ann.len();
            let crc = c.rc();
            annotate_seq(&crc, refdata, &mut ann, false, false, false);
            anns += ann.len();
            for i in ann {
                let h = &rheaders[i.2 as usize];
                if h.contains("TRAV")
                    || h.contains("TRBV")
                    || h.contains("IGHV")
                    || h.contains("IGKV")
                    || h.contains("IGLV")
                {
                    have_v = true;
                }
            }
        } else {
            let s = c.to_string();
            for (inner_primer, rprim) in zip(inner_primers, &rc_inner_primers) {
                let prim = stringme(inner_primer);
                if s.contains(&prim) || s.contains(rprim) {
                    anns += 1;
                }
            }
        }
        if anns == 0 {
            continue;
        }
        for &i in p {
            es.push(i);
            es.push(x.inv[i as usize] as i32);
        }
    }
    if !heur.free && !have_v {
        es.clear();
    }
    unique_sort(&mut es);
    if log_opts.strong_edges {
        fwriteln!(log, "strong edges = {:?}", es);
    }
    let mut goodtigs = Vec::<DnaString>::new();
    for &i in &es {
        goodtigs.push(x.h.g.edge_obj(i as u32).clone());
    }
    let mut xgood = Vec::<(Kmer20, i32, i32)>::new();
    make_kmer_lookup_20_single(&goodtigs, &mut xgood);
    let mut us = es
        .iter()
        .flat_map(|&e| x.ids[e as usize].iter().map(|&id| umi_id[id as usize]))
        .collect::<Vec<_>>();
    unique_sort(&mut us);
    let mut xucounts2 = Vec::<(i32, i32)>::new();
    let mut id = 0;
    let mut ux = 0;
    const MIN_GOODNESS: f64 = 0.5;
    while id < reads.len() {
        let id2 = next_diff(umi_id, id);
        let u = umi_id[id];
        if bin_member(&us, &u) {
            let mut n = (id2 - id) as i32;
            if !single_end {
                n /= 2;
            }
            if n > 1 || n50_n50_rpu == 2 {
                let (mut total, mut good) = (0, 0);
                for b in &reads[id..id2] {
                    if b.len() >= k {
                        for j in 0..b.len() - k + 1 {
                            let y: Kmer20 = b.get_kmer(j);
                            let low = lower_bound1_3(&xgood, &y) as usize;
                            total += 1;
                            if low < xgood.len() && y == xgood[low].0 {
                                good += 1;
                            }
                        }
                    }
                }
                let goodness = good as f64 / total as f64;
                if log_opts.survive {
                    fwriteln!(
                        log,
                        "umi {} has {} pairs and goodness {:.2}%",
                        ux,
                        n,
                        100_f64 * goodness
                    );
                }
                // ◼ Had to lower value from 0.75 to 0.5 for woofx.  Should figure
                // ◼ out why and change code.  A possible cause: kmers having low
                // ◼ quality scores.  One way to counteract this would be to count
                // ◼ only kmers whose quality score exceeded a given threshold.
                if goodness >= MIN_GOODNESS {
                    xucounts2.push((n, u));
                }
            }
        }
        ux += 1;
        id = id2;
    }
    xucounts2.sort_unstable();
    for i in &xucounts2 {
        barcode_data.xucounts.push(i.0);
        barcode_data.xuids.push(i.1);
    }
    if !log_opts.nucounts {
        fwriteln!(
            log,
            "surviving nonsolo ucounts = {}[{}]",
            barcode_data.xucounts.len(),
            abbrev_list(&barcode_data.xucounts)
        );
        if barcode_data.xuids.len() <= 10 {
            fwriteln!(log, "ids = {:?}", barcode_data.xuids);
        } else {
            let mut x = barcode_data.xuids.clone();
            x.truncate(10);
            fwriteln!(log, "ids = [{}, ...]", x.iter().format(", "));
        }
    }

    // Determine chain types so we can compute metrics for them.  First we use
    // strong paths to assign some of the UMIs.  Then we used a direct calculation
    // to try to assign the rest of them.  If we were instead to used a direct
    // calculation to try to assign all UMIs, total run time would increase a lot
    // (~30% as measured by running vdj_asm_demo on the standard samples).

    if !refdata.refs.is_empty() {
        let chain_types_list = CHAIN_TYPES.to_vec();
        let mut ann_memory = HashMap::<Vec<i32>, Vec<(i32, i32, i32, i32, i32)>>::new();
        for j in &strong {
            let p = &j.1;
            // Can't use entry.or_insert_with here because we don't want to
            // clone p unnecessarily.
            if !ann_memory.contains_key(p) {
                let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
                let c = x.cat(p);
                annotate_seq(&c, refdata_full, &mut ann, false, false, false);
                ann_memory.insert(p.clone(), ann);
            };
            for i in &ann_memory[p] {
                let h = &rheaders_full[i.2 as usize];
                let mut chain_types = Vec::<i8>::with_capacity(7);
                if h.contains("TRAV") || h.contains("TRAJ") {
                    chain_types.push(position(&chain_types_list, &"TRA") as i8);
                }
                if h.contains("TRBV") || h.contains("TRBJ") {
                    chain_types.push(position(&chain_types_list, &"TRB") as i8);
                }
                if h.contains("TRDV") || h.contains("TRDJ") {
                    chain_types.push(position(&chain_types_list, &"TRD") as i8);
                }
                if h.contains("TRGV") || h.contains("TRGJ") {
                    chain_types.push(position(&chain_types_list, &"TRG") as i8);
                }
                if h.contains("IGHV") || h.contains("IGHJ") {
                    chain_types.push(position(&chain_types_list, &"IGH") as i8);
                }
                if h.contains("IGKV") || h.contains("IGKJ") {
                    chain_types.push(position(&chain_types_list, &"IGK") as i8);
                }
                if h.contains("IGLV") || h.contains("IGLJ") {
                    chain_types.push(position(&chain_types_list, &"IGL") as i8);
                }
                unique_sort(&mut chain_types);
                if chain_types.len() == 1 {
                    barcode_data.umi_info[j.0 as usize].1 = chain_types[0];
                }
            }
        }
        let mut id = 0;
        let mut ucount = 0;
        let types = vec!["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
        while id < reads.len() {
            let id2 = next_diff(umi_id, id);
            if barcode_data.umi_info[ucount].1 == none_pos as i8 {
                let mut cts = Vec::<u8>::new();
                for b in &reads[id..id2] {
                    let ct = chain_type(b, rkmers_plus_full_20, &refdata_full.rtype);
                    if ct >= 0 {
                        cts.push(if ct >= 7 { ct - 7 } else { ct } as u8);
                    }
                }
                unique_sort(&mut cts);
                if cts.len() == 1 {
                    barcode_data.umi_info[ucount].1 =
                        position(&chain_types_list, &types[cts[0] as usize]) as i8;
                }
            }
            ucount += 1;
            id = id2;
        }
    }

    // Print final graph.

    log_opts.report_perf_stats(log, &t, "after printing final graph");
    if !log_opts.ngraph {
        fwriteln!(log, "\nFINAL GRAPH\n");
        fwriteln!(
            log,
            "edges = {}, checksum = {}",
            x.h.g.edge_count(),
            x.checksum()
        );
        print_with_annotations(
            heur.free,
            &x,
            reads,
            umi_id,
            refdata,
            log_opts.show_supp,
            log_opts.print_seq_edges,
            log,
        );
    }
    fwriteln!(
        log,
        "\n==========================================================\
         =========================="
    );
    *nedges = x.h.g.edge_count();

    // Make contigs.

    let mut con = Vec::<DnaString>::new(); // good contigs
    let mut con2 = Vec::<DnaString>::new(); // reject contigs
    log_opts.report_perf_stats(log, &t, "before making contigs");
    let mut jsupp = Vec::<(i32, i32)>::new(); // (umis, reads) supporting junction seqs of good contigs

    make_contigs(
        single_end,
        is_gd,
        inner_primers,
        reads,
        quals,
        &mut x,
        umi_id,
        uu,
        refdata,
        refdata_full,
        rkmers_plus_full_20,
        heur,
        &mut con,
        &mut jsupp,
        &mut con2,
        conxq,
        cumi,
        cids,
        barcode_data,
        validated_umis,
        non_validated_umis,
        invalidated_umis,
        log,
        log_opts,
    );
    log_opts.report_perf_stats(log, &t, "after making contigs");
    conx.extend(con.iter().cloned());
    conx.extend(con2.iter().cloned());
    productive.extend(vec![true; con.len()]); // good contigs == productive contigs (non-denovo mode)
    productive.extend(vec![false; con2.len()]); // reject contigs == unproductive contigs (non-denovo mode)
    junction_support.extend(jsupp.iter().map(|(umis, reads)| {
        Some(JunctionSupport {
            reads: *reads,
            umis: *umis,
        })
    }));
    junction_support.extend(vec![None; con2.len()]);
    barcode_data_brief.ncontigs = conx.len();
    barcode_data_brief.xucounts = barcode_data.xucounts.clone();
    barcode_data_brief.total_ucounts = barcode_data.total_ucounts;

    // Compute primer hits to good contigs.

    barcode_data.inner_hit_good_contigs = vec![0; inner_primers.len()];
    barcode_data.outer_hit_good_contigs = vec![0; outer_primers.len()];
    for (rc_inner, inner_hit) in zip(
        rc_inner_primers_bytes,
        &mut barcode_data.inner_hit_good_contigs,
    )
    .take(inner_primers.len())
    {
        for l in &con {
            if l.to_ascii_vec().ends_with(&rc_inner) {
                *inner_hit += 1;
            }
        }
    }
    for (rc_outer, outer_hit) in zip(
        rc_outer_primers_bytes,
        &mut barcode_data.outer_hit_good_contigs,
    )
    .take(outer_primers.len())
    {
        for l in &con {
            if l.to_ascii_vec().ends_with(&rc_outer) {
                *outer_hit += 1;
            }
        }
    }

    // Report good contigs.
    let mut ann_all = Vec::<Vec<(i32, i32, i32, i32, i32)>>::new();
    for b in &con {
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        if refdata.refs.is_empty() {
            ann_all.push(ann.clone());
        } else {
            // ◼ There must be a duplicate computation in make_contigs.
            annotate_seq_core(b, refdata, &mut ann, true, false, true, log, false);
            ann_all.push(ann.clone());
            barcode_data.good_refseqs.extend(ann.iter().map(|a| a.2));
        }
    }
    if !con.is_empty() && !log_opts.ngood {
        fwriteln!(log, "\nGOOD CONTIGS");
        for j in 0..con.len() {
            let tigname = format!("{}_contig_{}", barcode_data.barcode, j + 1);
            print_contig_info(
                &tigname,
                label,
                log,
                true,
                j,
                &con[j],
                &ann_all[j],
                &jsupp[j],
                &validated_umis[j],
                &non_validated_umis[j],
                &invalidated_umis[j],
                &conxq[j],
                &cumi[j],
                &cids[j],
                barcode_data,
                is_tcr,
                is_bcr,
                is_gd,
                refdata,
                refdatax,
                heur.free,
                log_opts,
            );

            // Align reads to the contig.

            /*
            // ◼ This code is off because it's insanely slow.
            // ◼ Also it crashes on snarl.
            let mut read_alignments = HashMap::<u32,AlignmentPacket>::new();
            align_reads_to_tig( &con[j], &cids[j], &reads, Some(MATCH_SCORE),
                Some(MISMATCH_SCORE), Some(GAP_OPEN), Some(GAP_EXTEND), Some(CLIP),
                None, &mut read_alignments );
            fwriteln!( log, "reads aligned: {}", read_alignments.len() );
            */
        }
    }

    // Report reject contigs.

    if !log_opts.nreject && !con2.is_empty() {
        fwriteln!(log, "\nREJECT CONTIGS");
        for j in 0..con2.len() {
            let tigname = format!("{}_{}", barcode_data.barcode, con.len() + j + 1);
            let js = (0, 0);
            let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
            annotate_seq_core(&con2[j], refdata, &mut ann, true, false, true, log, false);
            print_contig_info(
                &tigname,
                label,
                log,
                false,
                j,
                &con2[j],
                &ann,
                &js,
                &validated_umis[con.len() + j],
                &non_validated_umis[con.len() + j],
                &invalidated_umis[con.len() + j],
                &conxq[con.len() + j],
                &cumi[con.len() + j],
                &cids[con.len() + j],
                barcode_data,
                is_tcr,
                is_bcr,
                is_gd,
                refdata,
                refdatax,
                heur.free,
                log_opts,
            );

            // Check for long unannotated regions and note special case.

            if !heur.free {
                let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
                annotate_seq(&con2[j], refdata, &mut ann, true, false, true);
                // assert!(is_gd.unwrap_or(false)); -> is_gd == true
                is_valid(&con2[j], refdata, &ann, true, log, is_gd);
                annotate_seq(&con2[j], refdatax, &mut ann, true, false, true);
                if !ann.is_empty() && ann[0].0 >= 200 {
                    fwriteln!(log, "note long unannotated region");
                }
                if ann.len() == 1 {
                    let t = ann[0].2;
                    if refdatax.rheaders[t as usize].contains("segment preceding TRBD1 exon 1") {
                        let start = (ann[0].0 + ann[0].1) as usize;
                        let n = 20;
                        if start + n <= con2[j].len() {
                            let s = con2[j].slice(start, start + n).to_string();
                            fwriteln!(log, "TRBD1 signature: {}", s);
                        }
                    }
                }
            }
        }
    }

    // Finish up logging.

    log_opts.report_perf_stats(log, &t, "after reporting rejects");
}
