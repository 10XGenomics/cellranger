// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Find the optimal D segment, the runner up, and the delta between the scores.  This uses
// the donor V and J segments that are assigned to the clonotype.  Note that the optimal D
// segment may be null.  This is obvious from looking at data.

use crate::align_to_vdj_ref::{align_to_vdj_ref, match_bit_score, zero_one};
use enclone_proto::types::DonorReferenceItem;
use std::cmp::min;
use vdj_ann::annotate::Annotation;
use vdj_ann::refx::RefData;

pub fn vflank(_seq: &[u8], vref: &[u8]) -> usize {
    let mut flank = 13;
    if flank > vref.len() {
        flank = vref.len();
    }
    flank
}

pub fn jflank(seq: &[u8], jref: &[u8]) -> usize {
    let flank = 13;
    if flank > jref.len() {
        return jref.len();
    }

    // Let start be the first position on the J gene where there is a perfect match of length at
    // least five to the contig.

    const MATCHLEN: usize = 5;
    let mut start = 0;
    for i in 0..=jref.len() - MATCHLEN {
        let mut matchlen = 0;
        for j in 0..MATCHLEN {
            if seq[seq.len() - jref.len() + i + j] != jref[i + j] {
                break;
            }
            matchlen += 1;
        }
        if matchlen == MATCHLEN {
            start = i;
            break;
        }
    }

    // Add start to the flank, so long as that's possible.

    min(flank + start, jref.len())
}
#[allow(clippy::too_many_arguments)]
pub fn evaluate_d(
    tig: &[u8],
    vref: &[u8],
    seq_start: usize,
    ds: &[usize],
    jref: &[u8],
    refdata: &RefData,
    jscore_match: i32,
    jscore_mismatch: i32,
    jscore_gap_open: i32,
    jscore_gap_extend: i32,
    jscore_bits_multiplier: f64,
) -> (Vec<bio_edit::alignment::AlignmentOperation>, f64) {
    // Start to build reference concatenation.  First append the V segment.

    let mut concat = Vec::<u8>::new();
    let vstart = vref.len() - vflank(tig, vref);
    let vref = &vref[vstart..vref.len()];
    concat.extend(vref);

    // Append the D segment or segments.

    let mut dref = Vec::<u8>::new();
    let mut d2ref = Vec::<u8>::new();
    let mut drefname = String::new();
    for (j, &d) in ds.iter().enumerate() {
        if j == 0 {
            dref = refdata.refs[d].to_ascii_vec();
        } else if j == 1 {
            d2ref = refdata.refs[d].to_ascii_vec();
        }
        if j > 0 {
            drefname += ":";
        }
        drefname += refdata.name[d].as_str();
    }
    concat.extend(&dref);
    concat.extend(&d2ref);

    // Append the J segment.

    let jend = jflank(tig, jref);

    // Align the V..J sequence on the contig to the reference concatenation.

    let mut seq_end = tig.len() - (jref.len() - jend);
    if seq_end <= seq_start {
        seq_end = tig.len(); // bug fix for problem found by customer, couldn't reproduce internally
    }
    let seq = &tig[seq_start..seq_end];
    let jref = &jref[0..jend];
    concat.extend(jref);
    let (ops, count) = align_to_vdj_ref(
        seq,
        vref,
        &dref,
        &d2ref,
        jref,
        &drefname,
        true,
        jscore_match,
        jscore_mismatch,
        jscore_gap_open,
        jscore_gap_extend,
        jscore_bits_multiplier,
    );
    (ops, count)
}
#[allow(clippy::too_many_arguments)]
pub fn opt_d(
    v_ref_id: usize,     // ex.share[mid].v_ref_id
    j_ref_id: usize,     // ex.share[mid].j_ref_id
    tig: &[u8],          // ex.share[mid].seq_del
    annv: &[Annotation], // ex.share[mid].annv
    cdr3_aa: &str,       // ex.share[mid].cdr3_aa
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    dsx: &mut Vec<Vec<usize>>,
    jscore_match: i32,
    jscore_mismatch: i32,
    jscore_gap_open: i32,
    jscore_gap_extend: i32,
    jscore_bits_multiplier: f64,
    v_alt: Option<usize>,
) {
    let mut comp = 1000000.0;

    // Go through every D segment, or possibly every concatenation of D segments.

    let mut todo = vec![vec![]];
    for i in &refdata.ds {
        todo.push(vec![*i]);
    }
    let mut ds = Vec::<Vec<usize>>::new();
    let mut counts = Vec::<f64>::new();
    let mut good_d = Vec::<usize>::new();
    let mut vref = refdata.refs[v_ref_id].to_ascii_vec();
    if let Some(v_alt) = v_alt {
        vref.clone_from(&dref[v_alt].nt_sequence);
    }
    let vstart = vref.len() - vflank(tig, &vref);
    let mut seq_start = vstart as isize;
    // probably not exactly right
    if annv.len() > 1 {
        let q1 = annv[0].tig_start + annv[0].match_len;
        let q2 = annv[1].tig_start;

        seq_start += q1 as isize - q2 as isize;
    }
    let jref = refdata.refs[j_ref_id].to_ascii_vec();
    const MIN_BITS_FOR_D2: f64 = 14.0;
    for di in &todo {
        let (ops, count) = evaluate_d(
            tig,
            &vref,
            seq_start as usize,
            di,
            &jref,
            refdata,
            jscore_match,
            jscore_mismatch,
            jscore_gap_open,
            jscore_gap_extend,
            jscore_bits_multiplier,
        );
        counts.push(count);
        if !di.is_empty() {
            let drefx = refdata.refs[di[0]].to_ascii_vec();
            let vstart = vref.len() - vflank(tig, &vref);
            let vref = vref[vstart..vref.len()].to_vec();
            let zos = zero_one(&ops, vref.len(), vref.len() + drefx.len());
            let bits = match_bit_score(&zos);
            if bits >= MIN_BITS_FOR_D2 {
                good_d.push(di[0]);
            }
        }
        ds.push(di.clone());
        if count > comp {
            comp = count;
        }
    }
    if cdr3_aa.len() >= 20 {
        todo.clear();
        for &i1 in &good_d {
            for &i2 in &good_d {
                todo.push(vec![i1, i2]);
            }
        }
        for di in &todo {
            let (_ops, count) = evaluate_d(
                tig,
                &vref,
                seq_start as usize,
                di,
                &jref,
                refdata,
                jscore_match,
                jscore_mismatch,
                jscore_gap_open,
                jscore_gap_extend,
                jscore_bits_multiplier,
            );
            counts.push(count);
            ds.push(di.clone());
            if count > comp {
                comp = count;
            }
        }
    }

    // Reverse sort sync (counts, ds).

    let mut counts_ds = counts
        .iter()
        .zip(ds.iter())
        .map(|(&c, d)| (c, d.clone()))
        .collect::<Vec<_>>();
    counts_ds.sort_by(|a, b| b.partial_cmp(a).unwrap()); // reverse sort
    counts.clear();
    ds.clear();
    for count in counts_ds {
        counts.push(count.0);
        ds.push(count.1);
    }

    *dsx = ds;
}
