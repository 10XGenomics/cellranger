// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use super::align_to_vdj_ref::align_to_vdj_ref;
use super::defs::{ExactClonotype, Junction};
use super::opt_d::{jflank, opt_d};
use bio::alignment::AlignmentOperation::{Del, Ins, Subst};
use enclone_proto::types::DonorReferenceItem;
use rayon::prelude::*;
use vdj_ann::refx::RefData;

// This is largely copied from align_n.

pub(crate) fn heavy_complexity(
    refdata: &RefData,
    exact_clonotypes: &[ExactClonotype],
    dref: &[DonorReferenceItem],
) -> Vec<Junction> {
    let mut results = Vec::<(usize, Junction)>::new();
    for i in 0..exact_clonotypes.len() {
        results.push((i, Junction::default()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let ex = &exact_clonotypes[i];
        for r in 0..ex.share.len() {
            if ex.share[r].left && ex.share.len() == 2 && !ex.share[1 - r].left {
                let seq = ex.share[r].seq_del.as_ref();
                let mut vref = refdata.refs[ex.share[r].v_ref_id].to_ascii_vec();
                if ex.share[r].v_ref_id_donor.is_some() {
                    vref.clone_from(&dref[ex.share[r].v_ref_id_donor.unwrap()].nt_sequence);
                }
                let mut vstart = ex.share[r].cdr3_start - 2;

                // Compensate for indel.  Code here and next work imperfectly and
                // there would be value in investigating the error cases.

                if !ex.share[r].ins.is_empty() {
                    vstart -= ex.share[r].ins[0].1.len();
                } else if ex.share[r].seq.len() < ex.share[r].seq_del.len() {
                    vstart += ex.share[r].seq_del.len() - ex.share[r].seq.len();
                }

                // Prevent crash (working around bug).

                if vstart > vref.len() {
                    vstart = vref.len();
                }

                // Keep going.

                let vref = &vref[vstart..vref.len()];
                let mut drefx = Vec::<u8>::new();
                let mut d2ref = Vec::<u8>::new();
                let mut drefname = String::new();
                let mut ds = Vec::<Vec<usize>>::new();
                opt_d(
                    ex.share[r].v_ref_id,
                    ex.share[r].j_ref_id,
                    &ex.share[r].seq_del,
                    &ex.share[r].annv,
                    &ex.share[r].cdr3_aa,
                    refdata,
                    dref,
                    &mut ds,
                    JSCORE_MATCH,
                    JSCORE_MISMATCH,
                    JSCORE_GAP_OPEN,
                    JSCORE_GAP_EXTEND,
                    JSCORE_BITS_MULTIPLIER,
                    ex.share[r].v_ref_id_donor,
                );
                let mut opt = Vec::new();
                if !ds.is_empty() {
                    opt.clone_from(&ds[0]);
                }
                for (j, d) in opt.into_iter().enumerate() {
                    if j == 0 {
                        drefx = refdata.refs[d].to_ascii_vec();
                    } else {
                        d2ref = refdata.refs[d].to_ascii_vec();
                    }
                    if j > 0 {
                        drefname += ":";
                    }
                    drefname += &mut refdata.name[d].clone();
                }
                let jref = refdata.refs[ex.share[r].j_ref_id].to_ascii_vec();
                let jend = jflank(seq, &jref);
                let mut seq_start = vstart as isize;
                // probably not exactly right
                if ex.share[r].annv.len() > 1 {
                    let q1 = ex.share[r].annv[0].tig_start + ex.share[r].annv[0].match_len;
                    let q2 = ex.share[r].annv[1].tig_start;
                    seq_start += q2 as isize - q1 as isize;
                }
                let mut seq_end = seq.len() - (jref.len() - jend);
                // very flaky bug workaround
                // asserted on BCR=180030 CDR3=CARERDLIWFGPW JALIGN1
                if seq_start as usize > seq_end {
                    seq_start = vstart as isize;
                }
                if seq_end <= seq_start as usize {
                    seq_end = seq.len(); // bug fix for problem found by customer,
                    // couldn't reproduce internally
                }
                let seq = &seq[seq_start as usize..seq_end];
                let jref = &jref[..jend];
                let (ops, _score) = align_to_vdj_ref(
                    seq,
                    vref,
                    &drefx,
                    &d2ref,
                    jref,
                    &drefname,
                    true,
                    JSCORE_MATCH,
                    JSCORE_MISMATCH,
                    JSCORE_GAP_OPEN,
                    JSCORE_GAP_EXTEND,
                    JSCORE_BITS_MULTIPLIER,
                );
                let mut hcomp = 0;
                for i in 0..ops.len() {
                    if (ops[i] == Subst)
                        || (ops[i] == Ins)
                        || (ops[i] == Del && (i == 0 || ops[i - 1] != Del))
                    {
                        hcomp += 1;
                    }
                }
                res.1 = Junction { hcomp };
            }
        }
    });
    results.into_iter().map(|ri| ri.1).collect()
}

const JSCORE_MATCH: i32 = 20;
const JSCORE_MISMATCH: i32 = -20;
const JSCORE_BITS_MULTIPLIER: f64 = 2.2;
const JSCORE_GAP_OPEN: i32 = -120;
const JSCORE_GAP_EXTEND: i32 = -20;
