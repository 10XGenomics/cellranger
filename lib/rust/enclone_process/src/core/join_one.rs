// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use super::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use debruijn::Mer;
use debruijn::dna_string::ndiffs;
use qd::{Double, dd};
use std::cmp::min;
use std::collections::HashMap;
use string_utils::TextUtils;
use vdj_ann::refx::RefData;
// use stirling_numbers::p_at_most_m_distinct_in_sample_of_x_from_n;
use vector_utils::{meet, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is a copy of p_at_most_m_distinct_in_sample_of_x_from_n from the stirling_numbers crate,
// that has been modified to use higher precision internal math.  This should go into that crate
// (along with the stirling numbers ratio table code) when and if the qr crate is published.

pub(crate) fn p_at_most_m_distinct_in_sample_of_x_from_n_double(
    m: usize,
    x: usize,
    n: usize,
    sr: &[Vec<Double>],
) -> f64 {
    let mut p = dd![1.0];
    for u in m + 1..=x {
        let mut z = dd![sr[x][u]];
        for _ in 0..x {
            z *= dd![u as f64] / dd![n as f64];
        }
        for v in 1..=u {
            z *= dd![(n - v + 1) as f64] / dd![(u - v + 1) as f64];
        }
        p -= z;
    }
    if p < dd![0.0] {
        p = dd![0.0];
    }
    f64::from(p)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
#[allow(clippy::too_many_arguments)]
pub(crate) fn join_one<'a>(
    k1: usize,
    k2: usize,
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    to_bc: &'a HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    pot: &mut Vec<PotentialJoin>,
    refdata: &RefData,
) -> bool {
    let is_bcr = ctl.is_bcr();
    // Do not merge onesies or foursies with anything.  Deferred until later.
    // Note that perhaps some foursies should be declared doublets and deleted.
    // Note onesies merging above is turned off so this appears to be moot.

    let (info1, info2) = (&info[k1], &info[k2]);
    let (clono1, clono2) = (info1.clonotype_index, info2.clonotype_index);
    let chains1 = exact_clonotypes[clono1].share.len();
    let chains2 = exact_clonotypes[clono2].share.len();
    if !(2..=3).contains(&chains1) || !(2..=3).contains(&chains2) {
        return false;
    }
    // NEED FOR THIS SEEMS LIKE A BUG:
    if info1.vs.len() == 1 || info2.vs.len() == 4 {
        return false;
    }
    if info1.vs.len() > 2 {
        return false;
    }

    // Require that CDR3s have the same length.  Ugly.
    // First part should be a tautology.

    let (x1, x2) = (&info[k1].cdr3s, &info[k2].cdr3s);
    if x1.len() != x2.len() {
        return false;
    }
    for i in 0..x1.len() {
        if x1[i].len() != x2[i].len() {
            return false;
        }
    }

    // Put identity filter on CDR3s for BCR.

    if is_bcr {
        let (x1, x2) = (&info[k1].cdr3s, &info[k2].cdr3s);
        let mut cd = 0;
        let mut total = 0;
        for z in 0..2 {
            if x1[z].len() != x2[z].len() {
                return false;
            }
            for m in 0..x1[z].len() {
                if x1[z].as_bytes()[m] != x2[z].as_bytes()[m] {
                    cd += 1;
                }
            }
            total += x1[z].len();
        }
        if cd as f64 / total as f64 > 1.0 - JOIN_CDR3_IDENT / 100.0 {
            return false;
        }
    }

    // Compute number of differences.  The default behavior is that this is applied only to TCR.

    let (x1, x2) = (&info[k1].cdr3s, &info[k2].cdr3s);
    if !is_bcr || MAX_DIFFS < 1_000_000 {
        let mut diffs = 0_usize;
        for x in 0..info[k1].lens.len() {
            if !info[k1].has_del[x] && !info[k2].has_del[x] {
                // A great deal of time is spent in the call to ndiffs.  Notes on this:
                // 1. It is slower than if the computation is done outside
                //    the ndiffs function.  This is mysterious but must have something to
                //    do with the storage of the 256-byte lookup table.
                // 2. Adding #[inline(always)] in front of the ndiffs function definition
                //    doesn't help.
                // 3. Adding a bounds test for diffs > ctl.heur.max_diffs inside the ndiffs
                //    function doesn't help, whether placed in the inner loop or the other
                //    loop.
                diffs += ndiffs(&info[k1].tigsp[x], &info[k2].tigsp[x]);
            } else {
                for j in 0..info[k1].tigs[x].len() {
                    if info[k1].tigs[x][j] != info[k2].tigs[x][j] {
                        diffs += 1;
                    }
                }
            }
        }
        if diffs > MAX_DIFFS {
            return false;
        }
        if !is_bcr && diffs > 5 {
            return false;
        }
    }

    // Compute junction diffs.

    let mut cd = 0_isize;
    for l in 0..x1.len() {
        for m in 0..x1[l].len() {
            if x1[l].as_bytes()[m] != x2[l].as_bytes()[m] {
                cd += 1;
            }
        }
    }

    // Cap CDR3 diffs for TCR or as requested.

    if (MAX_CDR3_DIFFS < 1000 || !is_bcr) && (cd > MAX_CDR3_DIFFS as isize || (!is_bcr && cd > 0)) {
        return false;
    }

    // Unless MIX_DONORS specified, do not join across donors.
    // And test for error.
    //
    // WARNING!  There are actually two cases: where an individual exact subclonotype
    // itself crosses donors, and where we cross donors in making a join.  Note that
    // the former case is most improbable, unless there is cross-sample contamination.
    // And if that did happen, the output would be confusing and might have a greatly
    // exaggerated number of fails.

    let (mut donors1, mut donors2) = (Vec::<usize>::new(), Vec::<usize>::new());
    let ex1 = &exact_clonotypes[info[k1].clonotype_index];
    let ex2 = &exact_clonotypes[info[k2].clonotype_index];
    for j in 0..ex1.clones.len() {
        if ex1.clones[j][0].donor_index.is_some() {
            donors1.push(ex1.clones[j][0].donor_index.unwrap());
        }
    }
    for j in 0..ex2.clones.len() {
        if ex2.clones[j][0].donor_index.is_some() {
            donors2.push(ex2.clones[j][0].donor_index.unwrap());
        }
    }
    unique_sort(&mut donors1);
    unique_sort(&mut donors2);
    if !ctl.cr_opt.mix_donors && !donors1.is_empty() && !donors2.is_empty() && donors1 != donors2 {
        return false;
    }
    let err = donors1 != donors2 || donors1.len() != 1 || donors2.len() != 1;

    // Analyze the two clonotypes versus the reference.  First traverse the reference
    // sequences.  Either we use the references for k1 or the references for k2, but
    // these are nearly always the same.

    let mut nrefs = 1;
    for m in 0..2 {
        if info[k1].vs[m] != info[k2].vs[m] || info[k1].js[m] != info[k2].js[m] {
            nrefs = 2;
        }
    }
    let mut shares = vec![0; nrefs]; // shared mutations from reference
    let mut shares1 = vec![0; nrefs];
    let mut shares2 = vec![0; nrefs];
    let mut indeps = vec![0; nrefs]; // independent mutations from reference
    let mut total = vec![[0usize; 2]; nrefs]; // total differences from reference
    for u in 0..nrefs {
        let k = if u == 0 { k1 } else { k2 };

        // Traverse the chains in the clonotype.

        let nchains = info[k1].lens.len();
        for m in 0..nchains {
            let (tig1, tig2) = (&info[k1].tigs[m], &info[k2].tigs[m]);

            // Traverse the two segments (V and J).

            for si in 0..2 {
                let seg = if si == 0 {
                    &info[k].vs[m]
                } else {
                    &info[k].js[m]
                };
                let ref_trim = if si == 1 { REF_J_TRIM } else { REF_V_TRIM };
                for p in 0..seg.len() - ref_trim {
                    let (t1, t2);
                    let r;
                    if si == 0 {
                        // Ugly bailout arising very rarely if the two reference
                        // sequences have different lengths.
                        if p >= tig1.len() || p >= tig2.len() {
                            return false;
                        }
                        t1 = tig1[p];
                        t2 = tig2[p];
                        // r = seg.get(p);
                        let rx = seg.get(p);
                        if rx == 0 {
                            r = b'A';
                        } else if rx == 1 {
                            r = b'C';
                        } else if rx == 2 {
                            r = b'G';
                        } else {
                            r = b'T';
                        }
                    } else {
                        t1 = tig1[tig1.len() - p - 1];
                        t2 = tig2[tig2.len() - p - 1];
                        // r = seg.get( seg.len() - p - 1 );
                        let rx = seg.get(seg.len() - p - 1);
                        if rx == 0 {
                            r = b'A';
                        } else if rx == 1 {
                            r = b'C';
                        } else if rx == 2 {
                            r = b'G';
                        } else {
                            r = b'T';
                        }
                    }
                    if t1 == t2 && t1 != r {
                        shares[u] += 1;
                        if m == 1 {
                            shares1[u] += 1;
                        } else {
                            shares2[u] += 1;
                        }
                    } else if (t1 == r && t2 != r) || (t2 == r && t1 != r) {
                        indeps[u] += 1;
                    } else if t1 != r && t2 != r {
                        indeps[u] += 2;
                    }
                    if t1 != r {
                        total[u][0] += 1;
                    }
                    if t2 != r {
                        total[u][1] += 1;
                    }
                }
            }
        }
    }

    // Don't allow different references if one is strongly favored.
    // (not documented)

    if nrefs == 2 {
        for m in 0..2 {
            if total[0][m].abs_diff(total[1][m]) > MAX_DEGREDATION {
                return false;
            }
        }
    }

    // Another test for acceptable join.  (not fully documented)

    let min_shares = shares.iter().min().unwrap();
    let min_indeps = indeps.iter().min().unwrap();

    // Reject if barcode overlap. (not documented)

    let (mut bcs1, mut bcs2) = (Vec::<&'a str>::new(), Vec::<&'a str>::new());
    for origin in &info[k1].origin {
        bcs1.extend(
            to_bc[&(*origin, info[k1].clonotype_index)]
                .iter()
                .map(String::as_str),
        );
    }
    for origin in &info[k2].origin {
        bcs2.extend(
            to_bc[&(*origin, info[k2].clonotype_index)]
                .iter()
                .map(String::as_str),
        );
    }
    unique_sort(&mut bcs1);
    unique_sort(&mut bcs2);
    if meet(&bcs1, &bcs2) {
        return false;
    }

    // Test for concentration of SHM in the junction regions.
    let cdr3_mult = 5.0;
    if cd as f64 >= cdr3_mult * std::cmp::max(1, *min_indeps) as f64 {
        return false;
    }

    // Do not merge cells if they were assigned different light chain constant regions.
    // Unless cd = 0.
    for i in 0..info[k1].cdr3s.len() {
        let (j1, j2) = (info[k1].exact_cols[i], info[k2].exact_cols[i]);
        if !ex1.share[j1].left
            && ex1.share[j1].c_ref_id.is_some()
            && ex2.share[j2].c_ref_id.is_some()
            && ex1.share[j1].c_ref_id.unwrap() != ex2.share[j2].c_ref_id.unwrap()
            && cd > 0
        {
            return false;
        }
    }

    // Estimate the probability p1 that drawing k = min_indeps + 2 * min_shares
    // objects from n = 3 * (sum of VJ contig lengths) yields d = min_shares or
    // more duplicates.

    let n = 3 * (info[k1].tigs[0].len() + info[k1].tigs[1].len());
    let k = *min_indeps + 2 * *min_shares;
    let d = *min_shares;
    let p1 = p_at_most_m_distinct_in_sample_of_x_from_n_double((k - d) as usize, k as usize, n, sr);
    assert!(!p1.is_infinite()); // TODO: IS THIS SAFE?

    // Multiply by 80^cd.  This is sum( choose(3cn, m), m = 0..=cd ).
    // Changed to take into account CDR3 length.

    let mult = {
        let mut cd1 = 0;
        let n1 = x1[0].len();
        for m in 0..x1[0].len() {
            if x1[0].as_bytes()[m] != x2[0].as_bytes()[m] {
                cd1 += 1;
            }
        }
        let mut cd2 = 0;
        let n2 = x1[1].len();
        for m in 0..x1[1].len() {
            if x1[1].as_bytes()[m] != x2[1].as_bytes()[m] {
                cd2 += 1;
            }
        }
        let cdx = CDR3_NORMAL_LEN;
        let mult_pow: f64 = 80.0;
        mult_pow.powf(cdx as f64 * cd1 as f64 / n1 as f64)
            * mult_pow.powf(cdx as f64 * cd2 as f64 / n2 as f64)
    };

    // Compute score.

    let score = p1 * mult;

    // Apply JUN_SHARE.

    let mut accept = false;
    if COMP_FILT < 1_000_000
        && score > MAX_SCORE
        && *min_shares < AUTO_SHARE as isize
        && (COMP_FILT_BOUND == 0 || *min_indeps as usize <= COMP_FILT_BOUND)
        && ex1.share.len() == 2
        && ex2.share.len() == 2
        && ex1.share[0].left != ex1.share[1].left
    {
        let h1 = info[k1].exact_cols[0];
        let h2 = info[k2].exact_cols[0];
        let comp = min(ex1.share[h1].jun.hcomp, ex2.share[h2].jun.hcomp);
        if comp as isize - cd >= COMP_FILT as isize {
            /*
            println!("\nwould accept");
            println!("cdr3: {}", ex1.share[h1].cdr3_aa);
            println!("cdr3: {}", ex2.share[h2].cdr3_aa);
            */
            accept = true;
        }
    }

    // Threshold on score.

    if !accept && score > MAX_SCORE && *min_shares < AUTO_SHARE as isize {
        return false;
    }

    // If V gene names are different (after removing trailing *...), and either
    // • V gene reference sequences are different, after truncation on right to the same length
    // • or 5' UTR reference sequences are different, after truncation on left to the same length,
    // then the join is rejected.

    for i in 0..info[k1].cdr3s.len() {
        let (j1, j2) = (info[k1].exact_cols[i], info[k2].exact_cols[i]);
        let (x1, x2) = (&ex1.share[j1], &ex2.share[j2]);
        let (v1, v2) = (x1.v_ref_id, x2.v_ref_id);
        let (mut n1, mut n2) = (refdata.name[v1].clone(), refdata.name[v2].clone());
        if n1.contains('*') {
            n1 = n1.before("*").to_string();
        }
        if n2.contains('*') {
            n2 = n2.before("*").to_string();
        }
        if n1 != n2 {
            let (y1, y2) = (&refdata.refs[v1], &refdata.refs[v2]);
            if y1.len() == y2.len() {
                if y1 != y2 {
                    return false;
                }
            } else {
                let n = min(y1.len(), y2.len());
                for m in 0..n {
                    if y1.get(m) != y2.get(m) {
                        return false;
                    }
                }
            }
            let (u1, u2) = (x1.u_ref_id, x2.u_ref_id);
            if let (Some(u1), Some(u2)) = (u1, u2) {
                let (x1, x2) = (&refdata.refs[u1], &refdata.refs[u2]);
                let n = min(x1.len(), x2.len());
                for m in 0..n {
                    if x1.get(x1.len() - 1 - m) != x2.get(x2.len() - 1 - m) {
                        return false;
                    }
                }
            }
        }
    }

    // Require
    // percent heavy chain nuke identity on FWR1
    // minus
    // percent heavy chain nuke identity on CDR12
    // is less than 20.

    let nchains = info[k1].lens.len();
    let (mut fwr1_len, mut cdr1_len, mut cdr2_len) = (0, 0, 0);
    let (mut fwr1_diffs, mut cdr1_diffs, mut cdr2_diffs) = (0, 0, 0);
    for m in 0..nchains {
        let (j1, j2) = (info[k1].exact_cols[m], info[k2].exact_cols[m]);
        let (x1, x2) = (&ex1.share[j1], &ex2.share[j2]);
        if x1.left {
            if x1.cdr1_start.is_some() && x2.cdr1_start.is_some() {
                let fr1_start1 = x1.fr1_start;
                let fr1_stop1 = x1.cdr1_start.unwrap();
                let fr1_start2 = x2.fr1_start;
                let fr1_stop2 = x2.cdr1_start.unwrap();
                let len = fr1_stop1 - fr1_start1;
                if fr1_stop2 - fr1_start2 == len {
                    let mut diffs = 0;
                    for p in 0..len {
                        if x1.seq_del_amino[p + fr1_start1] != x2.seq_del_amino[p + fr1_start2] {
                            diffs += 1;
                        }
                    }
                    fwr1_len = len;
                    fwr1_diffs = diffs;
                }
            }
            if x1.cdr1_start.is_some()
                && x1.fr2_start.is_some()
                && x2.cdr1_start.is_some()
                && x2.fr2_start.is_some()
            {
                let cdr1_start1 = x1.cdr1_start.unwrap();
                let cdr1_stop1 = x1.fr2_start.unwrap();
                let cdr1_start2 = x2.cdr1_start.unwrap();
                let cdr1_stop2 = x2.fr2_start.unwrap();
                // there are some problematic annotations for which this is not upheld
                if cdr1_start1 <= cdr1_stop1 && cdr1_start2 <= cdr1_stop2 {
                    let len = cdr1_stop1 - cdr1_start1;
                    if cdr1_stop2 - cdr1_start2 == len {
                        let mut diffs = 0;
                        for p in 0..len {
                            if x1.seq_del_amino[p + cdr1_start1]
                                != x2.seq_del_amino[p + cdr1_start2]
                            {
                                diffs += 1;
                            }
                        }
                        cdr1_len = len;
                        cdr1_diffs = diffs;
                    }
                }
            }
            if x1.cdr2_start.is_some()
                && x1.fr3_start.is_some()
                && x2.cdr2_start.is_some()
                && x2.fr3_start.is_some()
            {
                let cdr2_start1 = x1.cdr2_start.unwrap();
                let cdr2_stop1 = x1.fr3_start.unwrap();
                let cdr2_start2 = x2.cdr2_start.unwrap();
                let cdr2_stop2 = x2.fr3_start.unwrap();
                // this was violated once when using IMGT reference
                if cdr2_start1 <= cdr2_stop1 && cdr2_start2 <= cdr2_stop2 {
                    let len = cdr2_stop1 - cdr2_start1;
                    if cdr2_stop2 - cdr2_start2 == len {
                        let mut diffs = 0;
                        for p in 0..len {
                            if x1.seq_del_amino[p + cdr2_start1]
                                != x2.seq_del_amino[p + cdr2_start2]
                            {
                                diffs += 1;
                            }
                        }
                        cdr2_len = len;
                        cdr2_diffs = diffs;
                    }
                }
            }
        }
    }
    if fwr1_len > 0 && cdr1_len > 0 && cdr2_len > 0 {
        let len = fwr1_len;
        let diffs = fwr1_diffs;
        let fwr1_identity = 100.0 * (len - diffs) as f64 / len as f64;
        let len = cdr1_len + cdr2_len;
        let diffs = cdr1_diffs + cdr2_diffs;
        let cdr12_identity = 100.0 * (len - diffs) as f64 / len as f64;
        if fwr1_identity - cdr12_identity >= FWR1_CDR12_DELTA {
            return false;
        }
    }

    // Save potential joins.  Note that this jacks up memory usage significantly,
    // so it would likely be more efficient to duplicate some of the computations
    // during the analysis phase.

    pot.push(PotentialJoin {
        k1,
        k2,
        cd,
        shares,
        err,
    });
    true
}

// Join algorithm parameters.
const MAX_SCORE: f64 = 100_000.0;
const MAX_CDR3_DIFFS: usize = 1000;
const JOIN_CDR3_IDENT: f64 = 85.0;
const FWR1_CDR12_DELTA: f64 = 20.0;
const CDR3_NORMAL_LEN: usize = 42;
const AUTO_SHARE: usize = 15;
pub(crate) const COMP_FILT: usize = 8;
const COMP_FILT_BOUND: usize = 80;

// Clonotyping heuristics.
const MAX_DIFFS: usize = 1_000_000;
const MAX_DEGREDATION: usize = 2;
pub(crate) const REF_V_TRIM: usize = 15;
pub(crate) const REF_J_TRIM: usize = 15;
