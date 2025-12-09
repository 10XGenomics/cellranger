// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use crate::core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use crate::core::enclone_structs::CandidateClonotype;
use crate::process::define_mat::{define_mat, setup_define_mat};
use equiv::EquivRel;
use qd::Double;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::{bin_position, unique_sort};

/// Check for disjoint candidate clonotypes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn split_candidate_clonotypes(
    candidate_clonotypes: &mut Vec<CandidateClonotype>,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    refdata: &RefData,
) {
    let mut candidate_clonotypes2 = Vec::<Vec<i32>>::new();
    for cc in candidate_clonotypes.iter() {
        let (od, exacts) = setup_define_mat(cc, info);
        let mat = define_mat(
            to_bc,
            sr,
            ctl,
            exact_clonotypes,
            &exacts,
            &od,
            info,
            raw_joins,
            refdata,
        );
        let cols = mat.len();

        // Define map of indices into exacts.

        let nexacts = exacts.len();
        let mut to_exacts = HashMap::<usize, usize>::with_capacity(nexacts);
        for (u, &e) in exacts.iter().enumerate() {
            to_exacts.insert(e, u);
        }

        // Get the info indices corresponding to this clonotype.

        let mut infos = Vec::<usize>::new();
        for &oi in cc {
            infos.push(oi as usize);
        }

        // Define map of exacts to infos.

        let mut to_infos = vec![Vec::<usize>::new(); nexacts];
        for (i, &infoi) in infos.iter().enumerate() {
            let u = to_exacts[&info[infoi].clonotype_index];
            to_infos[u].push(i);
        }

        // Determine which columns are "left", meaning IGH or TRB.

        let mut left = vec![false; cols];
        for m in 0..cols {
            for u in 0..mat[0].len() {
                if mat[m][u].is_some() {
                    let c = mat[m][u].unwrap();
                    let ex = &exact_clonotypes[exacts[u]];
                    if ex.share[c].left {
                        left[m] = true;
                    }
                    break;
                }
            }
        }

        // Determine which pairs of configurations share both chain types, and if so, call
        // them joined.

        let mut matu = Vec::<Vec<Option<usize>>>::with_capacity(nexacts);
        for u in 0..nexacts {
            let mut m = Vec::<Option<usize>>::with_capacity(cols);
            for mm in mat.iter().take(cols) {
                m.push(mm[u]);
            }
            matu.push(m);
        }
        unique_sort(&mut matu);
        let mut eqm = vec![vec![false; matu.len()]; matu.len()];
        for (mj1, eqm) in matu.iter().zip(eqm.iter_mut()) {
            for (mj2, eqm) in matu.iter().zip(eqm.iter_mut()) {
                let (mut l, mut r) = (false, false);
                for ((&mm1, &mm2), &ll) in mj1.iter().zip(mj2.iter()).zip(left.iter()).take(cols) {
                    if mm1.is_some() && mm2.is_some() {
                        if ll {
                            l = true;
                        } else {
                            r = true;
                        }
                    }
                }
                if l && r {
                    *eqm = true;
                }
            }
        }

        // Propagate this to an equivalence relation on the candidate clonotype elements.

        let mut eqx = EquivRel::new(cc.len() as u32);
        let mut lists = vec![Vec::<usize>::new(); matu.len()];
        for u in 0..nexacts {
            let mut m = Vec::<Option<usize>>::with_capacity(cols);
            for mat in mat.iter().take(cols) {
                m.push(mat[u]);
            }
            lists[bin_position(&matu, &m) as usize].push(u);
        }
        for (l1, eqm) in lists.iter().zip(eqm.into_iter()) {
            for (l2, eqm) in lists.iter().zip(eqm.into_iter()) {
                if eqm {
                    let u1 = l1[0];
                    for &u2 in l2 {
                        for &i1 in &to_infos[u1] {
                            for &i2 in &to_infos[u2] {
                                eqx.join(i1, i2);
                            }
                        }
                    }
                }
            }
        }

        // Join onesies where possible.  This should probably be more efficient.

        for (&e1, info1) in exacts.iter().zip(to_infos.iter()).take(nexacts) {
            let ex1 = &exact_clonotypes[e1];
            if ex1.share.len() == 1 {
                let mut is = Vec::<usize>::new();
                for (&e2, info2) in exacts.iter().take(nexacts).zip(to_infos.iter()) {
                    let ex2 = &exact_clonotypes[e2];
                    if ex2.share.len() == 1 {
                        if ex1.share[0].seq == ex2.share[0].seq {
                            eqx.join(info1[0], info2[0]);
                        }
                    } else {
                        for j in 0..ex2.share.len() {
                            if ex2.share[j].seq == ex1.share[0].seq {
                                is.push(info2[0]);
                            }
                        }
                    }
                }
                let mut rs = Vec::<usize>::new();
                for &ij in &is {
                    rs.push(eqx.set_id(ij));
                }
                unique_sort(&mut rs);
                if rs.len() == 1 {
                    eqx.join(info1[0], is[0]);
                }
            }
        }

        // Divide the candidate clonotype if needed.

        if eqx.n_sets() == 1 {
            candidate_clonotypes2.push(cc.clone());
        } else {
            for ox_iter in eqx.all_sets() {
                let mut o2 = Vec::<i32>::new();
                for ko in ox_iter {
                    o2.push(cc[ko]);
                }
                candidate_clonotypes2.push(o2);
            }
        }
    }
    *candidate_clonotypes = candidate_clonotypes2;
}
