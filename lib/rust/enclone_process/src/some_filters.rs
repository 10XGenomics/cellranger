// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use crate::core::enclone_structs::{BarcodeFates, CandidateClonotype};
use crate::core::join_one::{REF_J_TRIM, REF_V_TRIM};
use crate::process::define_column_info::define_column_info;
use crate::process::define_mat::{define_mat, setup_define_mat};
use qd::Double;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use vdj_ann::refx::RefData;
use vector_utils::{erase_if, next_diff1_2, unique_sort};

/// Given a signature s having at least two chains, if the total cells in the two-chain
/// signatures that are different from it but share a chain with it is at least 20 times
/// greater, delete s.
#[allow(clippy::too_many_arguments)]
pub(crate) fn signature_filter(
    candidate_clonotypes: &mut Vec<CandidateClonotype>,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    fate: &mut [BarcodeFates],
    refdata: &RefData,
) {
    // Note duplication of calls to define_mat with other code.  This is expensive.
    const SIG_MULT: usize = 20;
    let mut results = Vec::<(usize, Vec<(usize, String, BarcodeFate)>, Vec<usize>)>::new();
    for i in 0..candidate_clonotypes.len() {
        results.push((i, Vec::new(), Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let cc = candidate_clonotypes[i].clone();
        let (od, exacts) = setup_define_mat(&cc, info);
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

        // Find all the signatures and cell counts associated to each.

        let mut freq = Vec::<(usize, Vec<usize>)>::new();
        {
            let mut types = Vec::<(Vec<usize>, usize)>::new();
            for (u, &e) in exacts.iter().enumerate() {
                let mut t = Vec::<usize>::new();
                for (col, m) in mat.iter().enumerate() {
                    if m[u].is_some() {
                        t.push(col);
                    }
                }
                if t.len() >= 2 {
                    types.push((t, exact_clonotypes[e].ncells()));
                }
            }
            types.sort();
            let mut i = 0;
            while i < types.len() {
                let j = next_diff1_2(&types, i);
                let mut mult = 0;
                for t in &types[i..j] {
                    mult += t.1;
                }
                freq.push((mult, types[i].0.clone()));
                i = j;
            }
        }
        /*
        let mut msg = "\nfrequencies:\n".to_string(); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        use itertools::Itertools; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        for i in 0..freq.len() { // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            msg += &mut format!("{} ==> {}\n", freq[i].0, freq[i].1.iter().format(",")); // XXX
        } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        */

        // Decide which signatures to delete.

        let mut dels = HashSet::<Vec<usize>>::new();
        for i in 0..freq.len() {
            let mut n2 = 0;
            for j in 0..freq.len() {
                if j != i && freq[j].1.len() == 2 {
                    let mut share = false;
                    for x in &freq[j].1 {
                        if freq[i].1.contains(x) {
                            share = true;
                        }
                    }
                    if share {
                        n2 += freq[j].0;
                    }
                }
            }
            if n2 > SIG_MULT * freq[i].0 {
                dels.insert(freq[i].1.clone());
                /*
                msg += &mut format!("delete {}\n", freq[i].1.iter().format(",")); // XXXXXXXXXX
                */
            }
        }
        /*
        if dels.len() > 0 { println!("{}", msg); } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        */
        for u in 0..exacts.len() {
            let mut t = Vec::<usize>::new();
            for (col, m) in mat.iter().enumerate() {
                if m[u].is_some() {
                    t.push(col);
                }
            }
            if dels.contains(&t) {
                res.2.push(exacts[u]);
                let ex = &exact_clonotypes[exacts[u]];
                for i in 0..ex.ncells() {
                    res.1.push((
                        ex.clones[i][0].dataset_index,
                        ex.clones[i][0].barcode.clone(),
                        BarcodeFate::Signature,
                    ));
                }
            }
        }
    });
    let mut to_delete = vec![false; exact_clonotypes.len()];
    for i in 0..results.len() {
        for j in 0..results[i].1.len() {
            fate[results[i].1[j].0].insert(results[i].1[j].1.clone(), results[i].1[j].2);
        }
        for j in 0..results[i].2.len() {
            to_delete[results[i].2[j]] = true;
        }
    }
    let mut candidate_clonotypes2 = Vec::<Vec<i32>>::new();
    for cc in candidate_clonotypes.iter() {
        let mut cc = cc.clone();
        let mut del = vec![false; cc.len()];
        for j in 0..cc.len() {
            let id = info[cc[j] as usize].clonotype_index;
            if to_delete[id] {
                del[j] = true;
            }
        }
        erase_if(&mut cc, &del);
        candidate_clonotypes2.push(cc);
    }
    *candidate_clonotypes = candidate_clonotypes2;
}

/// Find and mark for deletion exact subclonotypes having a variant base in V..J that,
/// accounting for all the cells in all the exact subclonotypes, never occurs as Q60
/// doesn't occur as Q40 twice, and disagrees with the reference.
#[allow(clippy::too_many_arguments)]
pub(crate) fn qual_filter(
    candidate_clonotypes: &mut [CandidateClonotype],
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    fate: &mut [BarcodeFates],
    refdata: &RefData,
) {
    let mut results = Vec::<(usize, Vec<(usize, String, BarcodeFate)>, Vec<usize>)>::new();
    for i in 0..candidate_clonotypes.len() {
        results.push((i, Vec::new(), Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let cc = candidate_clonotypes[i].clone();
        let (od, exacts) = setup_define_mat(&cc, info);
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
        let rsi = define_column_info(&exacts, exact_clonotypes, mat);

        // Create vars, copied from vars_and_shares.

        let mut vars = Vec::<Vec<usize>>::new();
        for cx in 0..rsi.mat.len() {
            let mut n = 0;
            for z in 0..rsi.seqss[cx].len() {
                n = max(n, rsi.seqss[cx][z].len());
            }
            let mut v = Vec::<usize>::new();
            for p in 0..n {
                let mut bases = Vec::<u8>::new();
                for s in 0..rsi.seqss[cx].len() {
                    if p >= rsi.seqss[cx][s].len() {
                        continue;
                    }
                    bases.push(rsi.seqss[cx][s][p]);
                }
                unique_sort(&mut bases);
                if bases.len() > 1 {
                    v.push(p);
                }
            }
            vars.push(v);
        }

        // Proceed.

        // (column, pos, base, qual, row)
        let mut vquals = Vec::<(usize, usize, u8, u8, usize)>::new();
        for u in 0..exacts.len() {
            let clonotype_id = exacts[u];
            let ex = &exact_clonotypes[clonotype_id];
            #[allow(clippy::needless_range_loop)]
            for col in 0..rsi.mat.len() {
                if let Some(m) = rsi.mat[col][u] {
                    if ex.share[m].annv.len() > 1 {
                        continue;
                    }
                    let n = ex.share[m].seq_del.len();
                    let vref = &exact_clonotypes[exacts[u]].share[m].vs.to_ascii_vec();
                    let jref = &exact_clonotypes[exacts[u]].share[m].js.to_ascii_vec();
                    for z in 0..vars[col].len() {
                        let p = vars[col][z];
                        let b = ex.share[m].seq_del[p];
                        let mut refdiff = false;
                        if p < vref.len() - REF_V_TRIM && b != vref[p] {
                            refdiff = true;
                        }
                        if p >= n - (jref.len() - REF_J_TRIM) && b != jref[jref.len() - (n - p)] {
                            refdiff = true;
                        }
                        if refdiff {
                            for j in 0..ex.clones.len() {
                                let qual = ex.clones[j][m].quals[p];
                                vquals.push((col, p, b, qual, u));
                            }
                        }
                    }
                }
            }
        }
        vquals.sort_unstable();
        let mut j = 0;
        while j < vquals.len() {
            let mut k = j + 1;
            while k < vquals.len() {
                if vquals[k].0 != vquals[j].0
                    || vquals[k].1 != vquals[j].1
                    || vquals[k].2 != vquals[j].2
                {
                    break;
                }
                k += 1;
            }
            let mut q60 = false;
            let mut q40 = 0;
            for v in &vquals[j..k] {
                if v.3 >= 60 {
                    q60 = true;
                } else if v.3 >= 40 {
                    q40 += 1;
                }
            }
            if !q60 && q40 < 2 {
                let u = vquals[j].4;
                res.2.push(exacts[u]);

                let ex = &exact_clonotypes[exacts[u]];
                for i in 0..ex.ncells() {
                    res.1.push((
                        ex.clones[i][0].dataset_index,
                        ex.clones[i][0].barcode.clone(),
                        BarcodeFate::Qual,
                    ));
                }
            }
            j = k;
        }
    });
    let mut to_delete = vec![false; exact_clonotypes.len()];
    let mut dels = Vec::<i32>::new();
    for i in 0..results.len() {
        for j in 0..results[i].1.len() {
            fate[results[i].1[j].0].insert(results[i].1[j].1.clone(), results[i].1[j].2);
        }
        for x in &results[i].2 {
            to_delete[*x] = true;
        }
    }
    dels.sort_unstable();
    for cc in candidate_clonotypes.iter_mut() {
        let mut del = vec![false; cc.len()];
        for (&oj, d) in cc.iter().zip(del.iter_mut()) {
            let id = info[oj as usize].clonotype_index;
            if to_delete[id] {
                *d = true;
            }
        }
        erase_if(cc, &del);
    }
}
