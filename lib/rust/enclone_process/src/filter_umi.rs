//! Filter B cells based on UMI counts.
// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use crate::core::enclone_structs::{BarcodeFates, CandidateClonotype};
use std::cmp::max;
use vector_utils::{erase_if, next_diff1_5, reverse_sort};

/// Return the binomial CDF.
/// binomial_sum( n, k, p ): return sum_{i=0..k} choose(n,i) * p^i * (1-p)^(n-i)
/// No attempt has been made to make this efficient or to pay attention to
/// accuracy or overflow problems.
pub(crate) fn binomial_cdf(p: f64, n: usize, k: usize) -> f64 {
    assert!(n >= 1);
    assert!(k <= n);
    let mut sum = 0.0;
    let mut choose = 1.0;
    for _ in 0..n {
        choose *= 1.0 - p;
    }
    let q = p / (1.0 - p);
    for i in 0..=k {
        sum += choose;
        choose *= (n - i) as f64;
        choose /= (i + 1) as f64;
        choose *= q;
    }
    sum
}

/// For B cells, filter based on UMI counts.  More details in heuristics.html.
/// Find all clonotypes having one cell which has two chains,
/// one heavy and one light.  Get the sum of the chain UMI counts for this cell.
///
/// For each cell, let umish be the umi count for its heavy chain having the most umis, and
/// similarly define umisl.  Let umitot = umish + umisl.
///
/// If every cell in a clonotype would have been deleted, first find the exact subclonotype for
/// which the sum of its umitot values is greatest, and then in it, find the cell having
/// highest umitot value.  Protect this cell, so long as it has at least two chains.
pub(crate) fn filter_umi(
    ctl: &EncloneControl,
    exact_clonotypes: &mut [ExactClonotype],
    info: &[CloneInfo],
    candidate_clonotypes: Vec<CandidateClonotype>,
    fate: &mut [BarcodeFates],
) -> Vec<Vec<i32>> {
    let mut umis = vec![Vec::<usize>::new(); ctl.origin_info.len()];
    for cc in &candidate_clonotypes {
        if cc.len() == 1 {
            let x: &CloneInfo = &info[cc[0] as usize];
            let ex = &exact_clonotypes[x.clonotype_index];
            if ex.ncells() == 1 && ex.share.len() == 2 && ex.share[0].left != ex.share[1].left {
                umis[ex.clones[0][0].dataset_index]
                    .push(ex.clones[0][0].umi_count + ex.clones[0][1].umi_count);
            }
        }
    }
    let mut nu = vec![0; ctl.origin_info.len()];
    let mut umin = vec![0.0; ctl.origin_info.len()];
    for l in 0..ctl.origin_info.len() {
        umis[l].sort_unstable();
        nu[l] = umis[l].len();
        if nu[l] > 0 {
            let n10 = umis[l][nu[l] / 10] as f64;
            let n50 = umis[l][nu[l] / 2] as f64;
            umin[l] = n10.min(n50 - (4.0 * n50.sqrt()));
        }
    }

    // if ctl.clono_filt_opt_def.umi_filt || ctl.clono_filt_opt_def.umi_filt_mark {
    const MIN_BASELINE_CELLS: usize = 20;
    let candidate_clonotypes: Vec<_> = candidate_clonotypes
        .into_iter()
        .filter_map(|mut cc| {
            let ncells = cc
                .iter()
                .map(|x| exact_clonotypes[info[*x as usize].clonotype_index].ncells())
                .sum();

            if ncells > 1 {
                let (mut best_ex, mut best_ex_sum) = (0, 0);
                let (mut best_cell, mut best_cell_count) = (0, 0);
                let mut baselined = true;
                let mut protected = false;
                {
                    let mut nbads = 0;

                    for j in 0..cc.len() {
                        let x: &CloneInfo = &info[cc[j] as usize];
                        let ex = &mut exact_clonotypes[x.clonotype_index];
                        let mut ex_sum = 0;
                        for clone in &ex.clones {
                            let li = clone[0].dataset_index;
                            if nu[li] >= MIN_BASELINE_CELLS {
                                let (mut umish, mut umisl) = (0, 0);
                                for (s, c) in ex.share.iter().zip(clone.iter()) {
                                    if s.left {
                                        umish = max(umish, c.umi_count);
                                    } else {
                                        umisl = max(umish, c.umi_count);
                                    }
                                }
                                let umitot = umish + umisl;
                                ex_sum += umitot;
                                if (umitot as f64) < umin[li] {
                                    nbads += 1;
                                }
                            } else {
                                baselined = false;
                            }
                        }
                        if ex_sum > best_ex_sum {
                            best_ex = j;
                            best_ex_sum = ex_sum;
                        }
                    }

                    if nbads == 0 {
                        protected = true;
                    } else {
                        let p = 0.1;
                        let bound = 0.01;

                        // Find probability of observing nbads or more events of probability
                        // p in a sample of size ncells, and if that is at least bound,
                        // don't delete any cells (except onesies).
                        if binomial_cdf(1.0 - p, ncells, ncells - nbads) >= bound {
                            protected = true;
                        }
                    }
                }
                {
                    for j in 0..cc.len() {
                        let x: &CloneInfo = &info[cc[j] as usize];
                        let ex = &exact_clonotypes[x.clonotype_index];
                        for (k, clone) in ex.clones.iter().enumerate() {
                            let li = clone[0].dataset_index;
                            if nu[li] >= MIN_BASELINE_CELLS {
                                let (mut umish, mut umisl) = (0, 0);
                                for (s, c) in ex.share.iter().zip(clone.iter()) {
                                    if s.left {
                                        umish = max(umish, c.umi_count);
                                    } else {
                                        umisl = max(umish, c.umi_count);
                                    }
                                }
                                let umitot = umish + umisl;
                                if j == best_ex && umitot > best_cell_count && ex.share.len() > 1 {
                                    best_cell = k;
                                    best_cell_count = umitot;
                                }
                            } else {
                                baselined = false;
                            }
                        }
                    }
                }
                {
                    for j in 0..cc.len() {
                        let x: &CloneInfo = &info[cc[j] as usize];
                        let ex = &mut exact_clonotypes[x.clonotype_index];
                        let mut to_delete = vec![false; ex.ncells()];
                        for (k, (clone, d)) in
                            ex.clones.iter_mut().zip(to_delete.iter_mut()).enumerate()
                        {
                            let li = clone[0].dataset_index;
                            if nu[li] >= MIN_BASELINE_CELLS {
                                let (mut umish, mut umisl) = (0, 0);
                                for (s, c) in ex.share.iter().zip(clone.iter()) {
                                    if s.left {
                                        umish = max(umish, c.umi_count);
                                    } else {
                                        umisl = max(umish, c.umi_count);
                                    }
                                }
                                let umitot = umish + umisl;
                                if (umitot as f64) < umin[li] {
                                    if protected {
                                        if ex.share.len() == 1 {
                                            *d = true;
                                        }
                                    } else if !baselined
                                        || (best_ex, best_cell) != (j, k)
                                        || ex.share.len() == 1
                                    {
                                        *d = true;
                                    }
                                }
                            } else {
                                baselined = false;
                            }
                        }
                        if ctl.cr_opt.filter.umi_count {
                            for i in 0..ex.clones.len() {
                                if to_delete[i] {
                                    fate[ex.clones[i][0].dataset_index]
                                        .insert(ex.clones[i][0].barcode.clone(), BarcodeFate::Umi);
                                }
                            }
                            erase_if(&mut ex.clones, &to_delete);
                        }
                    }
                }
                cc.retain(|x| exact_clonotypes[info[*x as usize].clonotype_index].ncells() > 0);
            }
            if !cc.is_empty() { Some(cc) } else { None }
        })
        .collect();

    // Filter B cells based on UMI count ratios.  This assumes V..J identity to filter.
    const MIN_UMI_RATIO: usize = 500;
    candidate_clonotypes
        .into_iter()
        .filter_map(|cc| {
            let mut ncells = 0;
            for j in 0..cc.len() {
                let x: &CloneInfo = &info[cc[j] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                ncells += ex.ncells();
            }
            let mut nbads = 0;

            let mut z = Vec::<(Vec<u8>, usize, usize, usize, usize)>::new();
            let mut to_delete = Vec::<Vec<bool>>::new();
            for j in 0..cc.len() {
                let x: &CloneInfo = &info[cc[j] as usize];
                let ex = &exact_clonotypes[x.clonotype_index];
                to_delete.push(vec![false; ex.ncells()]);
                for k in 0..ex.ncells() {
                    let mut tot = 0;
                    for m in 0..ex.clones[k].len() {
                        tot += ex.clones[k][m].umi_count;
                    }
                    for m in 0..ex.clones[k].len() {
                        z.push((
                            ex.share[m].seq.clone(),
                            ex.clones[k][m].umi_count,
                            j,
                            k,
                            tot,
                        ));
                    }
                }
            }
            reverse_sort(&mut z);
            let mut j = 0;
            while j < z.len() {
                let k = next_diff1_5(&z, j);
                for l in j..k {
                    if z[j].1 >= MIN_UMI_RATIO * z[l].4 {
                        to_delete[z[l].2][z[l].3] = true;
                    }
                }
                j = k;
            }
            for j in 0..cc.len() {
                let x: &CloneInfo = &info[cc[j] as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                for l in 0..ex.ncells() {
                    if to_delete[j][l] {
                        nbads += 1;
                    }
                }
            }

            if nbads == 0 {
                return Some(cc);
            }
            let p = 0.1;
            let bound = 0.01;

            // Find probability of observing nbads or more events of probability
            // p in a sample of size ncells, and if that is at least bound,
            // don't delete any cells.
            if binomial_cdf(1.0 - p, ncells, ncells - nbads) >= bound {
                return Some(cc);
            }

            let mut to_deletex = vec![false; cc.len()];
            let mut z = Vec::<(Vec<u8>, usize, usize, usize, usize)>::new();
            let mut to_delete = Vec::<Vec<bool>>::new();
            for j in 0..cc.len() {
                let x: &CloneInfo = &info[cc[j] as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                to_delete.push(vec![false; ex.ncells()]);
                for k in 0..ex.ncells() {
                    let mut tot = 0;
                    for m in 0..ex.clones[k].len() {
                        tot += ex.clones[k][m].umi_count;
                    }
                    for m in 0..ex.clones[k].len() {
                        z.push((
                            ex.share[m].seq.clone(),
                            ex.clones[k][m].umi_count,
                            j,
                            k,
                            tot,
                        ));
                    }
                }
            }
            reverse_sort(&mut z);
            let mut j = 0;
            while j < z.len() {
                let k = next_diff1_5(&z, j);
                for l in j..k {
                    if z[j].1 >= MIN_UMI_RATIO * z[l].4 {
                        to_delete[z[l].2][z[l].3] = true;
                    }
                }
                j = k;
            }
            for j in 0..cc.len() {
                let x: &CloneInfo = &info[cc[j] as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                for l in 0..ex.ncells() {
                    if to_delete[j][l] {
                        nbads += 1;
                    }
                }

                if ctl.cr_opt.filter.umi_ratio {
                    for i in 0..ex.clones.len() {
                        if to_delete[j][i] {
                            fate[ex.clones[i][0].dataset_index]
                                .insert(ex.clones[i][0].barcode.clone(), BarcodeFate::UmiRatio);
                        }
                    }
                    erase_if(&mut ex.clones, &to_delete[j]);
                    if ex.ncells() == 0 {
                        to_deletex[j] = true;
                    }
                }
            }

            let mut cc = cc;
            if ctl.cr_opt.filter.umi_ratio {
                erase_if(&mut cc, &to_deletex);
            }
            if !cc.is_empty() { Some(cc) } else { None }
        })
        .collect()
}
