// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Merge onesies where totally unambiguous.  Possibly inefficient and should optimize.

use crate::core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use crate::core::enclone_structs::CandidateClonotype;
use equiv::EquivRel;
use vector_utils::{lower_bound1_2, unique_sort, upper_bound1_2};

pub(crate) fn merge_onesies(
    candidate_clonotypes: &mut Vec<CandidateClonotype>,
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    eq: &EquivRel,
    disintegrated: &[bool],
) {
    {
        // ctl.join_alg_opt.merge_onesies is always true
        let mut eqo = EquivRel::new(candidate_clonotypes.len() as u32);
        let mut to_candidate_clonotype = vec![None; info.len()];
        for i in 0..candidate_clonotypes.len() {
            for j in 0..candidate_clonotypes[i].len() {
                to_candidate_clonotype[candidate_clonotypes[i][j] as usize] = Some(i);
            }
        }
        let ncells_total = exact_clonotypes.iter().map(ExactClonotype::ncells).sum();
        let mut onesies = Vec::<usize>::new();
        for i in 0..info.len() {
            if to_candidate_clonotype[i].is_some()
                && info[i].tigs.len() == 1
                && !disintegrated[info[i].clonotype_index]
            {
                onesies.push(i);
            }
        }
        let mut alltigs2 = Vec::<(Vec<u8>, usize)>::new();
        for i in 0..info.len() {
            if to_candidate_clonotype[i].is_some() && info[i].tigs.len() >= 2 {
                for j in 0..info[i].tigs.len() {
                    alltigs2.push((info[i].tigs[j].clone(), i));
                }
            }
        }
        alltigs2.sort();
        for x in &onesies {
            let low = lower_bound1_2(&alltigs2, &info[*x].tigs[0]);
            let high = upper_bound1_2(&alltigs2, &info[*x].tigs[0]);
            let mut ms = Vec::<usize>::new();
            for m in low..high {
                if alltigs2[m as usize].0 == info[*x].tigs[0] {
                    ms.push(m as usize);
                }
            }
            let mut ok = !ms.is_empty();
            let mut exacts = Vec::<usize>::new();
            for j in 0..ms.len() {
                if eq.set_id(alltigs2[ms[j]].1) != eq.set_id(alltigs2[ms[0]].1) {
                    ok = false;
                }
                for z in eq.members_of_set_containing(alltigs2[ms[j]].1) {
                    exacts.push(info[z].clonotype_index);
                }
            }
            unique_sort(&mut exacts);

            let ncells0 = exact_clonotypes[info[*x].clonotype_index].ncells();
            if ncells0 * 10000 < ncells_total {
                ok = false;
            }

            if ok {
                let orb1 = to_candidate_clonotype[*x].unwrap();
                let orb2 = to_candidate_clonotype[alltigs2[ms[0]].1].unwrap();

                // Test for donor mixing.

                if !ctl.cr_opt.mix_donors {
                    let mut donors = vec![Vec::<Option<usize>>::new(); 2];
                    let orbs = [&orb1, &orb2];
                    for (pass, orb) in orbs.iter().enumerate() {
                        for id in &candidate_clonotypes[**orb] {
                            let ex = &exact_clonotypes[info[*id as usize].clonotype_index];
                            for i in 0..ex.clones.len() {
                                donors[pass].push(ex.clones[i][0].donor_index);
                            }
                        }
                        unique_sort(&mut donors[pass]);
                    }
                    if donors[0] != donors[1] {
                        continue;
                    }
                }

                // Make join.

                eqo.join(orb1, orb2);
            }
        }
        let mut candidate_clonotypes2 = Vec::<Vec<i32>>::new();
        for cc_iter in eqo.all_sets() {
            let mut cc = Vec::<i32>::new();
            for oj in cc_iter {
                cc.extend(&candidate_clonotypes[oj]);
            }
            candidate_clonotypes2.push(cc);
        }
        *candidate_clonotypes = candidate_clonotypes2;
    }
}
