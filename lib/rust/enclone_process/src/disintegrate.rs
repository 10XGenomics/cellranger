// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]
//
// If NWEAK_ONESIES is not specified, disintegrate certain onesie clonotypes into single cell
// clonotypes.  This requires editing of exact_clonotypes, info, eq, join_info and raw_joins.

use crate::core::defs::{CloneInfo, ExactClonotype};
use equiv::EquivRel;
use std::collections::HashMap;
use vector_utils::unique_sort;

pub(crate) struct DisintegrateOnesiesResult {
    pub(crate) disintegrated: Vec<bool>,
    pub(crate) eq: EquivRel,
    pub(crate) exact_clonotypes: Vec<ExactClonotype>,
    pub(crate) info: Vec<CloneInfo>,
    pub(crate) raw_joins: Vec<(i32, i32)>,
}

pub(crate) fn disintegrate_onesies(
    eq: EquivRel,
    exact_clonotypes: Vec<ExactClonotype>,
    info: Vec<CloneInfo>,
    raw_joins: Vec<(i32, i32)>,
) -> DisintegrateOnesiesResult {
    let ncells_total = exact_clonotypes.iter().map(ExactClonotype::ncells).sum();
    let mut to_info = HashMap::new();
    let mut exacts2 = Vec::new();
    for (i, inf) in info.iter().enumerate() {
        to_info.insert(inf.clonotype_index, i);
    }
    let mut to_exact_new = Vec::new();
    let mut disintegrated = Vec::new();
    for i in 0..exact_clonotypes.len() {
        let ex = &exact_clonotypes[i];
        let mut enew = Vec::<usize>::new();
        if ex.share.len() == 1
            && ex.ncells() > 1
            && ex.ncells() * 1000 < ncells_total
            && to_info.contains_key(&i)
            && eq.size_of_set_containing(to_info[&i]) == 1
        {
            for j in 0..ex.clones.len() {
                enew.push(exacts2.len());
                exacts2.push(ExactClonotype {
                    share: ex.share.clone(),
                    clones: vec![ex.clones[j].clone()],
                });
                disintegrated.push(true);
            }
        } else {
            enew.push(exacts2.len());
            exacts2.push(exact_clonotypes[i].clone());
            disintegrated.push(false);
        }
        to_exact_new.push(enew);
    }

    let exact_clonotypes = exacts2;
    let mut info2 = Vec::<CloneInfo>::new();
    let mut to_info2 = Vec::<Vec<usize>>::new();
    for clone_info in info {
        let mut x = Vec::<usize>::new();
        for clonotype_index in &to_exact_new[clone_info.clonotype_index] {
            let mut origins = Vec::<usize>::new();
            let ex = &exact_clonotypes[*clonotype_index];
            for i in 0..ex.clones.len() {
                origins.push(ex.clones[i][0].dataset_index);
            }
            unique_sort(&mut origins);
            let origin = origins;
            x.push(info2.len());
            info2.push(CloneInfo {
                clonotype_index: *clonotype_index,
                origin,
                ..clone_info.clone()
            });
        }
        to_info2.push(x);
    }

    let info = info2;
    let raw_joins = raw_joins
        .into_iter()
        .map(|raw_join| {
            (
                to_info2[raw_join.0 as usize][0] as i32,
                to_info2[raw_join.1 as usize][0] as i32,
            )
        })
        .collect();
    let mut eq2 = EquivRel::new(info.len() as u32);
    for cc in eq.all_sets().filter(|cc| cc.size() > 1) {
        for (o1, o2) in cc.duplicate().zip(cc.skip(1)) {
            eq2.join(to_info2[o1][0], to_info2[o2][0]);
        }
    }
    DisintegrateOnesiesResult {
        disintegrated,
        eq: eq2,
        exact_clonotypes,
        info,
        raw_joins,
    }
}
