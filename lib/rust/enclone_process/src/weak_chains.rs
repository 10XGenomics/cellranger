// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//! Based on the number of cells in each column, decide which exact subclonotypes
//! look like junk.  Preliminary heuristic.
#![deny(missing_docs)]

use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use crate::core::enclone_structs::{BarcodeFates, CandidateClonotype};
use crate::process::define_mat::{define_mat, setup_define_mat};
use qd::Double;
use rayon::prelude::*;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::erase_if;

#[derive(Default)]
struct Result {
    f1: Vec<(usize, String, BarcodeFate)>,
    f2: Vec<usize>,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn weak_chains(
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
    // Note mat calculation duplicated with print_clonotypes and also doublet detection.
    let results: Vec<_> = candidate_clonotypes
        .par_iter()
        .map(|cc| {
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
            let mut res = Result::default();
            if cols > 2 {
                let nexacts = exacts.len();
                let mut ncells = vec![0; cols];
                let mut col_entries = vec![Vec::<usize>::new(); cols];
                for (u, &clonotype_id) in exacts.iter().enumerate().take(nexacts) {
                    for ((nc, ce), m) in ncells
                        .iter_mut()
                        .zip(col_entries.iter_mut())
                        .zip(mat.iter().take(cols))
                    {
                        let mid = m[u];
                        if mid.is_some() {
                            ce.push(u);
                            *nc += exact_clonotypes[clonotype_id].clones.len();
                        }
                    }
                }
                let mut total_cells = 0;
                for j in 0..exacts.len() {
                    total_cells += exact_clonotypes[exacts[j]].ncells();
                }
                for j in 0..cols {
                    if ncells[j] <= 20 && 8 * ncells[j] < total_cells {
                        for d in &col_entries[j] {
                            res.f2.push(exacts[*d]);
                            let ex = &exact_clonotypes[exacts[*d]];
                            for i in 0..ex.ncells() {
                                res.f1.push((
                                    ex.clones[i][0].dataset_index,
                                    ex.clones[i][0].barcode.clone(),
                                    BarcodeFate::WeakChains,
                                ));
                            }
                        }
                    }
                }
            }
            res
        })
        .collect();
    let mut to_delete = vec![false; exact_clonotypes.len()];
    let mut dels = Vec::<i32>::new();
    for i in 0..results.len() {
        for j in 0..results[i].f1.len() {
            fate[results[i].f1[j].0].insert(results[i].f1[j].1.clone(), results[i].f1[j].2);
        }
        for x in &results[i].f2 {
            to_delete[*x] = true;
        }
    }
    dels.sort_unstable();
    for cc in candidate_clonotypes.iter_mut() {
        let del = cc
            .iter()
            .map(|&oj| {
                let id = info[oj as usize].clonotype_index;
                to_delete[id]
            })
            .collect::<Vec<_>>();
        erase_if(cc, &del);
    }
}
