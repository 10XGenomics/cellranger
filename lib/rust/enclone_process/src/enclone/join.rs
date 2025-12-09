// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// This file provides the single function join_exacts, which computes the equivalence relation
// on exact subclonotypes.
//
// Note that in principle the specificity of joining might be increased by using nonproductive
// contigs that represent the sequence of the "other" allele.  This does not look easy to
// execute.

use self::refx::RefData;
use super::join_core::join_core;
use super::join2::{JoinResult, finish_join};
use crate::core::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use equiv::EquivRel;
use qd::Double;
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashMap;
use vdj_ann::refx;
use vector_utils::erase_if;

pub(crate) fn join_exacts(
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    refdata: &RefData,
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    sr: &[Vec<Double>],
    report_whitelist_contamination: bool,
) -> (EquivRel, Vec<(i32, i32)>) {
    // Find potential joins.

    let mut i = 0;
    let mut results = Vec::<JoinResult>::new();
    while i < info.len() {
        let mut j = i + 1;
        while j < info.len() {
            // Note that the organization of the loop here separates info entries by their
            // contig lengths.  One could rejigger the code to also separate by CDR3 lengths,
            // but surprisingly this doesn't help much if any.  It does perturb results very
            // slightly.
            if info[j].lens != info[i].lens {
                break;
            }
            j += 1;
        }
        results.push(JoinResult::new(i, j));
        i = j;
    }

    let joinf = |r: &mut JoinResult| {
        let (i, j) = (r.i, r.j);
        let joins = &mut r.joins;
        let errors = &mut r.errors;
        let mut pot = Vec::<PotentialJoin>::new();

        // Main join logic.  If you change par_iter_mut to iter_mut above, and run happening,
        // a lot of time shows up on the following line.  If further you manually inline join_core
        // here, then instead the time shows up on the results.iter_mut line.  It's not clear
        // what this means.

        join_core(
            i,
            j,
            ctl,
            exact_clonotypes,
            info,
            to_bc,
            sr,
            &mut pot,
            refdata,
        );

        // Run two passes.

        for _ in 0..2 {
            // Form the equivalence relation implied by the potential joins.

            let mut eq: EquivRel = EquivRel::new((j - i) as u32);
            for pot in &pot {
                let (k1, k2) = (pot.k1, pot.k2);
                eq.join(k1 - i, k2 - i);
            }

            // Impose a higher bar on joins that involve only two cells. (not documented)

            let mut to_pot = vec![Vec::<usize>::new(); j - i];
            for (pj, pot) in pot.iter().enumerate() {
                let k1 = pot.k1;
                to_pot[k1 - i].push(pj);
            }
            let mut to_delete = vec![false; pot.len()];
            for mut cc_iter in eq.all_sets() {
                // Count the number of cells in the candidate clonotype.

                let mut ncells = 0;
                for t in cc_iter.duplicate() {
                    let k = t + i;
                    let mult = exact_clonotypes[info[k].clonotype_index].ncells();
                    ncells += mult;
                }

                // Impose more stringent conditions if number of cells is two.

                if ncells == 2 && cc_iter.size() == 2 {
                    let x0 = cc_iter.next().unwrap();
                    let x1 = cc_iter.next().unwrap();
                    let (k1, k2) = (x0 + i, x1 + i);
                    let k = min(k1, k2);
                    for pj in &to_pot[k - i] {
                        let cd = pot[*pj].cd;
                        let shares = &pot[*pj].shares;

                        // Impose higher bar.

                        let min_shares = shares.iter().min().unwrap();
                        if cd > *min_shares / 2 {
                            to_delete[*pj] = true;
                        }
                    }
                }
            }
            erase_if(&mut pot, &to_delete);
        }

        // Analyze potential joins.

        let mut eq: EquivRel = EquivRel::new((j - i) as u32);
        for pj in pot {
            let k1 = pj.k1;
            let k2 = pj.k2;
            let err = pj.err;

            // Do nothing if join could have no effect on equivalence relation.

            if eq.set_id(k1 - i) == eq.set_id(k2 - i) {
                continue;
            }

            // Save join and tally stats.

            r.join_list.push((k1, k2));
            *joins += 1;
            if err {
                *errors += 1;
            }
            eq.join(k1 - i, k2 - i);
        }
    };

    results.par_iter_mut().for_each(joinf);

    let mut raw_joins = Vec::new();
    for r in &results {
        for &j in &r.join_list {
            raw_joins.push((j.0 as i32, j.1 as i32));
        }
    }
    let eq = finish_join(info, results, report_whitelist_contamination);
    (eq, raw_joins)
}
