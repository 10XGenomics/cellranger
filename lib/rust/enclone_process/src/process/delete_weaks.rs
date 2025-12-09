// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use crate::core::defs::ExactClonotype;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub(crate) fn delete_weaks(
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    mat: &[Vec<Option<usize>>],
    bads: &mut [bool],
) {
    let nexacts = exacts.len();
    for u in 0..nexacts {
        // NOTE: these checks are vestigial, related to filter options that
        // have been removed. They are likely unreachable.
        if exact_clonotypes[exacts[u]].ncells() == 0 {
            bads[u] = true;
        }
        if exact_clonotypes[exacts[u]].share.is_empty() {
            bads[u] = true;
        }
    }

    // Remove onesies that do not have an exact match.
    let cols = mat.len();
    if cols > 1 {
        for u1 in 0..nexacts {
            let ex1 = &exact_clonotypes[exacts[u1]];
            if ex1.share.len() == 1 && !bads[u1] {
                let mut perf = false;
                'u2: for u2 in 0..nexacts {
                    let ex2 = &exact_clonotypes[exacts[u2]];
                    if ex2.share.len() > 1 && !bads[u2] {
                        for i in 0..ex2.share.len() {
                            if ex1.share[0].seq == ex2.share[i].seq {
                                perf = true;
                                break 'u2;
                            }
                        }
                    }
                }
                if !perf {
                    bads[u1] = true;
                }
            }
        }
    }
}
