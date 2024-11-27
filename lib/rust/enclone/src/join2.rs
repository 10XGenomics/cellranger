// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file provides the tail end code for join.rs, plus a small function used there.

use enclone_core::defs::CloneInfo;
use equiv::EquivRel;
use stats_utils::percent_ratio;
use vector_utils::next_diff1_2;

pub struct JoinResult {
    pub i: usize,
    pub j: usize,
    pub joins: usize,
    pub errors: usize,
    pub join_list: Vec<(usize, usize)>,
}

impl JoinResult {
    pub fn new(i: usize, j: usize) -> Self {
        Self {
            i,
            j,
            joins: 0,
            errors: 0,
            join_list: Default::default(),
        }
    }
}

pub fn finish_join(
    info: &[CloneInfo],
    results: Vec<JoinResult>,
    report_whitelist_contamination: bool,
) -> EquivRel {
    // Tally results.
    // Make equivalence relation.
    let (mut joins, mut errors) = (0, 0);
    let mut eq: EquivRel = EquivRel::new(info.len() as u32);

    for r in results {
        joins += r.joins;
        errors += r.errors;
        for j in &r.join_list {
            eq.join(j.0, j.1);
        }
    }
    // Report whitelist contamination.
    // WARNING: THIS ONLY WORKS IF YOU RUN WITH CLONES=1 AND NO OTHER FILTERS.
    // TODO: we should actually make an assertion that this is true.
    if report_whitelist_contamination {
        let bad_rate = percent_ratio(joins, errors);
        println!("whitelist contamination rate = {bad_rate:.2}%");
    }

    // Join candidate clonotypes that cross subclones of a clone.  This arose because we split up multi-chain
    // clonotypes into two-chain clonotypes.

    let mut ox = Vec::new();
    for (i, x) in info.iter().enumerate() {
        ox.push((x.clonotype_index, eq.set_id(i)));
    }
    ox.sort_unstable();
    let mut i = 0;
    while i < ox.len() {
        let j = next_diff1_2(&ox, i);
        for k in i..j - 1 {
            eq.join(ox[k].1, ox[k + 1].1);
        }
        i = j;
    }

    eq
}
