// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use enclone_core::join_one::join_one;
use equiv::EquivRel;
use qd::Double;
use std::collections::HashMap;
use vdj_ann::refx::RefData;

#[allow(clippy::too_many_arguments)]
pub fn join_core(
    i: usize,
    j: usize,
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    pot: &mut Vec<PotentialJoin>,
    refdata: &RefData,
) {
    let mut eq: EquivRel = EquivRel::new((j - i) as u32);
    for k1 in i..j {
        for k2 in k1 + 1..j {
            // Do nothing if join could have no effect on equivalence relation.
            // For certain samples, this hugely reduces run time.  That is the purpose of
            // having the equivalence relation.  Observed on MALT samples including 83808.
            // MALT is a B cell cancer in which j-i is very large and in fact the number of
            // exact subclonotypes in one clonotype is very large.

            if eq.set_id(k1 - i) == eq.set_id(k2 - i) {
                continue;
            }
            if join_one(k1, k2, ctl, exact_clonotypes, info, to_bc, sr, pot, refdata) {
                eq.join(k1 - i, k2 - i);
            }
        }
    }
}
