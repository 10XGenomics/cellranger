// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype, PotentialJoin};
use enclone_core::join_one::join_one;
use equiv::EquivRel;
use qd::Double;
use std::cmp::max;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::{bin_position, next_diff12_3, next_diff1_3, unique_sort};

// Define an equivalence relation on the chains, introducing connections defined by the
// raw joins.  Also join where there are identical V..J sequences.

fn joiner(
    infos: &[usize],
    info: &[CloneInfo],
    to_exacts: &HashMap<usize, usize>,
    raw_joinsx: &[Vec<usize>],
    chains: &[(usize, usize)],
    seq_chains: &[(Vec<u8>, usize, usize)],
) -> EquivRel {
    let mut e = EquivRel::new(chains.len() as u32);
    for i1 in 0..infos.len() {
        let j1 = infos[i1];
        let u1 = info[j1].clonotype_index;
        let v1 = to_exacts[&u1];
        let m1s = &info[j1].exact_cols;
        for i2 in &raw_joinsx[i1] {
            let j2 = infos[*i2];
            let u2 = info[j2].clonotype_index;
            let v2 = to_exacts[&u2];
            let m2s = &info[j2].exact_cols;
            for j in 0..2 {
                let z1 = bin_position(chains, &(v1, m1s[j]));
                let z2 = bin_position(chains, &(v2, m2s[j]));
                e.join(z1 as usize, z2 as usize);
            }
        }
    }
    let mut i = 0;
    while i < seq_chains.len() {
        let j = next_diff1_3(seq_chains, i);
        for k in i + 1..j {
            let (x1, x2) = (&seq_chains[i], &seq_chains[k]);
            let z1 = bin_position(chains, &(x1.1, x1.2));
            let z2 = bin_position(chains, &(x2.1, x2.2));
            e.join(z1 as usize, z2 as usize);
        }
        i = j;
    }
    e
}

// TOOD: refactor this into a struct
pub type Od = (Vec<usize>, usize, i32);

pub fn setup_define_mat(candidate_clonotype: &[i32], info: &[CloneInfo]) -> (Vec<Od>, Vec<usize>) {
    let mut od = Vec::<Od>::new();
    for id in candidate_clonotype {
        let x: &CloneInfo = &info[*id as usize];
        od.push((x.origin.clone(), x.clonotype_index, *id));
    }
    od.sort();
    let mut exacts = Vec::<usize>::new();
    let mut j = 0;
    while j < od.len() {
        let k = next_diff12_3(&od, j);
        exacts.push(od[j].1);
        j = k;
    }
    (od, exacts)
}

// This generates a cols x nexacts matrices for a given clonotype, where cols is defined by the
// algorithm, and is the number of columns (chains) in the clonotype table.
#[allow(clippy::too_many_arguments)]
pub fn define_mat(
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    exacts: &[usize],
    od: &[Od],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    refdata: &RefData,
) -> Vec<Vec<Option<usize>>> {
    // Define map of indices into exacts.

    let nexacts = exacts.len();
    let to_exacts: HashMap<usize, usize> =
        exacts.iter().enumerate().map(|(u, &x)| (x, u)).collect();

    // Get the info indices corresponding to this clonotype.

    let mut infos: Vec<usize> = od
        .iter()
        .filter_map(|i| {
            let x = i.2 as usize;
            if to_exacts.contains_key(&info[x].clonotype_index) {
                Some(x)
            } else {
                None
            }
        })
        .collect();
    infos.sort_unstable();

    // Define map of exacts to infos.

    let mut to_infos = vec![Vec::<usize>::new(); nexacts];
    for (i, &inf) in infos.iter().enumerate() {
        let u = to_exacts[&info[inf].clonotype_index];
        to_infos[u].push(i);
    }

    // Form the set of all chains that appear in an exact subclonotypes of this clonotype, and
    // also track the V..J sequences for the chains.

    let mut chains = Vec::<(usize, usize)>::new();
    let mut seq_chains = Vec::<(Vec<u8>, usize, usize)>::new();
    for (u, &exu) in exacts.iter().enumerate() {
        let ex = &exact_clonotypes[exu];
        for (m, share) in ex.share.iter().enumerate() {
            chains.push((u, m));
            seq_chains.push((share.seq.clone(), u, m));
        }
    }
    seq_chains.sort();

    // Gather the raw joins.

    let mut raw_joinsx = vec![Vec::<usize>::new(); infos.len()];
    for (&j1, raw_i1) in infos.iter().zip(raw_joinsx.iter_mut()) {
        for x in &raw_joins[j1] {
            let i2 = bin_position(&infos, x);
            if i2 >= 0 {
                raw_i1.push(i2 as usize);
            }
        }
    }

    // Look for additional raw joins.  In the join stage, for efficiency reasons, we don't try to
    // make a join if two info elements are already in the same equivalence class.  This causes
    // us to miss raw joins for two reasons:
    // • One reason is that where we have two exact subclonotypes, and at least one one has more
    //   than two chains, then we may miss a raw join that is needed here.
    // • Another reason is that exact subclonotypes may have been deleted since the original
    //   equivalence relation was built.  This we address partially.

    let mut extras = Vec::<(usize, usize)>::new();
    for (i1, (raw_i1, &j1)) in raw_joinsx.iter().zip(infos.iter()).enumerate() {
        for &i2 in raw_i1 {
            let j2 = infos[i2];
            let (u1, u2) = (info[j1].clonotype_index, info[j2].clonotype_index);
            let (ex1, ex2) = (&exact_clonotypes[u1], &exact_clonotypes[u2]);
            let (v1, v2) = (to_exacts[&u1], to_exacts[&u2]);
            if ex1.share.len() > 2 || ex2.share.len() > 2 {
                let (s1, s2) = (&to_infos[v1], &to_infos[v2]);
                for k1 in s1 {
                    for k2 in s2 {
                        let (k1, k2) = (*k1, *k2);
                        if (k1 == i1 && k2 == i2) || (k1 == i2 && k2 == i1) {
                            continue;
                        }
                        let (l1, l2) = (infos[k1], infos[k2]);
                        if info[l1].lens == info[l2].lens {
                            let mut pot = Vec::<PotentialJoin>::new();
                            if join_one(
                                l1,
                                l2,
                                ctl,
                                exact_clonotypes,
                                info,
                                to_bc,
                                sr,
                                &mut pot,
                                refdata,
                            ) {
                                extras.push((k1, k2));
                            }
                        }
                    }
                }
            }
        }
    }
    for x in &extras {
        raw_joinsx[x.0].push(x.1);
    }

    // Define an initial equivalence relation on the chains, and get set representatives.

    let mut e = joiner(&infos, info, &to_exacts, &raw_joinsx, &chains, &seq_chains);
    let r = e.set_reps();

    // First for each pair of chain sets with one "heavy" and one "light", pick some info
    // entries, if there are any.  This is effectively at random.  A parameter governs how
    // much we pick.

    let mut rxi = Vec::<(usize, usize, usize)>::new(); // (heavy set, light set, infos index)
    for (i, &inf_i) in infos.iter().enumerate() {
        let z = &info[inf_i];
        let u = z.clonotype_index;
        let v = to_exacts[&u];
        if z.exact_cols.len() != 2 {
            continue;
        }
        let (m1, m2) = (z.exact_cols[0], z.exact_cols[1]);
        let ex = &exact_clonotypes[u];
        if !ex.share[m1].left || ex.share[m2].left {
            continue; // maybe never happens
        }
        let p1 = e.set_id(bin_position(&chains, &(v, m1)) as usize);
        let p2 = e.set_id(bin_position(&chains, &(v, m2)) as usize);
        let q1 = bin_position(&r, &p1) as usize;
        let q2 = bin_position(&r, &p2) as usize;
        rxi.push((q1, q2, i));
    }
    rxi.sort_unstable();
    const MAX_USE: usize = 5; // knob set empirically
    let mut rxir = Vec::<(usize, usize, usize)>::new(); // (heavy set, light set, info index)
    let mut i = 0;
    while i < rxi.len() {
        let j = next_diff12_3(&rxi, i);
        rxir.extend(&rxi[i..j.min(i + MAX_USE)]);
        i = j;
    }

    // Now for each pair of these, if they are not effectively joined, attempt to join them.
    // This partially addresses the "second reason" described above.  It is partial because we
    // picked an info entry above at random, rather than trying them all.

    for f1 in &rxir {
        for f2 in &rxir {
            if f1.0 != f2.0 || f1.1 != f2.1 {
                let (i1, i2) = (infos[f1.2], infos[f2.2]);
                if info[i1].lens != info[i2].lens {
                    continue;
                }
                let mut pot = Vec::<PotentialJoin>::new();
                if join_one(
                    i1,
                    i2,
                    ctl,
                    exact_clonotypes,
                    info,
                    to_bc,
                    sr,
                    &mut pot,
                    refdata,
                ) {
                    e.join(r[f1.0], r[f2.0]);
                    e.join(r[f1.1], r[f2.1]);
                }
            }
        }
    }

    // Find the exact subclonotypes having three chains and list their sets, allowing for order
    // to vary.

    let r = e.set_reps();
    let mut threes = Vec::<Vec<usize>>::new();
    let mut threesp = HashMap::<Vec<usize>, usize>::new();
    for u in 0..nexacts {
        let ex = &exact_clonotypes[exacts[u]];
        if ex.share.len() == 3 {
            let zs = [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ];
            for z in &zs {
                if ex.share[z[0]].left && !ex.share[z[2]].left {
                    let p1 = e.set_id(bin_position(&chains, &(u, z[0])) as usize);
                    let p2 = e.set_id(bin_position(&chains, &(u, z[1])) as usize);
                    let p3 = e.set_id(bin_position(&chains, &(u, z[2])) as usize);
                    let q1 = bin_position(&r, &p1) as usize;
                    let q2 = bin_position(&r, &p2) as usize;
                    let q3 = bin_position(&r, &p3) as usize;
                    threes.push(vec![q1, q2, q3]);
                    threesp.insert(vec![q1, q2, q3], u);
                }
            }
        }
    }
    unique_sort(&mut threes);

    // There is one more case to deal with.  This is where we have two exact subclonotypes, each
    // with three chains, and we joined two of their chains, but not the third.  And where the
    // join algorithm would not have joined the third.  In this case, if the third chains are
    // "close enough", we join them anyway.  As before, we only test representatives.

    for t1 in &threes {
        't2_loop: for t2 in &threes {
            if t1 == t2 {
                continue;
            }
            let (mut matches, mut mismatch) = (0, 0);
            for i in 0..3 {
                if t1[i] == t2[i] {
                    matches += 1;
                } else {
                    mismatch = i;
                }
            }
            if matches != 2 {
                continue;
            }
            let (u1, u2) = (threesp[t1], threesp[t2]);
            let (ex1, ex2) = (&exact_clonotypes[exacts[u1]], &exact_clonotypes[exacts[u2]]);
            for (m1, ex1_sm1) in ex1.share.iter().enumerate() {
                let p1 = bin_position(&chains, &(u1, m1)) as usize;
                let q1 = bin_position(&r, &p1) as usize;
                if q1 == t1[mismatch] {
                    for (m2, ex2_sm2) in ex2.share.iter().enumerate() {
                        let p2 = bin_position(&chains, &(u2, m2)) as usize;
                        let q2 = bin_position(&r, &p2) as usize;
                        if q2 == t2[mismatch] {
                            let (seq1, seq2) = (&ex1_sm1.seq, &ex2_sm2.seq);
                            if seq1.len() == seq2.len() {
                                const MAX_DIFFS: usize = 10;
                                let diffs = seq1
                                    .iter()
                                    .zip(seq2.iter())
                                    .filter(|(&s1, &s2)| s1 != s2)
                                    .count();
                                if diffs <= MAX_DIFFS {
                                    e.join(p1, p2);
                                    break 't2_loop;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Get representatives for the chain sets.
    let mut r = e.set_reps();

    // Reorder the chains.  This is done to get the heavy chains before the light chains and also
    // to mimic the behavior of the previous version of this algorithm, to minimiize churn.  Then
    // update the representatives.

    let mut chainsp = Vec::<(String, usize, usize, usize)>::with_capacity(exacts.len());
    for (u, &exu) in exacts.iter().enumerate() {
        let ex = &exact_clonotypes[exu];
        for (m, share_m) in ex.share.iter().enumerate() {
            let mut c = share_m.chain_type.clone();
            if c.starts_with("TRB") {
                c = c.replacen("TRB", "TRX", 1);
            } else if c.starts_with("TRA") {
                c = c.replacen("TRA", "TRY", 1);
            }
            chainsp.push((format!("{c}:{}", share_m.cdr3_aa), share_m.seq.len(), u, m));
        }
    }
    chainsp.sort();
    let mut chainso = Vec::<(usize, usize)>::new();
    let mut chainsox = Vec::<(usize, usize, usize)>::new();
    for (i, c) in chainsp.into_iter().enumerate() {
        chainso.push((c.2, c.3));
        chainsox.push((c.2, c.3, i));
    }
    chainsox.sort_unstable();
    for ri in &mut r {
        *ri = chainsox[*ri].2;
    }
    r.sort_unstable();

    // Create rmap, that sends
    // (index into exact subclonotypes for this clonotype,
    //  index into chains for one of these exact subclonotypes)
    // to an index into the set reps for the chains.

    let mut rpos = HashMap::<(usize, usize), usize>::new();
    for (i, chain) in chains.into_iter().enumerate() {
        let c = e.set_id(i);
        let f = chainsox[c as usize].2;
        let q = bin_position(&r, &f) as usize;
        rpos.insert((chain.0, chain.1), q);
    }

    // Find the maximum multiplicity of each set, and the number of columns.

    let mut mm = vec![0; r.len()];
    for (u, &exu) in exacts.iter().enumerate() {
        let ex = &exact_clonotypes[exu];
        let mut mm0 = vec![0; r.len()];
        for m in 0..ex.share.len() {
            mm0[rpos[&(u, m)]] += 1;
        }
        for i in 0..r.len() {
            mm[i] = max(mm[i], mm0[i]);
        }
    }
    let cols = mm.iter().sum();

    // Define a matrix mat[col][ex] which is the column of the exact subclonotype ex corresponding
    // to the given column col of the clonotype, which may or may not be defined.

    let mut mat = vec![vec![None; nexacts]; cols];
    for (cx, cc) in mat.iter_mut().enumerate() {
        // for every column
        'exact: for (u, &exu) in exacts.iter().enumerate() {
            // for every exact subclonotype
            let ex = &exact_clonotypes[exu];
            let mut mm0 = vec![0; r.len()];
            for m in 0..ex.share.len() {
                // for every chain in the exact subclonotype:
                let q = rpos[&(u, m)];
                let mut col = mm0[q];
                col += mm.iter().take(q).sum::<usize>();
                mm0[q] += 1;
                if col == cx {
                    cc[u] = Some(m);
                    continue 'exact;
                }
            }
        }
    }
    mat
}
