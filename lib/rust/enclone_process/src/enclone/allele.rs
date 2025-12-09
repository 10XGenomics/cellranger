// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// This file provides functions to find alternate alleles and substitute them into references.

use self::refx::RefData;
use crate::core::defs::{AltRef, CloneInfo, ExactClonotype};
use crate::core::join_one::REF_V_TRIM;
use debruijn::Mer;
use rayon::prelude::*;
use std::cmp::{PartialOrd, max};
use std::iter::zip;
use vdj_ann::refx;
use vector_utils::{erase_if, next_diff, next_diff1_2, next_diff1_3, reverse_sort, unique_sort};

/// Derive consensus sequences for alternate alleles of V segments.
///
/// The priority of this algorithm is to reduce the likelihood of false positive joins.  It
/// has been tuned and optimized for this purpose, possibly at the expense of generating some
/// 'fake' alternate alleles.  However we do not actually know that this happens.  If we
/// wanted to test the true accuracy of allele finding, we would need whole-genome and VDJ
/// datasets from the same donors.
///
/// This calculation has to be made separately for each donor, which means that in a certain
/// sense the algorithm is not blinded to the truth data.  However, separating this calculation
/// out per donor is the right thing to do.
///
/// Alternate alleles might correspond to duplicated segments, which is fine, as
/// for purposes of this code that's functionally equivalent to bona fide alternate alleles.
///
/// We do not attempt to calculate the last 15 bases of an alternate allele.  These
/// bases are just copied.  If we want to really know these bases we may need to
/// have actual genomic sequences which could also provide a control for these calculations.
///
/// Attempts to also do this for J segments were unsuccessful.
///
/// Limitations and to do items:
/// 1. Hypothetically we could make a library of alternate alleles and use that
///    as a supplement or in place of calculating on the fly.
/// 2. Document as part of algorithm.
/// 3. Make alt_refs into a more efficient data structure.
/// 4. Speed up.
///
/// Organize data by reference ID.  Note that we ignore exact subclonotypes having four chains.
pub(crate) fn find_alleles(refdata: &RefData, exact_clonotypes: &[ExactClonotype]) -> Vec<AltRef> {
    let mut allxy =
        vec![Vec::<(usize, Vec<u8>, Vec<usize>, usize, usize, String)>::new(); refdata.refs.len()];
    for (m, x) in exact_clonotypes.iter().enumerate() {
        if x.share.len() >= 2 && x.share.len() <= 3 {
            for j in 0..x.share.len() {
                let y = &x.share[j];
                let id = y.v_ref_id;

                // Find the partner chains.

                let mut partner = Vec::<usize>::new();
                for ja in 0..x.share.len() {
                    if x.share[ja].left != x.share[j].left {
                        partner.push(x.share[ja].v_ref_id);
                    }
                }
                if !partner.is_empty() && y.seq_del.len() >= refdata.refs[id].len() - REF_V_TRIM {
                    for clone in &x.clones {
                        let donor = clone[j].donor_index;
                        if let Some(donor) = donor {
                            allxy[id].push((
                                donor,
                                y.seq_del.clone(),
                                partner.clone(),
                                m,
                                clone[j].dataset_index,
                                clone[0].barcode.clone(),
                            ));
                        }
                    }
                }
            }
        }
    }

    // Process each reference ID.

    let mut vs = Vec::<usize>::new();
    for id in 0..refdata.refs.len() {
        if refdata.is_v(id) {
            vs.push(id);
        }
    }
    let mut results = Vec::<(usize, Vec<AltRef>)>::new();
    for v in &vs {
        results.push((*v, Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let id = res.0;
        let mut allx = allxy[id].clone();

        // Divide by donor.

        allx.sort();
        let mut alls = Vec::<Vec<(usize, Vec<u8>, Vec<usize>, usize, usize, String)>>::new();
        let mut i = 0;
        while i < allx.len() {
            // let j = next_diff1_6(&allx, i);
            let mut j = i + 1;
            while j < allx.len() {
                if allx[j].0 != allx[i].0 {
                    break;
                }
                j += 1;
            }
            alls.push(allx[i..j].to_owned());
            i = j;
        }

        // Process donor by donor.

        for mut all in alls {
            // Data here are given by "all", the relevant entries of which are:
            // 1: V..J sequence for one chain of a given info entry
            // 2: the reference ID(s) of the partner chain(s) -- possibly not used
            // 3: the index in exact_clonotypes
            // 4: the dataset ID.

            let donor_id = all[0].0;

            // If two entries have
            // * the same length CDR3 sequence and
            // * the same partner chain reference ids for both V and J and
            // * the same length partner chain CDR3 sequence,
            // arbitrarily delete one of them.
            //
            // What we're trying to do here is reduce the incidence wherein multiple entries from
            // the same clonotype (which we have not yet computed) contribute to the count for a
            // putative alternate allele, resulting in alternate alleles that are not true
            // germline changes, but rather arise from SHM.
            //
            // In a better system we might first compute the clonotypes, then the alternate
            // alleles, however knowing the alternate alleles is in fact prerequisite to computing
            // the clonotypes.
            //
            // Note that a vulnerability of the algorithm is that if there is a very large
            // clonotype, then artifactual pairs arising from it could provide enough "evidence"
            // to create an alternate allele, which in fact should not exist.  This has not been
            // observed, but we haven't looked carefully. Biological cases such as lymphomas
            // could provide helpful test cases in this area.

            let mut to_delete = vec![false; all.len()];
            {
                let mut trace = Vec::<((usize, usize, usize, usize), usize)>::new();
                for (i, item) in all.iter().enumerate() {
                    let u = item.3;
                    let ex = &exact_clonotypes[u];
                    for j1 in 0..ex.share.len() {
                        if ex.share[j1].seq_del == item.1 {
                            for j2 in 0..ex.share.len() {
                                let (s1, s2) = (&ex.share[j1], &ex.share[j2]);
                                if s2.left != s1.left {
                                    trace.push((
                                        (
                                            s1.cdr3_dna.len(),
                                            s2.cdr3_dna.len(),
                                            s2.v_ref_id,
                                            s2.j_ref_id,
                                        ),
                                        i,
                                    ));
                                }
                            }
                        }
                    }
                }
                trace.sort_unstable();
                let mut i = 0;
                while i < trace.len() {
                    let j = next_diff1_2(&trace, i);
                    for k in i..j {
                        if !to_delete[trace[k].1] {
                            for l in i..j {
                                if l != k {
                                    to_delete[trace[l].1] = true;
                                }
                            }
                            break;
                        }
                    }
                    i = j;
                }
                erase_if(&mut all, &to_delete);
            }

            // Traverse the positions in the reference V segment.

            let mut ps = Vec::<usize>::new();
            for p in 0..refdata.refs[id].len() - REF_V_TRIM {
                // Let bases = {(contig base, contig index, index of the other ref V segment)}.
                // The other ref V segments are there only for diagnostic purposes.

                let mut bases = Vec::<(u8, usize, Vec<usize>)>::new();
                for (i, x) in all.iter().enumerate() {
                    bases.push((x.1[p], i, x.2.clone()));
                }
                bases.sort();

                // Let freqs = {( number of contig bases,
                //                {contig indices}, {other ref V segment indices}, base )}.

                let mut freqs = Vec::<(usize, Vec<usize>, Vec<Vec<usize>>, u8)>::new();
                let mut i = 0;
                while i < bases.len() {
                    let j = next_diff1_3(&bases, i);
                    let mut x = Vec::<usize>::new();
                    let mut y = Vec::<Vec<usize>>::new();
                    for base in &bases[i..j] {
                        x.push(base.1);
                        y.push(base.2.clone());
                    }
                    freqs.push((j - i, x, y, bases[i].0));
                    i = j;
                }
                reverse_sort(&mut freqs);

                // If frequency of second most frequent base is high enough, save
                // the position in "ps".  Likewise if frequency
                // of first most frequent base is high enough and it's non-reference.
                //
                // Note that this doesn't allow for the case where the third most frequent base
                // is frequent enough.

                let x = refdata.refs[id].get(p);
                let c;
                if x == 0 {
                    c = b'A';
                } else if x == 1 {
                    c = b'C';
                } else if x == 2 {
                    c = b'G';
                } else {
                    c = b'T';
                }
                if (!freqs.is_empty() && freqs[0].0 >= MIN_ALT && freqs[0].3 != c)
                    || (freqs.len() > 1
                        && MIN_MULT * freqs[1].0 >= bases.len()
                        && freqs[1].0 >= MIN_ALT)
                {
                    ps.push(p);
                }
            }
            if ps.is_empty() {
                continue;
            }

            // Define types = { ( (contig base at each position in ps), index of contig ) },
            // and sort.

            let mut types = Vec::<(Vec<u8>, usize)>::new();
            for (loc, item) in all.iter().enumerate() {
                let mut t = Vec::<u8>::new();
                for &p in &ps {
                    let base = item.1[p];
                    t.push(base);
                }
                types.push((t, loc));
            }
            types.sort();

            // Traverse the types, grouping contigs that have an identical footprint at
            // the positions in ps.

            let mut keep = Vec::<(Vec<u8>, usize, f64, bool, Vec<String>)>::new();
            let mut i = 0;
            let mut have_ref = false;
            while i < types.len() {
                let j = next_diff1_2(&types, i);

                // Determine if the contigs equal reference at the positions in ps.

                let mut is_ref = true;
                for (z, p) in ps.iter().enumerate() {
                    let x = refdata.refs[id].get(*p);
                    let c;
                    if x == 0 {
                        c = b'A';
                    } else if x == 1 {
                        c = b'C';
                    } else if x == 2 {
                        c = b'G';
                    } else {
                        c = b'T';
                    }
                    if c != types[i].0[z] {
                        is_ref = false;
                    }
                }
                if (j - i >= MIN_ALT && j - i >= types.len() / MIN_MULT) || is_ref {
                    let mut q = Vec::<Vec<usize>>::new();
                    let mut barcodes = Vec::<String>::new();
                    for t in &types[i..j] {
                        let m = t.1;
                        q.push(all[m].2.clone());
                        barcodes.push(all[m].5.clone());
                    }
                    q.sort();
                    let (mut m, mut r) = (0, 0);
                    while r < q.len() {
                        let s = next_diff(&q, r);
                        m = max(m, s - r);
                        r = s;
                    }
                    let purity = 100.0 * (1.0 - m as f64 / q.len() as f64);
                    keep.push((types[i].0.clone(), j - i, purity, is_ref, barcodes));
                    if is_ref {
                        have_ref = true;
                    }
                }
                i = j;
            }

            if keep.len() > 1 || (!keep.is_empty() && !have_ref) {
                // Remove columns that are pure reference.

                let mut to_delete = vec![false; keep[0].0.len()];
                for (i, (&p, d)) in ps
                    .iter()
                    .take(keep[0].0.len())
                    .zip(to_delete.iter_mut())
                    .enumerate()
                {
                    let mut is_ref = true;
                    for j in &keep {
                        let x = refdata.refs[id].get(p);
                        let c;
                        if x == 0 {
                            c = b'A';
                        } else if x == 1 {
                            c = b'C';
                        } else if x == 2 {
                            c = b'G';
                        } else {
                            c = b'T';
                        }
                        if c != j.0[i] {
                            is_ref = false;
                        }
                    }
                    if is_ref {
                        *d = true;
                    }
                }
                erase_if(&mut ps, &to_delete);
                for j in &mut keep {
                    erase_if(&mut j.0, &to_delete);
                }
                let mut keep0 = Vec::<Vec<u8>>::new();
                for i in &keep {
                    keep0.push(i.0.clone());
                }
                keep0.sort();
                keep.sort_by(|a, b| a.partial_cmp(b).unwrap());

                // Save alternate references.

                for x in &keep {
                    if !x.3 {
                        let mut b = refdata.refs[id].clone();
                        for (&x0, &ps) in x.0.iter().zip(ps.iter()) {
                            let c;
                            if x0 == b'A' {
                                c = 0;
                            } else if x0 == b'C' {
                                c = 1;
                            } else if x0 == b'G' {
                                c = 2;
                            } else {
                                c = 3;
                            }
                            b.set_mut(ps, c);
                        }
                        res.1.push(AltRef {
                            donor: donor_id,
                            ref_id: id,
                            alt_seq: b,
                            support: x.1,
                            is_ref: x.3,
                        });
                    }
                }
            }
        }
    });
    let mut alt_refs = Vec::new();
    for mut r in results {
        alt_refs.append(&mut r.1);
    }
    alt_refs.sort();
    alt_refs
}

/// Update reference sequences for V segments by substituting in alt alleles if better.
/// Computational performance dubious because of full alt_refs traversal.
#[allow(clippy::needless_range_loop)]
pub(crate) fn sub_alts(
    refdata: &RefData,
    alt_refs: &[AltRef],
    clone_info: &mut [CloneInfo],
    exact_clonotypes: &mut [ExactClonotype],
) {
    for info in clone_info {
        for ((v_seg, tig), vsid) in zip(zip(&mut info.vs, &info.tigs), &mut info.vsids) {
            if v_seg.len() - REF_V_TRIM <= tig.len() {
                let mut errs = 0;
                for l in 0..v_seg.len() - REF_V_TRIM {
                    let x = v_seg.get(l);
                    let c;
                    if x == 0 {
                        c = b'A';
                    } else if x == 1 {
                        c = b'C';
                    } else if x == 2 {
                        c = b'G';
                    } else {
                        c = b'T';
                    }
                    if tig[l] != c {
                        errs += 1;
                    }
                }
                let mut donors = Vec::<usize>::new();
                let ex = &exact_clonotypes[info.clonotype_index];
                for m in 0..ex.clones.len() {
                    if ex.clones[m][0].donor_index.is_some() {
                        donors.push(ex.clones[m][0].donor_index.unwrap());
                    }
                }
                unique_sort(&mut donors);
                for donor in donors {
                    for m in 0..alt_refs.len() {
                        if alt_refs[m].donor == donor
                            && refdata.name[alt_refs[m].ref_id] == refdata.name[*vsid]
                            && alt_refs[m].alt_seq.len() - REF_V_TRIM <= tig.len()
                        {
                            let mut alt_errs = 0;
                            for l in 0..alt_refs[m].alt_seq.len() - REF_V_TRIM {
                                let x = alt_refs[m].alt_seq.get(l);
                                let c;
                                if x == 0 {
                                    c = b'A';
                                } else if x == 1 {
                                    c = b'C';
                                } else if x == 2 {
                                    c = b'G';
                                } else {
                                    c = b'T';
                                }
                                if tig[l] != c {
                                    alt_errs += 1;
                                }
                            }
                            if alt_errs < errs {
                                *v_seg = alt_refs[m].alt_seq.clone();
                                *vsid = alt_refs[m].ref_id;
                                let ex = &mut exact_clonotypes[info.clonotype_index];
                                for shared_tig_data in &mut ex.share {
                                    if shared_tig_data.seq == *tig {
                                        shared_tig_data.v_ref_id = alt_refs[m].ref_id;
                                        shared_tig_data.v_ref_id_donor = Some(m);
                                        shared_tig_data.v_ref_id_donor_donor = Some(donor);
                                        let mut alts = 0;
                                        let mut mm = m;
                                        while mm >= 1 {
                                            mm -= 1;
                                            if alt_refs[mm].donor == donor
                                                && alt_refs[mm].ref_id == alt_refs[m].ref_id
                                            {
                                                alts += 1;
                                            }
                                        }
                                        shared_tig_data.v_ref_id_donor_alt_id = Some(alts);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

const MIN_MULT: usize = 4;
const MIN_ALT: usize = 4;
