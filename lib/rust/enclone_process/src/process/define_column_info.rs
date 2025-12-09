// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use crate::core::defs::{ColInfo, ExactClonotype};
use vector_utils::{erase_if, make_freq};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define reference segment identifiers, one per chain.  For the V segment, there is
// also an optional donor reference sequence alignment.  (For now.  We might do the
// same thing later for J segments.)
// ◼ 1. We vote using the number of cells, whereas a better way would be to test all
// ◼    the alternatives to find the best match.
// ◼ 2. Defining a constant region identifier for the entire clonotype is
// ◼    biologically dubious.
// ◼ 3. Maybe we only need to do this for pass 2.

pub(crate) fn define_column_info(
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    mat: Vec<Vec<Option<usize>>>,
) -> ColInfo {
    let cols = mat.len();

    // Compute reference info.

    let mut uids = vec![None; cols];
    let mut vids = vec![0; cols];
    let mut vpids = vec![None; cols];
    let mut vpids_d = vec![None; cols];
    let mut vpids_a = vec![None; cols];
    let mut dids = vec![None; cols];
    let mut jids = vec![0; cols];
    let mut cids = vec![None; cols];
    for (m, ((uids, (vids, (vpids, (vpids_d, vpids_a)))), (dids, (jids, cids)))) in mat
        .iter()
        .zip(
            uids.iter_mut()
                .zip(
                    vids.iter_mut().zip(
                        vpids
                            .iter_mut()
                            .zip(vpids_d.iter_mut().zip(vpids_a.iter_mut())),
                    ),
                )
                .zip(dids.iter_mut().zip(jids.iter_mut().zip(cids.iter_mut()))),
        )
        .take(cols)
    {
        let mut u = Vec::<usize>::new();
        let mut v = Vec::<usize>::new();
        let mut vp = Vec::<(usize, Option<usize>, Option<usize>, Option<usize>)>::new();
        let mut d = Vec::<usize>::new();
        let mut j = Vec::<usize>::new();
        let mut c = Vec::<usize>::new();
        for (&clonotype_id, &m) in exacts.iter().zip(m.iter()) {
            let ex = &exact_clonotypes[clonotype_id];
            if let Some(m) = m {
                let x = &ex.share[m];
                let ncells = ex.ncells();
                if let Some(u_ref_id) = x.u_ref_id {
                    u.resize(u.len() + ncells, u_ref_id);
                }
                // This is not actually correct.  It copies the consensus V gene assignment
                // for an exact subclonotype, rather than fetch the per cell entries.  However
                // it would be very rare for this to make a difference.
                v.resize(v.len() + ncells, x.v_ref_id);
                vp.resize(
                    vp.len() + ncells,
                    (
                        x.v_ref_id,
                        x.v_ref_id_donor,
                        x.v_ref_id_donor_donor,
                        x.v_ref_id_donor_alt_id,
                    ),
                );
                if let Some(d_ref_id) = x.d_ref_id {
                    d.resize(d.len() + ncells, d_ref_id);
                }
                j.resize(j.len() + ncells, x.j_ref_id);
                if let Some(c_ref_id) = x.c_ref_id {
                    c.resize(c.len() + ncells, c_ref_id);
                }
            }
        }
        u.sort_unstable();
        v.sort_unstable();
        vp.sort();
        d.sort_unstable();
        j.sort_unstable();
        c.sort_unstable();
        let mut uf = Vec::<(u32, usize)>::new();
        make_freq(&u, &mut uf);
        if !uf.is_empty() {
            *uids = Some(uf[0].1);
        }
        let mut vf = Vec::<(u32, usize)>::new();
        make_freq(&v, &mut vf);
        *vids = vf[0].1;
        let mut to_delete = vec![false; vp.len()];
        for i in 0..vp.len() {
            if vp[i].0 != *vids {
                to_delete[i] = true;
            }
        }
        erase_if(&mut vp, &to_delete);
        let mut vpf = Vec::<(u32, (usize, Option<usize>, Option<usize>, Option<usize>))>::new();
        make_freq(&vp, &mut vpf);
        *vpids = (vpf[0].1).1;
        *vpids_d = (vpf[0].1).2;
        *vpids_a = (vpf[0].1).3;
        let mut df = Vec::<(u32, usize)>::new();
        make_freq(&d, &mut df);
        if !df.is_empty() {
            *dids = Some(df[0].1);
        }
        let mut jf = Vec::<(u32, usize)>::new();
        make_freq(&j, &mut jf);
        *jids = jf[0].1;
        let mut cf = Vec::<(u32, usize)>::new();
        make_freq(&c, &mut cf);
        if !cf.is_empty() {
            *cids = Some(cf[0].1);
        }
    }

    // Compute seqss and seqss_amino.

    let mut seqss = Vec::<Vec<Vec<u8>>>::new();
    let nexacts = exacts.len();
    for cx in mat.iter().take(cols) {
        let mut seqs = Vec::<Vec<u8>>::new();
        for (&m, &e) in cx.iter().zip(exacts.iter()).take(nexacts) {
            if let Some(m) = m {
                seqs.push(exact_clonotypes[e].share[m].seq_del.clone());
            } else {
                seqs.push(Vec::<u8>::new());
            }
        }
        seqss.push(seqs.clone());
    }

    ColInfo {
        uids,
        vids,
        vpids,
        dids,
        jids,
        cids,
        seqss,
        mat,
    }
}
