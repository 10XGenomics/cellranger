// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use crate::innate::mark_innate;
use amino::nucleotide_to_aminoacid_sequence;
use debruijn::dna_string::DnaString;
use enclone_core::defs::{
    EncloneControl, ExactClonotype, Junction, OriginInfo, TigData, TigData0, TigData1,
};
use rayon::prelude::*;
use std::cmp::min;
use std::fmt::Write as _;
use vdj_ann::refx::RefData;
use vector_utils::{next_diff, next_diff1_2, next_diff1_3, reverse_sort};

pub fn create_exact_subclonotype_core(
    // inputs:
    tig_bc: &[Vec<TigData>],
    r: usize,
    s: usize,
) -> (Vec<TigData1>, Vec<Vec<TigData0>>) {
    let mut share = Vec::new();
    for m in 0..tig_bc[r].len() {
        // Form 5'-UTR consensus sequence.

        let mut utr = Vec::<u8>::new();
        let mut pos = 0;
        let mut last_calls = 0;
        loop {
            let mut calls = Vec::<(u8, u8)>::new(); // (base,qual)
            for this_tig_bc in &tig_bc[r..s] {
                if this_tig_bc[m].v_start > pos {
                    let p = this_tig_bc[m].v_start - pos - 1;
                    calls.push((this_tig_bc[m].full_seq[p], this_tig_bc[m].full_quals[p]));
                }
            }
            if calls.is_empty() || calls.len() < last_calls / 10 {
                break;
            }
            last_calls = calls.len();
            calls.sort_unstable();
            let mut callsx = Vec::<(usize, u8)>::new(); // (qual,base)
            let mut i = 0;
            while i < calls.len() {
                let j = next_diff1_2(&calls, i);
                let mut q = 0;
                for c in &calls[i..j] {
                    q += c.1 as usize;
                }
                callsx.push((q, calls[i].0));
                i = j;
            }
            reverse_sort(&mut callsx);
            utr.push(callsx[0].1);
            pos += 1;
        }
        utr.reverse();

        // Form constant region consensus sequence.

        let mut constx = Vec::<u8>::new();
        let mut pos = 0;
        let mut last_calls = 0;
        loop {
            let mut calls = Vec::<(u8, u8)>::new(); // (base,qual)
            for this_tig_bc in &tig_bc[r..s] {
                if this_tig_bc[m].j_stop + pos < this_tig_bc[m].full_seq.len() {
                    let p = this_tig_bc[m].j_stop + pos;
                    calls.push((this_tig_bc[m].full_seq[p], this_tig_bc[m].full_quals[p]));
                }
            }
            if calls.is_empty() || calls.len() < last_calls / 10 {
                break;
            }
            last_calls = calls.len();
            calls.sort_unstable();
            let mut callsx = Vec::<(usize, u8)>::new(); // (qual,base)
            let mut i = 0;
            while i < calls.len() {
                let j = next_diff1_2(&calls, i);
                let mut q = 0;
                for c in &calls[i..j] {
                    q += c.1 as usize;
                }
                callsx.push((q, calls[i].0));
                i = j;
            }
            reverse_sort(&mut callsx);
            constx.push(callsx[0].1);
            pos += 1;
        }

        // Form full sequence.

        let mut full = utr.clone();
        let mut z = tig_bc[r][m].seq().to_vec();
        full.append(&mut z);
        full.append(&mut constx);

        // Note that here we are taking the first entry (r), sort of assuming that all the entries
        // are the same, which in principle they should be, but this is not actually always true.
        // However this is hard to fix.

        let aa = nucleotide_to_aminoacid_sequence(tig_bc[r][m].seq(), 0);
        let mut d_start = None;
        if tig_bc[r][m].d_start.is_some() {
            d_start = Some(tig_bc[r][m].d_start.unwrap() + utr.len() - tig_bc[r][m].v_start);
        }
        share.push(TigData1 {
            cdr3_dna: tig_bc[r][m].cdr3_dna.clone(),
            seq: tig_bc[r][m].seq().to_vec(),
            seq_del: tig_bc[r][m].seq().to_vec(), // may get changed later
            seq_del_amino: tig_bc[r][m].seq().to_vec(), // may get changed later
            ins: Vec::new(),                      // may get changed later
            aa_mod_indel: aa,                     // may get changed later
            full_seq: full,
            v_start: utr.len(),
            v_stop: tig_bc[r][m].v_stop + utr.len() - tig_bc[r][m].v_start,
            v_stop_ref: tig_bc[r][m].v_stop_ref,
            d_start,
            j_start: tig_bc[r][m].j_start + utr.len() - tig_bc[r][m].v_start,
            j_start_ref: tig_bc[r][m].j_start_ref,
            j_stop: tig_bc[r][m].j_stop + utr.len() - tig_bc[r][m].v_start,
            u_ref_id: tig_bc[r][m].u_ref_id,
            v_ref_id: tig_bc[r][m].v_ref_id,
            v_ref_id_donor: None,
            v_ref_id_donor_alt_id: None,
            v_ref_id_donor_donor: None,
            d_ref_id: tig_bc[r][m].d_ref_id,
            j_ref_id: tig_bc[r][m].j_ref_id,
            c_ref_id: tig_bc[r][m].c_ref_id,
            fr1_start: tig_bc[r][m].fr1_start,
            fr2_start: tig_bc[r][m].fr2_start,
            fr3_start: tig_bc[r][m].fr3_start,
            cdr1_start: tig_bc[r][m].cdr1_start,
            cdr2_start: tig_bc[r][m].cdr2_start,
            cdr3_aa: tig_bc[r][m].cdr3_aa.clone(),
            cdr3_start: tig_bc[r][m].cdr3_start,
            left: tig_bc[r][m].left,
            chain_type: tig_bc[r][m].chain_type.clone(),
            annv: tig_bc[r][m].annv.clone(),
            // these get set when making CloneInfo objects:
            vs: DnaString::new(),
            js: DnaString::new(),
            // iNKT and MAIT annotations
            inkt_alpha_chain_gene_match: false,
            inkt_alpha_chain_junction_match: false,
            inkt_beta_chain_gene_match: false,
            inkt_beta_chain_junction_match: false,
            mait_alpha_chain_gene_match: false,
            mait_alpha_chain_junction_match: false,
            mait_beta_chain_gene_match: false,
            mait_beta_chain_junction_match: false,
            jun: Junction::default(),
        });
    }
    let mut clones = Vec::new();
    for this_tig_bc in &tig_bc[r..s] {
        let mut x = Vec::<TigData0>::new();
        for tbc in this_tig_bc {
            x.push(TigData0 {
                quals: tbc.quals.clone(),
                v_start: tbc.v_start,
                j_stop: tbc.j_stop,
                c_start: tbc.c_start,
                full_seq: tbc.full_seq.clone(),
                barcode: tbc.barcode.clone(),
                tigname: tbc.tigname.clone(),
                dataset_index: tbc.dataset_index,
                donor_index: tbc.donor_index,
                umi_count: tbc.umi_count,
                read_count: tbc.read_count,
                v_ref_id: tbc.v_ref_id,
            });
        }
        clones.push(x);
    }
    (share, clones)
}

/// Find exact subclonotypes.
pub fn find_exact_subclonotypes(
    ctl: &EncloneControl,
    tig_bc: &[Vec<TigData>],
    refdata: &RefData,
) -> Vec<ExactClonotype> {
    let mut r: usize = 0;
    let mut groups = Vec::<(usize, usize)>::new();
    while r < tig_bc.len() {
        let mut s = r + 1;
        while s < tig_bc.len() {
            let mut ok = true;
            if tig_bc[s].len() != tig_bc[r].len() {
                break;
            }
            for m in 0..tig_bc[r].len() {
                let (cid1, cid2) = (tig_bc[r][m].c_ref_id, tig_bc[s][m].c_ref_id);
                if tig_bc[s][m].cdr3_dna != tig_bc[r][m].cdr3_dna
                    || tig_bc[s][m].seq() != tig_bc[r][m].seq()

                    // Working around a bug here.  See above for explanation.

                    // || cid1 != cid2 {

                    || ( cid1.is_none() && cid2.is_some() ) || ( cid1.is_some() && cid2.is_none() )
                    || ( cid1.is_some() && cid2.is_some()
                        && refdata.name[cid1.unwrap()] != refdata.name[cid2.unwrap()] )

                    || ( cid1.is_some() && cid2.is_some()
                        && tig_bc[r][m].c_start.unwrap() + tig_bc[s][m].j_stop
                        < tig_bc[s][m].c_start.unwrap() + tig_bc[r][m].j_stop )

                    // Check for different donors if MIX_DONORS specified on command line.
                    // Note funky redundancy in checking each chain

                    || ( !ctl.cr_opt.mix_donors
                        && tig_bc[r][m].donor_index != tig_bc[s][m].donor_index )
                {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }
            s += 1;
        }
        groups.push((r, s));
        r = s;
    }

    let mut exact_clonotypes: Vec<_> = groups
        .into_par_iter()
        .filter_map(|(r, s)| {
            let mut bc = (r..s)
                .map(|t| (tig_bc[t][0].barcode.as_str(), t))
                .collect::<Vec<_>>();
            bc.sort_unstable();

            // Create the exact subclonotype.
            let (share, clones) = create_exact_subclonotype_core(tig_bc, r, s);

            // Save exact subclonotype.

            if !share.is_empty() && !clones.is_empty() {
                return Some(ExactClonotype { share, clones });
            }
            None
        })
        .collect();

    // Fill in iNKT and MAIT annotations.

    mark_innate(refdata, &mut exact_clonotypes);

    exact_clonotypes
}

/// Look for barcode reuse.  The primary purpose of this is to detect instances where two
/// datasets were obtained from the same cDNA (from the same GEM well).
pub fn check_for_barcode_reuse(
    origin_info: &OriginInfo,
    tig_bc: &[Vec<TigData>],
) -> Result<(), String> {
    {
        const MIN_REUSE_FRAC_TO_SHOW: f64 = 0.25;
        let mut all = Vec::<(&str, usize, usize)>::new();
        let mut total = vec![0; origin_info.dataset_id.len()];
        for (i, tig_i) in tig_bc.iter().enumerate() {
            all.push((tig_i[0].barcode.as_str(), tig_i[0].dataset_index, i));
            total[tig_i[0].dataset_index] += 1;
        }
        all.par_sort();
        let mut reuse = Vec::<(usize, usize)>::new();
        let mut i = 0;
        while i < all.len() {
            let j = next_diff1_3(&all, i);
            for k1 in i..j {
                for k2 in k1 + 1..j {
                    // We require identity on one cdr3_aa.  That seems to be about the right amount
                    // of concurrence that should be required.  If two datasets arise from the same
                    // barcode in the same cDNA (from the same GEM well), there can still be
                    // assembly differences.

                    let mut ok = false;
                    let (u1, u2) = (all[k1].2, all[k2].2);
                    for v1 in 0..tig_bc[u1].len() {
                        for v2 in 0..tig_bc[u2].len() {
                            if tig_bc[u1][v1].cdr3_aa == tig_bc[u2][v2].cdr3_aa {
                                ok = true;
                            }
                        }
                    }
                    if ok {
                        reuse.push((all[k1].1, all[k2].1));
                    }
                }
            }
            i = j;
        }
        reuse.sort_unstable();
        let mut i = 0;
        let mut found = false;
        let mut msg = String::new();
        while i < reuse.len() {
            let j = next_diff(&reuse, i);
            let n = j - i;
            let (l1, l2) = (reuse[i].0, reuse[i].1);
            let (n1, n2) = (total[l1], total[l2]);
            let frac = n as f64 / min(n1, n2) as f64;
            if frac >= MIN_REUSE_FRAC_TO_SHOW {
                if !found {
                    found = true;
                    msg += "\nSignificant barcode reuse detected.  If at least 25% of the barcodes \
                        in one dataset\nare present in another dataset, is is likely that two datasets \
                        arising from the\nsame library were included as input to enclone.  Since this \
                        would normally occur\nonly by accident, enclone exits.  \
                        If you wish to override this behavior,\nplease rerun with the argument \
                        ACCEPT_REUSE.\n\nHere are the instances of reuse that were observed:\n\n";
                }
                writeln!(
                    msg,
                    "{}, {} ==> {} of {}, {} barcodes ({:.1}%)",
                    origin_info.dataset_id[l1],
                    origin_info.dataset_id[l2],
                    n,
                    n1,
                    n2,
                    100.0 * frac
                )
                .unwrap();
            }
            i = j;
        }
        if found {
            return Err(msg);
        }
    }
    Ok(())
}
