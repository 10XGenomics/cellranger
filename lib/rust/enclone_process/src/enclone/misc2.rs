// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Miscellaneous functions.

use super::innate::mark_innate;
use crate::core::defs::{
    BarcodeContigs, CellContig, EncloneControl, ExactClonotype, Junction, SharedContig,
};
use amino::nucleotide_to_aminoacid_sequence;
use debruijn::dna_string::DnaString;
use rayon::prelude::*;
use vdj_ann::refx::RefData;
use vector_utils::{next_diff1_2, reverse_sort};

pub(crate) fn create_exact_subclonotype_core(
    // inputs:
    tig_bc: &[BarcodeContigs],
    r: usize,
    s: usize,
) -> (Vec<SharedContig>, Vec<Vec<CellContig>>) {
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
        share.push(SharedContig {
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
        let mut x = Vec::<CellContig>::new();
        for tbc in this_tig_bc {
            x.push(CellContig {
                quals: tbc.quals.clone(),
                barcode: tbc.barcode.clone(),
                tigname: tbc.tigname.clone(),
                dataset_index: tbc.dataset_index,
                donor_index: tbc.donor_index,
                umi_count: tbc.umi_count,
                read_count: tbc.read_count,
            });
        }
        clones.push(x);
    }
    (share, clones)
}

/// Find exact subclonotypes.
pub(crate) fn find_exact_subclonotypes(
    ctl: &EncloneControl,
    tig_bc: &[BarcodeContigs],
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
