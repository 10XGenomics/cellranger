// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Miscellaneous functions.

use enclone_core::defs::TigData;

use std::cmp::Ordering;

use vdj_ann::refx::RefData;

pub fn sort_tig_bc(tig_bc: &mut [Vec<TigData>], refdata: &RefData, mix_donors: bool) {
    tig_bc.sort_by(|x, y| -> Ordering {
        for i in 0..x.len() {
            // Order by number of chains.

            if i >= y.len() {
                return Ordering::Greater;
            }

            // Order by cdr3_dna.

            if x[i].cdr3_dna < y[i].cdr3_dna {
                return Ordering::Less;
            } else if x[i].cdr3_dna > y[i].cdr3_dna {
                return Ordering::Greater;

            // Order by chain length.
            } else if x[i].len < y[i].len {
                return Ordering::Less;
            } else if x[i].len > y[i].len {
                return Ordering::Greater;

            // Order by chain sequence.
            } else if x[i].seq() < y[i].seq() {
                return Ordering::Less;
            } else if x[i].seq() > y[i].seq() {
                return Ordering::Greater;
            }

            // Working around a bug here and below.  For TCR, there are two TRBC1 records in the,
            // reference, having a SNP after our primer, and we appear to pick one of
            // them at random.  Also not sure this fully respects the sort order.
            // And of course a customer could have the same feature in their reference.

            let (cid1, cid2) = (x[i].c_ref_id, y[i].c_ref_id);
            if cid1.is_none() && cid2.is_some() {
                return Ordering::Less;
            } else if cid2.is_none() && cid1.is_some() {
                return Ordering::Greater;

            // Order by constant region name.
            } else if cid1.is_some()
                && cid2.is_some()
                && refdata.name[cid1.unwrap()] < refdata.name[cid2.unwrap()]
            {
                return Ordering::Less;
            } else if cid1.is_some()
                && cid2.is_some()
                && refdata.name[cid1.unwrap()] > refdata.name[cid2.unwrap()]
            {
                return Ordering::Greater;

            // Order by JC delta.
            } else if x[i].c_start.is_some()
                && y[i].c_start.is_some()
                && x[i].c_start.unwrap() + y[i].j_stop < y[i].c_start.unwrap() + x[i].j_stop
            {
                return Ordering::Less;
            } else if x[i].c_start.is_some()
                && y[i].c_start.is_some()
                && x[i].c_start.unwrap() + y[i].j_stop > y[i].c_start.unwrap() + x[i].j_stop
            {
                return Ordering::Greater;

            // Order by donor if MIX_DONORS option used.
            } else if !mix_donors && x[i].donor_index < y[i].donor_index {
                return Ordering::Less;
            } else if !mix_donors && x[i].donor_index > y[i].donor_index {
                return Ordering::Greater;
            }
        }
        if x.len() < y.len() {
            return Ordering::Less;
        }
        Ordering::Equal
    });
}
