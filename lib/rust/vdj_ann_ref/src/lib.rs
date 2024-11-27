// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.

// This file contains code to make reference data.
//
// ◼ Experiment with building a reference from scratch.  This would be a better
// ◼ solution than ad hoc editing of a flawed reference.
//
// ◼ Document reference sequence requirements so that a customer who wishes to
// ◼ create a reference for a new species will know the conventions used by the
// ◼ code.

use io_utils::read_to_string_safe;
use vdj_ann::refx::{make_vdj_ref_data_core, RefData};

pub fn human_ref() -> &'static str {
    include_str!["../vdj_refs/human/fasta/regions.fa"]
}

pub fn human_ref_old() -> &'static str {
    include_str!["../vdj_refs_old/human/fasta/regions.fa"]
}

pub fn human_ref_2_0() -> &'static str {
    include_str!["../vdj_refs_2.0/human/fasta/regions.fa"]
}

pub fn human_ref_3_1() -> &'static str {
    include_str!["../vdj_refs_3.1/human/fasta/regions.fa"]
}

pub fn human_ref_4_0() -> &'static str {
    include_str!["../vdj_refs_4.0/human/fasta/regions.fa"]
}

pub fn human_ref_5_0() -> &'static str {
    include_str!["../vdj_refs_5.0/human/fasta/regions.fa"]
}

pub fn human_ref_7_0() -> &'static str {
    include_str!["../vdj_refs_7.0/human/fasta/regions.fa"]
}

pub fn mouse_ref() -> &'static str {
    include_str!["../vdj_refs/mouse/fasta/regions.fa"]
}

pub fn mouse_ref_old() -> &'static str {
    include_str!["../vdj_refs_old/mouse/fasta/regions.fa"]
}

pub fn mouse_ref_3_1() -> &'static str {
    include_str!["../vdj_refs_3.1/mouse/fasta/regions.fa"]
}

pub fn mouse_ref_4_0() -> &'static str {
    include_str!["../vdj_refs_4.0/mouse/fasta/regions.fa"]
}

pub fn mouse_ref_5_0() -> &'static str {
    include_str!["../vdj_refs_5.0/mouse/fasta/regions.fa"]
}

pub fn mouse_ref_7_0() -> &'static str {
    include_str!["../vdj_refs_7.0/mouse/fasta/regions.fa"]
}

// ids_to_use_opt: Optional hashSet of ids. If specified only reference
// entries with id in the HashSet is used to construct RefData
#[allow(clippy::too_many_arguments)]
pub fn make_vdj_ref_data(
    refdata: &mut RefData,
    imgt: bool,
    species: &str,
    extended: bool,
    is_tcr: bool,
    is_bcr: bool,
    human_supp_ref: &str,
    mouse_supp_ref: &str,
) {
    // Necessary for lifetime management of results from read_to_string_safe
    let x: String;
    let refx = match (imgt, species) {
        (false, "human") => human_ref(),
        (false, "mouse") => mouse_ref(),
        (true, "human") => {
            x = read_to_string_safe("vdj_IMGT_20170916-2.1.0/fasta/regions.fa");
            x.as_str()
        }
        (true, "mouse") => {
            x = read_to_string_safe("vdj_IMGT_mouse_20180723-2.2.0/fasta/regions.fa");
            x.as_str()
        }
        _ => panic!("Invalid species {}.", &species),
    };
    let ext_refx = if extended && !imgt {
        match species {
            "human" => human_supp_ref,
            "mouse" => mouse_supp_ref,
            _ => unreachable!("Invalid species {}.", &species),
        }
    } else {
        ""
    };
    assert!(!refx.is_empty(), "Reference file has zero length.");
    make_vdj_ref_data_core(refdata, refx, ext_refx, is_tcr, is_bcr, None);
}

#[cfg(test)]
mod tests {
    use super::*;

    // The following test checks for alignment of a D region.  This example was fixed by code
    // changes in March 2020.
    #[test]
    fn test_d_region_alignment() {
        use debruijn::dna_string::DnaString;
        use vdj_ann::annotate::annotate_seq;
        use vdj_ann::refx::make_vdj_ref_data_core;
        let seq = DnaString::from_acgt_bytes(
            b"GGAGGTGCGAATGACTCTGCTCTCTGTCCTGTCTCCTCATCTGCAAAATTAGGAAGCCTGTCTTGATTATCTCCAGGAA\
            CCTCCCACCTCTTCATTCCAGCCTCTGACAAACTCTGCACATTAGGCCAGGAGAAGCCCCCGAGCCAAGTCTCTTTTCTCATTCTC\
            TTCCAACAAGTGCTTGGAGCTCCAAGAAGGCCCCCTTTGCACTATGAGCAACCAGGTGCTCTGCTGTGTGGTCCTTTGTCTCCTGG\
            GAGCAAACACCGTGGATGGTGGAATCACTCAGTCCCCAAAGTACCTGTTCAGAAAGGAAGGACAGAATGTGACCCTGAGTTGTGAA\
            CAGAATTTGAACCACGATGCCATGTACTGGTACCGACAGGACCCAGGGCAAGGGCTGAGATTGATCTACTACTCACAGATAGTAAA\
            TGACTTTCAGAAAGGAGATATAGCTGAAGGGTACAGCGTCTCTCGGGAGAAGAAGGAATCCTTTCCTCTCACTGTGACATCGGCCC\
            AAAAGAACCCGACAGCTTTCTATCTCTGTGCCAGTAGTATTTTTCTTGCCGGGACAGGGGGCTGGAGCGGCACTGAAGCTTTCTTT\
            GGACAAGGCACCAGACTCACAGTTGTAGAGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGA",
        );
        let (refx, ext_refx) = (human_ref(), String::new());
        let (is_tcr, is_bcr) = (true, false);
        let mut refdata = RefData::new();
        make_vdj_ref_data_core(&mut refdata, refx, &ext_refx, is_tcr, is_bcr, None);
        let ann = annotate_seq(&seq, &refdata, true, false, true);
        let mut have_d = false;
        for a in ann {
            if refdata.is_d(a.ref_id as usize) {
                have_d = true;
            }
        }
        assert!(have_d, "\nFailed to find alignment of D region.\n");
    }
}
