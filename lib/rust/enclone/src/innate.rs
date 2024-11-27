// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Functions relating to the identification if iNKT and MAIT cells.

// species: return "human" or "mouse" or "unknown", based on a 60-base perfect match between
// the TRAC or IGHM sequence in the provided reference sequences, and internally provided reference
// sequences for human and mouse, for those regions.

use enclone_core::defs::ExactClonotype;
use string_utils::TextUtils;
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, reverse_sort, unique_sort};

pub fn species(refdata: &RefData) -> &'static str {
    let mut my_trac = Vec::<Vec<u8>>::new();
    for i in 0..refdata.refs.len() {
        if refdata.name[i].starts_with("TRAC") || refdata.name[i].starts_with("IGHM") {
            my_trac.push(refdata.refs[i].to_ascii_vec());
        }
    }
    const K: usize = 60;
    let mut kmers = Vec::<&[u8]>::new();
    for tr in &my_trac {
        if tr.len() >= K {
            for j in 0..=tr.len() - K {
                kmers.push(&tr[j..j + K]);
            }
        }
    }
    unique_sort(&mut kmers);
    let mut counts = Vec::<(usize, &'static str)>::new();
    for pass in 1..=2 {
        let mut count = 0;
        let species = if pass == 1 { "human" } else { "mouse" };

        // Build trac.  This is the concatenation, with single space separation, of the all
        // the human (pass = 1) or mouse (pass = 2) reference sequences that contain
        // |TRAC or |IGHM, for particular versions of these reference sequences (and probably
        // that choice doesn't matter much).

        let trac = if pass == 1 {
            b"GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGGCCG\
            TTGGCTGCCTCGCACAGGACTTCCTTCCCGACTCCATCACTTTCTCCTGGAAATACAAGAACAACTCTGACATCA\
            GCAGCACCCGGGGCTTCCCATCAGTCCTGAGAGGGGGCAAGTACGCAGCCACCTCACAGGTGCTGCTGCCTTCCA\
            AGGACGTCATGCAGGGCACAGACGAACACGTGGTGTGCAAAGTCCAGCACCCCAACGGCAACAAAGAAAAGAACG\
            TGCCTCTTCCAGTGATTGCTGAGCTGCCTCCCAAAGTGAGCGTCTTCGTCCCACCCCGCGACGGCTTCTTCGGCA\
            ACCCCCGCAAGTCCAAGCTCATCTGCCAGGCCACGGGTTTCAGTCCCCGGCAGATTCAGGTGTCCTGGCTGCGCG\
            AGGGGAAGCAGGTGGGGTCTGGCGTCACCACGGACCAGGTGCAGGCTGAGGCCAAAGAGTCTGGGCCCACGACCT\
            ACAAGGTGACCAGCACACTGACCATCAAAGAGAGCGACTGGCTCGGCCAGAGCATGTTCACCTGCCGCGTGGATC\
            ACAGGGGCCTGACCTTCCAGCAGAATGCGTCCTCCATGTGTGTCCCCGATCAAGACACAGCCATCCGGGTCTTCG\
            CCATCCCCCCATCCTTTGCCAGCATCTTCCTCACCAAGTCCACCAAGTTGACCTGCCTGGTCACAGACCTGACCA\
            CCTATGACAGCGTGACCATCTCCTGGACCCGCCAGAATGGCGAAGCTGTGAAAACCCACACCAACATCTCCGAGA\
            GCCACCCCAATGCCACTTTCAGCGCCGTGGGTGAGGCCAGCATCTGCGAGGATGACTGGAATTCCGGGGAGAGGT\
            TCACGTGCACCGTGACCCACACAGACCTGCCCTCGCCACTGAAGCAGACCATCTCCCGGCCCAAGGGGGTGGCCC\
            TGCACAGGCCCGATGTCTACTTGCTGCCACCAGCCCGGGAGCAGCTGAACCTGCGGGAGTCGGCCACCATCACGT\
            GCCTGGTGACGGGCTTCTCTCCCGCGGACGTCTTCGTGCAGTGGATGCAGAGGGGGCAGCCCTTGTCCCCGGAGA\
            AGTATGTGACCAGCGCCCCAATGCCTGAGCCCCAGGCCCCAGGCCGGTACTTCGCCCACAGCATCCTGACCGTGT\
            CCGAAGAGGAATGGAACACGGGGGAGACCTACACCTGCGTGGTGGCCCATGAGGCCCTGCCCAACAGGGTCACCG\
            AGAGGACCGTGGACAAGTCCACCGGTAAACCCACCCTGTACAACGTGTCCCTGGTCATGTCCGACACAGCTGGCA\
            CCTGCTAC GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCA\
            GCGTGGCCGTTGGCTGCCTCGCACAGGACTTCCTTCCCGACTCCATCACTTTCTCCTGGAAATACAAGAACAACT\
            CTGACATCAGCAGCACCCGGGGCTTCCCATCAGTCCTGAGAGGGGGCAAGTACGCAGCCACCTCACAGGTGCTGC\
            TGCCTTCCAAGGACGTCATGCAGGGCACAGACGAACACGTGGTGTGCAAAGTCCAGCACCCCAACGGCAACAAAG\
            AAAAGAACGTGCCTCTTCCAGTGATTGCTGAGCTGCCTCCCAAAGTGAGCGTCTTCGTCCCACCCCGCGACGGCT\
            TCTTCGGCAACCCCCGCAAGTCCAAGCTCATCTGCCAGGCCACGGGTTTCAGTCCCCGGCAGATTCAGGTGTCCT\
            GGCTGCGCGAGGGGAAGCAGGTGGGGTCTGGCGTCACCACGGACCAGGTGCAGGCTGAGGCCAAAGAGTCTGGGC\
            CCACGACCTACAAGGTGACCAGCACACTGACCATCAAAGAGAGCGACTGGCTCGGCCAGAGCATGTTCACCTGCC\
            GCGTGGATCACAGGGGCCTGACCTTCCAGCAGAATGCGTCCTCCATGTGTGTCCCCGATCAAGACACAGCCATCC\
            GGGTCTTCGCCATCCCCCCATCCTTTGCCAGCATCTTCCTCACCAAGTCCACCAAGTTGACCTGCCTGGTCACAG\
            ACCTGACCACCTATGACAGCGTGACCATCTCCTGGACCCGCCAGAATGGCGAAGCTGTGAAAACCCACACCAACA\
            TCTCCGAGAGCCACCCCAATGCCACTTTCAGCGCCGTGGGTGAGGCCAGCATCTGCGAGGATGACTGGAATTCCG\
            GGGAGAGGTTCACGTGCACCGTGACCCACACAGACCTGCCCTCGCCACTGAAGCAGACCATCTCCCGGCCCAAGG\
            GGGTGGCCCTGCACAGGCCCGATGTCTACTTGCTGCCACCAGCCCGGGAGCAGCTGAACCTGCGGGAGTCGGCCA\
            CCATCACGTGCCTGGTGACGGGCTTCTCTCCCGCGGACGTCTTCGTGCAGTGGATGCAGAGGGGGCAGCCCTTGT\
            CCCCGGAGAAGTATGTGACCAGCGCCCCAATGCCTGAGCCCCAGGCCCCAGGCCGGTACTTCGCCCACAGCATCC\
            TGACCGTGTCCGAAGAGGAATGGAACACGGGGGAGACCTACACCTGCGTGGTGGCCCATGAGGCCCTGCCCAACA\
            GGGTCACCGAGAGGACCGTGGACAAGTCCACCGAGGGGGAGGTGAGCGCCGACGAGGAGGGCTTTGAGAACCTGT\
            GGGCCACCGCCTCCACCTTCATCGTCCTCTTCCTCCTGAGCCTCTTCTACAGTACCACCGTCACCTTGTTCAAGG\
            TGAAA ATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCC\
            TATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTG\
            TGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTG\
            CAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGC\
            TGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCC\
            TCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGC"
        } else {
            b"AGAGTCAGTCCTTCCCAAATGTCTTCCCCCTCGTCTCCTGCGAGAGCCCCCTGTCTGATAAGAATCTGGTGGCCA\
            TGGGCTGCCTGGCCCGGGACTTCCTGCCCAGCACCATTTCCTTCACCTGGAACTACCAGAACAACACTGAAGTCA\
            TCCAGGGTATCAGAACCTTCCCAACACTGAGGACAGGGGGCAAGTACCTAGCCACCTCGCAGGTGTTGCTGTCTC\
            CCAAGAGCATCCTTGAAGGTTCAGATGAATACCTGGTATGCAAAATCCACTACGGAGGCAAAAACAAAGATCTGC\
            ATGTGCCCATTCCAGCTGTCGCAGAGATGAACCCCAATGTAAATGTGTTCGTCCCACCACGGGATGGCTTCTCTG\
            GCCCTGCACCACGCAAGTCTAAACTCATCTGCGAGGCCACGAACTTCACTCCAAAACCGATCACAGTATCCTGGC\
            TAAAGGATGGGAAGCTCGTGGAATCTGGCTTCACCACAGATCCGGTGACCATCGAGAACAAAGGATCCACACCCC\
            AAACCTACAAGGTCATAAGCACACTTACCATCTCTGAAATCGACTGGCTGAACCTGAATGTGTACACCTGCCGTG\
            TGGATCACAGGGGTCTCACCTTCTTGAAGAACGTGTCCTCCACATGTGCTGCCAGTCCCTCCACAGACATCCTAA\
            CCTTCACCATCCCCCCCTCCTTTGCCGACATCTTCCTCAGCAAGTCCGCTAACCTGACCTGTCTGGTCTCAAACC\
            TGGCAACCTATGAAACCCTGAATATCTCCTGGGCTTCTCAAAGTGGTGAACCACTGGAAACCAAAATTAAAATCA\
            TGGAAAGCCATCCCAATGGCACCTTCAGTGCTAAGGGTGTGGCTAGTGTTTGTGTGGAAGACTGGAATAACAGGA\
            AGGAATTTGTGTGTACTGTGACTCACAGGGATCTGCCTTCACCACAGAAGAAATTCATCTCAAAACCCAATGAGG\
            TGCACAAACATCCACCTGCTGTGTACCTGCTGCCACCAGCTCGTGAGCAACTGAACCTGAGGGAGTCAGCCACAG\
            TCACCTGCCTGGTGAAGGGCTTCTCTCCTGCAGACATCAGTGTGCAGTGGCTTCAGAGAGGGCAACTCTTGCCCC\
            AAGAGAAGTATGTGACCAGTGCCCCGATGCCAGAGCCTGGGGCCCCAGGCTTCTACTTTACCCACAGCATCCTGA\
            CTGTGACAGAGGAGGAATGGAACTCCGGAGAGACCTATACCTGTGTTGTAGGCCACGAGGCCCTGCCACACCTGG\
            TGACCGAGAGGACCGTGGACAAGTCCACTGGTAAACCCACACTGTACAATGTCTCCCTGATCATGTCTGACACAG\
            GCGGCACCTGCTAT AGAGTCAGTCCTTCCCAAATGTCTTCCCCCTCGTCTCCTGCGAGAGCCCCCTGTCTGATA\
            AGAATCTGGTGGCCATGGGCTGCCTGGCCCGGGACTTCCTGCCCAGCACCATTTCCTTCACCTGGAACTACCAGA\
            ACAACACTGAAGTCATCCAGGGTATCAGAACCTTCCCAACACTGAGGACAGGGGGCAAGTACCTAGCCACCTCGC\
            AGGTGTTGCTGTCTCCCAAGAGCATCCTTGAAGGTTCAGATGAATACCTGGTATGCAAAATCCACTACGGAGGCA\
            AAAACAAAGATCTGCATGTGCCCATTCCAGCTGTCGCAGAGATGAACCCCAATGTAAATGTGTTCGTCCCACCAC\
            GGGATGGCTTCTCTGGCCCTGCACCACGCAAGTCTAAACTCATCTGCGAGGCCACGAACTTCACTCCAAAACCGA\
            TCACAGTATCCTGGCTAAAGGATGGGAAGCTCGTGGAATCTGGCTTCACCACAGATCCGGTGACCATCGAGAACA\
            AAGGATCCACACCCCAAACCTACAAGGTCATAAGCACACTTACCATCTCTGAAATCGACTGGCTGAACCTGAATG\
            TGTACACCTGCCGTGTGGATCACAGGGGTCTCACCTTCTTGAAGAACGTGTCCTCCACATGTGCTGCCAGTCCCT\
            CCACAGACATCCTAACCTTCACCATCCCCCCCTCCTTTGCCGACATCTTCCTCAGCAAGTCCGCTAACCTGACCT\
            GTCTGGTCTCAAACCTGGCAACCTATGAAACCCTGAATATCTCCTGGGCTTCTCAAAGTGGTGAACCACTGGAAA\
            CCAAAATTAAAATCATGGAAAGCCATCCCAATGGCACCTTCAGTGCTAAGGGTGTGGCTAGTGTTTGTGTGGAAG\
            ACTGGAATAACAGGAAGGAATTTGTGTGTACTGTGACTCACAGGGATCTGCCTTCACCACAGAAGAAATTCATCT\
            CAAAACCCAATGAGGTGCACAAACATCCACCTGCTGTGTACCTGCTGCCACCAGCTCGTGAGCAACTGAACCTGA\
            GGGAGTCAGCCACAGTCACCTGCCTGGTGAAGGGCTTCTCTCCTGCAGACATCAGTGTGCAGTGGCTTCAGAGAG\
            GGCAACTCTTGCCCCAAGAGAAGTATGTGACCAGTGCCCCGATGCCAGAGCCTGGGGCCCCAGGCTTCTACTTTA\
            CCCACAGCATCCTGACTGTGACAGAGGAGGAATGGAACTCCGGAGAGACCTATACCTGTGTTGTAGGCCACGAGG\
            CCCTGCCACACCTGGTGACCGAGAGGACCGTGGACAAGTCCACTGAGGGGGAGGTGAATGCTGAGGAGGAAGGCT\
            TTGAGAACCTGTGGACCACTGCCTCCACCTTCATCGTCCTCTTCCTCCTGAGCCTCTTCTACAGCACCACCGTCA\
            CCCTGTTCAAGGTGAAA ACATCCAGAACCCAGAACCTGCTGTGTACCAGTTAAAAGATCCTCGGTCTCAGGACA\
            GCACCCTCTGCCTGTTCACCGACTTTGACTCCCAAATCAATGTGCCGAAAACCATGGAATCTGGAACGTTCATCA\
            CTGACAAAACTGTGCTGGACATGAAAGCTATGGATTCCAAGAGCAATGGGGCCATTGCCTGGAGCAACCAGACAA\
            GCTTCACCTGCCAAGATATCTTCAAAGAGACCAACGCCACCTACCCCAGTTCAGACGTTCCCTGTGATGCCACGT\
            TGACTGAGAAAAGCTTTGAAACAGATATGAACCTAAACTTTCAAAACCTGTCAGTTATGGGACTCCGAATCCTCC\
            TGCTGAAAGTAGCCGGATTTAACCTGCTCATGACGCTGAGGCTGTGGTCCAGT"
        };

        // Test the kmers.

        for i in 0..=trac.len() - K {
            let kmer = &trac[i..i + K];
            if bin_member(&kmers, &kmer) {
                count += 1;
            }
        }
        counts.push((count, species));
    }
    reverse_sort(&mut counts);
    if counts[0].0 == counts[1].0 {
        "unknown"
    } else {
        counts[0].1
    }
}

// innate_cdr3: for the given species and given class (iNKT or MAIT), return the list of CDR3_AA
// sequences that are known to occur for that class.  These are defined in files in this directory.

pub fn innate_cdr3(species: &str, class: &str) -> Vec<String> {
    assert!(class == "iNKT" || class == "MAIT");
    let mut json = String::new();
    if species == "human" && class == "iNKT" {
        json = include_str!["human_iNKT_CDR3.json"].to_string();
    } else if species == "human" && class == "MAIT" {
        json = include_str!["human_MAIT_CDR3.json"].to_string();
    }
    let mut cdr3 = Vec::<String>::new();
    for line in json.lines() {
        if line.contains("\"cdr3\": ") {
            cdr3.push(line.after("\"cdr3\": ").between("\"", "\"").to_string());
        }
    }
    unique_sort(&mut cdr3);
    cdr3
}

// mark_innate: for each exact subclonotype, fill in iNKT and MAIT fields.

pub fn mark_innate(refdata: &RefData, ex: &mut Vec<ExactClonotype>) {
    let species = species(refdata);
    let inkt_cdr3 = innate_cdr3(species, "iNKT");
    let mait_cdr3 = innate_cdr3(species, "MAIT");
    for e in ex {
        let (mut have_mait_tra, mut have_mait_trb) = (false, false);
        let (mut have_mait_tra_cdr3, mut have_mait_trb_cdr3) = (false, false);
        let (mut have_inkt_tra, mut have_inkt_trb) = (false, false);
        let (mut have_inkt_tra_cdr3, mut have_inkt_trb_cdr3) = (false, false);
        for share in &e.share {
            let mut vname = refdata.name[share.v_ref_id].as_str();
            if vname.contains('*') {
                vname = vname.before("*");
            }
            let mut jname = refdata.name[share.j_ref_id].as_str();
            if jname.contains('*') {
                jname = jname.before("*");
            }
            if species == "human" {
                if vname == "TRAV10" && jname == "TRAJ18" {
                    have_inkt_tra = true;
                }
                if vname == "TRBV25-1" {
                    have_inkt_trb = true;
                }
                if vname == "TRAV1-2"
                    && (jname == "TRAJ33" || jname == "TRAJ20" || jname == "TRAJ12")
                {
                    have_mait_tra = true;
                }
                if vname.starts_with("TRBV20") || vname.starts_with("TRBV6") {
                    have_mait_trb = true;
                }
            } else if species == "mouse" {
                if vname == "TRAV1" && jname == "TRAJ33" {
                    have_mait_tra = true;
                }
                if vname == "TRBV19"
                    || vname == "TRBV13-1"
                    || vname == "TRBV13-2"
                    || vname == "TRBV13-3"
                {
                    have_mait_trb = true;
                }
                if (vname == "TRAV11" || vname == "TRAV11D") && jname == "TRAJ18" {
                    have_inkt_tra = true;
                }
                if vname == "TRBV13-2" || vname == "TRBV1" || vname == "TRBV29" {
                    have_inkt_trb = true;
                }
            }
            if share.left {
                if bin_member(&inkt_cdr3, &share.cdr3_aa) {
                    have_inkt_trb_cdr3 = true;
                }
                if bin_member(&mait_cdr3, &share.cdr3_aa) {
                    have_mait_trb_cdr3 = true;
                }
            } else {
                if bin_member(&inkt_cdr3, &share.cdr3_aa) {
                    have_inkt_tra_cdr3 = true;
                }
                if bin_member(&mait_cdr3, &share.cdr3_aa) {
                    have_mait_tra_cdr3 = true;
                }
            }
        }
        for share in &mut e.share {
            share.inkt_alpha_chain_gene_match = have_inkt_tra;
            share.inkt_alpha_chain_junction_match = have_inkt_tra_cdr3;
            share.inkt_beta_chain_gene_match = have_inkt_trb;
            share.inkt_beta_chain_junction_match = have_inkt_trb_cdr3;
            share.mait_alpha_chain_gene_match = have_mait_tra;
            share.mait_alpha_chain_junction_match = have_mait_tra_cdr3;
            share.mait_beta_chain_gene_match = have_mait_trb;
            share.mait_beta_chain_junction_match = have_mait_trb_cdr3;
        }
    }
}
