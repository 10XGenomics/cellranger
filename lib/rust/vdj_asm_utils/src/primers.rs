// Specify inner primers.
//
// Usage: inner_primers( species, class )
//        - where species is either "human" or "mouse"
//        - and   class   is either "tcr"   or "bcr".
//
// Note that if we ever change the primers, we should append the new primers to
// this file (under each of the four categories), and not delete the old primers.
// In this way we would allow for the cases where the customer had data from the
// old primers, or data from the new primers, or a mixture of both (although that
// would correspond to a weird experimental design).
//
// These primers are found at
// https://support.10xgenomics.com/single-cell-vdj/library-prep/doc/technical-note-assay-scheme-and-configuration-of-chromium-single-cell-vdj-libraries
// This file has since been deleted but the primers can be found in another document.
//
// ◼ Hypothetically we could infer primers from the data rather than specify them
// ◼ at all.
//
// TODO: These primers are also defined in cellranger/vdj/constants.py
// for use in preflight checks and needs to be in sync with the list below.
// It would be better to maintain just one copy of this.

use crate::constants::PRIMER_EXT_LEN;
use tenkit2::pack_dna::reverse_complement;
use vdj_ann::refx::RefData;
use vector_utils::{contains_at, unique_sort};

pub fn inner_primers(species: &str, class: &str) -> Vec<Vec<u8>> {
    match (species, class) {
        ("human", "tcr") => {
            vec![
                b"AGTCTCTCAGCTGGTACACG".to_vec(),
                b"TCTGATGGCTCAAACACAGC".to_vec(),
            ]
        }
        ("human", "bcr") => {
            vec![
                b"GGGAAGTTTCTGGCGGTCA".to_vec(),
                b"GGTGGTACCCAGTTATCAAGCAT".to_vec(),
                b"GTGTCCCAGGTCACCATCAC".to_vec(),
                b"TCCTGAGGACTGTAGGACAGC".to_vec(),
                b"CACGCTGCTCGTATCCGA".to_vec(),
                b"TAGCTGCTGGCCGC".to_vec(),
                b"GCGTTATCCACCTTCCACTGT".to_vec(),
            ]
        }
        ("mouse", "tcr") => {
            vec![
                b"AGTCAAAGTCGGTGAACAGGCA".to_vec(),
                b"GGCCAAGCACACGAGGGTA".to_vec(),
            ]
        }
        ("mouse", "bcr") => {
            vec![
                b"TACACACCAGTGTGGCCTT".to_vec(),
                b"CAGGCCACTGTCACACCACT".to_vec(),
                b"CAGGTCACATTCATCGTGCCG".to_vec(),
                b"GAGGCCAGCACAGTGACCT".to_vec(),
                b"GCAGGGAAGTTCACAGTGCT".to_vec(),
                b"CTGTTTGAGATCAGTTTGCCATCCT".to_vec(),
                b"TGCGAGGTGGCTAGGTACTTG".to_vec(),
                b"CCCTTGACCAGGCATCC".to_vec(),
                b"AGGTCACGGAGGAACCAGTTG".to_vec(),
                b"GGCATCCCAGTGTCACCGA".to_vec(),
                b"AGAAGATCCACTTCACCTTGAAC".to_vec(),
                b"GAAGCACACGACTGAGGCAC".to_vec(),
            ]
        }
        _ => {
            panic!("Illegal arguments to inner_primers.");
        }
    }
}

pub fn outer_primers(species: &str, class: &str) -> Vec<Vec<u8>> {
    match (species, class) {
        ("human", "tcr") => {
            vec![
                b"TGAAGGCGTTTGCACATGCA".to_vec(),
                b"TCAGGCAGTATCTGGAGTCATTGAG".to_vec(),
            ]
        }
        ("human", "bcr") => {
            vec![
                b"CAGGGCACAGTCACATCCT".to_vec(),
                b"TGCTGGACCACGCATTTGTA".to_vec(),
                b"GGTTTTGTTGTCGACCCAGTCT".to_vec(),
                b"TTGTCCACCTTGGTGTTGCT".to_vec(),
                b"CATGACGTCCTTGGAAGGCA".to_vec(),
                b"TGTGGGACTTCCACTG".to_vec(),
                b"TTCTCGTAGTCTGCTTTGCTCAG".to_vec(),
            ]
        }
        ("mouse", "tcr") => {
            vec![
                b"CTGGTTGCTCCAGGCAATGG".to_vec(),
                b"TGTAGGCCTGAGGGTCCGT".to_vec(),
            ]
        }
        ("mouse", "bcr") => {
            vec![
                b"TCAGCACGGGACAAACTCTTCT".to_vec(),
                b"GCAGGAGACAGACTCTTCTCCA".to_vec(),
                b"AACTGGCTGCTCATGGTGT".to_vec(),
                b"TGGTGCAAGTGTGGTTGAGGT".to_vec(),
                b"TGGTCACTTGGCTGGTGGTG".to_vec(),
                b"CACTTGGCAGGTGAACTGTTTTCT".to_vec(),
                b"AACCTTCAAGGATGCTCTTGGGA".to_vec(),
                b"GGACAGGGATCCAGAGTTCCA".to_vec(),
                b"AGGTGACGGTCTGACTTGGC".to_vec(),
                b"GCTGGACAGGGCTCCATAGTT".to_vec(),
                b"GGCACCTTGTCCAATCATGTTCC".to_vec(),
                b"ATGTCGTTCATACTCGTCCTTGGT".to_vec(),
            ]
        }
        _ => {
            panic!("Illegal arguments to outer_primers.");
        }
    }
}

// For each enrichment primer, find each 40-mer on a constant region that ends with
// the reverse complement of the given primer.  You can call this on the inner or
// outer enrichment primers.  The 40-mers are represented as Vec<u8>s.

pub fn get_primer_exts(list_of_primers: &[Vec<u8>], refdata: &RefData) -> Vec<Vec<Vec<u8>>> {
    let mut exts = Vec::<Vec<Vec<u8>>>::with_capacity(list_of_primers.len());
    for primer in list_of_primers {
        let mut primer = primer.clone();
        reverse_complement(&mut primer);
        let mut x = Vec::<Vec<u8>>::new();
        for ref_idx in 0..refdata.refs.len() {
            if refdata.is_c(ref_idx) {
                let constant_region = refdata.refs[ref_idx].to_ascii_vec();
                for pos in 0..constant_region.len() {
                    if contains_at(&constant_region, &primer, pos)
                        && pos + primer.len() >= PRIMER_EXT_LEN
                    {
                        x.push(
                            constant_region
                                [pos + primer.len() - PRIMER_EXT_LEN..pos + primer.len()]
                                .to_vec(),
                        );
                    }
                }
            }
        }
        unique_sort(&mut x);
        exts.push(x);
    }
    exts
}
