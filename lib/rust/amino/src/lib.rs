//! amino
// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use debruijn::Mer;
use debruijn::dna_string::DnaString;

/// Test to see if a given DnaString has a start or stop codon at a given position.
pub fn have_start(b: &DnaString, j: usize) -> bool {
    let (a, g, t) = (0u8, 2u8, 3u8);
    if b.get(j) == a && b.get(j + 1) == t && b.get(j + 2) == g {
        return true;
    }
    false
}

/// Return whether position `j` is a stop codon
pub fn have_stop(b: &DnaString, j: usize) -> bool {
    let (a, g, t) = (0u8, 2u8, 3u8);
    if b.get(j) == t && b.get(j + 1) == a && b.get(j + 2) == g {
        return true;
    }
    if b.get(j) == t && b.get(j + 1) == a && b.get(j + 2) == a {
        return true;
    }
    if b.get(j) == t && b.get(j + 1) == g && b.get(j + 2) == a {
        return true;
    }
    false
}

/// Convert a codon to an amino acid
pub fn codon_to_aa(codon: &[u8]) -> u8 {
    assert!(codon.len() == 3);

    match codon {
        b"GGT" => b'G',
        b"GGC" => b'G',
        b"GGA" => b'G',
        b"GGG" => b'G',
        b"TGG" => b'W',
        b"TGT" => b'C',
        b"TGC" => b'C',
        b"TTT" => b'F',
        b"TTC" => b'F',
        b"TTA" => b'L',
        b"TTG" => b'L',
        b"CTT" => b'L',
        b"CTC" => b'L',
        b"CTA" => b'L',
        b"CTG" => b'L',
        b"ATT" => b'I',
        b"ATC" => b'I',
        b"ATA" => b'I',
        b"GTT" => b'V',
        b"GTC" => b'V',
        b"GTA" => b'V',
        b"GTG" => b'V',
        b"TCT" => b'S',
        b"TCC" => b'S',
        b"TCA" => b'S',
        b"TCG" => b'S',
        b"AGT" => b'S',
        b"AGC" => b'S',
        b"CCT" => b'P',
        b"CCC" => b'P',
        b"CCA" => b'P',
        b"CCG" => b'P',
        b"ACT" => b'T',
        b"ACC" => b'T',
        b"ACA" => b'T',
        b"ACG" => b'T',
        b"GCT" => b'A',
        b"GCC" => b'A',
        b"GCA" => b'A',
        b"GCG" => b'A',
        b"TAT" => b'Y',
        b"TAC" => b'Y',
        b"CAT" => b'H',
        b"CAC" => b'H',
        b"CAA" => b'Q',
        b"CAG" => b'Q',
        b"AAT" => b'N',
        b"AAC" => b'N',
        b"AAA" => b'K',
        b"AAG" => b'K',
        b"GAT" => b'D',
        b"GAC" => b'D',
        b"GAA" => b'E',
        b"GAG" => b'E',
        b"CGT" => b'R',
        b"CGC" => b'R',
        b"CGA" => b'R',
        b"CGG" => b'R',
        b"AGA" => b'R',
        b"AGG" => b'R',
        b"ATG" => b'M',
        b"TAG" => b'*',
        b"TAA" => b'*',
        b"TGA" => b'*',
        _ => panic!("Unexpected codon {}.", std::str::from_utf8(codon).unwrap()),
    }
}

/// Convert a given DNA sequence to amino acids, starting at a given position.
pub fn nucleotide_to_aminoacid_sequence(dna_seq: &[u8], start: usize) -> Vec<u8> {
    let mut a = Vec::<u8>::new();
    if dna_seq.len() >= 3 {
        for pos in (start..dna_seq.len() - 3 + 1).step_by(3) {
            if dna_seq[pos] == b'-' && dna_seq[pos + 1] == b'-' && dna_seq[pos + 2] == b'-' {
                a.push(b'-');
            } else {
                a.push(codon_to_aa(&dna_seq[pos..pos + 3]));
            }
        }
    }
    a
}
