//! Crate for dealing with dna sequences represented as a byte vector
#![deny(missing_docs)]

use std::borrow::Borrow;

const COMPLEMENT: [u8; 256] = make_complement_map();

const fn make_complement_map() -> [u8; 256] {
    let mut comp = [0; 256];
    // For some reason for isn't allowed in const contexts, but while is.
    let mut i = 0;
    while i < comp.len() {
        comp[i] = i as u8;
        i += 1;
    }
    const C1: &[u8; 15] = b"AGCTYRWSKMDVHBN";
    const C2: &[u8; 15] = b"TCGARYWSMKHBDVN";
    i = 0;
    while i < C1.len() {
        let (a, b) = (C1[i], C2[i]);
        comp[a as usize] = b;
        comp[a as usize + 32] = b + 32; // lowercase variants
        i += 1;
    }
    comp
}

/// Return complement of given DNA alphabet character (IUPAC alphabet supported).
const fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// An iterator over the reverse complement of given text (IUPAC alphabet supported).
pub fn revcomp_iter<C, T>(text: T) -> impl Iterator<Item = u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter().rev().map(|a| complement(*a.borrow()))
}

/// Calculate reverse complement of given text (IUPAC alphabet supported).
pub fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}
