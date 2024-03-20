//! Crate for dealing with dna sequences represented as a byte vector

use lazy_static::lazy_static;
use std::borrow::Borrow;
use std::iter::zip;

// Code copied from rust-bio
lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in zip(b"AGCTYRWSKMDVHBN", b"TCGARYWSMKHBDVN") {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

/// Return complement of given DNA alphabet character (IUPAC alphabet supported).
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

pub fn maybe_revcomp_iter<'a, C, T>(text: T, rev: bool) -> Box<dyn Iterator<Item = u8> + 'a>
where
    C: Borrow<u8> + 'a,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator + 'a,
{
    if rev {
        Box::new(text.into_iter().rev().map(|a| complement(*a.borrow())))
    } else {
        Box::new(text.into_iter().map(|a| *a.borrow()))
    }
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
