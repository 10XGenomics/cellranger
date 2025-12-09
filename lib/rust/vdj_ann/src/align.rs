#![expect(missing_docs)]
use bio::alignment::Alignment;
use bio::alignment::pairwise::Aligner;
use debruijn::dna_string::DnaString;

/// Return a "standard" affine alignment of x to y.  This is intended to be
/// applied to the case where x is to be fully aligned to part of y.
pub fn affine_align(x: &DnaString, y: &DnaString) -> Alignment {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let mut aligner = Aligner::new(-6, -1, &score);
    aligner.semiglobal(&x.to_ascii_vec(), &y.to_ascii_vec())
}
