//! Align homopolymer sequences.

use std::cmp::max;
use std::ops::Range;

/// Construct the score vector.
fn construct_score_vector(nt: u8, seq: &[u8]) -> Vec<i64> {
    const MATCH_SCORE: i64 = 1;
    const MISMATCH_SCORE: i64 = -9;
    let scoring_function = |x| {
        if x == nt {
            MATCH_SCORE
        } else {
            MISMATCH_SCORE
        }
    };

    // Construct the score vector.
    seq.iter()
        .map(|&x| scoring_function(x))
        .fold(Vec::with_capacity(seq.len()), |mut acc, x| {
            acc.push(x + max(0, *acc.last().unwrap_or(&0)));
            acc
        })
}

/// Find the best scoring 3' homopolymer run.
pub fn align_homopolymer_threeprime(nt: u8, seq: &[u8]) -> Range<usize> {
    if seq.is_empty() {
        return 0..0;
    }

    let score_vec = construct_score_vector(nt, seq);

    if score_vec.is_empty() || (*score_vec.last().unwrap() <= 0) {
        return 0..0;
    }

    // Backtrack
    let start = score_vec.iter().rposition(|&x| x < 0).map_or(0, |x| 1 + x);
    start..seq.len()
}

/// Count the number of 3' homopolymer matches, allowing for mismatches.
pub fn count_homopolymer_matches_threeprime(nt: u8, seq: &[u8]) -> usize {
    let range = align_homopolymer_threeprime(nt, seq);
    let subseq = &seq[range];
    subseq.iter().filter(|&&x| x == nt).count()
}

/// Find the best scoring homopolymer run.
pub fn align_homopolymer(nt: u8, seq: &[u8]) -> Range<usize> {
    if seq.is_empty() {
        return 0..0;
    }

    let score_vec = construct_score_vector(nt, seq);
    let max_score = *score_vec
        .iter()
        .max()
        .expect("The score vec did not have a max.");

    if max_score > 0 {
        // Backtrack
        // TODO: Return the longest match when there are multiple matches with the best score.
        let end = 1 + score_vec.iter().rposition(|&x| x == max_score).unwrap();
        let start = score_vec[0..end]
            .iter()
            .rposition(|&x| x < 0)
            .map_or(0, |x| 1 + x);
        start..end
    } else {
        0..0
    }
}

/// Count the number of matches in the longest homopolymer run, allowing for mismatches.
pub fn count_homopolymer_matches(nt: u8, seq: &[u8]) -> usize {
    let range = align_homopolymer(nt, seq);
    let subseq = &seq[range];
    subseq.iter().filter(|&&x| x == nt).count()
}

#[test]
fn test_align_homopolymer() {
    let seq = "AAA".as_bytes();
    assert_eq!(align_homopolymer(b'T', seq), 0..0);
    assert_eq!(count_homopolymer_matches(b'T', seq), 0);

    let seq = "ATATTTTTTTTTATTTTTTTTTAATTTTTTTTA".as_bytes();
    assert_eq!(align_homopolymer(b'T', seq), 3..22);
    assert_eq!(count_homopolymer_matches(b'T', seq), 18);
}

#[test]
fn test_align_homopolymer_threeprime() {
    let seq = "ATATTTTTTTTTATTTTTTTTTAATTTTTTTTTA".as_bytes();
    assert_eq!(align_homopolymer_threeprime(b'T', seq), 0..0);
    assert_eq!(count_homopolymer_matches_threeprime(b'T', seq), 0);

    let seq = "ATATTTTTTTTTATTTTTTTTTAATTTTTTTTTTA".as_bytes();
    assert_eq!(align_homopolymer_threeprime(b'T', seq), 24..35);
    assert_eq!(count_homopolymer_matches_threeprime(b'T', seq), 10);

    let seq = "ATATTTTTTTTTAATTTTTTTTTTATTTTTTTTTA".as_bytes();
    assert_eq!(align_homopolymer_threeprime(b'T', seq), 14..35);
    assert_eq!(count_homopolymer_matches_threeprime(b'T', seq), 19);
}
