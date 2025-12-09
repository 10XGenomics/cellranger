// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

//! Sized, stack-allocated container for a short DNA sequence.

use crate::array::{ArrayContent, ByteArray};
use itertools::Itertools;
use std::iter::Iterator;
use std::str;

const UPPER_ACGTN: &[u8; 5] = b"ACGTN";
const N_BASE_INDEX: usize = 4;

#[derive(Clone, Copy, PartialOrd, Ord, PartialEq, Eq)]
pub struct SSeqContents;

impl ArrayContent for SSeqContents {
    /// Make sure that the input byte slice contains only
    /// "ACGTN" characters. Panics otherwise with an error
    /// message describing the position of the first character
    /// that is not an ACGTN.
    fn validate_bytes(seq: &[u8]) {
        for (i, &s) in seq.iter().enumerate() {
            assert!(
                UPPER_ACGTN.contains(&s),
                "Non ACGTN character {s} at position {i}"
            );
        }
    }
    fn expected_contents() -> &'static str {
        "An [ACGTN]* string"
    }
}

/// Fixed-sized container for a short DNA sequence, with capacity determined by type `N`.
/// Used as a convenient container for barcode or UMI sequences.
/// An `SSeqGen` is guaranteed to contain only "ACGTN" alphabets
pub type SSeqGen<const N: usize> = ByteArray<SSeqContents, N>;

/// Fixed-sized container for a short DNA sequence, up to 23bp in length.
/// Used as a convenient container for barcode or UMI sequences.
/// An `SSeq` is guaranteed to contain only "ACGTN" alphabets
pub type SSeq = SSeqGen<23>;

impl<const N: usize> SSeqGen<N> {
    /// Returns a byte slice of this sequence's contents.
    /// A synonym for as_bytes().
    pub fn seq(&self) -> &[u8] {
        self.as_bytes()
    }

    /// Returns a byte slice of this sequence's contents.
    /// A synonym for as_bytes().
    pub fn seq_mut(&mut self) -> &mut [u8] {
        self.as_mut_bytes()
    }

    /// Returns true if this sequence contains an N.
    pub fn has_n(&self) -> bool {
        self.iter().any(|&c| c == b'N' || c == b'n')
    }

    /// Returns true if this sequence is a homopolymer.
    pub fn is_homopolymer(&self) -> bool {
        assert!(!self.is_empty());
        self.iter().all(|&c| c == self.seq()[0])
    }

    /// Returns true if the last n characters of this sequence are the specified homopolymer.
    pub fn has_homopolymer_suffix(&self, c: u8, n: usize) -> bool {
        self.len() >= n && self.iter().rev().take(n).all(|&x| x == c)
    }

    /// Returns true if the last n characters of this sequence are T.
    pub fn has_polyt_suffix(&self, n: usize) -> bool {
        self.has_homopolymer_suffix(b'T', n)
    }

    /// Returns a 2-bit encoding of this sequence.
    pub fn encode_2bit_u32(self) -> u32 {
        let mut res: u32 = 0;
        assert!(self.len() <= 16);

        let seq = self.seq();
        for (bit_pos, str_pos) in (0..self.len()).rev().enumerate() {
            let byte: u32 = match seq[str_pos] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("non-ACGT sequence"),
            };

            let v = byte << (bit_pos * 2);

            res |= v;
        }

        res
    }

    /// Return the reverse-complement sequence.
    pub fn revcomp(&self) -> Self {
        self.iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'N' => b'N',
                other => unreachable!("non-ACGTN character with ASCII code {other} in {self}"),
            })
            .collect()
    }

    /// Iterator over sequences one substitution (hamming distance) away.
    pub fn one_hamming_iter(&self, opt: HammingIterOpt) -> SSeqOneHammingIter<N> {
        SSeqOneHammingIter::new(*self, opt)
    }

    /// Iterator over sequences one deletion away.
    pub fn one_deletion_iter(&self) -> impl Iterator<Item = Self> + '_ {
        (0..self.len()).map(move |i| {
            let mut seq = *self;
            seq.remove(i);
            seq
        })
    }

    /// Iterator over sequences one insertion away. The length of the sequence
    /// must be less than the capacity so that there is room for 1 insertion.
    pub fn one_insertion_iter(&self, opt: InsertionIterOpt) -> impl Iterator<Item = Self> + '_ {
        let last_index = N_BASE_INDEX
            + match opt {
                InsertionIterOpt::IncludeNBase => 1,
                InsertionIterOpt::ExcludeNBase => 0,
            };
        (0..=self.len()).flat_map(move |pos| {
            UPPER_ACGTN[0..last_index].iter().map(move |base| {
                let mut seq = *self;
                seq.insert_unchecked(pos, *base);
                seq
            })
        })
    }

    /// Iterator over sequences one edit distance away. The length of the sequence
    /// must be less than the capacity so that there is room for 1 insertion.
    pub fn one_edit_iter(
        &self,
        ham: HammingIterOpt,
        ins: InsertionIterOpt,
    ) -> impl Iterator<Item = Self> + '_ {
        self.one_hamming_iter(ham)
            .chain(self.one_deletion_iter())
            .chain(self.one_insertion_iter(ins))
    }
}

#[derive(Copy, Clone)]
pub enum InsertionIterOpt {
    /// Insert N bases in addition to ACGT
    IncludeNBase,
    /// Do not insert N bases. Only insert ACGT
    ExcludeNBase,
}

/// Helper struct used to represent a subset of ACGTN characters
#[derive(Clone, Copy)]
pub struct SSeqChars(SSeqGen<5>);

impl SSeqChars {
    pub fn new(chars: &[u8]) -> Self {
        SSeqChars(chars.iter().unique().collect())
    }
    pub fn actg() -> Self {
        SSeqChars(SSeqGen::from_bytes_unchecked(b"ACGT"))
    }
    pub fn actgn() -> Self {
        SSeqChars(SSeqGen::from_bytes_unchecked(b"ACGTN"))
    }
    pub fn n() -> Self {
        SSeqChars(SSeqGen::from_bytes_unchecked(b"N"))
    }
    pub fn none() -> Self {
        SSeqChars(SSeqGen::new())
    }
    pub fn contains(&self, char: &u8) -> bool {
        self.0.as_bytes().contains(char)
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

#[derive(Copy, Clone)]
pub enum HammingIterOpt {
    /// Mutate the N base. The allowed mutations at each base is ACGT
    MutateNBase,
    /// Don't mutate the N base. The allowed mutations for other bases are ACGT
    SkipNBase,
    Custom {
        /// Characters to skip.
        skip_chars: SSeqChars,
        /// Allowed mutations for each base.
        mutation_chars: SSeqChars,
    },
}

impl HammingIterOpt {
    fn skip_chars(self) -> SSeqChars {
        match self {
            HammingIterOpt::SkipNBase => SSeqChars::n(),
            HammingIterOpt::MutateNBase => SSeqChars::none(),
            HammingIterOpt::Custom { skip_chars, .. } => skip_chars,
        }
    }
    fn mutation_chars(self) -> SSeqChars {
        match self {
            HammingIterOpt::SkipNBase | HammingIterOpt::MutateNBase => SSeqChars::actg(),
            HammingIterOpt::Custom { mutation_chars, .. } => mutation_chars,
        }
    }
}

/// An iterator over sequences which are one hamming distance away
/// from an `SSeq`. `SSeq` is guaranteed to contain "ACGTN" alphabets.
/// Positions containing "N" or "n" are mutated or skipped
/// depending on the `HammingIterOpt`
pub struct SSeqOneHammingIter<const N: usize> {
    /// Original SSeq from which we need to generate values
    source: SSeqGen<N>,
    /// Index into SSeq where last base was mutated
    position: usize,
    /// Characters used for mutation
    chars: SSeqChars,
    /// The last base which was used
    chars_index: usize,
    /// Characters to skip
    skip_chars: SSeqChars,
}

impl<const N: usize> SSeqOneHammingIter<N> {
    fn new(sseq: SSeqGen<N>, opt: HammingIterOpt) -> Self {
        SSeqOneHammingIter {
            source: sseq,
            position: 0,
            chars: opt.mutation_chars(),
            chars_index: 0,
            skip_chars: opt.skip_chars(),
        }
    }
}

impl<const N: usize> Iterator for SSeqOneHammingIter<N> {
    type Item = SSeqGen<N>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.source.len() {
            return None;
        }
        let base_at_pos = self.source[self.position];
        let skip_pos = self.skip_chars.contains(&base_at_pos);
        if skip_pos || (self.chars_index >= self.chars.len()) {
            // Skip this position or we have exhausted all mutations for this position
            self.position += 1;
            self.chars_index = 0;
            self.next()
        } else if base_at_pos == self.chars.0[self.chars_index] {
            self.chars_index += 1;
            self.next()
        } else {
            let mut next_sseq = self.source;
            next_sseq[self.position] = self.chars.0[self.chars_index];
            self.chars_index += 1;
            Some(next_sseq)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let positions = self.source.len().saturating_sub(self.position);

        if !self.skip_chars.is_empty() {
            // The number of sequences remaining if none of the source characters
            // are in skip_chars = self.chars.len() mutations per position.
            // We assume that if the char_index is >= 0 we may have skipped one already for this
            // position.
            let no_skip_remaining = positions
                .saturating_mul(self.chars.len().min(4))
                .saturating_sub(self.chars_index.saturating_sub(1));
            (0, Some(no_skip_remaining))
        } else {
            // The lower bound is the number of positions times self.chars.len()-1, with the
            // assumption that the matching character hasn't yet been passed.
            // The upper bound is if they're all Ns so we get 4 for each
            // position.
            (
                positions
                    .saturating_mul(self.chars.len().saturating_sub(1))
                    .saturating_sub(self.chars_index),
                Some(
                    positions
                        .saturating_mul(self.chars.len().min(4))
                        .saturating_sub(self.chars_index.saturating_sub(1)),
                ),
            )
        }
    }
}

#[cfg(test)]
mod sseq_test {
    use super::*;
    use bincode;
    use itertools::{Itertools, assert_equal};
    use proptest::collection::vec;
    use proptest::{prop_assert, prop_assert_eq, proptest};
    use std::cmp::Ordering;
    use std::collections::HashSet;
    use std::iter::FromIterator;

    #[test]
    fn sort_test() {
        let s1: &[u8] = b"ACNGTA";
        let s2 = b"TAGTCGGC";
        let s3 = b"CATC";
        let s4 = b"TGTG";
        let s5 = b"";
        let s6 = b"A";
        let s7 = b"AACCATAGCCGGNATC";
        let s8 = b"GAACNAGNTGGA";

        let mut seqs = [s1, s2, s3, s4, s5, s6, s7, s8];
        let mut sseqs: Vec<SSeq> = seqs.iter().map(|x| SSeq::from_bytes(x)).collect();

        seqs.sort();
        sseqs.sort();

        for i in 0..seqs.len() {
            assert_eq!(seqs[i], sseqs[i].seq());
        }
    }

    #[test]
    fn size_hint_test() {
        let s1: &[u8] = b"ACNGTA";
        let s2 = b"CATC";
        let s3 = b"T";
        let s4 = b"";
        let s5 = b"A";
        let s6 = b"AA";
        let s7 = b"NNN";

        for seq in [s1, s2, s3, s4, s5, s6, s7] {
            for opt in [HammingIterOpt::MutateNBase, HammingIterOpt::SkipNBase] {
                let iter = SSeq::from_bytes(seq).one_hamming_iter(opt);
                let (lower, upper) = iter.size_hint();
                let actual = iter.count();
                assert!(lower <= actual);
                assert!(upper.unwrap() >= actual);
                let mut iter = SSeq::from_bytes(seq).one_hamming_iter(opt);
                if iter.next().is_some() {
                    let (lower, upper) = iter.size_hint();
                    let actual = iter.count();
                    assert!(lower <= actual);
                    assert!(upper.unwrap() >= actual);
                }
            }
        }
    }

    proptest! {
        #[test]
        fn prop_test_sort_sseq(
            ref seqs_str in vec("[ACGTN]{0, 23}", 0usize..=10usize),
        ) {
            let mut seqs = seqs_str.iter().map(|s| s.clone().into_bytes()).collect_vec();
            let mut sseqs: Vec<SSeq> = seqs.iter().map(|x| SSeq::from_bytes(x)).collect();

            seqs.sort();
            sseqs.sort();

            for i in 0..seqs.len() {
                assert_eq!(seqs[i], sseqs[i].seq());
            }
        }

        #[test]
        fn prop_test_into_iter(
            ref seq_str in "[ACGTN]{0, 23}",
        ) {
            let seq = SSeq::from_bytes(seq_str.as_bytes());
            assert_eq!(seq.len(), seq.into_iter().count());
        }
    }

    #[test]
    fn dna_encode() {
        let s1 = SSeq::from_bytes(b"AAAAA");
        assert_eq!(s1.encode_2bit_u32(), 0);

        let s1 = SSeq::from_bytes(b"AAAAT");
        assert_eq!(s1.encode_2bit_u32(), 3);

        let s1 = SSeq::from_bytes(b"AAACA");
        assert_eq!(s1.encode_2bit_u32(), 4);

        let s1 = SSeq::from_bytes(b"AACAA");
        assert_eq!(s1.encode_2bit_u32(), 16);

        let s1 = SSeq::from_bytes(b"AATA");
        assert_eq!(s1.encode_2bit_u32(), 12);
    }

    #[test]
    fn test_serde() {
        let seq = b"AGCTAGTCAGTCAGTA";
        let mut sseqs = Vec::new();
        for _ in 0..4 {
            let s = SSeq::from_bytes(seq);
            sseqs.push(s);
        }

        let mut buf = Vec::new();
        bincode::serialize_into(&mut buf, &sseqs).unwrap();
        let roundtrip: Vec<SSeq> = bincode::deserialize_from(&buf[..]).unwrap();
        assert_eq!(sseqs, roundtrip);
    }

    #[test]
    fn test_serde_json() {
        let seq = SSeq::from_bytes(b"AGCTAGTCAGTCAGTA");
        let json_str = serde_json::to_string(&seq).unwrap();
        assert_eq!(json_str, r#""AGCTAGTCAGTCAGTA""#);
    }

    proptest! {
        #[test]
        fn prop_test_serde_sseq(
            ref seq in "[ACGTN]{0, 23}",
        ) {
            let target = SSeq::from_bytes(seq.as_bytes());
            let encoded: Vec<u8> = bincode::serialize(&target).unwrap();
            let decoded: SSeq = bincode::deserialize(&encoded[..]).unwrap();
            prop_assert_eq!(target, decoded);
        }
        #[test]
        fn prop_test_serde_json_sseq(
            ref seq in "[ACGTN]{0, 23}",
        ) {
            let target = SSeq::from_bytes(seq.as_bytes());
            let encoded = serde_json::to_string_pretty(&target).unwrap();
            let decoded: SSeq = serde_json::from_str(&encoded).unwrap();
            prop_assert_eq!(target, decoded);
        }

        #[test]
        fn prop_test_sseq_push(
            ref seq1 in "[ACGTN]{0, 23}",
            ref seq2 in "[ACGTN]{0, 23}",
        ) {
            if seq1.len() + seq2.len() <= 23 {
                let mut s = SSeq::new();
                s.push(seq1.as_bytes());
                s.push(seq2.as_bytes());
                assert_eq!(s, seq1.as_bytes().iter().chain(seq2.as_bytes()).collect());
            }
        }
    }

    fn test_hamming_helper(seq: &str, opt: HammingIterOpt, n: u8) {
        let sseq = SSeq::from_bytes(seq.as_bytes());
        // Make sure that the hamming distance is 1 for all elements
        for neighbor in sseq.one_hamming_iter(opt) {
            assert_eq!(
                sseq.seq()
                    .iter()
                    .zip_eq(neighbor.seq().iter())
                    .filter(|(a, b)| a != b)
                    .count(),
                1
            );
        }
        // Make sure that the total number of elements is equal to what we expect.
        let non_n = sseq.seq().iter().filter(|&&s| s != n).count();
        let n_bases = sseq.len() - non_n;
        assert_eq!(
            sseq.one_hamming_iter(opt).collect::<HashSet<_>>().len(),
            match opt {
                HammingIterOpt::SkipNBase => 3 * non_n,
                HammingIterOpt::MutateNBase => 3 * non_n + 4 * n_bases,
                HammingIterOpt::Custom { .. } => unimplemented!(),
            }
        );
    }

    proptest! {
        #[test]
        fn prop_test_one_hamming_upper(
            seq in "[ACGTN]{0, 23}", // 0 and 23 are inclusive bounds
        ) {
            test_hamming_helper(&seq, HammingIterOpt::SkipNBase, b'N');
            test_hamming_helper(&seq, HammingIterOpt::MutateNBase, b'N');
        }

        #[test]
        fn prop_test_one_deletion_iter(
            seq in "[ACGTN]{0, 23}", // 0 and 23 are inclusive bounds
        ) {
            prop_assert_eq!(SSeq::from_bytes(seq.as_bytes()).one_deletion_iter().count(), seq.len());
            prop_assert!(SSeq::from_bytes(seq.as_bytes()).one_deletion_iter().all(|s| s.len() == seq.len()-1));
        }

        #[test]
        fn prop_test_insertion(
            seq in "[ACGTN]{0, 22}", // 0 and 22 are inclusive bounds
            pos in 0usize..22usize,
            base_index in 0usize..5usize,
        ) {
            let pos = if seq.is_empty() { 0 } else { pos % seq.len() };
            let mut sseq = SSeq::from_bytes(seq.as_bytes());
            let base = UPPER_ACGTN[base_index];
            sseq.insert_unchecked(pos, base);
            prop_assert_eq!(sseq[pos], base);
            prop_assert_eq!(sseq.len(), seq.len() + 1);

        }

        #[test]
        fn prop_test_one_insertion_iter(
            seq in "[ACGTN]{0, 22}", // 0 and 22 are inclusive bounds
        ) {
            prop_assert_eq!(
                SSeq::from_bytes(seq.as_bytes())
                    .one_insertion_iter(InsertionIterOpt::ExcludeNBase)
                    .count(),
                4 * (seq.len() + 1)
            );
            prop_assert_eq!(
                SSeq::from_bytes(seq.as_bytes())
                    .one_insertion_iter(InsertionIterOpt::IncludeNBase)
                    .count(),
                5 * (seq.len() + 1)
            );
        }
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_1() {
        let _ = SSeq::from_bytes(b"ASDF");
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_2() {
        let _ = SSeq::from_bytes(b"ag");
    }

    #[test]
    #[should_panic]
    fn test_sseq_too_long() {
        let _ = SSeq::from_bytes(b"GGGACCGTCGGTAAAGCTACAGTGAGGGATGTAGTGATGC");
    }

    #[test]
    fn test_as_bytes() {
        assert_eq!(SSeq::from_bytes(b"ACGT").as_bytes(), b"ACGT");
    }

    #[test]
    fn test_has_n() {
        assert!(SSeq::from_bytes(b"ACGTN").has_n());
        assert!(!SSeq::from_bytes(b"ACGT").has_n());
    }

    #[test]
    fn test_is_homopolymer() {
        assert!(SSeq::from_bytes(b"AAAA").is_homopolymer());
        assert!(!SSeq::from_bytes(b"ACGT").is_homopolymer());
    }

    #[test]
    fn test_has_homopolymer_suffix() {
        assert!(SSeq::from_bytes(b"ACGTAAAAA").has_homopolymer_suffix(b'A', 5));
        assert!(!SSeq::from_bytes(b"ACGTTAAAA").has_homopolymer_suffix(b'A', 5));
        assert!(SSeq::from_bytes(b"CCCCC").has_homopolymer_suffix(b'C', 5));
        assert!(!SSeq::from_bytes(b"GGGG").has_homopolymer_suffix(b'G', 5));
    }

    #[test]
    fn test_has_polyt_suffix() {
        assert!(SSeq::from_bytes(b"CGCGTTTTT").has_polyt_suffix(5));
        assert!(!SSeq::from_bytes(b"CGCGAAAAA").has_polyt_suffix(5));
    }

    #[test]
    fn test_one_hamming_simple() {
        assert_equal(
            SSeq::from_bytes(b"GAT")
                .one_hamming_iter(HammingIterOpt::SkipNBase)
                .collect_vec(),
            vec![
                b"AAT", b"CAT", b"TAT", b"GCT", b"GGT", b"GTT", b"GAA", b"GAC", b"GAG",
            ]
            .into_iter()
            .map(|x| SSeq::from_bytes(x)),
        );

        assert_equal(
            SSeq::from_bytes(b"GNT")
                .one_hamming_iter(HammingIterOpt::SkipNBase)
                .collect_vec(),
            vec![b"ANT", b"CNT", b"TNT", b"GNA", b"GNC", b"GNG"]
                .into_iter()
                .map(|x| SSeq::from_bytes(x)),
        );
    }

    #[test]
    fn test_remove() {
        let mut seq = SSeq::from_bytes(b"GCAT");
        assert_eq!(seq.remove(2), b'A');
        assert_eq!(seq.cmp(&SSeq::from_bytes(b"GCT")), Ordering::Equal);
        assert_eq!(seq, SSeq::from_bytes(b"GCT"));
    }

    #[test]
    fn test_insert() {
        let mut seq = SSeq::from_bytes(b"GCAT");
        seq.insert_unchecked(2, b'G');
        assert_eq!(seq, SSeq::from_bytes(b"GCGAT"));
    }

    #[test]
    #[should_panic]
    fn test_insert_panic_out_of_capacity() {
        let mut seq = SSeq::from_bytes(b"AGCTTTGCGTTAGGCAGGTTTAC");
        seq.insert_unchecked(2, b'G');
    }

    #[test]
    #[should_panic]
    fn test_insert_panic_bad_index() {
        let mut seq = SSeq::from_bytes(b"AG");
        seq.insert_unchecked(10, b'G');
    }

    #[test]
    fn test_one_deletion() {
        assert_equal(
            SSeq::from_bytes(b"GAT").one_deletion_iter().collect_vec(),
            vec![b"AT", b"GT", b"GA"]
                .into_iter()
                .map(|x| SSeq::from_bytes(x)),
        );
    }

    #[test]
    fn test_one_insertion() {
        assert_equal(
            SSeq::from_bytes(b"GA")
                .one_insertion_iter(InsertionIterOpt::ExcludeNBase)
                .collect_vec(),
            vec![
                b"AGA", b"CGA", b"GGA", b"TGA", b"GAA", b"GCA", b"GGA", b"GTA", b"GAA", b"GAC",
                b"GAG", b"GAT",
            ]
            .into_iter()
            .map(|x| SSeq::from_bytes(x)),
        );
    }

    #[test]
    fn test_from_iter() {
        let seq = SSeq::from_bytes(b"ACGT");
        let _ = SSeq::from_iter(seq.as_bytes());
        let seq_vec = seq.as_bytes().to_vec();
        let _ = SSeq::from_iter(seq_vec);
    }

    #[test]
    fn test_hamming_custom() {
        let seq = SSeq::from_bytes(b"AN");
        assert_equal(
            seq.one_hamming_iter(HammingIterOpt::Custom {
                skip_chars: SSeqChars::n(),
                mutation_chars: SSeqChars::actgn(),
            }),
            ["CN", "GN", "TN", "NN"]
                .into_iter()
                .map(|x| SSeq::from_bytes(x.as_bytes())),
        );

        assert_equal(
            seq.one_hamming_iter(HammingIterOpt::Custom {
                skip_chars: SSeqChars::new(b"A"),
                mutation_chars: SSeqChars::actg(),
            }),
            ["AA", "AC", "AG", "AT"]
                .into_iter()
                .map(|x| SSeq::from_bytes(x.as_bytes())),
        );
    }

    #[test]
    fn test_revcomp() {
        assert_eq!(
            SSeq::from_bytes(b"ACGTN").revcomp(),
            SSeq::from_bytes(b"NACGT")
        );
    }
}
