//!
//! Module that deals with translating a UMI sequence
//!
//! Typically in all our assays, the UMI sequence is a random k-mer with a
//! length of 10 or 12 bases.
//!

use arrayvec::ArrayVec;
use fastq_set::sseq::SSeqGen;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

type SplintSeq = SSeqGen<NUM_SPLINT_BASES>;
const BASE_OPTS: [u8; 4] = [b'A', b'C', b'G', b'T'];

const NUM_SPLINT_BASES: usize = 4;

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug)]
pub struct SplintToUmiTranslator {
    sequences: ArrayVec<[SplintSeq; 4]>,
    max_allowed_hamming_distance: u32,
}

impl SplintToUmiTranslator {
    pub fn new<'a>(
        seq_iter: impl IntoIterator<Item = &'a [u8]>,
        max_allowed_hamming_distance: u32,
    ) -> Self {
        let mut sequences = ArrayVec::new();
        for seq in seq_iter.into_iter().map(SplintSeq::from_bytes) {
            sequences
                .try_push(seq)
                .expect("Cannot handle more than 4 sequences in SplintToUmiTranslator");
        }
        assert!(
            !sequences.is_empty(),
            "Expected at least 1 sequence in the SplintToUmiTranslator"
        );
        assert!(
            sequences.iter().unique().count() == sequences.len(),
            "Duplicates found in the input sequences to SplintToUmiTranslator"
        );
        SplintToUmiTranslator {
            sequences,
            max_allowed_hamming_distance,
        }
    }

    pub fn translate(&self, query: &[u8]) -> u8 {
        let hamming_distances: ArrayVec<[u32; 4]> = self
            .sequences
            .iter()
            .map(|seq| triple_accel::hamming(seq.as_bytes(), query))
            .collect();
        let min_hd = *hamming_distances.iter().min().unwrap();
        if min_hd <= self.max_allowed_hamming_distance
            && hamming_distances.iter().filter(|&&hd| hd == min_hd).count() == 1
        // Unique min HD
        {
            BASE_OPTS[hamming_distances
                .iter()
                .find_position(|&&hd| hd == min_hd)
                .unwrap()
                .0]
        } else {
            b'N'
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_translate() {
        let translator = SplintToUmiTranslator::new(
            ["AATA", "CGCC", "GCGG", "TTAT"]
                .iter()
                .map(|x| x.as_bytes()),
            0,
        );

        // EXACT MATCHES
        assert_eq!(translator.translate(b"AATA"), b'A');
        assert_eq!(translator.translate(b"CGCC"), b'C');
        assert_eq!(translator.translate(b"GCGG"), b'G');
        assert_eq!(translator.translate(b"TTAT"), b'T');

        // 1-HD matches
        assert_eq!(translator.translate(b"ACTA"), b'N');

        let translator = SplintToUmiTranslator::new(
            ["AATA", "CGCC", "GCGG", "TTAT"]
                .iter()
                .map(|x| x.as_bytes()),
            1,
        );
        assert_eq!(translator.translate(b"ACTA"), b'A');
        assert_eq!(translator.translate(b"ANTA"), b'A');
        assert_eq!(translator.translate(b"AGCC"), b'C');
        assert_eq!(translator.translate(b"GCTG"), b'G');
        assert_eq!(translator.translate(b"TTAA"), b'T');
        assert_eq!(translator.translate(b"TAAA"), b'N');
    }

    #[should_panic]
    #[test]
    fn test_zero_seqs() {
        let _ = SplintToUmiTranslator::new([], 0);
    }

    #[should_panic]
    #[test]
    fn test_more_than_four_seqs() {
        let _ = SplintToUmiTranslator::new(
            ["AATA", "CGCC", "GCGG", "TTAT", "AACG"]
                .iter()
                .map(|x| x.as_bytes()),
            0,
        );
    }

    #[should_panic]
    #[test]
    fn test_duplicate_seqs() {
        let _ = SplintToUmiTranslator::new(
            ["AATA", "AATA", "TTAT", "AACG"]
                .iter()
                .map(|x| x.as_bytes()),
            0,
        );
    }
}
