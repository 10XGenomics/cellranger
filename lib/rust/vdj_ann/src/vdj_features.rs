// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

// This file contains functions that take as input an amino acid reference sequence for a V
// segment, along with its chain type (IGH, IGK, IGL, TRA or TRB), and attempt to find features in
// the sequence.
//
// These functions can fail, in particular for reference sequences that are damaged and possibly
// for other sequences.
//
// The same functions can also be applied to a protein sequence, however the sequence needs to be
// modified to add a fake leader sequence on the left (we use MXXXXXXXXXXXXXXXXXXXX), and to
// truncate on the right to trim a bit beyond the start of the CDR3.

use vector_utils::reverse_sort;

type Score = usize;

/// A scored set of sequence elements to look for at a particular offset.
struct ScoredElements {
    /// The offest into the sequence to search for these elements.
    offset: usize,
    /// The elements to match; each byte in the string is one candidate.
    elements: &'static str,
    /// The score for matching an element at offset.
    score: Score,
}

impl ScoredElements {
    /// Return the score if this matches the sequence.
    pub fn score(&self, seq: &[u8]) -> Option<Score> {
        let seq_element = seq[self.offset];
        for search_element in self.elements.bytes() {
            if seq_element == search_element {
                return Some(self.score);
            }
        }
        None
    }
}

struct Motif {
    elements: Vec<ScoredElements>,
    /// The length of this motif.
    len: usize,
}

impl Motif {
    pub fn from_items(elements: impl IntoIterator<Item = (usize, &'static str, Score)>) -> Self {
        let elements: Vec<_> = elements
            .into_iter()
            .map(|(offset, elements, score)| ScoredElements {
                offset,
                elements,
                score,
            })
            .collect();
        Self {
            len: elements.iter().map(|motif| motif.offset).max().unwrap() + 1,
            elements,
        }
    }

    /// Score this motif against the provided sequence.
    /// If the sequence is not long enough to be matched against the full motif,
    /// return None.
    fn score(&self, seq: &[u8]) -> Option<Score> {
        if seq.len() < self.len {
            return None;
        }
        let score: Score = self
            .elements
            .iter()
            .filter_map(|element| element.score(seq))
            .sum();
        if score == 0 { None } else { Some(score) }
    }

    /// Search for the best match in a sliding window through a sequence.
    /// The search starts matching the motif at the start element, and tries to
    /// match against all motif start positions in the range
    /// [start, start + width). If the sequence is not long
    /// enough to match the full motif against the sequence starting from
    /// the first start position, return TooShort.
    pub fn best_match(&self, seq: &[u8], start: usize, width: usize) -> MatchResult {
        if seq.len() < start + self.len {
            return Err(NoMatch::TooShort);
        }
        let end = start + width;

        match (start..usize::min(end, seq.len()))
            .filter_map(|motif_start| {
                self.score(&seq[motif_start..])
                    .map(|score| (score, motif_start))
            })
            .max_by_key(|(score, _)| *score)
        {
            Some((_, motif_start)) => Ok(motif_start),
            None => Err(NoMatch::NotFound),
        }
    }
}

/// Represent the result of matching a motif against a sequence window.
type MatchResult = Result<usize, NoMatch>;

/// Ways we can fail to match a motif against a sequence window.
#[derive(Debug, PartialEq, Eq)]
enum NoMatch {
    TooShort,
    NotFound,
}

trait WithOffset {
    fn with_offset(self, offset: impl FnOnce(usize) -> i64) -> Self;
}

impl WithOffset for MatchResult {
    /// Call the provided offset closure with the found start position.
    /// Add the resulting offset to the start position.
    /// Internal math is performed as signed integers, so the offset can
    /// safely be negative.
    fn with_offset(self, offset: impl FnOnce(usize) -> i64) -> Self {
        self.map(|p| {
            let off = offset(p);
            let with_offset = p as i64 + off;
            assert!(
                with_offset >= 0,
                "applying offset {off} resulting in a negative sequence index {with_offset}",
            );

            with_offset as usize
        })
    }
}

/// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
/// framework one region, which is immmediately after the signal peptide.
///
/// Scoring motif:
///
/// pos   weight  values
/// 1     150     Q, E, D, G, K
/// 2      50     V, I, Q, A
/// 4     100     L, V, M
/// 6     250     Q, E
/// 22    250     C
/// 23    250     C
///
/// If the starting amino acid is C, we add one to the start position.
pub fn fr1_start(aa: &[u8], chain_type: &str) -> Option<usize> {
    let motif = Motif::from_items([
        (0, "QEDGK", 150),
        (1, "VIQA", 50),
        (3, "LVM", 100),
        (5, "QE", 250),
        (21, "C", if chain_type == "IGH" { 500 } else { 250 }),
        (22, "C", 250),
    ]);

    motif
        .best_match(aa, 0, if chain_type == "IGL" { 26 } else { 41 })
        .with_offset(|p| if aa[p] == b'C' { 1 } else { 0 })
        .ok()
}

/// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
/// CDR1 region.  Note that there is more than one convention regarding the CDR1 start, and these
/// conventions appear to differ by fixed offsets.  The convention used here is for IMGT.
/// Chain type is one of IGH, IGK, IGL, TRA or TRB.
pub fn cdr1_start(aa: &[u8], chain_type: &str, _verbose: bool) -> Option<usize> {
    let motif = Motif::from_items([
        (0, "V", 50),
        (1, "T", 30),
        (2, "LIVM", 200),
        (3, "STR", 80),
        (4, "C", 250),
        (7, "SID", 100),
    ]);

    let (start, end) = {
        let start = 25;
        let fr1 = fr1_start(aa, chain_type)?;
        let end = fr1 + 20;
        if start <= end {
            (start, end)
        } else {
            (end, start)
        }
    };
    let width = end - start;

    motif
        .best_match(aa, start, width)
        .with_offset(|_| {
            if chain_type.starts_with("TR") || chain_type == "IGH" {
                8
            } else {
                5
            }
        })
        .ok()
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// FR2 region.  Chain type is one of IGH, IGK, IGL, TRA or TRB'

pub fn fr2_start(aa: &[u8], chain_type: &str, _verbose: bool) -> Option<usize> {
    let motif = Motif::from_items([
        (0, "LMVF", 50),
        (1, "Y", if chain_type == "IGL" { 250 } else { 0 }),
        (2, "W", 250),
        (3, "Y", 150),
        (4, "R", 100),
        (5, "Q", 250),
        (8, "G", 110),
        (9, "KQ", 60),
        (10, "GKA", 40),
    ]);

    motif
        .best_match(aa, 40, if chain_type == "IGH" { 23 } else { 34 })
        .with_offset(|_| {
            let mut offset = -1;
            if chain_type != "IGH" {
                offset += 1;
            }
            if chain_type == "IGK" || chain_type == "IGL" {
                offset += 2;
            }
            offset
        })
        .ok()
}

/// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
/// CDR2 region.  Chain type is one of IGH, IGK, IGL, TRA or TRB.
pub fn cdr2_start(aa: &[u8], chain_type: &str, _verbose: bool) -> Option<usize> {
    let fwr2_start = fr2_start(aa, chain_type, false)?;

    if chain_type == "IGH" {
        // Six amino acids preceeding the CDR2 start.
        let motif = Motif::from_items([
            (0, "L", 80),
            (1, "E", 80),
            (2, "W", 80),
            (3, "VMIL", 40),
            (4, "GSA", 40),
            // NOTE: this originally had a trailing sixth empty entry, not sure if it was doing anything
            // besides effectively asserting that the sequence needed to have at
            // least six elements.
        ]);

        motif
            .best_match(aa, fwr2_start + 8, 6)
            .with_offset(|_| 7)
            .ok()
    } else if chain_type == "TRA" {
        let motif = Motif::from_items([
            (0, "PL", 15),
            (1, "QVETI", 15),
            (2, "LF", 20),
            (3, "L", 35),
            (4, "LI", 15),
            // NOTE: this originally had a trailing sixth empty entry, not sure if it was doing anything
            // besides effectively asserting that the sequence needed to have at
            // least six elements.
        ]);

        motif
            .best_match(aa, fwr2_start + 10, 3)
            .with_offset(|_| 6)
            .ok()
    } else {
        Some(
            fwr2_start
                + (if chain_type == "IGK" || chain_type == "IGL" {
                    15
                } else {
                    17
                }),
        )
    }
}

pub fn cdr3_start(aa: &[u8], _chain_type: &str, _verbose: bool) -> usize {
    let motif = [b"LQPEDSAVYYC", b"VEASQTGTYFC", b"ATSGQASLYLC"];
    let nm = motif[0].len();
    let reach = 18;
    let mut scores = Vec::<(usize, usize)>::new();
    for j in aa.len().saturating_sub(nm + reach)..=aa.len().saturating_sub(nm) {
        let mut score = 0;
        for k in 0..nm {
            for m in motif {
                if aa[j + k] == m[k] {
                    score += 1;
                    if aa[j + k] == b'Q' {
                        break;
                    }
                }
            }
        }
        scores.push((score, j + nm));
    }
    reverse_sort(&mut scores);
    scores[0].1
}

pub fn cdr3_score(aa: &[u8], _chain_type: &str, _verbose: bool) -> usize {
    let motif = [b"LQPEDSAVYYC", b"VEASQTGTYFC", b"ATSGQASLYLC"];
    let nm = motif[0].len();
    const REACH: usize = 18;
    let mut scores = Vec::<(usize, usize)>::with_capacity(REACH + 1);
    for j in aa.len().saturating_sub(nm + REACH)..=aa.len().saturating_sub(nm) {
        let mut score = 0;
        for k in 0..nm {
            for m in motif {
                if aa[j + k] == m[k] {
                    score += 1;
                    if aa[j + k] == b'Q' {
                        break;
                    }
                }
            }
        }
        scores.push((score, j + nm));
    }
    reverse_sort(&mut scores);
    scores[0].0
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given the amino acid sequence for a V reference sequence, attempt to find the start of the
// FR3 region.

pub fn fr3_start(aa: &[u8], chain_type: &str, verbose: bool) -> Option<usize> {
    let cdr3_start = cdr3_start(aa, chain_type, verbose);

    if chain_type == "IGK" || chain_type == "IGL" {
        let motif = Motif::from_items([
            (0, "G", 100),
            (2, "P", 100),
            (4, "R", 100),
            (5, "F", 100),
            (7, "G", 100),
        ]);

        // Score positions.

        if cdr3_start < 35 {
            return None;
        }

        let start = cdr3_start - 35;
        let width = 8;

        motif.best_match(aa, start, width).ok()
    } else if chain_type == "IGH" {
        let motif = Motif::from_items([
            (0, "YN", 600),
            (1, "Y", 500),
            (2, "AN", 400),
            (5, "FL", 850),
            (6, "KQR", 800),
            (8, "RK", 1000),
            (9, "FVAL", 700),
        ]);

        // Score positions.
        if cdr3_start < 40 {
            return None;
        }

        motif
            .best_match(aa, cdr3_start - 40, 11)
            .with_offset(|_| -1)
            .ok()
    } else if chain_type == "TRA" {
        let motif = Motif::from_items([
            (0, "EVNK", 50),
            (1, "TKAE", 50),
            (2, "ES", 50),
            (3, "NDS", 50),
            (4, "GN", 50),
            (5, "RGM", 80),
            (6, "FYAI", 50),
            (7, "ST", 50),
            (8, "AV", 50),
            (9, "TE", 50),
            (11, "ND", 50),
        ]);

        // Score positions.
        if cdr3_start < 36 {
            return None;
        };

        motif
            .best_match(aa, cdr3_start - 36, 4)
            .with_offset(|_| 1)
            .ok()
    } else {
        // Do TRB.
        let motif = Motif::from_items([
            (2, "KED", 50),
            (3, "GSQ", 200),
            (4, "DEGS", 200),
            (5, "IVLM", 200),
            (6, "PS", 100),
        ]);

        // Score positions.
        if cdr3_start < 38 {
            return None;
        }

        motif.best_match(aa, cdr3_start - 38, 4).ok()
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn score_fwr3(aa: &[u8], r: usize, freqs: &[Vec<Vec<(u32, u8)>>]) -> f64 {
    let chain_type;
    if r == 0 {
        chain_type = "IGH";
    } else if r == 1 {
        chain_type = "IGK";
    } else if r == 2 {
        chain_type = "IGL";
    } else if r == 3 {
        chain_type = "TRA";
    } else {
        chain_type = "TRB";
    }
    let cdr3 = cdr3_start(aa, chain_type, false);
    let motif = freqs[0].len();
    let mut score = 0.0;
    for j in 0..motif {
        let x = aa[cdr3 - j - 1];
        let mut m = 0;
        let mut total = 0;
        for k in 0..freqs[r][j].len() {
            let count = freqs[r][j][k].0;
            let y = freqs[r][j][k].1;
            total += count;
            if y == x {
                m += count;
            }
        }
        score += m as f64 / total as f64;
    }
    score
}

pub fn score4(aa: &[u8], r: usize) -> usize {
    let chain_type;
    if r == 0 {
        chain_type = "IGH";
    } else if r == 1 {
        chain_type = "IGK";
    } else if r == 2 {
        chain_type = "IGL";
    } else if r == 3 {
        chain_type = "TRA";
    } else {
        chain_type = "TRB";
    }
    let cdr3 = cdr3_start(aa, chain_type, false);
    let n = aa.len();
    assert!(n >= 22);
    let mut score = 0;
    let x = aa[cdr3 - 4];
    if x == b'V' || x == b'T' || x == b'L' {
        score += 1;
    }
    let x = aa[cdr3 - 3];
    if x == b'Y' {
        score += 1;
    }
    let x = aa[cdr3 - 2];
    if x == b'Y' || x == b'F' || x == b'L' {
        score += 1;
    }
    let x = aa[cdr3 - 1];
    if x == b'C' {
        score += 3;
    }
    score
}

#[cfg(test)]
mod test {
    use super::{Motif, NoMatch, ScoredElements};

    #[test]
    fn test_scored_element() {
        let e = ScoredElements {
            offset: 1,
            elements: "ABC",
            score: 1,
        };
        assert_eq!(None, e.score("XX".as_bytes()));
        assert_eq!(Some(1), e.score("XA".as_bytes()));
        assert_eq!(Some(1), e.score("XC".as_bytes()));
        assert_eq!(None, e.score("XXA".as_bytes()));
    }

    #[test]
    fn test_motif_match() {
        use NoMatch::*;
        let motif = Motif::from_items([(0, "ABC", 10), (2, "D", 3)]);
        assert_eq!(TooShort, motif.best_match("".as_bytes(), 0, 1).unwrap_err());
        assert_eq!(
            TooShort,
            motif.best_match("A".as_bytes(), 0, 1).unwrap_err()
        );
        assert_eq!(
            TooShort,
            motif.best_match("AA".as_bytes(), 0, 1).unwrap_err()
        );
        assert_eq!(
            TooShort,
            motif.best_match("AAA".as_bytes(), 1, 1).unwrap_err()
        );
        assert_eq!(Ok(0), motif.best_match("AAA".as_bytes(), 0, 5));

        assert_eq!(Ok(0), motif.best_match("AXX".as_bytes(), 0, 1));
        assert_eq!(Ok(0), motif.best_match("AAD".as_bytes(), 0, 1));
        assert_eq!(Ok(1), motif.best_match("AADD".as_bytes(), 0, 2));

        assert_eq!(
            NotFound,
            motif.best_match("XXX".as_bytes(), 0, 1).unwrap_err()
        );
        // potential change to legacy behavior/implicit in current impl:
        // later matches end up being preferred for same score.
        assert_eq!(Ok(1), motif.best_match("AAAA".as_bytes(), 0, 2));
    }
}
