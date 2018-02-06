//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use debruijn::dna_string::DnaString;
use debruijn::Mer;

use std::cmp;
use utils;
use graph_read::Read;
use itertools::Itertools;
use constants::UmiType;
use bio;
use std::borrow::Borrow;

use bio::alignment::{AlignmentOperation, Alignment, AlignmentMode};
use bio::alignment::sparse::{HashMapFx, sdpkpp};
use bio::alignment::pairwise::banded;
pub type Scoring = bio::alignment::pairwise::Scoring<bio::alignment::pairwise::MatchParams>;

const MAX_OUT_QUAL : u8 = 60;
const MAX_READ_QUAL : u8 = 30;
const MIN_LOG_PROB : f64 = -100000.0;

lazy_static! {
    static ref MATCH_SCORE_VS_QUAL: [f64; (MAX_READ_QUAL + 1) as usize] = {
        let mut table: [f64; (MAX_READ_QUAL + 1) as usize] = [0.0f64; (MAX_READ_QUAL + 1) as usize];
        for q in 0..(MAX_READ_QUAL + 1) as usize {
            table[q] = (1.0 - 10.0_f64.powf(-(q as f64) / 10.0)).log10();
        }
        table
    };
    static ref MISMATCH_SCORE_VS_QUAL: [f64; (MAX_READ_QUAL + 1) as usize] = {
        let mut table: [f64; (MAX_READ_QUAL + 1) as usize] = [0.0f64; (MAX_READ_QUAL + 1) as usize];
        for q in 0..(MAX_READ_QUAL + 1) as usize {
            table[q] = -(q as f64) / 10.0 - 3.0_f64.log10();;
        }
        table
    };
}

/// Alignment against a reference
/// This struct pack the index of the reference and the alignment together.
#[derive(Debug, Clone)]
pub struct AlignmentPacket {
    pub ref_idx: usize, // index of reference
    pub alignment: Alignment,
}

impl AlignmentPacket {
    pub fn new(ref_idx:usize, alignment: Alignment) -> Self {
        AlignmentPacket {
            ref_idx: ref_idx,
            alignment: alignment,
        }
    }

    /// Compute the edit distance implied by the alignment.
    /// The clipped regions are not considered in the computation
    pub fn edit_distance(&self) -> usize {
        let mut d = 0;
        for &op in self.alignment.operations.iter() {
            d += match op {
                AlignmentOperation::Match => 0,
                AlignmentOperation::Xclip(_) => 0,
                AlignmentOperation::Yclip(_) => 0,
                _ => 1
            };
        }
        d
    }
}

pub struct AlignHelper<'a> {
    pub refs: Vec<DnaString>,
    pub kmers_hash: HashMapFx<&'a [u8], Vec<(usize, u32)>>,
    scoring: Scoring,
    pub k: usize,
    pub w: usize,
}

impl<'a> AlignHelper<'a> {
    /// Creates a new AlignHelper instance given the set of reference sequences,
    /// the scoring parameters, and the banding parameters (kmer length and window size)
    ///
    /// # Arguments
    ///
    /// * `refs` - vector of reference sequences (stored as String)
    /// * `scoring` - Scoring struct
    /// * `k` - kmer length for constructing the band (see bio::alignment::pairwise::banded)
    /// * `w` - window size for constructing the band (see bio::alignment::pairwise::banded)
    ///
    pub fn new(refs: &'a Vec<String>, scoring: Scoring, k: usize, w: usize) -> Self {

        // Make sure that the scoring implies that the alignment is local in the reference y
        assert!(scoring.yclip_prefix == 0);
        assert!(scoring.yclip_suffix == 0);

        let mut kmers_hash: HashMapFx<&'a [u8], Vec<(usize, u32)>> = HashMapFx::default();
        let mut dna_strings = Vec::new();

        for (ref_idx, ref_seq) in refs.iter().enumerate() {
            let dna_string = DnaString::from_dna_string(ref_seq);
            dna_strings.push(dna_string);
            let seq = ref_seq.as_bytes();
            for i in 0..(seq.len() + 1).saturating_sub(k) {
                kmers_hash.entry(&seq[i..i + k])
                    .or_insert_with(|| Vec::new())
                    .push((ref_idx as usize, i as u32));
            }
        }

        AlignHelper {
            refs: dna_strings,
            kmers_hash: kmers_hash,
            scoring: scoring,
            k: k,
            w: w
        }
    }

    pub fn set_scoring(&mut self, scoring: Scoring) {
        // Make sure that the scoring implies that the alignment is local in the reference y
        assert!(scoring.yclip_prefix == 0);
        assert!(scoring.yclip_suffix == 0);
        self.scoring = scoring;
    }

    /// Align the read with at least a minimum SW score. The read is aligned with the
    /// sequences in refs and the first acceptable alignment is returned. We use the 
    /// banded aligner from rust-bio (see bio::alignment::pairwise::banded). The function
    /// returns Some(AlignmentPacket) if we could find a good enough alignment, otherwise None.
    /// First we find all the kmer matches between the read and all the refs using the
    /// precomputed 'kmers_hash'. We sort the matches based on the number of kmer
    /// matches per reference and then call the aligner.
    ///
    /// # Arguments
    ///
    /// * `read` - Read object. 
    /// * `min_align_score` - Minimum SW score for accepting an alignment
    ///
    pub fn find_read_matches(&self, read: &Read, min_align_score: i32) -> Option<AlignmentPacket> {
        
        let mut all_kmer_matches: Vec<(i32, usize, usize, usize)> = Vec::new();
        let read_seq = read.seq.to_dna_string();
        let read_seq = read_seq.as_bytes();
        let mut ref_kmer_counts = vec![0i32; self.refs.len()];

        for i in 0..(read_seq.len() + 1).saturating_sub(self.k) {
            let slc = &read_seq[i..i + self.k];
            match self.kmers_hash.get(slc) {
                Some(matches) => {
                    for &ref_pos in matches {
                        all_kmer_matches.push((0, ref_pos.0, i, ref_pos.1 as usize));
                        ref_kmer_counts[ref_pos.0] += 1;
                    }
                },
                None => {},
            }
        }
        for i in 0..all_kmer_matches.len() {
            all_kmer_matches[i].0 = -ref_kmer_counts[all_kmer_matches[i].1];
        }
        all_kmer_matches.sort();
        
        let mut aligner = banded::Aligner::with_scoring(self.scoring.clone(), self.k, self.w);
        for (ref_idx, match_group) in &all_kmer_matches.into_iter().group_by(|x| x.1) {
            let matches = match_group.into_iter().map(|x| (x.2 as u32, x.3 as u32)).collect();
            let ref_seq = &self.refs[ref_idx].to_dna_string();
            let alignment = aligner.custom_with_expanded_matches(read_seq, ref_seq.as_bytes(), matches, Some(1), true);
            if alignment.score > min_align_score {
                let align_pack = AlignmentPacket::new(ref_idx, alignment);
                assert!(read.validate_alignment(&align_pack));
                return Some(align_pack);
            }
        };

        None
    }

    pub fn find_read_matches_sparse(&self, read: &Read, min_align_score: i32) -> Option<AlignmentPacket> {

        assert!(min_align_score >= 0);
        
        let mut all_kmer_matches: Vec<(i32, usize, usize, usize)> = Vec::new();
        let read_seq = read.seq.to_dna_string();
        let read_seq = read_seq.as_bytes();
        let mut ref_kmer_counts = vec![0i32; self.refs.len()];
        for i in 0..(read_seq.len() + 1).saturating_sub(self.k) {
            let slc = &read_seq[i..i + self.k];
            match self.kmers_hash.get(slc) {
                Some(matches) => {
                    for &ref_pos in matches {
                        all_kmer_matches.push((0, ref_pos.0, i, ref_pos.1 as usize));
                        ref_kmer_counts[ref_pos.0] += 1;
                    }
                },
                None => {},
            }
        }
        for i in 0..all_kmer_matches.len() {
            all_kmer_matches[i].0 = -ref_kmer_counts[all_kmer_matches[i].1];
        }
        all_kmer_matches.sort();
        
        for (ref_idx, match_group) in &all_kmer_matches.into_iter().group_by(|x| x.1) {
            let matches = match_group.into_iter().map(|x| (x.2 as u32, x.3 as u32)).collect();
            let sparse_al = sdpkpp(&matches, self.k, self.scoring.match_fn.match_score as u32, 
                self.scoring.gap_open, self.scoring.gap_extend);
            if sparse_al.score > min_align_score as u32 {
                let match_path: Vec<(u32,u32)> = sparse_al.path.iter().map(|i| matches[*i]).collect();
                let alignment = self.alignment_from_sparse_alignment(match_path, sparse_al.score as i32, read_seq.len(), self.refs[ref_idx].len());
                let align_pack = AlignmentPacket::new(ref_idx, alignment);
                assert!(read.validate_alignment(&align_pack));
                return Some(align_pack);
            }
        };
        None
    }

    fn alignment_from_sparse_alignment(&self, match_path: Vec<(u32,u32)>, score: i32,
                                     read_len: usize, ref_len: usize) -> Alignment {

        assert!(!match_path.is_empty());
        let mut steps = Vec::new();
        let first = match_path[0];
        if first.0 !=0 { steps.push(AlignmentOperation::Xclip(first.0 as usize)); }
        if first.1 !=0 { steps.push(AlignmentOperation::Yclip(first.1 as usize)); }
        let xstart = first.0 as usize;
        let ystart = first.1 as usize;
        // Each step suggests a kmer match.
        for _ in 0..self.k {
            steps.push(AlignmentOperation::Match);
        }
        for (idx, s) in match_path.iter().skip(1).enumerate() {
            // If previous step was (i, j) and this step is (i+1, j+1),
            // then we just advanced the kmer match by one position.
            // Otherwise, we skipped some bases in either one or both of
            // the sequences and then we added a new kmer match.
            let prev_step = match_path[idx];
            let skip_ref = s.1 - prev_step.1;
            let skip_read = s.0 - prev_step.0;

            if skip_ref > 1 || skip_read > 1 {
                for _ in 0..skip_read.saturating_sub(self.k as u32) {
                    steps.push(AlignmentOperation::Ins);
                }
                for _ in 0..skip_ref.saturating_sub(self.k as u32) {
                    steps.push(AlignmentOperation::Del);
                }
                // After the indel, there's a new kmer match.
                for _ in 0..self.k {
                    steps.push(AlignmentOperation::Match);
                }
            } else {
                steps.push(AlignmentOperation::Match);
            }
        }

        let last = match_path.last().unwrap();
        let xend = last.0 as usize + self.k;
        let yend = last.1 as usize + self.k;
        if read_len > xend { steps.push(AlignmentOperation::Xclip(read_len - xend)); }
        if ref_len > xend { steps.push(AlignmentOperation::Yclip(ref_len - yend)); }

        Alignment {
            score: score,
            ystart: ystart,
            xstart: xstart,
            yend: yend,
            xend: xend,
            ylen: ref_len,
            xlen: read_len,
            operations: steps,
            mode: AlignmentMode::Custom,
        }
    }

    /// Given a reference (specified by ref_idx) and a vector of (AlignmentPacket, Read) pairs,
    /// return a vector of length equal to that of the reference sequence. Each entry on the
    /// returned vector is a set of triplets (umi, base, quality) corresponding to the
    /// alignment at that position in the reference sorted by the umi
    pub fn pileup<T: Borrow<Read>>(&self, ref_idx: usize, alignments: &Vec<(AlignmentPacket, T)>) -> Vec<Vec<(u32, u8, u8)>> {
        let ref_len = self.refs[ref_idx].len();
        let mut pos_pileup = Vec::with_capacity(ref_len);
        for _ in 0..ref_len {
            pos_pileup.push(Vec::new());
        }

        for &(ref al, ref read_borrow) in alignments.iter() {
            let read = read_borrow.borrow();
            let ref alignment = al.alignment;
            let mut ref_idx = alignment.ystart;
            let mut query_idx = alignment.xstart;

            for &op in alignment.operations.iter() {
                match op {
                    AlignmentOperation::Ins => {
                        query_idx += 1
                    },
                    AlignmentOperation::Del => {
                        ref_idx += 1
                    },
                    AlignmentOperation::Match | AlignmentOperation::Subst => {
                        pos_pileup[ref_idx].push((read.umi,
                                                  read.seq.get(query_idx),
                                                  read.quals[query_idx]));
                        ref_idx += 1;
                        query_idx += 1;
                    },
                    _ => {}
                }
            }
        }

        for ref_idx in 0..ref_len {
            pos_pileup[ref_idx].sort_by_key(|x| x.0);
        }
        pos_pileup
    }


    pub fn base_quals(&self, ref_idx: usize, pileup: &Vec<Vec<(UmiType, u8, u8)>>, rt_err: f64) -> Vec<u8> {

        assert!(self.refs[ref_idx].len() == pileup.len());

        let mut quals : Vec<u8> = Vec::with_capacity(pileup.len());

        for (pos, base_pileup) in pileup.iter().enumerate() {

            if base_pileup.is_empty() {
                quals.push(0);
                continue;
            }

            let contig_base = self.refs[ref_idx].get(pos) as usize;
            let base_quals = pos_base_quals(base_pileup, rt_err, false);
            quals.push(*(base_quals.get(contig_base).unwrap()));
        }
        quals
    }
}

pub fn pos_base_quals_helper<F>(reads: &Vec<(UmiType, u8, u8)>, rt_err: f64, lse: F) -> Vec<u8> where F: Fn(&Vec<f64>, f64) -> f64 {
    // log-probability of the transcript reads at that position given
    // that the real base is A, C, G, or T.
    let mut probs = vec![0.0_f64; 4];
    let factor1 = (1.0_f64 - rt_err).log10();
    let factor2 = (rt_err / 3.0_f64).log10();

    for (_, umi_reads) in &reads.iter().group_by(|x| x.0) {

        // base_probs[r][b] is the log-probability of the observed reads given that
        // the real base is r and the transcript/UMI base is b.
        let mut base_probs = vec![vec![0.0_f64; 4]; 4];

        for &(_, base, qual) in umi_reads {
            let adj_qual = cmp::min(qual, MAX_READ_QUAL);
            let missmatch_score = MISMATCH_SCORE_VS_QUAL[adj_qual as usize];//-(adj_qual as f64) / 10.0 - 3.0_f64.log10();
            let match_score = MATCH_SCORE_VS_QUAL[adj_qual as usize];//(1.0 - 10.0_f64.powf(-(adj_qual as f64) / 10.0)).log10();
            for r in 0..4 {
                for b in 0..4 {
                    if base as usize == b {
                        base_probs[r][b] = base_probs[r][b] + match_score;
                    } else {
                        base_probs[r][b] = base_probs[r][b] + missmatch_score;
                    }
                }
            }
        }

        for r in 0..4 {
            for b in 0..4 {
                if base_probs[r][b] >= 0.0 {
                    base_probs[r][b] = MIN_LOG_PROB;
                } else if b == r {
                    base_probs[r][b] = base_probs[r][b] + factor1;
                } else {
                    base_probs[r][b] = base_probs[r][b] + factor2;
                }
            }
        }

        for r in 0..4 {
            let mut umi_prob = lse(&base_probs[r], 10.0f64);

            if umi_prob > 0.0 {
                umi_prob = 0.0;
            }
            if umi_prob < MIN_LOG_PROB {
                umi_prob = MIN_LOG_PROB;
            }
            probs[r] = probs[r] + umi_prob;
        }

    }

    for r in 0..4 {
        if probs[r] >= 0.0 {
            probs[r] = MIN_LOG_PROB;
        }
    }

    let mut int_quals = vec![0u8; 4];
    let denominator = lse(&probs, 10.0f64);

    for r in 0..4 {
        let mut other_probs = Vec::new();
        for k in 0..4 {
            if k!=r { other_probs.push(probs[k]); }
        }
        let numerator = lse(&other_probs, 10.0f64);
        let mut final_qual = -10.0 * (numerator - denominator);

        if final_qual < 0.0 {
            final_qual = 0.0;
        } else if final_qual > MAX_OUT_QUAL as f64 {
            final_qual = MAX_OUT_QUAL as f64;
        }
        int_quals[r] = (final_qual as f32) as u8;
    }
    int_quals
}

// log-probability of the transcript reads at that position given
// that the real base is A, C, G, or T.
// The 'approx' argument specifies whether to compute the log sum exp
// approximately (faster) or exactly (slower)
pub fn pos_base_quals(reads: &Vec<(UmiType, u8, u8)>, rt_err: f64, approx: bool) -> Vec<u8> {
    if approx {
        pos_base_quals_helper(reads, rt_err, utils::logaddexp_approx_arr)
    } else {
        pos_base_quals_helper(reads, rt_err, utils::logaddexp_arr)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use debruijn::dna_string::DnaString;
    use graph_read::Read;

    #[test]
    fn test_find_read_matches_1() {
        let refs = vec!["CCCCACGTACGTGGGGGGA".to_string(), "ACGTACGTACGTGGGGG".to_string()];
        
        let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
        let mut align_helper = AlignHelper::new(&refs, scoring, 3, 8);

        let quals = vec![30; 8];
        let mut reads = Vec::new();
        reads.push(Read::new(0, 0, "0".to_string(), DnaString::from_dna_string("ACGTACCC"), quals.clone()));
        reads.push(Read::new(1, 0, "1".to_string(), DnaString::from_dna_string("TTTTTTTT"), quals.clone()));
        reads.push(Read::new(2, 0, "2".to_string(), DnaString::from_dna_string("ACGTACGT"), quals.clone()));
        reads.push(Read::new(3, 0, "3".to_string(), DnaString::from_dna_string("ACGAACGA"), quals.clone()));
        reads.push(Read::new(4, 0, "4".to_string(), DnaString::from_dna_string("TGGGGGTT"), quals.clone()));
        reads.push(Read::new(5, 0, "5".to_string(), DnaString::from_dna_string("TTTGGGGG"), quals.clone()));

        for read in reads.iter() {
            let matches = align_helper.find_read_matches(&read, 100);
            assert!(matches.is_none());
        }

        let matches = align_helper.find_read_matches(&reads[2], 15);
        assert_eq!(matches.unwrap().alignment.score, 16);

        align_helper.set_scoring(Scoring::from_scores(-4, -1, 2, -2).xclip(-5).yclip(0));
        let cigars = vec!["6=2X", "", "8=", "3=1X3=1X", "6=2X", "2X6="];
        let dists = vec![2, 0, 0, 2, 2, 2, 0];
        for i in 0..reads.len() {
            let matches = align_helper.find_read_matches(&reads[i], 0);
            if i == 1 {
                assert!(matches.is_none());
            } else {
                let unwrapped = matches.unwrap();
                assert_eq!(unwrapped.alignment.cigar(false), cigars[i]);
                assert_eq!(unwrapped.edit_distance(), dists[i]);
            }
        }

        align_helper.set_scoring(Scoring::from_scores(-4, -2, 2, -2).xclip(-3).yclip(0));
        let cigars = vec!["6=2S", "", "8=", "3=1X3=1X", "6=2S", "2S6=", "1="];
        let dists = vec![0, 0, 0, 2, 0, 0, 0];
        for i in 0..reads.len() {
            let matches = align_helper.find_read_matches(&reads[i], 0);
            if i == 1 {
                assert!(matches.is_none());
            } else {
                let unwrapped = matches.unwrap();
                assert_eq!(unwrapped.alignment.cigar(false), cigars[i]);
                assert_eq!(unwrapped.edit_distance(), dists[i]);
            }
        }
        
    }

    #[test]
    fn test_find_read_matches_2() {
        let refs = vec!["CCCCACGTACGTGGGGGGTTTT".to_string()];
        let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
        let mut align_helper = AlignHelper::new(&refs, scoring, 3, 5);

        let quals = vec![30; 8];
        let mut reads = Vec::new();
        reads.push(Read::new(1, 0, "1".to_string(), DnaString::from_dna_string("CCACGTCC"), quals.clone()));
        reads.push(Read::new(2, 0, "2".to_string(), DnaString::from_dna_string("ACGTACGT"), quals.clone()));
        reads.push(Read::new(3, 0, "3".to_string(), DnaString::from_dna_string("ACGAACGA"), quals.clone()));
        reads.push(Read::new(4, 0, "4".to_string(), DnaString::from_dna_string("ACGTTTTT"), quals.clone()));
        reads.push(Read::new(5, 0, "5".to_string(), DnaString::from_dna_string("ACGTTT"), quals.clone()));


        for read in reads.iter() {
            let matches = align_helper.find_read_matches(&read, 100);
            assert!(matches.is_none());
        }

        let matches = align_helper.find_read_matches(&reads[1], 15);
        assert_eq!(matches.unwrap().alignment.score, 16);

        align_helper.set_scoring(Scoring::from_scores(0, 0, 2, -2).xclip(-1).yclip(0));
        let matches = align_helper.find_read_matches(&reads[3], 15);
        assert_eq!(matches.unwrap().alignment.cigar(true), "4=6D4=");
        
        align_helper.set_scoring(Scoring::from_scores(-2, -1, 2, -2).xclip(-3).yclip(0));
        let cigars = vec!["6=1X1=", "8=", "3=1X3=1X", "4=6D4=", "4=2S"];
        let scores = vec![12, 16, 8, 8, 5];
        for i in 0..5 {
            let matches = align_helper.find_read_matches(&reads[i], 0);
            let unwrapped = matches.unwrap();
            assert_eq!(unwrapped.alignment.cigar(false), cigars[i]);
            assert_eq!(unwrapped.alignment.score, scores[i]);

        }
    }
}
