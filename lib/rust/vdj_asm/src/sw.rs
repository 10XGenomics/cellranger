//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use tada::kmer;
use tada::bitenc;
use std::cmp;
use utils;
use graph_read::Read;
use itertools::Itertools;
use std::collections::{HashMap};
use constants::{EPSILON, UmiType};
use bio;
use rust_htslib::bam;
use time::PreciseTime;
use std::borrow::Borrow;

const MAX_OUT_QUAL : u8 = 60;
const MAX_READ_QUAL : u8 = 30;
const MIN_LOG_PROB : f64 = -100000.0;

/// Smith-Waterman alignment parameters.
/// All of them should be non-negative (algos add the negative sign for penalties)
pub struct AlignParams {
    pub match_score: f32,
    pub miss_score: f32,
    pub gap_open: f32,
    pub gap_extend: f32,
    pub clip: f32,
}

impl AlignParams {
    pub fn new(match_score: f32, miss_score: f32, gap_open: f32, gap_extend: f32, clip: f32) -> AlignParams {
        assert!(gap_open >= 0.0);
        assert!(gap_extend >= 0.0);
        assert!(miss_score >= 0.0);

        AlignParams{
            match_score: match_score,
            miss_score: miss_score,
            gap_open: gap_open,
            gap_extend: gap_extend,
            clip: clip,
        }
    }
}

/// Smith-Waterman step.
#[derive(Eq, PartialEq, Debug, Clone, Copy)]
pub enum AlignmentStep {
    Match,
    Mismatch,
    Clip,
    Del,
    Ins,
}

/// Alignment against a reference
#[derive(Debug, Clone)]
pub struct Alignment {
    pub score: f32, // alignment score
    pub ref_idx: usize, // index of reference
    pub start_pos: usize, // starting position on the reference
    pub alignment: Vec<AlignmentStep>,
}

impl Alignment {
    pub fn new(score: f32, ref_idx:usize, pos:usize, alignment: Vec<AlignmentStep>) -> Alignment {
        Alignment {
            score: score,
            ref_idx: ref_idx,
            start_pos: pos,
            alignment: alignment,
        }
    }

    pub fn edit_distance(&self) -> usize {
        let mut d = 0;
        for step in self.alignment.iter() {
            let new_d = ((*step != AlignmentStep::Match) && (*step != AlignmentStep::Clip)) as usize;
            d += new_d;
        }
        d
    }

    /// Last aligned position
    pub fn alignment_end(&self) -> usize {
        let mut pos = self.start_pos;
        for step in self.alignment.iter() {
            pos += match *step {
                AlignmentStep::Match | AlignmentStep::Mismatch | AlignmentStep::Del => 1,
                _ => 0,
            };
        }
        pos
    }

    pub fn cigar(&self) -> String {
        let mut cigar_str = "".to_owned();
        if self.alignment.len() == 0 {
            return cigar_str;
        }
        let mut last_step = self.alignment[0].clone();
        let mut last_len = 1;
        for step in self.alignment.iter().skip(1) {
            if *step != last_step {
                let new_cigar_str = match last_step {
                    AlignmentStep::Clip => "S",
                    AlignmentStep::Del => "D",
                    AlignmentStep::Ins => "I",
                    AlignmentStep::Mismatch => "X",
                    _ => "=",
                };
                cigar_str = format!("{}{}{}", cigar_str, last_len, new_cigar_str);
                last_len = 1;
                last_step = (*step).clone();
            } else {
                last_len += 1;
            }
        }

        let new_cigar_str = match last_step {
            AlignmentStep::Clip => "S",
            AlignmentStep::Del => "D",
            AlignmentStep::Ins => "I",
            AlignmentStep::Mismatch => "X",
            _ => "=",
        };
        cigar_str = format!("{}{}{}", cigar_str, last_len, new_cigar_str);

        cigar_str
    }

    pub fn bam_cigar(&self) -> Vec<bam::record::Cigar> {
        let mut cigar = Vec::new();
        if self.alignment.len() == 0 {
            return cigar;
        }
        let mut last_step = self.alignment[0].clone();
        let mut last_len = 1;
        for step in self.alignment.iter().skip(1) {
            if *step != last_step {
                cigar.push(match last_step {
                    AlignmentStep::Clip => bam::record::Cigar::SoftClip(last_len),
                    AlignmentStep::Del => bam::record::Cigar::Del(last_len),
                    AlignmentStep::Ins => bam::record::Cigar::Ins(last_len),
                    AlignmentStep::Mismatch => bam::record::Cigar::Diff(last_len),
                    _ => bam::record::Cigar::Equal(last_len),

                });
                last_len = 1;
                last_step = (*step).clone();
            } else {
                last_len += 1;
            }
        }

        cigar.push(match last_step {
            AlignmentStep::Clip => bam::record::Cigar::SoftClip(last_len),
            AlignmentStep::Del => bam::record::Cigar::Del(last_len),
            AlignmentStep::Ins => bam::record::Cigar::Ins(last_len),
            AlignmentStep::Mismatch => bam::record::Cigar::Diff(last_len),
            _ => bam::record::Cigar::Equal(last_len)
        });
        cigar
    }
}

fn steps_from_bio_steps(bio_alignment: &bio::alignment::Alignment, seq_len: usize) -> Vec<AlignmentStep> {
    let mut steps = Vec::new();
    for _ in 0..bio_alignment.xstart {
        steps.push(AlignmentStep::Clip);
    }
    for tmp_step in bio_alignment.operations.clone() {
        match tmp_step {
            bio::alignment::AlignmentOperation::Match => steps.push(AlignmentStep::Match),
            bio::alignment::AlignmentOperation::Subst => steps.push(AlignmentStep::Mismatch),
            bio::alignment::AlignmentOperation::Del => steps.push(AlignmentStep::Del),
            bio::alignment::AlignmentOperation::Ins => steps.push(AlignmentStep::Ins),
        }
    }
    for _ in 0..(seq_len - bio_alignment.xend) {
        steps.push(AlignmentStep::Clip);
    }
    steps
}

pub struct Aligner {
    pub refs: Vec<bitenc::BitEnc>,
    // From kmer to a list of reference sequences and positions where the kmer appears
    pub kmers: HashMap<kmer::Lmer, Vec<(usize, usize)>>,
    pub k: usize,
}

impl Aligner {
    /// seed_length: length of seed in base-pairs
    pub fn new(refs: &Vec<String>, seed_length: usize) -> Aligner {

        let mut kmers : HashMap<kmer::Lmer, Vec<(usize, usize)>> = HashMap::new();
        let mut bitenc_refs = Vec::new();

        let cc_start = PreciseTime::now();

        for (ref_idx, ref_seq) in refs.iter().enumerate() {
            let bitenc_ref = bitenc::BitEnc::from_dna_string(ref_seq);

            for &(kmer, pos) in bitenc_ref.lmers(seed_length).iter() {
                if kmers.contains_key(&kmer) {
                    let curr_val = kmers.get_mut(&kmer).unwrap();
                    (*curr_val).push((ref_idx, pos));
                } else {
                    kmers.insert(kmer, vec![(ref_idx, pos)]);
                }
            }
            bitenc_refs.push(bitenc_ref);
        }
        println!("Build index in {}sec: {} kmers", cc_start.to(PreciseTime::now()), kmers.len());

        Aligner {
            refs:bitenc_refs,
            kmers:kmers,
            k:seed_length,
        }
    }

    pub fn global_align(&self, read_bitenc: &bitenc::BitEnc, params: &AlignParams) -> Alignment {
        let read_seq = read_bitenc.to_dna_string();
        let score = |a: u8, b: u8| if a == b {params.match_score as i32} else {-params.miss_score as i32};

        let mut max_alignment : Option<bio::alignment::Alignment> = None;
        let mut max_ref_idx = 0;
        let mut max_score = 0;

        for (ref_idx, ref_bitenc) in self.refs.iter().enumerate() {
            let ref_seq = ref_bitenc.to_dna_string();
            let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(read_seq.len(), ref_seq.len(),
                                                                           -params.gap_open as i32, -params.gap_extend as i32, &score);
            let tmp_alignment = aligner.global(read_seq.as_bytes(), ref_seq.as_bytes());
            //println!("{:?}", tmp_alignment.pretty(read_seq.as_bytes(), ref_seq.as_bytes()));

            if max_alignment.is_none() || max_score < tmp_alignment.score {
                max_score = tmp_alignment.score;
                max_alignment = Some(tmp_alignment);
                max_ref_idx = ref_idx;
            }
        }
        let final_alignment = max_alignment.unwrap();
        assert_eq!(final_alignment.xstart, 0);
        assert_eq!(final_alignment.ystart, 0);

        let steps = steps_from_bio_steps(&final_alignment, read_seq.len());
        let alignment = Alignment::new(final_alignment.score as f32, max_ref_idx, final_alignment.ystart, steps);
        alignment
    }

    /// SW of a read in the vicinity of a kmer match.
    #[inline(never)]
    fn align_to_partial_ref(&self, ref_idx: usize, ref_pos: usize, extent: usize, 
                            read: &Read, params: &AlignParams, local:bool) -> Alignment {

        let mut align_end = ref_pos + extent;
        if align_end > self.refs[ref_idx].len() {
            align_end = self.refs[ref_idx].len()
        }
        let mut align_start = 0;
        if ref_pos > extent {
            align_start = ref_pos - extent;
        }

        let read_seq = read.seq.to_dna_string();
        let partial_ref = bitenc::BitEnc::from_bytes(&self.refs[ref_idx].iter().skip(align_start).take(align_end - align_start).collect());
        let ref_seq = partial_ref.to_dna_string();

        let score = |a: u8, b: u8| if a == b {params.match_score as i32} else {-params.miss_score as i32};
        let mut aligner = bio::alignment::pairwise::Aligner::with_capacity(read_seq.len(), ref_seq.len(),
                                                                           -params.gap_open as i32, -params.gap_extend as i32, &score);
        // read is globally aligned, partial_ref is local
        let tmp_alignment = match local {
            true => aligner.local(read_seq.as_bytes(), ref_seq.as_bytes()),
            false => aligner.semiglobal(read_seq.as_bytes(), ref_seq.as_bytes()),
        };

        assert!(tmp_alignment.xstart == 0 || local);
        assert!(tmp_alignment.xend == read.seq.len() || local);

        let steps = steps_from_bio_steps(&tmp_alignment, read_seq.len());

        let alignment = Alignment::new(tmp_alignment.score as f32, ref_idx, align_start + tmp_alignment.ystart, steps);
        alignment
    }

    /// Find matches of a read with at least a minimum SW score.
    ///
    /// # Args
    /// - read: Read object. We'll try to find a match with at least the input SW score.
    /// Note that this function DOES NOT try to find the best match. Just A match that exceeds
    /// the minimum requested score.
    /// - params: SW parameters
    ///
    /// # Return value
    /// Some(Alignment) if we could find a good enough alignment, otherwise None.
    pub fn find_read_matches(&self, read: &Read, params: &AlignParams, 
                             min_align_score: f32, 
                             local:bool) -> Option<Alignment> {

        /// Get maximum possible gap assuming no mismatches
        /// Total_score = (len - gap) * match - gap_open - gap_extend * (gap - 1)
        /// Total = len * match - gap * match - gap_open - gap_extend * gap - gap_extend
        /// gap * (match + gap_extend) = len * match - gap_open - gap_extend - total
        let numerator = (read.len() as f32) * params.match_score - params.gap_open - params.gap_extend - min_align_score;
        let mut max_gap_len =  numerator / (params.match_score + params.gap_extend);
        if max_gap_len < 0 as f32 {
            max_gap_len = 0 as f32;
        } else if max_gap_len > 2.0 * (read.len() as f32) {
            max_gap_len = 2.0 * (read.len() as f32); // Make sure you never extend too much.
        }

        for &(kmer, _) in read.seq.lmers(self.k).iter() {
            if self.kmers.contains_key(&kmer) {
                let kmer_poses = self.kmers.get(&kmer).unwrap();

                for &(idx, pos) in kmer_poses.iter() {

                    let alignment = self.align_to_partial_ref(idx, pos, 
                                                              read.len() + max_gap_len as usize, 
                                                              &read, &params, local);
                    if alignment.score > min_align_score as f32 {
                        assert!(read.validate_alignment(&alignment));
                        return Some(alignment);
                    }
                }
            }
        }
        None
    }

    /// Gets a match path of LCSk++ and computes the implied score and alignment.
    ///
    /// # Returns
    /// A tuple (score, alignment_steps).
    fn alignment_from_sparse_matches(&self, match_path: &Vec<(u32, u32)>,
                                     params: &AlignParams, read_len: usize) -> (f32, Vec<AlignmentStep>) {
        let mut steps = Vec::new();
        let mut score = 0.0;

        for (step_idx, s) in match_path.iter().enumerate() {
            if step_idx == 0 {
                for _ in 0..(s.0) {
                    score -= params.clip;
                    steps.push(AlignmentStep::Clip);
                }
                // Each step suggests a kmer match.
                for _ in 0..self.k {
                    score += params.match_score;
                    steps.push(AlignmentStep::Match);
                }
            } else {
                // If previous step was (i, j) and this step is (i+1, j+1),
                // then we just advanced the kmer match by one position.
                // Otherwise, we skipped some bases in either one or both of
                // the sequences and then we added a new kmer match.
                let prev_step = match_path[step_idx - 1];
                let skip_ref = s.1 as f32 - prev_step.1 as f32;
                let skip_read = s.0 as f32 - prev_step.0 as f32;

                if skip_ref > 1.0 || skip_read > 1.0 {
                    let mut indel_len = skip_ref - skip_read;
                    let step_type = match indel_len > 0.0 {
                        true => AlignmentStep::Del,
                        false => AlignmentStep::Ins,
                    };
                    if indel_len < 0.0 {
                        indel_len = -indel_len;
                    }
                    for _ in 0..(indel_len as usize) {
                        steps.push(step_type);
                    }
                    score -= params.gap_open + (indel_len - 1.0) * params.gap_extend;
                    // After the indel, there's a new kmer match.
                    for _ in 0..self.k {
                        score += params.match_score;
                        steps.push(AlignmentStep::Match);
                    }
                } else {
                    steps.push(AlignmentStep::Match);
                    score += params.match_score;
                }
            }
            if step_idx == match_path.len() - 1 {
                for _ in 0..(read_len - (s.0 as usize) - self.k) {
                    steps.push(AlignmentStep::Clip);
                    score -= params.clip;
                }
            }
        }
        (score, steps)
    }

    /// This is similar to find_read_matches, but uses a different algorithm which is much
    /// faster than doing a full local SW.
    pub fn find_read_matches_sparse(&self, read_seq: &bitenc::BitEnc, params: &AlignParams,
                                    min_align_score: f32) -> Option<Alignment> {

        let mut all_kmer_matches = Vec::new();

        for &(kmer, kmer_pos) in read_seq.lmers(self.k).iter() {
            match self.kmers.get(&kmer) {
                Some(ref_pos_list) => {
                    for ref_pos in ref_pos_list {
                        all_kmer_matches.push((ref_pos.0, kmer_pos, ref_pos.1));
                    }
                },
                None => {},
            }
        }

        all_kmer_matches.sort();

        for (ref_idx, matches) in all_kmer_matches.iter().group_by(|x| x.0) {
            let tmp_matches = matches.iter().map(|x| (x.1 as u32, x.2 as u32)).collect();
            let sparse_al = bio::alignment::sparse::lcskpp(&tmp_matches, self.k);
            let mut match_path: Vec<(u32,u32)> = sparse_al.path.iter().map(|i| tmp_matches[*i]).collect();
            if match_path.len() == 0 {
                match_path = vec![tmp_matches[0]];
            }
            let (score, steps) = self.alignment_from_sparse_matches(&match_path, &params, read_seq.len());
            if score > min_align_score {
                return Some(Alignment::new(score, ref_idx, match_path[0].1 as usize, steps));
            }
        }

        None
    }

    pub fn pileup<T: Borrow<Read>>(&self, ref_idx: usize, alignments: &Vec<(Alignment, T)>) -> Vec<Vec<(u32, u8, u8)>> {
        let ref_len = self.refs[ref_idx].len();
        let mut pos_pileup = Vec::with_capacity(ref_len);
        for _ in 0..ref_len {
            pos_pileup.push(Vec::new());
        }

        for &(ref alignment, ref read_borrow) in alignments.iter() {
            let read = read_borrow.borrow();
            let mut ref_idx = alignment.start_pos;
            let mut query_idx = 0;

            for step in alignment.alignment.iter() {
                match *step {
                    AlignmentStep::Clip => {
                        query_idx += 1
                    },
                    AlignmentStep::Ins => {
                        query_idx += 1
                    },
                    AlignmentStep::Del => {
                        ref_idx += 1
                    },
                    _ => {
                        pos_pileup[ref_idx].push((read.umi,
                                                  read.seq.get(query_idx).unwrap(),
                                                  read.quals[query_idx]));
                        ref_idx += 1;
                        query_idx += 1;
                    },
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

            let contig_base = self.refs[ref_idx].get(pos).unwrap() as usize;
            let base_quals = pos_base_quals(base_pileup, rt_err);
            quals.push(*(base_quals.get(contig_base).unwrap()));
        }
        quals
    }
}

pub fn pos_base_quals(reads: &Vec<(UmiType, u8, u8)>, rt_err: f64) -> Vec<u8> {
    // log-probability of the transcript reads at that position given
    // that the real base is A, C, G, or T.
    let mut probs = vec![0.0_f64; 4];

    for (_, umi_reads) in reads.iter().group_by(|x| x.0) {

        // base_probs[r][b] is the log-probability of the observed reads given that
        // the real base is r and the transcript/UMI base is b.
        let mut base_probs = vec![vec![0.0_f64; 4]; 4];

        for &(_, base, qual) in umi_reads {
            let adj_qual = cmp::min(qual, MAX_READ_QUAL);
            let missmatch_score = -(adj_qual as f64) / 10.0 - 3.0_f64.log10();
            let match_score = (1.0 - 10.0_f64.powf(-(adj_qual as f64) / 10.0)).log10();
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
                    base_probs[r][b] = base_probs[r][b] + (1.0_f64 - rt_err).log10();
                } else {
                    base_probs[r][b] = base_probs[r][b] + (rt_err / 3.0_f64).log10();
                }
            }
        }

        for r in 0..4 {
            let mut umi_prob = utils::logaddexp_arr(&base_probs[r], 10.0);

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

    let mut denominator = utils::logaddexp_arr(&probs, 10.0);

    for r in 0..4 {
        let qual = probs[r];

        // Deal with numerical issues.
        if denominator > 0.0 {
            denominator = 0.0;
        }
        if denominator <= qual {
            denominator = qual + EPSILON;
        }
        let mut final_qual = -10.0 * (1.0 - 10.0_f64.powf(qual - denominator)).log10();
        if final_qual < 0.0 {
            final_qual = 0.0;
        } else if final_qual > MAX_OUT_QUAL as f64 {
            final_qual = MAX_OUT_QUAL as f64;
        }
        int_quals[r] = (final_qual as f32) as u8;
    }

    int_quals
}


#[cfg(test)]
mod tests {
    use super::*;
    use tada::bitenc;
    use tada::kmer;
    use graph_read::Read;
    use bam;

    #[test]
    fn test_cigar() {
        {
            let steps = vec![AlignmentStep::Match, AlignmentStep::Mismatch, AlignmentStep::Mismatch,
            AlignmentStep::Ins, AlignmentStep::Del, AlignmentStep::Del, AlignmentStep::Clip];

            let al = Alignment::new(10.0, 0, 0, steps.clone());
            assert_eq!(al.cigar(), "1=2X1I2D1S");
            assert_eq!(al.bam_cigar(), vec![bam::record::Cigar::Equal(1),
                                            bam::record::Cigar::Diff(2),
                                            bam::record::Cigar::Ins(1),
                                            bam::record::Cigar::Del(2),
                                            bam::record::Cigar::SoftClip(1)]);
        }
        {
            let steps = vec![AlignmentStep::Mismatch, AlignmentStep::Match, AlignmentStep::Mismatch];

            let al = Alignment::new(10.0, 0, 0, steps.clone());
            assert_eq!(al.cigar(), "1X1=1X");
            assert_eq!(al.bam_cigar(), vec![bam::record::Cigar::Diff(1),
                                            bam::record::Cigar::Equal(1),
                                            bam::record::Cigar::Diff(1)]);
        }
        {
            let steps = vec![AlignmentStep::Match];

            let al = Alignment::new(10.0, 0, 0, steps.clone());
            assert_eq!(al.cigar(), "1=");
            assert_eq!(al.bam_cigar(), vec![bam::record::Cigar::Equal(1)]);
        }
    }

    #[test]
    fn test_alignment_end() {
        let mut steps = vec![AlignmentStep::Match, AlignmentStep::Match];

        let al = Alignment::new(10.0, 0, 0, steps.clone());
        assert_eq!(al.alignment_end(), 2);

        steps.push(AlignmentStep::Ins);
        let al = Alignment::new(10.0, 0, 0, steps.clone());
        assert_eq!(al.alignment_end(), 2);

        steps.push(AlignmentStep::Del);
        let al = Alignment::new(10.0, 0, 0, steps.clone());
        assert_eq!(al.alignment_end(), 3);

        steps.push(AlignmentStep::Clip);
        let al = Alignment::new(10.0, 0, 0, steps);
        assert_eq!(al.alignment_end(), 3);

    }

    #[test]
    fn test_aligner() {
        let refs = vec!["ACGTACGTACGTGGGGG".to_string(),
                        "CCCCACGTACGTGGGGGGA".to_string()];
        let aligner = Aligner::new(&refs, 4);

        let a = kmer::base_to_bits(b'A');
        let c = kmer::base_to_bits(b'C');
        let g = kmer::base_to_bits(b'G');
        let t = kmer::base_to_bits(b'T');
        let acgt = kmer::Lmer::new(vec![a, c, g, t].iter());

        assert!(aligner.kmers.contains_key(&acgt));
        assert!(*aligner.kmers.get(&acgt).unwrap() == vec![(0, 0), (0, 4), (0, 8), (1, 4), (1, 8)]);
        assert_eq!(aligner.kmers.len(), 13); // 13 distinct kmers

        let gggg = kmer::Lmer::new(vec![g, g, g, g].iter());
        assert!(*aligner.kmers.get(&gggg).unwrap() == vec![(0, 12), (0, 13), (1, 12), (1, 13), (1, 14)]);

        let aligner = Aligner::new(&refs, 1);
        assert_eq!(aligner.kmers.len(), 4);
    }

    #[test]
    fn test_global_align() {
        let refs = vec!["ACAAGTACGTAGTACTGAA".to_string(),
                        "AACCACGTACGTGGGGGGA".to_string()];

        let aligner = Aligner::new(&refs, 20);
        let ref_seq1 = bitenc::BitEnc::from_dna_string(&refs[0]);
        let ref_seq2 = bitenc::BitEnc::from_dna_string("CCACGTACGT");
        let ref_seq3 = bitenc::BitEnc::from_dna_string("AAGGGGGG");
        let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);

        let alignment1 = aligner.global_align(&ref_seq1, &params);
        assert_eq!(alignment1.cigar(), "19=");

        // AACCACGTACGTGGGGGGA
        //   CCACGTACGT
        let alignment2 = aligner.global_align(&ref_seq2, &params);
        assert_eq!(alignment2.ref_idx, 1);
        assert_eq!(alignment2.score, -2.0);
        assert_eq!(alignment2.cigar(), "2D10=7D");

        let params = AlignParams::new(2.0, 2.0, 400.0, 2.0, 1.0);
        let alignment3 = aligner.global_align(&ref_seq3, &params);
        assert_eq!(alignment3.ref_idx, 1);
    }

    #[test]
    fn test_find_read_matches() {
        let refs = vec!["CCCCACGTACGTGGGGGGA".to_string(), "ACGTACGTACGTGGGGG".to_string()];
        {
            let aligner = Aligner::new(&refs, 4);

            let quals = vec![30; 8];
            let mut reads = Vec::new();
            reads.push(Read::new(0, 0, "0".to_string(), bitenc::BitEnc::from_dna_string("ACGTACCC"), quals.clone()));
            reads.push(Read::new(1, 0, "1".to_string(), bitenc::BitEnc::from_dna_string("TTTTTTTT"), quals.clone())); // no kmer match
            reads.push(Read::new(2, 0, "2".to_string(), bitenc::BitEnc::from_dna_string("ACGTACGT"), quals.clone()));
            reads.push(Read::new(3, 0, "3".to_string(), bitenc::BitEnc::from_dna_string("ACGAACGA"), quals.clone())); // no kmer match
            reads.push(Read::new(4, 0, "4".to_string(), bitenc::BitEnc::from_dna_string("TGGGGGTT"), quals.clone()));
            reads.push(Read::new(5, 0, "5".to_string(), bitenc::BitEnc::from_dna_string("TTTGGGGG"), quals.clone()));
            reads.push(Read::new(6, 0, "6".to_string(), bitenc::BitEnc::from_dna_string("T"), vec![30;1]));

            let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
            for read in reads.iter() {
                let matches = aligner.find_read_matches(&read, &params, 100.0, false);
                assert!(matches.is_none());
            }

            let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
            let matches = aligner.find_read_matches(&reads[2], &params, 15.9, false);
            assert_eq!(matches.unwrap().score, 16.0);

            let params = AlignParams::new(2.0, 2.0, 4.0, 1.0, 1.0);
            let cigars = vec!["6=2X", "", "8=", "", "6=2X", "2X6=", ""];
            let dists = vec![2, 0, 0, 0, 2, 2, 0];
            for i in 0..reads.len() {
                let matches = aligner.find_read_matches(&reads[i], &params, 0.0, false);
                if i == 1 || i == 3 || i == 6 {
                    assert!(matches.is_none());
                } else {
                    let unwrapped = matches.unwrap();
                    assert_eq!(unwrapped.cigar(), cigars[i]);
                    assert_eq!(unwrapped.edit_distance(), dists[i]);
                }
            }

            let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
            let cigars = vec!["6=2S", "", "8=", "", "6=2S", "2S6=", ""];
            let dists = vec![0, 0, 0, 0, 0, 0, 0];
            for i in 0..reads.len() {
                let matches = aligner.find_read_matches(&reads[i], &params, 0.0, true);
                if i == 1 || i == 3 || i == 6 {
                    assert!(matches.is_none());
                } else {
                    let unwrapped = matches.unwrap();
                    assert_eq!(unwrapped.cigar(), cigars[i]);
                    assert_eq!(unwrapped.edit_distance(), dists[i]);
                }
            }
        }
        {
            let aligner = Aligner::new(&refs, 1);
            let reads = Read::new(5, 0, "1".to_string(), bitenc::BitEnc::from_dna_string("T"), vec![30;1]);
            let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
            let matches = aligner.find_read_matches(&reads, &params, 0.0, true);
            assert!(matches.unwrap().cigar() == "1=");
        }
    }

    #[test]
    fn test_find_read_matches_sparse() {
        let refs = vec!["CCCCACGTACGTGGGGGGTTTT".to_string()];
        {
            let aligner = Aligner::new(&refs, 4);

            let quals = vec![30; 8];
            let mut reads = Vec::new();
            reads.push(Read::new(1, 0, "1".to_string(), bitenc::BitEnc::from_dna_string("CCACGTCC"), quals.clone()));
            reads.push(Read::new(2, 0, "2".to_string(), bitenc::BitEnc::from_dna_string("ACGTACGT"), quals.clone()));
            reads.push(Read::new(3, 0, "3".to_string(), bitenc::BitEnc::from_dna_string("ACGAACGA"), quals.clone())); // no kmer match
            reads.push(Read::new(4, 0, "4".to_string(), bitenc::BitEnc::from_dna_string("ACGTTTTT"), quals.clone()));
            reads.push(Read::new(5, 0, "5".to_string(), bitenc::BitEnc::from_dna_string("ACGTCGT"), quals.clone())); // the last 3 can't match with k=4

            {
                let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
                for read in reads.iter() {
                    let matches = aligner.find_read_matches_sparse(&read.seq, &params, 100.0);
                    assert!(matches.is_none());
                }
            }
            {
                let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
                let matches = aligner.find_read_matches_sparse(&reads[1].seq, &params, 15.0);
                assert_eq!(matches.unwrap().score, 16.0);
            }
            {
                let params = AlignParams::new(2.0, 2.0, 0.0, 0.0, 1.0);
                let matches = aligner.find_read_matches_sparse(&reads[3].seq, &params, 15.0);
                assert_eq!(matches.unwrap().cigar(), "4=6D4=");
            }
            {
                let params = AlignParams::new(2.0, 2.0, 2.0, 1.0, 1.0);
                let cigars = vec!["6=2S", "8=", "", "4=6D4=", "4=3S"];
                let scores = vec![10.0, 16.0, 0.0, 9.0, 5.0];
                for i in 0..5 {
                    let matches = aligner.find_read_matches_sparse(&reads[i].seq, &params, 0.0);
                    if i == 2 {
                        assert!(matches.is_none());
                    } else {
                        let unwrapped = matches.unwrap();
                        assert_eq!(unwrapped.cigar(), cigars[i]);
                        assert_eq!(unwrapped.score, scores[i]);
                    }
                }
            }
        }
        {
            let aligner = Aligner::new(&refs, 1);
            let params = AlignParams::new(2.0, 2.0, 4.0, 2.0, 1.0);
            let matches = aligner.find_read_matches_sparse(&bitenc::BitEnc::from_dna_string("T"), &params, 0.0);
            assert_eq!(matches.unwrap().cigar(), "1=");
        }
    }
}
