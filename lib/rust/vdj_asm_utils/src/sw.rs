use crate::constants::UmiType;
use bio;
use bio::alignment::{Alignment, AlignmentOperation};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::cmp;
pub type Scoring = bio::alignment::pairwise::Scoring<bio::alignment::pairwise::MatchParams>;

const MAX_OUT_QUAL: u8 = 60;
const MAX_READ_QUAL: u8 = 30;
const MIN_LOG_PROB: f64 = -100000.0;

lazy_static! {
    static ref MATCH_SCORE_VS_QUAL: [f64; (MAX_READ_QUAL + 1) as usize] = {
        let mut table = [0.0f64; (MAX_READ_QUAL + 1) as usize];
        for (q, v) in table.iter_mut().enumerate() {
            *v = (1.0 - 10.0_f64.powf(-(q as f64) / 10.0)).log10();
        }
        table
    };
    static ref MISMATCH_SCORE_VS_QUAL: [f64; (MAX_READ_QUAL + 1) as usize] = {
        let mut table = [0.0f64; (MAX_READ_QUAL + 1) as usize];
        for (q, v) in table.iter_mut().enumerate() {
            *v = -(q as f64) / 10.0 - 3.0_f64.log10();
        }
        table
    };
}

/// Alignment against a reference
/// This struct pack the index of the reference and the alignment together.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentPacket {
    pub ref_idx: usize, // index of reference
    pub alignment: Alignment,
}

impl AlignmentPacket {
    pub fn new(ref_idx: usize, alignment: Alignment) -> Self {
        AlignmentPacket { ref_idx, alignment }
    }

    /// Compute the edit distance implied by the alignment.
    /// The clipped regions are not considered in the computation
    pub fn edit_distance(&self) -> usize {
        let mut d = 0;
        for &op in &self.alignment.operations {
            d += match op {
                AlignmentOperation::Match => 0,
                AlignmentOperation::Xclip(_) => 0,
                AlignmentOperation::Yclip(_) => 0,
                _ => 1,
            };
        }
        d
    }
}

pub fn pos_base_quals_helper<F>(reads: &[(UmiType, u8, u8)], rt_err: f64, lse: F) -> [u8; 4]
where
    F: Fn(&[f64], f64) -> f64,
{
    // log-probability of the transcript reads at that position given
    // that the real base is A, C, G, or T.
    let mut probs = [0.0_f64; 4];
    let factor1 = (1.0_f64 - rt_err).log10();
    let factor2 = (rt_err / 3.0_f64).log10();

    for (_, umi_reads) in &reads.iter().group_by(|x| x.0) {
        // base_probs[r][b] is the log-probability of the observed reads given that
        // the real base is r and the transcript/UMI base is b.
        let mut base_probs = [[0.0_f64; 4]; 4];

        for &(_, base, qual) in umi_reads {
            let adj_qual = cmp::min(qual, MAX_READ_QUAL);
            let missmatch_score = MISMATCH_SCORE_VS_QUAL[adj_qual as usize]; //-(adj_qual as f64) / 10.0 - 3.0_f64.log10();
            let match_score = MATCH_SCORE_VS_QUAL[adj_qual as usize]; //(1.0 - 10.0_f64.powf(-(adj_qual as f64) / 10.0)).log10();
            for probs in &mut base_probs {
                for (b, prob) in probs.iter_mut().enumerate() {
                    *prob += if base as usize == b {
                        match_score
                    } else {
                        missmatch_score
                    };
                }
            }
        }

        for (r, probs) in base_probs.iter_mut().enumerate() {
            for (b, prob) in probs.iter_mut().enumerate() {
                if *prob >= 0.0 {
                    *prob = MIN_LOG_PROB;
                } else if b == r {
                    *prob += factor1;
                } else {
                    *prob += factor2;
                }
            }
        }

        for r in 0..4 {
            let mut umi_prob = lse(&base_probs[r], 10.0f64);

            umi_prob = umi_prob.clamp(MIN_LOG_PROB, 0.0);
            probs[r] += umi_prob;
        }
    }

    for prob in &mut probs {
        if *prob >= 0.0 {
            *prob = MIN_LOG_PROB;
        }
    }

    let mut int_quals = [0u8; 4];
    let denominator = lse(&probs, 10.0f64);

    for (r, int_qual) in int_quals.iter_mut().enumerate() {
        let mut other_probs = [0_f64; 3];
        other_probs[..r].clone_from_slice(&probs[..r]);
        other_probs[r..].clone_from_slice(&probs[r + 1..]);
        let numerator = lse(&other_probs, 10.0f64);
        let mut final_qual = -10.0 * (numerator - denominator);

        final_qual = final_qual.clamp(0.0, MAX_OUT_QUAL as f64);
        *int_qual = (final_qual as f32) as u8;
    }
    int_quals
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_read::Read;
    use bio::alignment::pairwise::banded;
    use bio::alignment::sparse::HashMapFx;
    use debruijn::dna_string::DnaString;

    struct AlignHelper<'a> {
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
        pub fn new(refs: &'a [&str], scoring: Scoring, k: usize, w: usize) -> Self {
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
                    kmers_hash
                        .entry(&seq[i..i + k])
                        .or_default()
                        .push((ref_idx, i as u32));
                }
            }
            AlignHelper {
                refs: dna_strings,
                kmers_hash,
                scoring,
                k,
                w,
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
        pub fn find_read_matches(
            &self,
            read: &Read,
            min_align_score: i32,
        ) -> Option<AlignmentPacket> {
            let mut all_kmer_matches: Vec<(i32, usize, usize, usize)> = Vec::new();
            let read_seq = read.seq.to_string();
            let read_seq = read_seq.as_bytes();
            let mut ref_kmer_counts = vec![0i32; self.refs.len()];
            for i in 0..(read_seq.len() + 1).saturating_sub(self.k) {
                let slc = &read_seq[i..i + self.k];
                if let Some(matches) = self.kmers_hash.get(slc) {
                    for &ref_pos in matches {
                        all_kmer_matches.push((0, ref_pos.0, i, ref_pos.1 as usize));
                        ref_kmer_counts[ref_pos.0] += 1;
                    }
                }
            }
            for i in 0..all_kmer_matches.len() {
                all_kmer_matches[i].0 = -ref_kmer_counts[all_kmer_matches[i].1];
            }
            all_kmer_matches.sort_unstable();
            let mut aligner = banded::Aligner::with_scoring(self.scoring, self.k, self.w);
            for (ref_idx, match_group) in &all_kmer_matches.into_iter().group_by(|x| x.1) {
                let matches = match_group
                    .into_iter()
                    .map(|x| (x.2 as u32, x.3 as u32))
                    .collect();
                let ref_seq = &self.refs[ref_idx].to_string();
                let alignment = aligner.custom_with_expanded_matches(
                    read_seq,
                    ref_seq.as_bytes(),
                    matches,
                    Some(1),
                    true,
                );
                if alignment.score > min_align_score {
                    let align_pack = AlignmentPacket::new(ref_idx, alignment);
                    assert!(read.validate_alignment(&align_pack));
                    return Some(align_pack);
                }
            }
            None
        }
    }
    #[test]
    fn test_find_read_matches_1() {
        let refs = ["CCCCACGTACGTGGGGGGA", "ACGTACGTACGTGGGGG"];

        let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
        let mut align_helper = AlignHelper::new(&refs, scoring, 3, 8);

        let quals = vec![30; 8];
        let reads = vec![
            Read::new(
                0,
                0,
                "0".to_string(),
                DnaString::from_dna_string("ACGTACCC"),
                quals.clone(),
            ),
            Read::new(
                1,
                0,
                "1".to_string(),
                DnaString::from_dna_string("TTTTTTTT"),
                quals.clone(),
            ),
            Read::new(
                2,
                0,
                "2".to_string(),
                DnaString::from_dna_string("ACGTACGT"),
                quals.clone(),
            ),
            Read::new(
                3,
                0,
                "3".to_string(),
                DnaString::from_dna_string("ACGAACGA"),
                quals.clone(),
            ),
            Read::new(
                4,
                0,
                "4".to_string(),
                DnaString::from_dna_string("TGGGGGTT"),
                quals.clone(),
            ),
            Read::new(
                5,
                0,
                "5".to_string(),
                DnaString::from_dna_string("TTTGGGGG"),
                quals,
            ),
        ];

        for read in &reads {
            let matches = align_helper.find_read_matches(read, 100);
            assert!(matches.is_none());
        }

        let matches = align_helper.find_read_matches(&reads[2], 15);
        assert_eq!(matches.unwrap().alignment.score, 16);

        align_helper.set_scoring(Scoring::from_scores(-4, -1, 2, -2).xclip(-5).yclip(0));
        let cigars = ["6=2X", "", "8=", "3=1X3=1X", "6=2X", "2X6="];
        let dists = [2, 0, 0, 2, 2, 2, 0];
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
        let cigars = ["6=2S", "", "8=", "3=1X3=1X", "6=2S", "2S6=", "1="];
        let dists = [0, 0, 0, 2, 0, 0, 0];
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
        let refs = ["CCCCACGTACGTGGGGGGTTTT"];
        let scoring = Scoring::from_scores(-4, -2, 2, -2).xclip(-5).yclip(0);
        let mut align_helper = AlignHelper::new(&refs, scoring, 3, 5);

        let quals = vec![30; 8];
        let reads = [
            Read::new(
                1,
                0,
                "1".to_string(),
                DnaString::from_dna_string("CCACGTCC"),
                quals.clone(),
            ),
            Read::new(
                2,
                0,
                "2".to_string(),
                DnaString::from_dna_string("ACGTACGT"),
                quals.clone(),
            ),
            Read::new(
                3,
                0,
                "3".to_string(),
                DnaString::from_dna_string("ACGAACGA"),
                quals.clone(),
            ),
            Read::new(
                4,
                0,
                "4".to_string(),
                DnaString::from_dna_string("ACGTTTTT"),
                quals.clone(),
            ),
            Read::new(
                5,
                0,
                "5".to_string(),
                DnaString::from_dna_string("ACGTTT"),
                quals,
            ),
        ];

        for read in &reads {
            let matches = align_helper.find_read_matches(read, 100);
            assert!(matches.is_none());
        }

        let matches = align_helper.find_read_matches(&reads[1], 15);
        assert_eq!(matches.unwrap().alignment.score, 16);

        align_helper.set_scoring(Scoring::from_scores(0, 0, 2, -2).xclip(-1).yclip(0));
        let matches = align_helper.find_read_matches(&reads[3], 15);
        assert_eq!(matches.unwrap().alignment.cigar(true), "4=6D4=");

        align_helper.set_scoring(Scoring::from_scores(-2, -1, 2, -2).xclip(-3).yclip(0));
        let cigars = ["6=1X1=", "8=", "3=1X3=1X", "4=6D4=", "4=2S"];
        let scores = [12, 16, 8, 8, 5];
        for i in 0..5 {
            let matches = align_helper.find_read_matches(&reads[i], 0);
            let unwrapped = matches.unwrap();
            assert_eq!(unwrapped.alignment.cigar(false), cigars[i]);
            assert_eq!(unwrapped.alignment.score, scores[i]);
        }
    }
}
