//!
//! ALign reads to a contig

use bio::alignment::pairwise::{self, banded, MatchParams};
use bio::alignment::sparse::{find_kmer_matches_seq2_hashed, hash_kmers, HashMapFx};
use bio::alignment::Alignment;
pub type Scoring = pairwise::Scoring<MatchParams>;

pub struct ContigAligner<'a> {
    contig: &'a [u8], // Ascii encoded
    kmers_hash: HashMapFx<&'a [u8], Vec<u32>>,
    scoring: Scoring,
    k: usize,
    w: usize,
}

impl<'a> ContigAligner<'a> {
    /// Creates a new ContigAligner instance given the contig sequence,
    /// the scoring parameters, and the banded align parameters (kmer length and window size)
    ///
    /// # Arguments
    ///
    /// * `contig` - Reference sequences (stored as a byte vector)
    /// * `scoring` - Scoring struct
    /// * `k` - kmer length for constructing the band (see bio::alignment::pairwise::banded)
    /// * `w` - window size for constructing the band (see bio::alignment::pairwise::banded)
    ///
    pub fn new(contig: &'a [u8], scoring: Scoring, k: usize, w: usize) -> Self {
        // Make sure that the scoring implies that the alignment is local in the reference y
        assert!(scoring.yclip_prefix == 0);
        assert!(scoring.yclip_suffix == 0);

        let kmers_hash = hash_kmers(contig, k);

        ContigAligner {
            contig,
            kmers_hash,
            scoring,
            k,
            w,
        }
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
    /// * `read_seq` - Read sequence.
    /// * `min_align_score` - Minimum SW score for accepting an alignment
    ///
    pub fn align_read(&self, read_seq: &[u8], min_align_score: i32) -> Option<Alignment> {
        let kmer_matches = find_kmer_matches_seq2_hashed(read_seq, &self.kmers_hash, self.k);
        if kmer_matches.is_empty() {
            return None;
        }

        let mut aligner = banded::Aligner::with_scoring(self.scoring, self.k, self.w);
        let alignment = aligner.custom_with_expanded_matches(
            read_seq,
            self.contig,
            kmer_matches,
            Some(1),
            true,
        );

        if alignment.score > min_align_score {
            Some(alignment)
        } else {
            None
        }
    }
}
