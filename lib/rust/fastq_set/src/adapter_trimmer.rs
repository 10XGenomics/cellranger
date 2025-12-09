//! Trim adapters from reads using a combination of
//! k-mer matches, sparse alignment and banded alignment.
//! Inspired by the [cutadapt](https://cutadapt.readthedocs.io/) tool.
//!
//! # Features/Limitations
//! * Supports regular, anchored and non-internal, 3' and 5' adapter trimming
//! * Linked adapters are not supported as of now
//! * Allowed error rate for the adapter is 10%
//!
//! # Algorithm
//!
//! # Example
//!
//! # TODO
//! A 5' adapter can be treated as a 3' adapter with
//! the read reversed. Perhaps it makes sense to use
//! this symmetry in the algorithm?
#![expect(missing_docs)]
use crate::WhichEnd;
use crate::read_pair::WhichRead;
use bio::alignment::pairwise::{self, MatchParams, Scoring};
use bio::alignment::sparse;
use bio::alignment::sparse::HashMapFx;
use serde::{Deserialize, Serialize};
use std::cmp::{Ordering, max, min};
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::BuildHasher;
use std::ops::Range;

type Aligner = pairwise::banded::Aligner<MatchParams>;

//------------------ CONSTANTS -----------------------//
// kmer length used for computing kmer matches
// between the read and the adapter.
const KMER_LEN: usize = 5;

// Banded aligner window length
const WINDOW_LEN: usize = 2;

// Expected read length to preallocate memory
// Reads longer than this will trigger reallocation
const EXPECTED_READ_LEN: usize = 150;

// Maximum allowed error rate in the adapter alignment
const ALLOWED_ERROR_RATE: f64 = 0.1_f64;

const MATCH_SCORE: i32 = 1;
const EDIT_SCORE: i32 = -2;

// Minimum score threshold for sparse alignment
const MIN_PATH_SCORE: i32 = -5;
//----------------END OF CONSTANTS--------------------//

/// An `enum` to specify the location of the adapter.
///
/// # 3' adapter
/// An adapter is considered a 3' adapter if you want to trim everything after the
/// position of the adapter. (i.e <Useful sequence>ADAPTER<Useless sequence>)
/// The following table summarizes the various read layouts for a 3' adapters
/// and whether they are recognized by the different `AdapterLoc` variants.
/// (Similar to [this table](https://cutadapt.readthedocs.io/en/stable/guide.html#id3)
/// in cutadapt docs)
///
/// | 3' Adapter location                        | Read Layout         | Found by `Anywhere` | Found by `NonInternal` | Found by `Anchored` |
/// |--------------------------------------------|---------------------|---------------------|------------------------|---------------------|
/// | Full adapter sequence internal to the read | acgtacgtADAPTERacgt | Yes                 | No                     | No                  |
/// | Partial adapter sequence at 3’ end         | acgtacgtacgtADAP    | Yes                 | Yes                    | No                  |
/// | Full adapter sequence at 3’ end            | acgtacgtacgtADAPTER | Yes                 | Yes                    | Yes                 |
///
/// # 5' adapter
/// An adapter is considered a 5' adapter if you want to trim everything before the
/// position of the adapter. (i.e <Useless sequence>ADAPTER<Useful sequence>)
/// The following table summarizes the various read layouts for a 5' adapters
/// and whether they are recognized by the different `AdapterLoc` variants.
/// (Similar to [this table](https://cutadapt.readthedocs.io/en/stable/guide.html#id4)
/// in cutadapt docs)
///
/// | 5' Adapter location                        | Read Layout         | Found by `Anywhere` | Found by `NonInternal` | Found by `Anchored` |
/// |--------------------------------------------|---------------------|---------------------|------------------------|---------------------|
/// | Full adapter sequence internal to the read | acgtacgtADAPTERacgt | Yes                 | No                     | No                  |
/// | Partial adapter sequence at 5’ end         | PTERacgtacgtacgt    | Yes                 | Yes                    | No                  |
/// | Full adapter sequence at 5’ end            | ADAPTERacgtacgtacgt | Yes                 | Yes                    | Yes                 |
///
/// # Serde
/// The variants are renamed into snake case by serde:
/// * `Anywhere` - "anywhere"
/// * `NonInternal` - "non_internal"
/// * `Anchored` - "anchored"
///
/// Hence, in order to deserialize from a JSON file/string for example, make sure
/// you use the snake case
/// ```rust
/// use fastq_set::adapter_trimmer::AdapterLoc;
/// assert_eq!(serde_json::from_str::<AdapterLoc>(r#""anywhere""#).unwrap(), AdapterLoc::Anywhere);
/// assert_eq!(serde_json::from_str::<AdapterLoc>(r#""non_internal""#).unwrap(), AdapterLoc::NonInternal);
/// assert_eq!(serde_json::from_str::<AdapterLoc>(r#""anchored""#).unwrap(), AdapterLoc::Anchored);
/// ```
#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq, Eq)]
pub enum AdapterLoc {
    /// Look for the full adapter sequence anywhere in the read.
    /// The adapter could be completely inside the read, or could be a prefix.
    #[serde(rename = "anywhere")]
    Anywhere,
    /// Disallow adapters that are internal to the read.
    #[serde(rename = "non_internal")]
    NonInternal,
    /// The read starts with the full adapter sequence (5' adapters) or
    /// the read ends with the full adapter sequence (3' adapters).
    #[serde(rename = "anchored")]
    Anchored,
}

/// A `struct` to store data associated with an adapter
///
/// An adapter is defined by:
/// * `name`: Name of the adapter
/// * `end`: whether it's a `FivePrime` adapter or a `ThreePrime` adapter
/// * `location`: Specify the location of the adapter (See [`AdapterLoc`](enum.AdapterLoc.html))
/// * `seq`: The sequence of the adapter as a `String`. One could use a `Vec<u8>` here, but
///   choose a String for the ease of auto (de)serialization from a json file.
///
/// # Example
/// The example below shows how you can create an `Adapter` from a JSON string
/// ```rust
/// use fastq_set::adapter_trimmer::{Adapter, AdapterLoc};
/// use fastq_set::WhichEnd;
/// let adapter_json_str = r#"{
///     "name": "spacer",
///     "end": "five_prime",
///     "location": "anchored",
///     "seq": "TTTCTTATATGGG"
/// }"#;
/// let adapter: Adapter = serde_json::from_str(adapter_json_str).unwrap();
/// assert_eq!(adapter.name, "spacer".to_string());
/// assert_eq!(adapter.end, WhichEnd::FivePrime);
/// assert_eq!(adapter.location, AdapterLoc::Anchored);
/// assert_eq!(adapter.seq, "TTTCTTATATGGG".to_string());
/// ```
///
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    pub name: String,
    pub end: WhichEnd,
    pub location: AdapterLoc,
    pub seq: String,
}

impl Adapter {
    /// Create a new `Adapter` object
    ///
    /// # Arguments
    /// * `name`: Name of the adaper. Can be any type which implements `ToString`
    /// * `end`: Specify whether the adapter is `FivePrime` or `ThreePrime`
    /// * `location`: Specify the adapter location
    /// * `seq`: Adapter sequence. Can be any type which implements `ToString`
    ///
    /// # Example
    /// ```rust
    /// use fastq_set::adapter_trimmer::{Adapter, AdapterLoc};
    /// use fastq_set::WhichEnd;
    /// let adapter = Adapter::new("custom_primer", WhichEnd::ThreePrime, AdapterLoc::Anywhere, "ACCGGTAACCGTTTAGC");
    /// ```
    pub fn new(
        name: impl ToString,
        end: WhichEnd,
        location: AdapterLoc,
        seq: impl ToString,
    ) -> Self {
        let name = name.to_string();
        let seq = seq.to_string();
        Adapter {
            name,
            end,
            location,
            seq,
        }
    }
}

/// This `struct` is used for searching and trimming adapters in read sequences.
///
/// Internally it stores a reference to the `Adapter`, a hashmap containing the
/// kmers of the adapter sequences, an aligner and scoring details.
///
/// # Example
/// The following example demonstates how to trim an adapter from a single read
///
/// ```rust
/// use fastq_set::adapter_trimmer::{Adapter, AdapterLoc, AdapterTrimmer};
/// use fastq_set::WhichEnd;
/// // Define the adapter
/// let adapter = Adapter::new("my_primer",
///     WhichEnd::ThreePrime,
///     AdapterLoc::Anywhere,
///     "TGCTTAAACTGACTACGTGTCGAAG");
/// // Create a trimmer
/// let mut trimmer = AdapterTrimmer::new(&adapter);
/// // The read has adapter at position 25 with a single base mismatch.
/// // Since the allowed error rate is 10%, upto 2 edits are okay
/// //                                    <--------Adapter-------->
/// let read = b"ACAATCCGGAACTGAGCGGAGTATTTGCTTAAACTGACTCCGTGTCGAAGGTAGGACAGC";
///
/// // Find the adapter
/// let result = trimmer.find(read);
/// assert!(result.is_some()); // Make sure we found the adapter
/// let result = result.unwrap();
/// // The adapter starts at position 25
/// assert_eq!(result.adapter_range.start, 25);
/// assert_eq!(result.adapter_range.len(), adapter.seq.len());
/// // We need to trim from base 25 to the end of the read
/// assert_eq!(result.trim_range, 25..read.len());
/// // We need to retain the first 25 bases in the read
/// assert_eq!(result.retain_range, 0..25);
/// ```
///
/// # Algorithm
/// [TODO]
pub struct AdapterTrimmer<'a> {
    pub adapter: &'a Adapter,
    seq: &'a [u8],
    kmer_hash: HashMapFx<&'a [u8], Vec<u32>>,
    aligner: Aligner,
    cut_scores: CutScores,
    homopolymer: bool,
}

impl<'a> AdapterTrimmer<'a> {
    /// Create a new `AdapterTrimmer` from an `Adapter`.
    pub fn new(adapter: &'a Adapter) -> Self {
        assert!(
            adapter.seq.len() >= KMER_LEN,
            "Adapter sequence cannot be shorter than {KMER_LEN}"
        );

        let cut_scores = CutScores::new(adapter.end, adapter.location);
        let scoring = Scoring::from_scores(0, EDIT_SCORE, MATCH_SCORE, EDIT_SCORE)
            .xclip_prefix(cut_scores.read_prefix)
            .xclip_suffix(cut_scores.read_suffix)
            .yclip_prefix(cut_scores.adapter_prefix)
            .yclip_suffix(cut_scores.adapter_suffix);

        let aligner = Aligner::with_capacity_and_scoring(
            EXPECTED_READ_LEN,
            adapter.seq.len(),
            scoring,
            KMER_LEN,
            WINDOW_LEN,
        );
        let seq = adapter.seq.as_bytes();
        let homopolymer = seq.iter().all(|&x| x == seq[0]);
        AdapterTrimmer {
            adapter,
            seq,
            kmer_hash: sparse::hash_kmers(adapter.seq.as_bytes(), KMER_LEN),
            aligner,
            cut_scores,
            homopolymer,
        }
    }

    /// Search for the adapter in the read.
    ///
    /// Arguments:
    /// * `read`: Byte array of the read sequence
    ///
    /// Ouput:
    /// * `Option<TrimResult>`: `None` if the adapter is not found,
    ///   otherwise `Some(TrimResult)` (See [`TrimResult`](struct.TrimResult.html))
    pub fn find(&mut self, read: &[u8]) -> Option<TrimResult> {
        use bio::alignment::AlignmentOperation::{Del, Ins, Match, Subst};

        if self.homopolymer {
            return self.find_homopolymer(read);
        }

        let matches = sparse::find_kmer_matches_seq2_hashed(read, &self.kmer_hash, KMER_LEN);
        let matches = sparse::expand_kmer_matches(read, self.seq, KMER_LEN, &matches, 1);
        let band_path = compute_path(
            &matches,
            read.len(),
            self.seq.len(),
            KMER_LEN,
            &self.cut_scores,
        );
        match band_path {
            Some(p) => {
                let score = p.total_score();
                if score > MIN_PATH_SCORE {
                    // Set the clipping score correctly if the
                    // adapter location is `Anywhere`
                    if self.adapter.location == AdapterLoc::Anywhere {
                        let scoring = self.aligner.get_mut_scoring();
                        match self.adapter.end {
                            WhichEnd::ThreePrime => {
                                let dx = read.len() - p.last_match.0 as usize;
                                let dy = self.adapter.seq.len() - p.last_match.1 as usize;
                                // if (max(dx, dy) - min(dx, dy)) < 3 {
                                match dx.cmp(&dy) {
                                    Ordering::Equal => {
                                        scoring.xclip_suffix = pairwise::MIN_SCORE;
                                        scoring.yclip_suffix = pairwise::MIN_SCORE;
                                    }
                                    Ordering::Greater => {
                                        scoring.xclip_suffix = 0;
                                        scoring.yclip_suffix = pairwise::MIN_SCORE;
                                    }
                                    Ordering::Less => {
                                        scoring.xclip_suffix = pairwise::MIN_SCORE;
                                        scoring.yclip_suffix = 0;
                                    }
                                }
                            }
                            WhichEnd::FivePrime => {
                                let dx = p.first_match.0;
                                let dy = p.first_match.1;
                                match dx.cmp(&dy) {
                                    Ordering::Equal => {
                                        scoring.xclip_prefix = pairwise::MIN_SCORE;
                                        scoring.yclip_prefix = pairwise::MIN_SCORE;
                                    }
                                    Ordering::Greater => {
                                        scoring.xclip_prefix = 0;
                                        scoring.yclip_prefix = pairwise::MIN_SCORE;
                                    }
                                    Ordering::Less => {
                                        scoring.xclip_prefix = pairwise::MIN_SCORE;
                                        scoring.yclip_prefix = 0;
                                    }
                                }
                            }
                        }
                    }
                    let alignment = self
                        .aligner
                        .custom_with_match_path(read, self.seq, &matches, &p.path);

                    let num_ops = alignment
                        .operations
                        .iter()
                        .filter(|&op| *op == Match || *op == Subst || *op == Ins || *op == Del)
                        .count();
                    let max_score = min(alignment.ylen, num_ops);
                    let num_errs = alignment
                        .operations
                        .iter()
                        .filter(|&op| *op == Subst || *op == Ins || *op == Del)
                        .count();
                    let err_rate = (num_errs as f64) / (max_score as f64);
                    if err_rate <= ALLOWED_ERROR_RATE {
                        Some(TrimResult {
                            adapter_range: alignment.xstart..alignment.xend,
                            trim_range: match self.adapter.end {
                                WhichEnd::ThreePrime => alignment.xstart..alignment.xlen,
                                WhichEnd::FivePrime => 0..alignment.xend,
                            },
                            retain_range: match self.adapter.end {
                                WhichEnd::ThreePrime => 0..alignment.xstart,
                                WhichEnd::FivePrime => alignment.xend..alignment.xlen,
                            },
                            score: alignment.score,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            None => None,
        }
    }

    fn find_homopolymer(&mut self, read: &[u8]) -> Option<TrimResult> {
        assert!(self.homopolymer);
        let base = self.seq[0];
        let cumulative_matches: Vec<_> = vec![0]
            .into_iter()
            .chain(read.iter().map(|&x| (x == base) as usize))
            .scan(0, |state, x| {
                *state += x;
                Some(*state)
            })
            .collect();

        let get_counts_if_valid = |i: usize, j: usize| {
            let num_bases = j.saturating_sub(i);
            if num_bases < KMER_LEN {
                return None;
            }
            let num_matches = cumulative_matches[j] - cumulative_matches[i];
            let num_errs = num_bases - num_matches;
            let err_rate = (num_errs as f64) / (num_bases as f64);
            if err_rate <= ALLOWED_ERROR_RATE {
                Some((num_matches, num_errs))
            } else {
                None
            }
        };
        let mut best_num_matches = 0;
        let mut best_err = self.seq.len();
        let mut best_pos = None;
        match self.adapter.end {
            WhichEnd::ThreePrime => {
                let i_range = match self.adapter.location {
                    AdapterLoc::Anywhere => 0..read.len(),
                    AdapterLoc::NonInternal => {
                        read.len().saturating_sub(self.seq.len())..read.len()
                    }
                    AdapterLoc::Anchored => {
                        read.len().saturating_sub(self.seq.len())
                            ..read.len().saturating_sub(self.seq.len() - 1)
                    }
                };
                for i in i_range {
                    let j = min(read.len(), i + self.seq.len());
                    if let Some((num_matches, num_errs)) = get_counts_if_valid(i, j)
                        && ((num_matches > best_num_matches)
                            || (num_matches == best_num_matches && num_errs < best_err))
                    {
                        best_num_matches = num_matches;
                        best_err = num_errs;
                        best_pos = Some(i);
                    }
                }
                match best_pos {
                    Some(p) => Some(TrimResult {
                        adapter_range: p..min(p + self.seq.len(), read.len()),
                        trim_range: p..read.len(),
                        retain_range: 0..p,
                        score: MATCH_SCORE * best_num_matches as i32 + EDIT_SCORE * best_err as i32,
                    }),
                    None => None,
                }
            }
            WhichEnd::FivePrime => {
                let j_range = match self.adapter.location {
                    AdapterLoc::Anywhere => 0..read.len(),
                    AdapterLoc::NonInternal => 0..self.seq.len() + 1,
                    AdapterLoc::Anchored => self.seq.len()..self.seq.len() + 1,
                };
                for j in j_range {
                    let i = j.saturating_sub(self.seq.len());
                    if let Some((num_matches, num_errs)) = get_counts_if_valid(i, j)
                        && ((num_matches > best_num_matches)
                            || (num_matches == best_num_matches && num_errs < best_err))
                    {
                        best_num_matches = num_matches;
                        best_err = num_errs;
                        best_pos = Some(j);
                    }
                }
                match best_pos {
                    Some(p) => Some(TrimResult {
                        adapter_range: p.saturating_sub(self.seq.len())..p,
                        trim_range: 0..p,
                        retain_range: p..read.len(),
                        score: MATCH_SCORE * best_num_matches as i32 + EDIT_SCORE * best_err as i32,
                    }),
                    None => None,
                }
            }
        }
    }
}

/// The result of an adapter search, containing ranges to
/// locate the position of the adapter within the read as well
/// as ranges in the read that need to be trimmed and retained.
///
/// The following figure illustrates the three ranges for a 5' and
/// 3' adapter whose location can be `Anywhere`.
/// ![TrimResult](../../../../doc-media/fastq/adapter_trimming_result.png)
///
/// In the case of `NonInternal` or `Anchored` adapters, the `trim_range`
/// would be identical to the `adapter_range`.
#[derive(Debug)]
pub struct TrimResult {
    /// `Range` in the read sequence where the adapter is located
    pub adapter_range: Range<usize>,
    /// `Range` in the read sequence which needs to be trimmed
    pub trim_range: Range<usize>,
    /// `Range` in the read sequence after trimming the adapter
    pub retain_range: Range<usize>,
    /// The alignment score
    pub score: i32,
}

/// A function to compute the intersection of two `Range<usize>`. I would
/// have liked this to be part of the std library. Takes two ranges as input
/// and computes the intersection. This is a useful function when you want to
/// combine `TrimResult`s from multiple adapters to determine what range to
/// actually keep. The `Range` returned is guaranteed to have start<=end even
/// when the ranges do not overlap
///
/// # Test
/// * test_intersect_ranges() - Tests this function for a small number of inputs
pub fn intersect_ranges(a: &Range<usize>, b: &Range<usize>) -> Range<usize> {
    let end = min(a.end, b.end);
    let start = min(max(a.start, b.start), end);
    Range { start, end }
}

#[derive(Debug, Copy, Clone)]
struct CutScores {
    match_score: i32,
    edit_score: i32,
    read_prefix: i32,
    read_suffix: i32,
    adapter_prefix: i32,
    adapter_suffix: i32,
    right_cut: i32,
    left_cut: i32,
}

impl CutScores {
    fn new(end: WhichEnd, loc: AdapterLoc) -> Self {
        use self::AdapterLoc::{Anchored, Anywhere, NonInternal};
        use self::WhichEnd::{FivePrime, ThreePrime};
        let mut scores = CutScores {
            match_score: MATCH_SCORE,
            edit_score: EDIT_SCORE,
            read_prefix: pairwise::MIN_SCORE,
            read_suffix: pairwise::MIN_SCORE,
            adapter_prefix: pairwise::MIN_SCORE,
            adapter_suffix: pairwise::MIN_SCORE,
            right_cut: pairwise::MIN_SCORE,
            left_cut: pairwise::MIN_SCORE,
        };

        match (end, loc) {
            (ThreePrime, Anywhere) => {
                scores.read_prefix = 0;
                scores.right_cut = 0;
            }
            (ThreePrime, NonInternal) => {
                scores.read_prefix = 0;
                scores.adapter_suffix = 0;
            }
            (ThreePrime, Anchored) => {
                scores.read_prefix = 0;
            }
            (FivePrime, Anywhere) => {
                scores.read_suffix = 0;
                scores.left_cut = 0;
            }
            (FivePrime, NonInternal) => {
                scores.read_suffix = 0;
                scores.adapter_prefix = 0;
            }
            (FivePrime, Anchored) => {
                scores.read_suffix = 0;
            }
        }

        scores
    }
}

#[derive(Debug, Clone)]
struct PathTracker {
    parent: Option<usize>,
    index: usize,
    last_match: (u32, u32),
    left_score: i32,
    right_score: i32,
    internal_score: i32,
}
impl PathTracker {
    fn total_score(&self) -> i32 {
        self.left_score + self.internal_score + self.right_score
    }
}

#[derive(Debug, Clone)]
struct PathData {
    path: Vec<usize>,
    first_match: (u32, u32),
    last_match: (u32, u32),
    left_score: i32,
    right_score: i32,
    internal_score: i32,
}
impl PathData {
    fn total_score(&self) -> i32 {
        self.left_score + self.internal_score + self.right_score
    }
}

fn compute_path(
    matches: &[(u32, u32)],
    read_len: usize,
    adapter_len: usize,
    k: usize,
    cut_scores: &CutScores,
) -> Option<PathData> {
    if matches.is_empty() {
        return None;
    }

    let max_allowed_indel = (adapter_len as f64 * ALLOWED_ERROR_RATE) as u32;

    // incoming matches must be sorted.
    for i in 1..matches.len() {
        assert!(matches[i - 1] < matches[i]);
    }
    let mut best_path_ending_here: Vec<PathTracker> = Vec::with_capacity(matches.len());

    let l_score = |this_match: (u32, u32)| {
        let read_score = max(
            (this_match.0 as i32) * cut_scores.edit_score,
            cut_scores.read_prefix,
        );
        let adapter_score = max(
            (this_match.1 as i32) * cut_scores.edit_score,
            cut_scores.adapter_prefix,
        );
        max(
            read_score + adapter_score,
            cut_scores.left_cut + max(read_score, adapter_score),
        )
    };

    let r_score = |this_match: (u32, u32)| {
        let read_edit_score =
            ((read_len - k - (this_match.0 as usize)) as i32) * cut_scores.edit_score;
        let adapter_edit_score =
            ((adapter_len - k - (this_match.1 as usize)) as i32) * cut_scores.edit_score;
        let read_score = max(read_edit_score, cut_scores.read_suffix);
        let adapter_score = max(adapter_edit_score, cut_scores.adapter_suffix);
        max(
            read_score + adapter_score,
            cut_scores.right_cut + max(read_score, adapter_score),
        )
    };

    {
        let mut max_l_score = pairwise::MIN_SCORE;
        let mut max_r_score = pairwise::MIN_SCORE;
        for &this_match in matches {
            max_l_score = max(max_l_score, l_score(this_match));
            max_r_score = max(max_r_score, r_score(this_match));
        }
        if (max_l_score + max_r_score) <= (-(adapter_len as i32) * cut_scores.match_score) {
            return None;
        }
    }

    for (i, &this_match) in matches.iter().enumerate() {
        let left_score = l_score(this_match);
        let right_score = r_score(this_match);

        let mut best_path = PathTracker {
            parent: None,
            index: i,
            last_match: this_match,
            left_score,
            right_score,
            internal_score: (k as i32) * cut_scores.match_score,
        };

        for (j, p) in best_path_ending_here.iter().enumerate().take(i) {
            if (p.last_match.0 + 1 == this_match.0) && (p.last_match.1 + 1 == this_match.1) {
                let score = p.left_score + p.internal_score + cut_scores.match_score + right_score;
                if score > best_path.total_score() {
                    best_path.parent = Some(j);
                    best_path.left_score = p.left_score;
                    best_path.internal_score += cut_scores.match_score;
                }
            } else if ((p.last_match.0 + k as u32) <= this_match.0)
                && ((p.last_match.1 + k as u32) <= this_match.1)
            {
                let dx = this_match.0 - p.last_match.0 - k as u32;
                let dy = this_match.1 - p.last_match.1 - k as u32;
                let gap = max(dx, dy);
                let min_indel = gap - min(dx, dy);
                let internal_score = p.internal_score
                    + (gap as i32) * cut_scores.edit_score
                    + (k as i32) * cut_scores.match_score;
                let score = p.left_score + internal_score + right_score;
                if score > best_path.total_score() && min_indel <= max_allowed_indel {
                    best_path.parent = Some(j);
                    best_path.left_score = p.left_score;
                    best_path.internal_score = internal_score;
                }
            }
        }
        best_path_ending_here.push(best_path);
    }

    let mut best_end = best_path_ending_here
        .iter()
        .max_by_key(|x| (x.total_score(), -(x.last_match.0 as i32)))
        .unwrap();
    let best_ind = best_end.index;

    let mut path = Vec::with_capacity(matches.len());
    path.push(best_end.index);
    while let Some(p) = best_end.parent {
        path.push(p);
        best_end = &best_path_ending_here[p];
    }
    path.reverse();

    let first_match = matches[path[0]];
    Some(PathData {
        path,
        first_match,
        last_match: best_path_ending_here[best_ind].last_match,
        left_score: best_path_ending_here[best_ind].left_score,
        right_score: best_path_ending_here[best_ind].right_score,
        internal_score: best_path_ending_here[best_ind].internal_score,
    })
}

#[derive(Default)]
pub struct ReadAdapterCatalog<'a> {
    adapter_names: HashSet<&'a str>,
    read1_trimmers: Vec<AdapterTrimmer<'a>>,
    read2_trimmers: Vec<AdapterTrimmer<'a>>,
}

impl<'a> ReadAdapterCatalog<'a> {
    pub fn new() -> Self {
        ReadAdapterCatalog::default()
    }

    pub fn push_trimmer(&mut self, read: WhichRead, trimmer: AdapterTrimmer<'a>) {
        assert!(
            self.adapter_names.insert(&trimmer.adapter.name),
            "Adapter names need to be unique in ReadAdapterCatalog. '{}' already exists",
            trimmer.adapter.name
        );

        match read {
            WhichRead::R1 => self.read1_trimmers.push(trimmer),
            WhichRead::R2 => self.read2_trimmers.push(trimmer),
            _ => panic!("ReadAdapterCatalog push() can only accept R1/R2. Found {read:?}"),
        }
    }

    pub fn add_adapter(&mut self, read: WhichRead, adapter: &'a Adapter) {
        self.push_trimmer(read, AdapterTrimmer::new(adapter));
    }

    pub fn get_mut_trimmers(&mut self, read: WhichRead) -> &mut Vec<AdapterTrimmer<'a>> {
        match read {
            WhichRead::R1 => &mut self.read1_trimmers,
            WhichRead::R2 => &mut self.read2_trimmers,
            _ => {
                panic!("ReadAdapterCatalog get_mut_trimmer() can only accept R1/R2. Found {read:?}")
            }
        }
    }
}

impl<'a, S> From<&'a HashMap<WhichRead, Vec<Adapter>, S>> for ReadAdapterCatalog<'a>
where
    S: BuildHasher,
{
    fn from(adapter_map: &'a HashMap<WhichRead, Vec<Adapter>, S>) -> Self {
        let mut catalog = ReadAdapterCatalog::new();
        for (read, adapters) in adapter_map {
            for adapter in adapters {
                catalog.add_adapter(*read, adapter);
            }
        }
        catalog
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_anywhere_3p() {
        test_helper(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            WhichEnd::ThreePrime,
            AdapterLoc::Anywhere,
            "tests/input.fa",
            "tests/cutadapt_v1_18_anywhere_3p.fa",
        );
    }

    #[test]
    fn test_anywhere_3p_poly_a() {
        test_helper(
            "AAAAAAAAAAAAAAAAAAAA",
            WhichEnd::ThreePrime,
            AdapterLoc::Anywhere,
            "tests/input_polyA.fa",
            "tests/cutadapt_v1_18_anywhere_3p_polyA.fa",
        );
    }

    #[test]
    fn test_anywhere_5p() {
        test_helper(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            WhichEnd::FivePrime,
            AdapterLoc::Anywhere,
            "tests/input.fa",
            "tests/cutadapt_v1_18_anywhere_5p.fa",
        );
    }

    #[test]
    fn test_anywhere_5p_poly_a() {
        test_helper(
            "AAAAAAAAAAAAAAAAAAAA",
            WhichEnd::FivePrime,
            AdapterLoc::Anywhere,
            "tests/input_polyA.fa",
            "tests/cutadapt_v1_18_anywhere_5p_polyA.fa",
        );
    }

    #[test]
    fn test_anchored_3p() {
        test_helper(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            WhichEnd::ThreePrime,
            AdapterLoc::Anchored,
            "tests/input.fa",
            "tests/cutadapt_v1_18_anchored_3p.fa",
        );
    }

    #[test]
    fn test_anchored_3p_poly_a() {
        test_helper(
            "AAAAAAAAAAAAAAAAAAAA",
            WhichEnd::ThreePrime,
            AdapterLoc::Anchored,
            "tests/input_polyA.fa",
            "tests/cutadapt_v1_18_anchored_3p_polyA.fa",
        );
    }

    #[test]
    fn test_anchored_5p() {
        test_helper(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            WhichEnd::FivePrime,
            AdapterLoc::Anchored,
            "tests/input.fa",
            "tests/cutadapt_v1_18_anchored_5p.fa",
        );
    }

    #[test]
    fn test_anchored_5p_poly_a() {
        test_helper(
            "AAAAAAAAAAAAAAAAAAAA",
            WhichEnd::FivePrime,
            AdapterLoc::Anchored,
            "tests/input_polyA.fa",
            "tests/cutadapt_v1_18_anchored_5p_polyA.fa",
        );
    }

    #[test]
    fn test_non_internal_3p() {
        test_helper(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            WhichEnd::ThreePrime,
            AdapterLoc::NonInternal,
            "tests/input.fa",
            "tests/cutadapt_v1_18_non_internal_3p.fa",
        );
    }

    #[test]
    fn test_non_internal_3p_poly_a() {
        test_helper(
            "AAAAAAAAAAAAAAAAAAAA",
            WhichEnd::ThreePrime,
            AdapterLoc::NonInternal,
            "tests/input_polyA.fa",
            "tests/cutadapt_v1_18_non_internal_3p_polyA.fa",
        );
    }

    #[test]
    fn test_non_internal_5p() {
        test_helper(
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            WhichEnd::FivePrime,
            AdapterLoc::NonInternal,
            "tests/input.fa",
            "tests/cutadapt_v1_18_non_internal_5p.fa",
        );
    }

    #[test]
    fn test_non_internal_5p_poly_a() {
        test_helper(
            "AAAAAAAAAAAAAAAAAAAA",
            WhichEnd::FivePrime,
            AdapterLoc::NonInternal,
            "tests/input_polyA.fa",
            "tests/cutadapt_v1_18_non_internal_5p_polyA.fa",
        );
    }

    fn test_helper(
        adapter_seq: impl ToString,
        end: WhichEnd,
        loc: AdapterLoc,
        input_path: impl AsRef<Path> + Debug,
        cutadapt_output_path: impl AsRef<Path> + Debug,
    ) {
        use bio::io::fasta::Reader;

        let adapter = Adapter::new("primer", end, loc, adapter_seq);
        let mut trimmer = AdapterTrimmer::new(&adapter);

        let input = Reader::from_file(input_path).unwrap();
        let cutadapt_output = Reader::from_file(cutadapt_output_path).unwrap();

        let mut total = 0;
        let mut correct_untrimmed = 0;
        let mut self_trimmed_only = 0;
        let mut self_untrimmed_only = 0;
        let mut correct_trimmed = 0;
        let mut inconsistent_trimmed = 0;

        for (r_i, r_o) in input.records().zip(cutadapt_output.records()) {
            let r_i = r_i.unwrap();
            let r_o = r_o.unwrap();

            total += 1;

            let seq = r_i.seq();
            let seq_len = seq.len();
            let expected_trimmed_len = r_o.seq().len();

            let cutadapt_result = trimmer.find(seq);
            let trimmed_len = match cutadapt_result {
                Some(r) => r.retain_range.len(),
                None => seq.len(),
            };

            if seq_len == expected_trimmed_len && trimmed_len == seq_len {
                correct_untrimmed += 1;
            } else if seq_len == expected_trimmed_len {
                self_trimmed_only += 1;
            } else if trimmed_len == seq_len {
                self_untrimmed_only += 1;
            } else if (max(expected_trimmed_len, trimmed_len)
                - min(expected_trimmed_len, trimmed_len))
                <= 2
            {
                correct_trimmed += 1;
            } else {
                inconsistent_trimmed += 1;
            }
        }

        // Sensitivity = (Self & Ref) / Ref
        let sensitivity = (correct_trimmed as f64 + inconsistent_trimmed as f64)
            / (correct_trimmed as f64 + inconsistent_trimmed as f64 + self_untrimmed_only as f64);
        // PPV = (Self & Ref)/Self
        let ppv = (correct_trimmed as f64 + inconsistent_trimmed as f64)
            / (correct_trimmed as f64 + inconsistent_trimmed as f64 + self_trimmed_only as f64);
        // Concordance = (Self & Ref & Self==Ref) / (Self & Ref)
        let concordance =
            (correct_trimmed as f64) / (correct_trimmed as f64 + inconsistent_trimmed as f64);

        println!("{adapter:?}");
        println!("Reads not trimmed {correct_untrimmed}/{total}");
        println!("Sensitivity {:.2}%", 100.0 * sensitivity);
        println!("PPV         {:.2}%", 100.0 * ppv);
        println!("Concordance {:.2}%", 100.0 * concordance);

        assert!(sensitivity > 0.99_f64);
        assert!(ppv > 0.99_f64);
        assert!(concordance > 0.99_f64);
    }

    #[test]
    fn test_multiple_occurence_3p() {
        let adapter = Adapter::new(
            "primer",
            WhichEnd::ThreePrime,
            AdapterLoc::Anywhere,
            "AGATCGGAAGAGCACAC",
        );
        // (25 bases)(adapter)(25 bases)(adapter)(10 bases)
        let read = b"ACAAATAGGTCAGTCCAGTAGATTCAGATCGGAAGAGCACACGTGGTGGACTGTAGCGGACAAGAAGAGATCGGAAGAGCACACATTTTCCCCA";
        let mut trimmer = AdapterTrimmer::new(&adapter);
        let result = trimmer.find(read).unwrap();
        assert_eq!(result.retain_range, 0..25);
        assert_eq!(result.adapter_range, 25..42);
        assert_eq!(result.trim_range, 25..94);
    }

    #[test]
    fn test_multiple_occurence_5p() {
        let adapter = Adapter::new(
            "primer",
            WhichEnd::FivePrime,
            AdapterLoc::Anywhere,
            "AGATCGGAAGAGCACAC",
        );
        // (25 bases)(adapter)(25 bases)(adapter)(10 bases)
        let read = b"ACAAATAGGTCAGTCCAGTAGATTCAGATCGGAAGAGCACACGTGGTGGACTGTAGCGGACAAGAAGAGATCGGAAGAGCACACATTTTCCCCA";
        let mut trimmer = AdapterTrimmer::new(&adapter);
        let result = trimmer.find(read).unwrap();
        assert_eq!(result.retain_range, 42..94);
        assert_eq!(result.adapter_range, 25..42);
        assert_eq!(result.trim_range, 0..42);
    }

    #[test]
    fn test_multiple_occurence_3p_non_internal() {
        let adapter = Adapter::new(
            "primer",
            WhichEnd::ThreePrime,
            AdapterLoc::NonInternal,
            "AGATCGGAAGAGCACAC",
        );
        // (25 bases)(adapter)(25 bases)(partial adapter)
        let read =
            b"ACAAATAGGTCAGTCCAGTAGATTCAGATCGGAAGAGCACACGTGGTGGACTGTAGCGGACAAGAAGAGATCGGAAGA";
        let mut trimmer = AdapterTrimmer::new(&adapter);
        let result = trimmer.find(read).unwrap();
        assert_eq!(result.retain_range, 0..67);
        assert_eq!(result.adapter_range, 67..78);
        assert_eq!(result.trim_range, 67..78);
    }

    #[test]
    #[should_panic]
    fn test_empty_adapter_seq() {
        let adapter = Adapter::new("p", WhichEnd::ThreePrime, AdapterLoc::NonInternal, "");
        let _ = AdapterTrimmer::new(&adapter);
    }

    #[test]
    #[should_panic]
    fn test_too_short_adapter_seq() {
        let adapter = Adapter::new("p", WhichEnd::ThreePrime, AdapterLoc::NonInternal, "ACG");
        let _ = AdapterTrimmer::new(&adapter);
    }

    #[test]
    fn test_homopolymer() {
        let adapter = Adapter::new(
            "p",
            WhichEnd::ThreePrime,
            AdapterLoc::NonInternal,
            "AAAAAAAA",
        );
        let trimmer = AdapterTrimmer::new(&adapter);
        assert!(trimmer.homopolymer);

        let adapter = Adapter::new(
            "p",
            WhichEnd::ThreePrime,
            AdapterLoc::NonInternal,
            "GCATCCAATC",
        );
        let trimmer = AdapterTrimmer::new(&adapter);
        assert!(!trimmer.homopolymer);
    }

    #[test]
    #[should_panic]
    fn test_read_adapter_catalog_repeat() {
        let adapter_1 = Adapter::new(
            "polya",
            WhichEnd::FivePrime,
            AdapterLoc::Anywhere,
            "AAAAAAAAAA",
        );
        let adapter_2 = Adapter::new(
            "polya",
            WhichEnd::FivePrime,
            AdapterLoc::Anywhere,
            "AAAAAAAAAA",
        );
        let mut catalog = ReadAdapterCatalog::new();
        catalog.push_trimmer(WhichRead::R1, AdapterTrimmer::new(&adapter_1));
        catalog.push_trimmer(WhichRead::R2, AdapterTrimmer::new(&adapter_2));
    }

    #[test]
    fn test_intersect_ranges() {
        assert_eq!(intersect_ranges(&(0..10), &(5..15)), 5..10);
        assert_eq!(intersect_ranges(&(5..15), &(0..10)), 5..10);
        assert_eq!(intersect_ranges(&(0..15), &(5..10)), 5..10);
        assert_eq!(intersect_ranges(&(0..20), &(50..100)), 20..20);
    }
}
