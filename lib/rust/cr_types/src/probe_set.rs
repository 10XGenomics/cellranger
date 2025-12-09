#![expect(missing_docs)]

use crate::reference::probe_set_reference::{
    PROBE_SET_FILE_FORMAT, PROBE_SET_HEADER, TargetSetFile,
};
use crate::reference::reference_info::{MULTI_GENOME_SEPARATOR, ReferenceInfo};
use crate::types::{AlignerParam, GenomeName};
use anyhow::{Context, Result, anyhow, bail, ensure};
use bio::alignment::distance::simd::{bounded_levenshtein, levenshtein};
use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::{Alignment, AlignmentOperation};
use itertools::{Itertools, chain, zip_eq};
use lazy_static::lazy_static;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::{Ord, Ordering, PartialOrd, min};
use std::collections::HashMap;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, Write};
use std::iter::zip;
use std::ops::Deref;
use std::path::Path;
use std::str::FromStr;
use std::{fmt, iter};
use strum::{AsRefStr, Display, EnumString};
use tempfile::NamedTempFile;
use transcriptome::Gene;

/// One read half maps to a probe and the other half does not.
pub const MAPQ_HALF_MAPPED: u8 = 1;

/// Each read half maps to a different probe.
pub const MAPQ_SPLIT_MAPPED: u8 = 3;

/// Both LHS and RHS map, but the gap sequence is not within the max error of expected gap
pub const MAPQ_GAP_MAPPED_NOT_WITHIN_MAX_ERR: u8 = 5;

/// Maximum rate of error in the gap sequence. The value is chosen
/// based on the expectations that the errors in gap should be really small.
/// Will revisit this once we have more data.
pub const MAX_GAP_ERROR_RATE: f64 = 0.25;
/// Minmum number of basepairs difference allowed for gap to still be aligned
/// to the expected gap and alignment metrics to be reported.
pub const MIN_ALLOWED_GAP_ERROR: u32 = 2;

/// Return the intersection of two sorted iterators.
fn intersect_sorted_lists<I, T>(xs: I, ys: I) -> impl Iterator<Item = T>
where
    I: IntoIterator<Item = T>,
    T: Ord,
{
    let mut xs = xs.into_iter();
    let mut ys = ys.into_iter();
    iter::from_fn(move || {
        let mut opt_x = xs.next();
        let mut opt_y = ys.next();
        while let (Some(x), Some(y)) = (opt_x, opt_y) {
            match x.cmp(&y) {
                Ordering::Less => {
                    opt_x = xs.next();
                    opt_y = Some(y);
                    assert!(
                        opt_x.is_none() || &x < opt_x.as_ref().unwrap(),
                        "xs is not sorted"
                    );
                }
                Ordering::Greater => {
                    opt_x = Some(x);
                    opt_y = ys.next();
                    assert!(
                        opt_y.is_none() || &y < opt_y.as_ref().unwrap(),
                        "ys is not sorted"
                    );
                }
                Ordering::Equal => {
                    return Some(x);
                }
            }
        }
        None
    })
}

/// Deprecated probes are excluded from the filtered_feature_bc_matrix.h5.
const DEPRECATED_PROBE_ID_PREFIX: &str = "DEPRECATED_";

/// Return true if this probe ID is excluded based on its probe ID.
pub fn is_deprecated_probe(probe_id: &str) -> bool {
    probe_id.starts_with(DEPRECATED_PROBE_ID_PREFIX)
}

/// Strip deprecated prefix from the probe or gene ID.
pub(super) fn strip_deprecated_prefix(id: &str) -> &str {
    id.strip_prefix(DEPRECATED_PROBE_ID_PREFIX).unwrap_or(id)
}

// This was added in to help look for gDNA
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
pub enum ProbeRegion {
    // Overlaps a splice site on the gene, likely to only be available in processed mRNA
    Spliced,
    // Does not overlap splice sites, presumed to be available in genome/pre+post mRNA.
    Unspliced,
    // Undefined.
    Other,
}

impl ProbeRegion {
    pub fn new(kind: &str) -> ProbeRegion {
        match kind {
            "unspliced" => ProbeRegion::Unspliced,
            "spliced" => ProbeRegion::Spliced,
            _ => ProbeRegion::Other,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match *self {
            ProbeRegion::Unspliced => "unspliced",
            ProbeRegion::Spliced => "spliced",
            ProbeRegion::Other => "other",
        }
    }
}

impl fmt::Display for ProbeRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

// do we even need the probe type?? discard it before merge if it ends up unnecessary
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Deserialize,
    Serialize,
    Copy,
    EnumString,
    AsRefStr,
    Default,
    Display,
)]
pub enum ProbeType {
    #[default]
    RTL,
    PairedGapAlign,
    UnpairedGapAlign,
}

/// A probe index and its gene.
#[derive(Clone, Debug, Eq, Deserialize, Serialize)]
pub struct Probe {
    /// The probe ID.
    pub probe_id: String,

    /// The probe sequencd.
    pub probe_seq: String,

    /// The probe gene ID and name.
    pub gene: Gene,

    /// An optional boolean to indicate if the probe is included
    pub included: bool,

    /// The region of the probe
    pub region: Option<ProbeRegion>,

    // The type of the probe
    pub probe_type: ProbeType,

    // reference sequence (chromosome) name where probe is expected to map
    pub ref_sequence_name: String,

    // reference sequence (chromosome) position where probe is expected to map
    pub ref_sequence_pos: Option<usize>,

    // CIGAR string describing expected alignment
    pub cigar_string: String,

    // Reference genome name
    pub genome: GenomeName,
}

impl Probe {
    /// Return whether this probe ID is excluded based on its probe ID or `included` is false.
    pub fn is_excluded_probe(&self) -> bool {
        !self.included || is_deprecated_probe(&self.probe_id)
    }

    /// Return a record for csv.Writer::write_record.
    pub fn as_record(&self) -> [&str; 6] {
        [
            &self.gene.id,
            &self.probe_seq,
            &self.probe_id,
            if self.included { "TRUE" } else { "FALSE" },
            self.region.as_ref().map_or("", ProbeRegion::as_str),
            &self.gene.name,
        ]
    }
}

impl PartialEq for Probe {
    /// Compare two probe IDs for equality.
    fn eq(&self, other: &Self) -> bool {
        self.probe_id == other.probe_id
    }
}

impl PartialOrd for Probe {
    /// Compare two probe IDs.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Probe {
    /// Compare two probe IDs.
    fn cmp(&self, other: &Self) -> Ordering {
        self.probe_id.cmp(&other.probe_id)
    }
}

impl Hash for Probe {
    /// Return the hash value of this probe.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.probe_id.hash(state);
    }
}

impl fmt::Display for Probe {
    // Output the probe ID.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", &self.probe_id)
    }
}

/// A pair of left and right half-probe sequences.
pub struct ProbeSequence {
    /// LHS half-probe sequence.
    pub lhs: Vec<u8>,

    /// RHS half-probe sequence.
    pub rhs: Vec<u8>,

    /// Expected sequence of the optional gap.
    gap: Vec<u8>,
}

impl FromStr for ProbeSequence {
    type Err = anyhow::Error;

    fn from_str(probe_seq_str: &str) -> Result<Self> {
        let probe_seq = probe_seq_str.as_bytes();
        Ok(match probe_seq.iter().filter(|&&x| x == b'-').count() {
            0 => {
                // Skip the middle base if the sequence has odd length.
                let lhs_end = probe_seq.len() / 2;
                let rhs_start = probe_seq.len().div_ceil(2);
                ProbeSequence {
                    lhs: probe_seq[..lhs_end].to_vec(),
                    rhs: probe_seq[rhs_start..].to_vec(),
                    gap: Vec::new(),
                }
            }
            1 => {
                let (lhs, rhs) = probe_seq.split(|&x| x == b'-').collect_tuple().unwrap();
                ProbeSequence {
                    lhs: lhs.to_vec(),
                    rhs: rhs.to_vec(),
                    gap: Vec::new(),
                }
            }
            2 => {
                let (lhs, gap, rhs) = probe_seq.split(|&x| x == b'-').collect_tuple().unwrap();
                ProbeSequence {
                    lhs: lhs.to_vec(),
                    rhs: rhs.to_vec(),
                    gap: gap.to_vec(),
                }
            }
            3.. => bail!("Too many hyphens in probe sequence: {probe_seq_str}"),
        })
    }
}

/// A map of probe IDs to probe sequences.
type ProbeToSeqMap = TxHashMap<Probe, ProbeSequence>;

/// A map of probe IDs to integer probe indices.
type ProbeIDToIndex = TxHashMap<String, usize>;

/// A map of half-probe sequences to a set of probes.
type SeqToProbeMap = TxHashMap<Vec<u8>, Vec<Probe>>;

/// The metadata associated with a probe set reference.
pub struct ProbeSetReferenceMetadata(HashMap<String, String>);

impl Deref for ProbeSetReferenceMetadata {
    type Target = HashMap<String, String>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl ProbeSetReferenceMetadata {
    /// Read probe set reference metadata from a file.
    pub fn load_from(path: &TargetSetFile) -> Result<Self> {
        Ok(Self(read_csv_comment_metadata(path)?))
    }

    /// Return true if the this metadata contains probe_set_file_format.
    pub fn is_probe_set_metadata(&self) -> bool {
        self.0.contains_key("probe_set_file_format")
    }

    /// Return the reference_genome as seen in the metadata.
    pub fn reference_genome(&self) -> &str {
        &self["reference_genome"]
    }

    /// Return the reference_genome split by MULTI_GENOME_SEPARATOR collected into a Vec.
    pub fn reference_genomes(&self) -> Vec<GenomeName> {
        self["reference_genome"]
            .split(MULTI_GENOME_SEPARATOR)
            .map(GenomeName::from)
            .collect()
    }

    /// Return the reference_version as seen in the metadata.
    pub fn reference_version(&self) -> &str {
        &self["reference_version"]
    }

    /// Return the panel_name as seen in the metadata.
    pub fn panel_name(&self) -> &str {
        &self["panel_name"]
    }

    /// Return the panel_type as seen in the metadata.
    pub fn panel_type(&self) -> &str {
        &self["panel_type"]
    }

    fn write_to<W: Write>(&self, mut writer: W) -> Result<()> {
        for (key, value) in self.0.iter().sorted_by(|(a, _), (b, _)| a.cmp(b)) {
            writeln!(writer, "#{key}={value}")?;
        }
        Ok(())
    }
}

/// A probe set reference.
/// The CSV comments contain metadata.
/// ```csv
/// #probe_set_file_format=1.0
/// #panel_name=Visium Human Transcriptome Probe Set
/// #panel_type=predesigned
/// #reference_genome=GRCh38
/// #reference_version=2020-A
/// ```
pub struct ProbeSetReference {
    /// key=value metadata headers
    pub metadata: ProbeSetReferenceMetadata,

    /// The minimum alignment score threshold.
    transcriptome_min_score: usize,

    /// The reference genome names
    pub genomes: Vec<GenomeName>,

    /// The length of a probe's left-hand side half.
    /// All LHS probe half sequences must have the same length.
    lhs_length: usize,

    /// The length of a probe's right-hand side half.
    /// All RHS probe half sequences must have the same length.
    rhs_length: usize,

    /// The map of probes to probe sequences.
    probe_to_seq: ProbeToSeqMap,

    /// The map of probes to integer indices, sorted by probe ID.
    probe_id_to_index: ProbeIDToIndex,

    /// The map of left probe sequences to probes.
    lhs_seq_to_probe: SeqToProbeMap,

    /// The map of right probe sequences to probes.
    rhs_seq_to_probe: SeqToProbeMap,

    /// True if this probe set contains one or more gap aligned probes.
    pub has_gap_probes: bool,

    /// True if the probe set contains one or more unpaired gap probes.
    pub has_unpaired_gap_probes: bool,

    /// True if the probe set does not have any paired gap probes
    pub has_no_paired_gap_probes: bool,
}

pub const LHS_START: usize = 0;

impl ProbeSetReference {
    /// Align a half read to the probe set reference, allowing up to one mismatch.
    /// Return no match if there is more than one possible match.
    pub fn align_half_read<'a>(
        &self,
        seq_to_probe: &'a SeqToProbeMap,
        seq: &[u8],
    ) -> Option<(&'a [Probe], i32)> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

        if seq.len() < min(self.lhs_length, self.rhs_length) {
            // Half read sequence is shorter than half probe sequence.
            return None;
        }

        if let Some(probes) = seq_to_probe.get(seq) {
            // Return the perfect match.
            return Some((probes, seq.len() as i32));
        }

        let mut first_match = None;
        let mut mutant_seq = seq.to_vec();
        for (pos, &orig_nt) in seq.iter().enumerate() {
            for mutant_nt in BASES {
                if mutant_nt != orig_nt {
                    mutant_seq[pos] = mutant_nt;
                    match (first_match, seq_to_probe.get(&mutant_seq)) {
                        (_, None) => (),
                        (None, Some(probes)) => first_match = Some(probes.as_slice()),
                        (Some(_), Some(_)) => return None, // Multiple matches, so return none.
                    }
                }
            }
            mutant_seq[pos] = orig_nt;
        }
        if let Some(probes) = first_match {
            let mismatches = 1;
            let matches = seq.len() as i32 - mismatches;
            let score = matches - mismatches;
            Some((probes, score))
        } else {
            None
        }
    }

    /// Align a read to the probe set reference.
    pub fn align_probe_read(&self, seq: &[u8]) -> MappedProbe {
        if seq.len() < min(self.lhs_length, self.rhs_length) {
            // Read sequence is shorter than half the probe sequence.
            return MappedProbe::new(None, None);
        }

        // intial alignment of the lhs probe half
        let lhs_seq = &seq[LHS_START..self.lhs_length];
        let (lhs_probes, lhs_score) = self
            .align_half_read(&self.lhs_seq_to_probe, lhs_seq)
            .unwrap_or((&[], 0));

        // initial alignment of rhs probe half
        let (rhs_probes, rhs_score, rhs_start) = {
            // return first encountered kmer with RHS (if any)
            let mut read_kmers = get_rhs_read_kmers(
                self.has_gap_probes,
                self.lhs_length,
                lhs_probes,
                self.rhs_length,
                seq,
            );
            read_kmers.find_map(|(read_pos, rhs_seq)| {
                self.align_half_read(&self.rhs_seq_to_probe, rhs_seq)
                    .map(|(probes, score)| (probes, score, read_pos))
            })
        }
        .unwrap_or((&[], 0, self.lhs_length));

        let mut mapped_probe = match (lhs_probes, rhs_probes) {
            ([], []) => MappedProbe::new(None, None),
            ([..], []) => {
                // return first encountered kmer with RHS (if any)
                let mut read_kmers = get_rhs_read_kmers(
                    self.has_gap_probes,
                    self.lhs_length,
                    lhs_probes,
                    self.rhs_length,
                    seq,
                );
                read_kmers
                    .find_map(|(read_pos, rhs_seq)| {
                        let res =
                            self.rescue_rhs(LHS_START, lhs_probes, lhs_score, read_pos, rhs_seq);
                        if res.rhs.is_some() { Some(res) } else { None }
                    })
                    .unwrap_or_else(|| {
                        MappedProbe::new(
                            Some(MappedProbeHalf::new(
                                lhs_probes[0].clone(),
                                lhs_score,
                                LHS_START,
                                LHS_START + self.lhs_length,
                            )),
                            None,
                        )
                    })
            }
            ([], [..]) => self.rescue_lhs(LHS_START, lhs_seq, rhs_start, rhs_probes, rhs_score),
            ([lhs], [rhs]) => MappedProbe::new(
                Some(MappedProbeHalf::new(
                    lhs.clone(),
                    self.lhs_length as i32,
                    LHS_START,
                    LHS_START + self.lhs_length,
                )),
                Some(MappedProbeHalf::new(
                    rhs.clone(),
                    self.rhs_length as i32,
                    rhs_start,
                    rhs_start + self.rhs_length,
                )),
            ),
            ([..], [..]) => {
                // Confident matches are the intersection of lhs_probes and rhs_probes.
                // Multiple matches are caused by multiple probes with identical sequence,
                // in which case use the lexicographically minimal probe ID.
                if let Some(probe) = intersect_sorted_lists(lhs_probes, rhs_probes).next() {
                    MappedProbe::new(
                        Some(MappedProbeHalf::new(
                            probe.clone(),
                            self.lhs_length as i32,
                            LHS_START,
                            LHS_START + self.lhs_length,
                        )),
                        Some(MappedProbeHalf::new(
                            probe.clone(),
                            self.rhs_length as i32,
                            rhs_start,
                            rhs_start + self.rhs_length,
                        )),
                    )
                } else {
                    MappedProbe::new(
                        Some(MappedProbeHalf::new(
                            lhs_probes[0].clone(),
                            self.lhs_length as i32,
                            LHS_START,
                            LHS_START + self.lhs_length,
                        )),
                        Some(MappedProbeHalf::new(
                            rhs_probes[0].clone(),
                            self.rhs_length as i32,
                            rhs_start,
                            rhs_start + self.rhs_length,
                        )),
                    )
                }
            }
        };
        if mapped_probe.is_paired_gap_align() {
            let probe = mapped_probe.lhs_probe().unwrap();
            let expected_gap_seq = &self.probe_to_seq[probe].gap;
            assert!(!expected_gap_seq.is_empty());

            let gap_seq = mapped_probe.gap_info().unwrap().gap_seq(seq);

            mapped_probe.gap = Some(MappedGap {
                gap_seq: std::str::from_utf8(gap_seq).unwrap().to_string(),
                expected_gap_seq: std::str::from_utf8(expected_gap_seq).unwrap().to_string(),
            });
        }
        mapped_probe
    }

    /// Rescue an unaligned half read. Return the mapping and alignment score.
    /// When one half of a read aligns to the probe sequence lookup table and the
    /// other half does not, attempt to rescue the unaligned half by aligning its
    /// read sequence to its corresponding partner probe sequence.
    fn align_unaligned_half(
        &self,
        mapped_probe: &Probe,
        mapped_score: i32,
        probe_seq: &[u8],
        read_pos: usize,
        read_seq: &[u8],
    ) -> Option<MappedProbeHalf> {
        assert!(read_seq.len() <= probe_seq.len());
        let mismatches = if read_seq.len() == probe_seq.len() {
            bio::alignment::distance::simd::hamming(probe_seq, read_seq) as i32
        } else {
            zip(read_seq, probe_seq).filter(|(a, b)| a != b).count() as i32
        };
        let matches = read_seq.len() as i32 - mismatches;
        let score = matches - mismatches;
        // The score of the half-probe alignment must be positive, and the sum
        // of both half-probe score must be at least transcriptome_min_score.
        if score > 0 && mapped_score + score >= self.transcriptome_min_score as i32 {
            Some(MappedProbeHalf::new(
                mapped_probe.clone(),
                score,
                read_pos,
                read_pos + probe_seq.len(),
            ))
        } else {
            None
        }
    }

    /// Rescue an unaligned half-probe sequence.
    #[allow(clippy::too_many_arguments)]
    fn rescue_unaligned_half(
        &self,
        mapped_probes: &[Probe],
        mapped_score: i32,
        mapped_probe_pos: usize,
        mapped_probe_length: usize,
        read_pos: usize,
        read_seq: &[u8],
        rescue_which_half: impl Fn(&ProbeSequence) -> &[u8],
    ) -> (MappedProbeHalf, Option<MappedProbeHalf>) {
        let best = mapped_probes
            .iter()
            .rev() // Return the first probe with max score.
            .filter_map(|probe| {
                let probe_seq = rescue_which_half(&self.probe_to_seq[probe]);
                if !probe_seq.is_empty() {
                    self.align_unaligned_half(probe, mapped_score, probe_seq, read_pos, read_seq)
                } else {
                    None
                }
            })
            .max_by_key(|x| x.score);
        if let Some(ref x) = best {
            (
                MappedProbeHalf::new(
                    x.probe.clone(),
                    mapped_score,
                    mapped_probe_pos,
                    mapped_probe_pos + mapped_probe_length,
                ),
                best,
            )
        } else {
            (
                MappedProbeHalf::new(
                    mapped_probes[0].clone(),
                    mapped_score,
                    mapped_probe_pos,
                    mapped_probe_pos + mapped_probe_length,
                ),
                None,
            )
        }
    }

    /// Rescue an unaligned LHS sequence.
    fn rescue_lhs(
        &self,
        lhs_read_pos: usize,
        lhs_read_seq: &[u8],
        rhs_read_pos: usize,
        rhs_probes: &[Probe],
        rhs_score: i32,
    ) -> MappedProbe {
        let (rhs, opt_lhs) = self.rescue_unaligned_half(
            rhs_probes,
            rhs_score,
            rhs_read_pos,
            self.lhs_length,
            lhs_read_pos,
            lhs_read_seq,
            |x| &x.lhs,
        );
        MappedProbe::new(opt_lhs, Some(rhs))
    }

    /// Rescue an unaligned RHS sequence.
    fn rescue_rhs(
        &self,
        lhs_read_pos: usize,
        lhs_probes: &[Probe],
        lhs_score: i32,
        rhs_read_pos: usize,
        rhs_read_seq: &[u8],
    ) -> MappedProbe {
        let (lhs, opt_rhs) = self.rescue_unaligned_half(
            lhs_probes,
            lhs_score,
            lhs_read_pos,
            self.rhs_length,
            rhs_read_pos,
            rhs_read_seq,
            |x| &x.rhs,
        );
        MappedProbe::new(Some(lhs), opt_rhs)
    }

    /// Construct a probe set refrence from a probe set CSV file.
    /// CSV header must include #probe_set_file_format.
    /// CSV header row must be gene_id,probe_seq,probe_id,included,region.
    pub fn from_path(
        target_set: &TargetSetFile,
        transcriptome_reference_path: Option<&Path>,
        transcriptome_min_score: usize,
    ) -> Result<Self> {
        let metadata = ProbeSetReferenceMetadata::load_from(target_set)?;
        assert!(metadata.is_probe_set_metadata());

        let genomes = if let Some(transcriptome_reference_path) = transcriptome_reference_path {
            ReferenceInfo::from_reference_path(transcriptome_reference_path)?.genomes
        } else {
            ReferenceInfo::from_probe_set_csv(target_set)?.genomes
        };

        // Read the probe set CSV file.
        let mut reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .from_path(target_set)
            .with_context(|| target_set.display().to_string())?;

        // Parse the probe set CSV file.
        let mut has_gap_probes: bool = false;
        let mut has_unpaired_gap_probes: bool = false;
        let mut has_no_paired_gap_probes: bool = true;
        let mut lhs_length = None;
        let mut rhs_length = None;
        let probe_to_seq: ProbeToSeqMap = zip_eq(
            TargetSetFile::read(target_set, transcriptome_reference_path)?,
            reader.records(),
        )
        .map(|(probe, records)| {
            has_gap_probes |= match probe.probe_type {
                ProbeType::UnpairedGapAlign => true,
                ProbeType::RTL => false,
                ProbeType::PairedGapAlign => true,
            };
            has_unpaired_gap_probes |= probe.probe_type == ProbeType::UnpairedGapAlign;
            has_no_paired_gap_probes &= probe.probe_type != ProbeType::PairedGapAlign;

            // Ensure that the length of the LHS and RHS sequences is fixed.
            let probe_seq = ProbeSequence::from_str(&records?[1])?;
            if !probe_seq.lhs.is_empty() {
                if let Some(lhs_length) = lhs_length {
                    assert_eq!(probe_seq.lhs.len(), lhs_length);
                } else {
                    lhs_length = Some(probe_seq.lhs.len());
                }
            }
            if !probe_seq.rhs.is_empty() {
                if let Some(rhs_length) = rhs_length {
                    assert_eq!(probe_seq.rhs.len(), rhs_length);
                } else {
                    rhs_length = Some(probe_seq.rhs.len());
                }
            }

            anyhow::Ok((probe, probe_seq))
        })
        .try_collect()?;

        let Some(lhs_length) = lhs_length else {
            bail!("probe_set CSV has no LHS records: {}", target_set.display())
        };
        let Some(rhs_length) = rhs_length else {
            bail!("probe_set CSV has no RHS records: {}", target_set.display())
        };

        // Construct a map of probes to integer indices, sorted by probe ID.
        let probe_id_to_index: TxHashMap<_, _> = probe_to_seq
            .keys()
            .sorted()
            .enumerate()
            .map(|(i, probe)| (probe.probe_id.clone(), i))
            .collect();

        // Construct a map of probe sequences to probes.
        let mut lhs_seq_to_probe: SeqToProbeMap = probe_to_seq
            .iter()
            .filter(|(_, probe_sequence)| !probe_sequence.lhs.is_empty())
            .fold(Default::default(), |mut map, (probe, seq)| {
                map.entry(seq.lhs.clone()).or_default().push(probe.clone());
                map
            });
        let mut rhs_seq_to_probe: SeqToProbeMap = probe_to_seq
            .iter()
            .filter(|(_, probe_sequence)| !probe_sequence.rhs.is_empty())
            .fold(Default::default(), |mut map, (probe, seq)| {
                map.entry(seq.rhs.clone()).or_default().push(probe.clone());
                map
            });

        // Sort probes by probe ID.
        for probes in lhs_seq_to_probe.values_mut() {
            probes.sort();
        }
        for probes in rhs_seq_to_probe.values_mut() {
            probes.sort();
        }

        Ok(Self {
            metadata,
            transcriptome_min_score,
            genomes,
            lhs_length,
            rhs_length,
            probe_to_seq,
            probe_id_to_index,
            lhs_seq_to_probe,
            rhs_seq_to_probe,
            has_gap_probes,
            has_unpaired_gap_probes,
            has_no_paired_gap_probes,
        })
    }

    /// Return an iterator over the probes in an arbitrary order.
    pub fn unsorted_probes(&self) -> impl Iterator<Item = &Probe> {
        self.probe_to_seq.keys()
    }

    /// Return a vector of the probes sorted by probe ID.
    pub fn sorted_probes(&self) -> Vec<&Probe> {
        self.probe_to_seq.keys().sorted().collect()
    }

    /// Return the number of probes.
    pub fn number_of_probes(&self) -> usize {
        self.probe_to_seq.len()
    }

    /// Return the integer index of this probe in the probe set sorted by probe ID.
    pub fn probe_id_to_index(&self, probe_id: &str) -> usize {
        self.probe_id_to_index[probe_id]
    }

    pub fn probe_sequence_from_id(&self, probe_id: &str) -> &ProbeSequence {
        self.probe_to_seq
            .get(
                self.sorted_probes()
                    .get(self.probe_id_to_index(probe_id))
                    .unwrap(),
            )
            .unwrap()
    }
}

// Helper function returning an iterator over read kmers to try match to the RHS probe halves in the probe set reference
fn get_rhs_read_kmers<'a>(
    has_gap_probes: bool,
    lhs_length: usize,
    lhs_probes: &'a [Probe],
    rhs_length: usize,
    seq: &'a [u8],
) -> Box<dyn Iterator<Item = (usize, &'a [u8])> + 'a> {
    // If there is room for the RHS seq, set up the kmers to screen
    if lhs_length <= seq.len() - rhs_length {
        if has_gap_probes {
            // To do gap alignment scan full read in reverse direction to return the LAST match, i.e. longest gap (maybe change this to return best match later?)
            // For rescue_rhs during gap alignment this returns the LAST rhs match in the reads sequence that results in a score >= self.transcriptome_min_score, i.e. longest gap
            Box::new(
                seq.windows(rhs_length)
                    .enumerate()
                    .skip(if lhs_probes.is_empty() { 0 } else { lhs_length })
                    .rev(),
            )
        } else {
            // For non gapped alignment only check pos adjecent to LHS
            Box::new(vec![(lhs_length, &seq[lhs_length..(lhs_length + rhs_length)])].into_iter())
        }
    } else {
        // if the the read is too short, do nothing by returning an empty iterator
        Box::new(iter::empty::<(usize, &[u8])>())
    }
}

/// Read key/value metadata from the beginning of a CSV file.
/// Each line should be of the form #key=value
pub fn read_csv_comment_metadata(path: &Path) -> Result<HashMap<String, String>> {
    let mut metadata = HashMap::new();
    for line in
        BufReader::new(File::open(path).with_context(|| path.display().to_string())?).lines()
    {
        let line = line?;
        if !line.starts_with('#') {
            break;
        }
        let (key, value) = line[1..]
            .splitn(2, '=')
            .collect_tuple()
            .ok_or_else(|| anyhow!("CSV comment metadata line missing value: \"{line}\""))?;
        metadata.insert(key.to_string(), value.to_string());
    }
    Ok(metadata)
}

/// The result of mapping one half of a read to a probe set reference.
#[derive(Deserialize, Serialize)]
pub struct MappedProbeHalf {
    /// The mapped probe.
    pub probe: Probe,

    /// The alignment score.
    pub score: i32,

    /// The alignment start position.
    pub start_pos: usize,

    /// The alignment end position.
    pub end_pos: usize,
}

impl MappedProbeHalf {
    /// Return a new MappedProbeHalf.
    pub fn new(probe: Probe, score: i32, start_pos: usize, end_pos: usize) -> Self {
        assert!(score > 0);
        MappedProbeHalf {
            probe,
            score,
            start_pos,
            end_pos,
        }
    }
}

const MAPPED_GAP_ALIGNMENT_MATCH_SCORE: i32 = 1;
const MAPPED_GAP_ALIGNMENT_MISMATCH_SCORE: i32 = -1;
const MAPPED_GAP_ALIGNMENT_GAP_OPEN_SCORE: i32 = -1;
const MAPPED_GAP_ALIGNMENT_GAP_EXTEND_SCORE: i32 = -1;
const MAPPED_GAP_ALIGNMENT_CLIP_SCORE: i32 = -3;

/// The gap part of the probe mapping result
/// Mapped gap is only created for reads that LHS and RHS map to the same probe
#[derive(Deserialize, Serialize, Clone)]
pub struct MappedGap {
    /// Gap sequence extracted from the read
    pub gap_seq: String,

    /// Expected gap sequence from probe set
    pub expected_gap_seq: String,
}

impl MappedGap {
    pub fn get_expected_gap_seq(&self) -> String {
        self.expected_gap_seq.clone()
    }

    pub fn get_alignment(&self) -> Alignment {
        // custom scoring to enable clipping in mapped and expected gap sequence
        let scoring = Scoring::from_scores(
            MAPPED_GAP_ALIGNMENT_GAP_OPEN_SCORE,
            MAPPED_GAP_ALIGNMENT_GAP_EXTEND_SCORE,
            MAPPED_GAP_ALIGNMENT_MATCH_SCORE,
            MAPPED_GAP_ALIGNMENT_MISMATCH_SCORE,
        )
        .xclip(MAPPED_GAP_ALIGNMENT_CLIP_SCORE)
        .yclip(MAPPED_GAP_ALIGNMENT_CLIP_SCORE);

        let mut aligner = Aligner::with_capacity_and_scoring(
            self.gap_seq.len(),
            self.expected_gap_seq.len(),
            scoring,
        );
        aligner.custom(self.gap_seq.as_bytes(), self.expected_gap_seq.as_bytes())
    }
    pub fn get_max_gap_error(&self) -> u32 {
        std::cmp::max(
            MIN_ALLOWED_GAP_ERROR,
            (self.expected_gap_seq.len() as f64 * MAX_GAP_ERROR_RATE).ceil() as u32,
        )
    }
    pub fn gap_within_max_error(&self) -> bool {
        bounded_levenshtein(
            self.expected_gap_seq.as_bytes(),
            self.gap_seq.as_bytes(),
            self.get_max_gap_error(),
        )
        .is_some()
    }
    pub fn get_gap_levenshtein_distance(&self) -> u32 {
        levenshtein(self.expected_gap_seq.as_bytes(), self.gap_seq.as_bytes())
    }
}

pub struct MappedGapAlignmentInfo {
    pub gap_seq: String,
    pub expected_gap_seq: String,
    pub alignment_operations: Vec<AlignmentOperation>,
}

impl MappedGapAlignmentInfo {
    pub fn gap_exactly_matches_expected(&self) -> bool {
        self.gap_seq == self.expected_gap_seq
    }

    pub fn get_num_matches(&self) -> usize {
        self.alignment_operations
            .iter()
            .filter(|&a| *a == AlignmentOperation::Match)
            .count()
    }
    pub fn get_num_mismatches(&self) -> usize {
        self.alignment_operations
            .iter()
            .filter(|&a| *a == AlignmentOperation::Subst)
            .count()
    }

    pub fn ends_with_insertion(&self) -> bool {
        matches!(
            self.alignment_operations.last(),
            Some(&AlignmentOperation::Xclip(_) | &AlignmentOperation::Ins)
        )
    }

    pub fn starts_with_insertion(&self) -> bool {
        matches!(
            self.alignment_operations.first(),
            Some(&AlignmentOperation::Xclip(_) | &AlignmentOperation::Ins)
        )
    }

    pub fn ends_with_deletion(&self) -> bool {
        matches!(
            self.alignment_operations.last(),
            Some(&AlignmentOperation::Yclip(_) | &AlignmentOperation::Del)
        )
    }

    pub fn starts_with_deletion(&self) -> bool {
        matches!(
            self.alignment_operations.first(),
            Some(&AlignmentOperation::Yclip(_) | &AlignmentOperation::Del)
        )
    }

    /// Number of insertions in gap compared to expected gap (including soft clips before or after)
    pub fn get_num_insertions(&self) -> usize {
        self.alignment_operations
            .iter()
            .map(|op| match op {
                AlignmentOperation::Xclip(l) => *l,
                AlignmentOperation::Ins => 1,
                _ => 0,
            })
            .sum()
    }

    /// Number of deletions in gap compared to expected gap (including soft clips before or after)
    pub fn get_num_deletions(&self) -> usize {
        self.alignment_operations
            .iter()
            .map(|op| match op {
                AlignmentOperation::Yclip(l) => *l,
                AlignmentOperation::Del => 1,
                _ => 0,
            })
            .sum()
    }
}

/// The result of mapping a read to a probe set reference.
#[derive(Deserialize, Serialize)]
pub struct MappedProbe {
    /// The result of mapping the left sequence.
    lhs: Option<MappedProbeHalf>,

    /// The result of mapping the right sequence.
    rhs: Option<MappedProbeHalf>,

    /// The read maps confidently to a probe but multimapped with STAR.
    is_rescued: bool,

    /// The result of mapping the gap sequence.
    /// Only populated when the LHS/RHS sequences are mapped to the same probe
    /// and the probe is a gap align probe.
    gap: Option<MappedGap>,
}

impl MappedProbe {
    /// Create a new mapping result.
    pub fn new(lhs: Option<MappedProbeHalf>, rhs: Option<MappedProbeHalf>) -> Self {
        MappedProbe {
            lhs,
            rhs,
            is_rescued: false,
            gap: None,
        }
    }

    /// Return the result of mapping the left probe sequence.
    pub fn lhs_probe(&self) -> Option<&Probe> {
        self.lhs.as_ref().map(|x| &x.probe)
    }

    /// Return the result of mapping the right probe sequence.
    pub fn rhs_probe(&self) -> Option<&Probe> {
        self.rhs.as_ref().map(|x| &x.probe)
    }

    /// Return the alignment score of the left half.
    pub fn lhs_score(&self) -> i32 {
        self.lhs.as_ref().map_or(0, |x| x.score)
    }

    /// Return the alignment score of the right half.
    pub fn rhs_score(&self) -> i32 {
        self.rhs.as_ref().map_or(0, |x| x.score)
    }

    /// Return the alignment score.
    pub fn score(&self) -> i32 {
        self.lhs_score() + self.rhs_score()
    }

    /// Return gap
    pub fn gap(&self) -> Option<MappedGap> {
        self.gap.clone()
    }

    /// Return an iterator over the genomes.
    /// Return one genome when both probes map to the same genome.
    /// Return two genomes when the two probes map to different genomes.
    /// Return one genome when one probe maps and the other does not.
    /// Return no genomes when no probes map.
    pub fn genomes(&self) -> impl Iterator<Item = &GenomeName> {
        let lhs = self.lhs_probe().map(|x| &x.genome);
        let rhs = self.rhs_probe().map(|x| &x.genome);
        chain(lhs, if lhs == rhs { None } else { rhs })
    }

    /// Return whether either left or right sequences mapped to a probe.
    pub fn is_mapped(&self) -> bool {
        self.lhs.is_some() || self.rhs.is_some()
    }

    pub fn is_lhs_rhs_mapped_to_same_probe(&self) -> bool {
        match (self.lhs_probe(), self.rhs_probe()) {
            (Some(lhs), Some(rhs)) => lhs == rhs,
            _ => false,
        }
    }

    /// Return whether both left and right sequences mapped to the same probe.
    pub fn is_conf_mapped(&self) -> bool {
        self.is_lhs_rhs_mapped_to_same_probe()
            && self
                .gap
                .as_ref()
                .is_none_or(MappedGap::gap_within_max_error)
    }

    /// Return the confidently mapped gene.
    pub fn conf_gene(&self) -> Option<&Gene> {
        if self.is_conf_mapped() {
            Some(&self.lhs_probe().unwrap().gene)
        } else {
            None
        }
    }

    /// Return an iterator over the mapped genes.
    /// Return one gene when both probes map to the same gene.
    /// Return two genes when the two probes map to different genes.
    /// Return one gene when one probe maps and the other does not.
    /// Return no genes when no probes map.
    pub fn genes(&self) -> impl Iterator<Item = &Gene> {
        let lhs = self.lhs_probe().map(|x| &x.gene);
        let rhs = self.rhs_probe().map(|x| &x.gene);
        chain(lhs, if lhs == rhs { None } else { rhs })
    }

    /// Return an iterator over the mapped (genome, gene)'s.
    /// Return one (genome, gene) when both probes map to the same (genome, gene).
    /// Return two (genome, gene)'s when the two probes map to different (genome, gene)'s.
    /// Return one (genome, gene) when one probe maps and the other does not.
    /// Return no (genome, gene)'s when no probes map.
    pub fn mapped_genes(&self) -> impl Iterator<Item = (&GenomeName, &Gene)> {
        let lhs = self.lhs_probe().map(|x| (&x.genome, &x.gene));
        let rhs = self.rhs_probe().map(|x| (&x.genome, &x.gene));
        chain(lhs, if lhs == rhs { None } else { rhs })
    }

    /// Return an iterator over the mapped probes.
    /// Return one probe when both halves of the read map to the same probe.
    /// Return two probes when the two halves of the read map to different probes or either probe does not map.
    /// Return no probes when neither half of the read maps.
    /// A string of "NA" indicates that half of the read does not map.
    pub fn probes(&self) -> impl Iterator<Item = &Probe> {
        lazy_static! {
            /// A static reference to an unmapped probe.
            static ref PROBE_NA: Probe = Probe {
                probe_id: "NA".to_string(),
                probe_seq: String::new(),
                gene: Gene {
                    id: String::new(),
                    name: String::new(),
                },
                included: true,
                region: None,
                probe_type: ProbeType::default(),
                ref_sequence_name: String::new(),
                ref_sequence_pos: None,
                cigar_string: String::new(),
                genome: GenomeName::default(),
            };
        }

        match (self.lhs_probe(), self.rhs_probe()) {
            (None, None) => chain(None, None),
            (lhs, rhs) if lhs == rhs => chain(lhs, None),
            (lhs, rhs) => chain(lhs.or(Some(&PROBE_NA)), rhs.or(Some(&PROBE_NA))),
        }
    }

    /// Return the mapping quality of this alignment.
    /// 255 when both halves map to the same probe.
    /// 3 when the two halves map to different probes.
    /// 1 when one half maps to a probe and the other half does not.
    /// 0 when neither half maps to a probe.
    pub fn mapq(&self) -> u8 {
        match (&self.lhs, &self.rhs) {
            (Some(_), Some(_)) => {
                if self.is_conf_mapped() {
                    AlignerParam::Hurtle.high_conf_mapq()
                } else if self.is_lhs_rhs_mapped_to_same_probe() {
                    assert!(self.gap.is_some());
                    MAPQ_GAP_MAPPED_NOT_WITHIN_MAX_ERR
                } else {
                    MAPQ_SPLIT_MAPPED
                }
            }
            (Some(_), None) | (None, Some(_)) => MAPQ_HALF_MAPPED,
            (None, None) => 0,
        }
    }

    /// The read maps confidently to a probe but multimapped with STAR.
    pub fn is_rescued(&self) -> bool {
        self.is_rescued
    }

    /// Set the flag returned by is_rescued.
    pub fn set_rescued(&mut self, is_rescued: bool) {
        self.is_rescued = is_rescued;
    }

    /// Whether both halves are LHS and RHS of a PairedGapAlign probes
    pub fn is_paired_gap_align(&self) -> bool {
        match (self.lhs_probe(), self.rhs_probe()) {
            (Some(lhs), Some(rhs)) => lhs == rhs && rhs.probe_type == ProbeType::PairedGapAlign,
            _ => false,
        }
    }

    pub fn is_unpaired_gap_align(&self) -> bool {
        match (self.lhs_probe(), self.rhs_probe()) {
            (Some(lhs), Some(rhs)) => {
                lhs.probe_type == ProbeType::UnpairedGapAlign
                    && rhs.probe_type == ProbeType::UnpairedGapAlign
            }
            _ => false,
        }
    }

    pub fn is_gap_align(&self) -> bool {
        self.is_unpaired_gap_align() || self.is_paired_gap_align()
    }

    /// Return the gap information if available
    /// Option tuple of (lhs end position, rhs start position)
    /// if probe type gapped and both probes aligned else None
    pub fn gap_info(&self) -> Option<GapInfo> {
        if !self.is_gap_align() {
            return None;
        }
        match (&self.lhs, &self.rhs) {
            (Some(lhs), Some(rhs)) => Some(GapInfo {
                probe_type: lhs.probe.probe_type, // == rhs.probe.probe_type
                lhs_end: lhs.end_pos,
                rhs_start: rhs.start_pos,
            }),
            _ => unreachable!(), // is_gap_align() would return false
        }
    }
}

pub struct GapInfo {
    pub probe_type: ProbeType,
    pub lhs_end: usize,
    pub rhs_start: usize,
}

impl GapInfo {
    pub fn gap_seq<'a>(&self, read: &'a [u8]) -> &'a [u8] {
        let start = self.lhs_end;
        let end = self.rhs_start.max(start);
        &read[start..end]
    }
}

/// Merge probe set CSVs, and return its combined probe set name.
///
/// Remove duplicate probes from the resulting merged set, allowing the input
/// probe sets to contain some overlapping probes.
pub fn merge_probe_set_csvs<W: Write>(
    probe_sets: &[TargetSetFile],
    mut writer: W,
    transcriptome_reference_path: Option<&Path>,
) -> Result<String> {
    assert!(probe_sets.len() >= 2);
    let metadatas: Vec<_> = probe_sets
        .iter()
        .map(ProbeSetReferenceMetadata::load_from)
        .try_collect()?;

    let panel_name = metadatas
        .iter()
        .map(ProbeSetReferenceMetadata::panel_name)
        .join("|");
    let panel_type = metadatas
        .iter()
        .map(ProbeSetReferenceMetadata::panel_type)
        .unique()
        .join("|");

    for (genome, versions) in &metadatas
        .iter()
        .map(|x| (x.reference_genome(), x.reference_version()))
        .into_group_map()
    {
        ensure!(
            versions.iter().unique().count() == 1,
            "All reference_version of reference_genome {genome} must be identical: {}",
            versions.join(", ")
        );
    }

    #[expect(unstable_name_collisions)]
    let (reference_genome, reference_version): (String, String) = metadatas
        .iter()
        .map(|x| (x.reference_genome(), x.reference_version()))
        .unique()
        .intersperse((MULTI_GENOME_SEPARATOR, MULTI_GENOME_SEPARATOR))
        .collect();

    let combined_metadata = ProbeSetReferenceMetadata(HashMap::from([
        (
            "probe_set_file_format".to_string(),
            PROBE_SET_FILE_FORMAT.to_string(),
        ),
        ("reference_genome".to_string(), reference_genome),
        ("reference_version".to_string(), reference_version),
        ("panel_name".to_string(), panel_name),
        ("panel_type".to_string(), panel_type),
    ]));

    let unique_genome_count = metadatas
        .iter()
        .map(ProbeSetReferenceMetadata::reference_genome)
        .unique()
        .count();
    let max_genome_len = metadatas
        .iter()
        .map(|x| x.reference_genome().len())
        .max()
        .unwrap();

    fn infix_genome_to_gene_id_or_gene_name(s: &str, padded_genome: &str) -> String {
        if let Some(depr_gid) = s.strip_prefix(DEPRECATED_PROBE_ID_PREFIX) {
            format!("{DEPRECATED_PROBE_ID_PREFIX}{padded_genome}_{depr_gid}")
        } else {
            format!("{padded_genome}_{s}")
        }
    }

    // Write the combined probe set CSV to a temp file first transcriptome_reference_path = None.
    // The temp file will have the genome prefixed gene_id, gene_name, and probe_id.
    // Now that the temp file has the correct gene_id, we parse it with transcriptome_reference_path
    // which allows us to extract the gene_name from the transcriptome reference for the probes.

    let mut tmp_file = NamedTempFile::with_suffix(".csv")?;
    let tmp_file_path = tmp_file.path().to_owned();
    combined_metadata.write_to(&mut tmp_file)?;
    let mut tmp_csv_writer = csv::WriterBuilder::new().from_writer(tmp_file);
    tmp_csv_writer.write_record(PROBE_SET_HEADER)?;

    for (probe_set, metadata) in zip_eq(probe_sets, metadatas) {
        let reference_genome = metadata.reference_genome();
        let padded_genome = format!("{reference_genome:_<max_genome_len$}");
        for probe in probe_set.read(None)? {
            // Add the genome name to the gene_id, gene_name, and probe_id.
            tmp_csv_writer.write_record(
                if unique_genome_count >= 2 {
                    Probe {
                        gene: Gene {
                            id: infix_genome_to_gene_id_or_gene_name(
                                &probe.gene.id,
                                &padded_genome,
                            ),
                            name: infix_genome_to_gene_id_or_gene_name(
                                &probe.gene.name,
                                &padded_genome,
                            ),
                        },
                        probe_id: infix_genome_to_gene_id_or_gene_name(
                            &probe.probe_id,
                            &padded_genome,
                        ),
                        ..probe
                    }
                } else {
                    probe
                }
                .as_record(),
            )?;
        }
    }
    tmp_csv_writer.flush()?;

    // Write the combined probe set CSV to the final file.
    combined_metadata.write_to(&mut writer)?;
    let mut csv_writer = csv::Writer::from_writer(writer);
    csv_writer.write_record(PROBE_SET_HEADER)?;
    let temp_probe_set_file = TargetSetFile::from(tmp_file_path);

    // Filter out duplicate probes in case our inputs had some of the same probes.
    for probe in temp_probe_set_file
        .read(transcriptome_reference_path)?
        .into_iter()
        .unique()
    {
        csv_writer.write_record(probe.as_record())?;
    }
    csv_writer.flush()?;
    Ok(combined_metadata.panel_name().to_string())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_intersect_sorted_lists() {
        assert!(
            intersect_sorted_lists(&[0, 2, 4, 6, 8], &[1, 3, 5, 7, 9])
                .next()
                .is_none()
        );
        assert_eq!(
            intersect_sorted_lists(&[1, 2, 3, 4, 5], &[1, 3, 5, 7, 9]).collect::<Vec<_>>(),
            [&1, &3, &5]
        );
        assert_eq!(
            intersect_sorted_lists([1, 2, 3, 4, 5].iter(), [5].iter()).collect::<Vec<_>>(),
            [&5]
        );
        assert_eq!(
            intersect_sorted_lists(
                [1, 5, 10, 12, 13, 14, 15, 16, 18, 19].iter().by_ref(),
                [1, 2, 4, 5, 7, 8, 11, 14, 17, 20].iter().by_ref()
            )
            .collect::<Vec<_>>(),
            [&1, &5, &14]
        );
    }

    #[test]
    fn test_load_metadata() -> Result<()> {
        assert!(
            !ProbeSetReferenceMetadata::load_from(&TargetSetFile::from(
                "../cr_lib/test/target_panels/Immunology_targeting_hybrid.csv"
            ))?
            .is_probe_set_metadata()
        );
        assert!(
            ProbeSetReferenceMetadata::load_from(&TargetSetFile::from(
                "../cr_lib/test/target_panels/Immunology_targeting_templated_ligation.csv"
            ))?
            .is_probe_set_metadata()
        );
        Ok(())
    }

    #[test]
    fn test_from_path() -> Result<()> {
        ProbeSetReference::from_path(
            &TargetSetFile::from(
                "../cr_lib/test/target_panels/Immunology_targeting_templated_ligation.csv",
            ),
            None,
            0,
        )?;
        Ok(())
    }

    #[test]
    fn test_file_format_3() -> Result<()> {
        let probe_set = ProbeSetReference::from_path(
            &TargetSetFile::from("../cr_lib/test/target_panels/mock_probe_set_file_format_3.csv"),
            None,
            0,
        )?;

        let probes = probe_set.sorted_probes();
        assert_eq!(probes.len(), 87);

        let first_probe = probes.first().unwrap();
        assert_eq!(first_probe.region, Some(ProbeRegion::Spliced));
        assert_eq!(first_probe.probe_type, ProbeType::RTL);
        assert_eq!(first_probe.ref_sequence_name, "not_a_real_ref_2");
        assert_eq!(first_probe.ref_sequence_pos, Some(5));
        assert_eq!(first_probe.cigar_string, "50M");
        let mapped_probe = probe_set.align_probe_read(b"TCATACTCCTGCTTGCTGATCCACATCTGCTGGAAGGTGGACAGCGAGGCGTCAGCGACTACGTACGTACGTAGCTGGGCATGCGATCG");
        assert_eq!(mapped_probe.lhs_probe(), mapped_probe.rhs_probe());
        assert_eq!(&mapped_probe.probes().next().unwrap(), first_probe);

        let last_probe = probes.last().unwrap();
        assert_eq!(last_probe.region, Some(ProbeRegion::Unspliced));
        assert_eq!(last_probe.probe_type, ProbeType::PairedGapAlign);
        assert_eq!(last_probe.ref_sequence_name, "not_a_real_ref");
        assert_eq!(last_probe.ref_sequence_pos, Some(77));
        assert_eq!(last_probe.cigar_string, "200M");
        let mapped_probe = probe_set.align_probe_read(b"TGGCCATCGGGCAGCTCGTAGCTCTNNNNNNNTCTCCAGAGAAGAGGAGGATNNNNNAGTCGGTCAGGTCCCGGCCAGCCAGNNNNNNNGCGGCGGTGGCCATCTCCTGCTCGAAGTCCANNNNNAGGATCTTCATGAGGT");
        assert_eq!(mapped_probe.lhs_probe(), mapped_probe.rhs_probe());
        assert_eq!(&mapped_probe.probes().next().unwrap(), last_probe);
        assert_eq!(mapped_probe.gap.clone().unwrap().expected_gap_seq, "TCTCCAGAGAAGAGGAGGATGCGGCGGTGGCCATCTCCTGCTCGAAGTCCAGGGCGACGTAGCACAGCTTCTCCTTGATGTCGCGCACGATTTCCCGCTCGGCCGTGGTGGTGAAGCTGTAGCCTCGCTCAGTGAGGATCTTCATGAGGT".to_string());
        assert_eq!(
            mapped_probe.gap.clone().unwrap().gap_seq,
            "NNNNNNNTCTCCAGAGAAGAGGAGGATNNNNN".to_string()
        );
        assert!(!mapped_probe.gap.clone().unwrap().gap_within_max_error());
        Ok(())
    }

    #[test]
    fn test_mapped_gap_alignment_info() -> Result<()> {
        // Alignment:
        // gap:    TGGCTTACACTTTCAACTTG
        //         |||||\|||||\||||||
        // exp: AACTGGCTAACACTCTCAACT
        // Alignment Ops:
        // [Yclip(3), Match, Match, Match, Match, Match, Subst, Match, Match,
        //   Match, Match, Match, Subst, Match, Match, Match, Match, Match, Match, Xclip(2)]
        //
        let expected_gap_seq = "AACTGGCTAACACTCTCAACT".to_string();
        let gap_seq = "TGGCTTACACTTTCAACTGT".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"TGGCTTACACTTTCAACTGT", b"AACTGGCTAACACTCTCAACT", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        assert_eq!(mg.get_max_gap_error(), 6);
        assert!(!mg.gap_within_max_error());
        assert_eq!(mg.get_gap_levenshtein_distance(), 7);
        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 16);
        assert_eq!(mgi.get_num_mismatches(), 2);
        assert_eq!(mgi.get_num_insertions(), 2);
        assert_eq!(mgi.get_num_deletions(), 3);
        assert!(!mgi.ends_with_deletion());
        assert!(mgi.ends_with_insertion());
        assert!(mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // Alignment:
        // gap: CATTTTCTT---CCG
        //      |||||||||xxx|||
        // exp: CATTTTCTTCCACCG

        // Alignment ops:
        // [Match, Match, Match, Match, Match, Match, Match, Match, Match, Del, Del, Del, Match, Match, Match]
        let expected_gap_seq = "CATTTTCTTCCACCG".to_string();
        let gap_seq = "CATTTTCTTCCG".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"CATTTTCTTCCG", b"CATTTTCTTCCACCG", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        assert_eq!(mg.get_max_gap_error(), 4);
        assert!(mg.gap_within_max_error());
        assert_eq!(mg.get_gap_levenshtein_distance(), 3);
        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 12);
        assert_eq!(mgi.get_num_mismatches(), 0);
        assert_eq!(mgi.get_num_insertions(), 0);
        assert_eq!(mgi.get_num_deletions(), 3);
        assert!(!mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // Alignment:
        // gap: TAGTTTCCCCTTCA-
        //      ||||||||\|||||x
        // exp: TAGTTTCCACTTCAT

        // Alignment ops:
        // [Match, Match, Match, Match, Match, Match, Match, Match, Subst, Match, Match, Match, Match, Match, Del]
        let expected_gap_seq = "TAGTTTCCACTTCAT".to_string();
        let gap_seq = "TAGTTTCCCCTTCA".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"TAGTTTCCCCTTCA", b"TAGTTTCCACTTCAT", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 13);
        assert_eq!(mgi.get_num_mismatches(), 1);
        assert_eq!(mgi.get_num_insertions(), 0);
        assert_eq!(mgi.get_num_deletions(), 1);
        assert!(mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // Alignment:
        // gap: CTTCCTTC-CTT
        //      ||||||||x|||
        // exp: CTTCCTTCGCTTCTT

        // Alignment ops:
        // [Match, Match, Match, Match, Match, Match, Match, Match, Del, Match,
        //    Match, Match, Yclip(3)]
        let expected_gap_seq = "CTTCCTTCGCTTTTT".to_string();
        let gap_seq = "CTTCCTTCCTT".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"CTTCCTTCCTT", b"CTTCCTTCGCTTTTT", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 11);
        assert_eq!(mgi.get_num_mismatches(), 0);
        assert_eq!(mgi.get_num_insertions(), 0);
        assert_eq!(mgi.get_num_deletions(), 4);
        assert!(mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // Alignment:
        // gap: ACTGCTCAG---TCA
        //      |||||||||xxx\||
        // exp: ACTGCTCAGACCACA

        // Alignment ops:
        // [Match, Match, Match, Match, Match, Match, Match, Match, Match, Del,
        //    Del, Del, Subst, Match, Match]
        let expected_gap_seq = "ACTGCTCAGACCACA".to_string();
        let gap_seq = "ACTGCTCAGTCA".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"ACTGCTCAGTCA", b"ACTGCTCAGACCACA", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 11);
        assert_eq!(mgi.get_num_mismatches(), 1);
        assert_eq!(mgi.get_num_insertions(), 0);
        assert_eq!(mgi.get_num_deletions(), 3);
        assert!(!mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // Alignment:
        // gap: TCATATCTGCCTCAAA
        //      ||\||\||||\|+|||
        // exp: TCTTAGCTGCAT-AAA

        // Alignment ops:
        // [Match, Match, Subst, Match, Match, Subst, Match, Match, Match, Match,
        //   Subst, Match, Ins, Match, Match, Match]
        let expected_gap_seq = "TCTTAGCTGCATAAA".to_string();
        let gap_seq = "TCATATCTGCCTCAAA".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"TCATATCTGCCTCAAA", b"TCTTAGCTGCATAAA", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 12);
        assert_eq!(mgi.get_num_mismatches(), 3);
        assert_eq!(mgi.get_num_insertions(), 1);
        assert_eq!(mgi.get_num_deletions(), 0);
        assert!(!mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // Alignment:
        // gap: ATTCTCTTGAGACGTT
        //      \||||||||\|\||+|
        // exp: CTTCTCTTGGGCCG-T

        // Alignment ops:
        // [Subst, Match, Match, Match, Match, Match, Match, Match, Match, Subst,
        //   Match, Subst, Match, Match, Ins, Match]
        let expected_gap_seq = "CTTCTCTTGGGCCGT".to_string();
        let gap_seq = "ATTCTCTTGAGACGTT".to_string();
        let mg = MappedGap {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
        };
        println!(
            "Alignment: \n{}",
            mg.get_alignment()
                .pretty(b"ATTCTCTTGAGACGTT", b"CTTCTCTTGGGCCGT", 80),
        );
        println!("Alignment ops: \n{:?}", mg.get_alignment().operations);

        let mgi = MappedGapAlignmentInfo {
            gap_seq: gap_seq.clone(),
            expected_gap_seq: expected_gap_seq.clone(),
            alignment_operations: mg.get_alignment().operations.clone(),
        };
        assert_eq!(mgi.get_num_matches(), 12);
        assert_eq!(mgi.get_num_mismatches(), 3);
        assert_eq!(mgi.get_num_insertions(), 1);
        assert_eq!(mgi.get_num_deletions(), 0);
        assert!(!mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // assert!(false);
        Ok(())
    }

    #[test]
    fn test_gap_alignment() -> Result<()> {
        let probe_set = ProbeSetReference::from_path(
            &TargetSetFile::from("../cr_lib/test/target_panels/mock_probe_set_file_format_3.csv"),
            None,
            0,
        )?;

        let probes = probe_set.sorted_probes();
        assert_eq!(probes.len(), 87);

        let gapfill_probe = probes[85];
        let lhs_seq = "CTTCTCCAGGGAGGAGCTGGAAGCA";
        let rhs_seq = "TCTCTTGCTCGAAGTCCAGGGCGAC";
        let expected_gap_seq = "GCCGTGGCCA";
        assert_eq!(gapfill_probe.region, Some(ProbeRegion::Unspliced));
        assert_eq!(gapfill_probe.probe_type, ProbeType::PairedGapAlign);
        assert_eq!(gapfill_probe.ref_sequence_name, "not_a_real_ref");
        assert_eq!(gapfill_probe.ref_sequence_pos, Some(13));
        assert_eq!(gapfill_probe.cigar_string, "60M");
        let gapfill_probe_seq = probe_set.probe_to_seq.get(gapfill_probe).unwrap();
        assert_eq!(
            std::str::from_utf8(gapfill_probe_seq.lhs.as_ref()).unwrap(),
            lhs_seq
        );
        assert_eq!(
            std::str::from_utf8(gapfill_probe_seq.rhs.as_ref()).unwrap(),
            rhs_seq
        );
        assert_eq!(
            std::str::from_utf8(gapfill_probe_seq.gap.as_ref()).unwrap(),
            expected_gap_seq
        );

        let perfect_read = lhs_seq.to_string() + expected_gap_seq + rhs_seq + "AAACTGGCTGACTGAC";
        let mapped_probe = probe_set.align_probe_read(perfect_read.as_bytes());
        assert_eq!(mapped_probe.lhs_probe(), mapped_probe.rhs_probe());
        assert_eq!(mapped_probe.probes().next().unwrap(), gapfill_probe);
        assert_eq!(
            mapped_probe.gap.clone().unwrap().expected_gap_seq,
            expected_gap_seq.to_string()
        );
        assert_eq!(
            mapped_probe.gap.clone().unwrap().gap_seq,
            expected_gap_seq.to_string()
        );
        assert!(mapped_probe.gap.clone().unwrap().gap_within_max_error());

        // gap: GCCGTGGCC-
        //      |||||||||x
        // exp: GCCGTGGCCG
        // Alignment ops:
        //  [Match, Match, Match, Match, Match, Match, Match, Match, Match, Del]
        let read_with_del_end_of_gap = lhs_seq.to_string()
            + &expected_gap_seq[0..expected_gap_seq.len() - 1]
            + rhs_seq
            + "AAACTGGCTGACTGAC";
        let mapped_probe = probe_set.align_probe_read(read_with_del_end_of_gap.as_bytes());
        assert_eq!(mapped_probe.lhs_probe(), mapped_probe.rhs_probe());
        assert_eq!(mapped_probe.probes().next().unwrap(), gapfill_probe);
        assert_eq!(
            mapped_probe.gap.clone().unwrap().expected_gap_seq,
            expected_gap_seq.to_string()
        );
        assert_eq!(
            mapped_probe.gap.clone().unwrap().gap_seq,
            expected_gap_seq[0..expected_gap_seq.len() - 1].to_string()
        );
        assert!(mapped_probe.gap.clone().unwrap().gap_within_max_error()); // 1 base del is within 90% so it has expected gap
        println!(
            "Alignment: \n{}",
            mapped_probe.gap.clone().unwrap().get_alignment().pretty(
                b"GCCGTGGCC",
                b"GCCGTGGCCG",
                80
            ),
        );
        println!(
            "Alignment ops: \n{:?}",
            mapped_probe.gap.clone().unwrap().get_alignment().operations
        );
        let mgi = MappedGapAlignmentInfo {
            gap_seq: mapped_probe.gap.clone().unwrap().gap_seq,
            expected_gap_seq: mapped_probe.gap.clone().unwrap().expected_gap_seq,
            alignment_operations: mapped_probe
                .gap
                .clone()
                .unwrap()
                .get_alignment()
                .operations
                .clone(),
        };
        assert_eq!(mgi.get_num_matches(), 9);
        assert_eq!(mgi.get_num_mismatches(), 0);
        assert_eq!(mgi.get_num_insertions(), 0);
        assert_eq!(mgi.get_num_deletions(), 1);
        assert!(mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // gap:     TTGACCATT
        //          \||\|||
        // exp:  GCCGTGGCCA
        // Alignment ops:
        //  [Yclip(3), Subst, Match, Match, Subst, Match, Match, Match, Xclip(2)]
        let read_with_complex_alignment =
            lhs_seq.to_string() + "TTGACCAAC" + rhs_seq + "AAACTGGCTGACTGAC";
        let mapped_probe = probe_set.align_probe_read(read_with_complex_alignment.as_bytes());
        assert_eq!(mapped_probe.lhs_probe(), mapped_probe.rhs_probe());
        assert_eq!(mapped_probe.probes().next().unwrap(), gapfill_probe);
        assert_eq!(
            mapped_probe.gap.clone().unwrap().expected_gap_seq,
            expected_gap_seq.to_string()
        );
        assert_eq!(
            mapped_probe.gap.clone().unwrap().gap_seq,
            "TTGACCAAC".to_string()
        );
        assert!(!mapped_probe.gap.clone().unwrap().gap_within_max_error());
        println!(
            "Alignment: \n{}",
            mapped_probe.gap.clone().unwrap().get_alignment().pretty(
                b"TTGACCAAC",
                b"GCCGTGGCCA",
                80
            ),
        );
        println!(
            "Alignment ops: \n{:?}",
            mapped_probe.gap.clone().unwrap().get_alignment().operations
        );

        let mgi = MappedGapAlignmentInfo {
            gap_seq: mapped_probe.gap.clone().unwrap().gap_seq,
            expected_gap_seq: mapped_probe.gap.clone().unwrap().expected_gap_seq,
            alignment_operations: mapped_probe
                .gap
                .clone()
                .unwrap()
                .get_alignment()
                .operations
                .clone(),
        };

        assert_eq!(mgi.get_num_matches(), 5);
        assert_eq!(mgi.get_num_mismatches(), 2);
        assert_eq!(mgi.get_num_insertions(), 2);
        assert_eq!(mgi.get_num_deletions(), 3);
        assert!(!mgi.ends_with_deletion());
        assert!(mgi.ends_with_insertion());
        assert!(mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        // gap: GCCGTGG-CA
        //      |||||||x||
        // exp: GCCGTGGCCA
        // Alignment ops:
        //  [Match, Match, Match, Match, Match, Match, Match, Del, Match, Match]
        let read_with_complex_alignment =
            lhs_seq.to_string() + "GCCGTGGCA" + rhs_seq + "AAACTGGCTGACTGAC";
        let mapped_probe = probe_set.align_probe_read(read_with_complex_alignment.as_bytes());
        assert_eq!(mapped_probe.lhs_probe(), mapped_probe.rhs_probe());
        assert_eq!(mapped_probe.probes().next().unwrap(), gapfill_probe);
        assert_eq!(
            mapped_probe.gap.clone().unwrap().expected_gap_seq,
            expected_gap_seq.to_string()
        );
        assert_eq!(
            mapped_probe.gap.clone().unwrap().gap_seq,
            "GCCGTGGCA".to_string()
        );
        assert!(mapped_probe.gap.clone().unwrap().gap_within_max_error());
        println!(
            "Alignment: \n{}",
            mapped_probe.gap.clone().unwrap().get_alignment().pretty(
                b"GCCGTGGCA",
                b"GCCGTGGCCA",
                80
            ),
        );
        println!(
            "Alignment ops: \n{:?}",
            mapped_probe.gap.clone().unwrap().get_alignment().operations
        );

        let mgi = MappedGapAlignmentInfo {
            gap_seq: mapped_probe.gap.clone().unwrap().gap_seq,
            expected_gap_seq: mapped_probe.gap.clone().unwrap().expected_gap_seq,
            alignment_operations: mapped_probe
                .gap
                .clone()
                .unwrap()
                .get_alignment()
                .operations
                .clone(),
        };

        assert_eq!(mgi.get_num_matches(), 9);
        assert_eq!(mgi.get_num_mismatches(), 0);
        assert_eq!(mgi.get_num_insertions(), 0);
        assert_eq!(mgi.get_num_deletions(), 1);
        assert!(!mgi.ends_with_deletion());
        assert!(!mgi.ends_with_insertion());
        assert!(!mgi.starts_with_deletion());
        assert!(!mgi.starts_with_insertion());

        Ok(())
    }

    #[test]
    fn test_probe_set_merge_one_genome() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_a.csv"),
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_b.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_preserve_genome_input_order() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/mm10-fmt3-refv20_a.csv"),
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_a.csv"),
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_b.csv"),
                TargetSetFile::from("test/probe_set_merge/mm10-fmt3-refv20_b.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38-and-mm10_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge_two_genomes_same_version() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv20_a.csv"),
                TargetSetFile::from("test/probe_set_merge/mm10-fmt3-refv20_a.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38-and-mm10_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_preserve_record_input_order() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_a.csv"),
                TargetSetFile::from("test/probe_set_merge/mm10-fmt3-refv20_a.csv"),
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_b.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38-and-mm10_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge_two_genomes_same_version_format2() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt2-refv20_a.csv"),
                TargetSetFile::from("test/probe_set_merge/mm10-fmt2-refv20_a.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38-and-mm10_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge_two_genomes_different_version() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_a.csv"),
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_b.csv"),
                TargetSetFile::from("test/probe_set_merge/mm10-fmt3-refv20_a.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38-and-mm10_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge_diff_format() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt2-refv20_a.csv"),
                TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv20_a.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/set1.csv"),
                TargetSetFile::from("test/probe_set_merge/set2.csv"),
            ],
            &mut writer,
            None,
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge_different_headers() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/set1.csv"),
                TargetSetFile::from("test/probe_set_merge/additional_column.csv"),
            ],
            &mut writer,
            Some(Path::new("test/reference/GRCh38_ref_tiny")),
        )?;
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }

    #[test]
    fn test_probe_set_merge_error_same_genome_different_version() {
        let mut writer = Vec::new();
        assert!(
            merge_probe_set_csvs(
                &[
                    TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv20_a.csv"),
                    TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_a.csv"),
                    TargetSetFile::from("test/probe_set_merge/GRCh38-fmt3-refv24_b.csv"),
                ],
                &mut writer,
                None,
            )
            .is_err()
        );
    }

    #[test]
    fn test_probe_set_merge_deduplicate() -> Result<()> {
        let mut writer = Vec::new();
        merge_probe_set_csvs(
            &[
                TargetSetFile::from("test/probe_set_merge/set1.csv"),
                TargetSetFile::from("test/probe_set_merge/set1.csv"),
            ],
            &mut writer,
            None,
        )?;
        // Merging the same file with itself should result in the original file
        // due to probe deduplication, assuming that the original file was valid
        // and contained no duplicate probes itself.
        insta::assert_snapshot!(String::from_utf8(writer)?);
        Ok(())
    }
}
