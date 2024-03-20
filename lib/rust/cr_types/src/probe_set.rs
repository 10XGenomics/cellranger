use crate::rna_read::HIGH_CONF_MAPQ;
use crate::types::GenomeName;
use anyhow::{anyhow, Context, Result};
use itertools::{chain, Itertools};
use lazy_static::lazy_static;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::{min, Ord, Ordering, PartialOrd};
use std::collections::HashMap;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader};
use std::iter::zip;
use std::ops::Deref;
use std::path::Path;
use std::{fmt, iter};
use transcriptome::{Gene, Transcriptome};

/// One read half maps to a probe and the other half does not.
pub const MAPQ_HALF_MAPPED: u8 = 1;

/// Each read half maps to a different probe.
pub const MAPQ_SPLIT_MAPPED: u8 = 3;

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

// TODO: Handle constants better to make sure both python and
// rust code has access to them.
/// List of probe id prefixes to be ignored
const EXCLUDED_PROBE_ID_PREFIXES: [&str; 8] = [
    "DEPRECATED",
    "Hum-",
    "IGNORE",
    "INTERGENIC",
    "IR",
    "NC",
    "VAR",
    "VDJ",
];

/// Return true if this probe ID is excluded based on its probe ID.
pub fn is_deprecated_probe(probe_id: &str) -> bool {
    EXCLUDED_PROBE_ID_PREFIXES
        .iter()
        .any(|z| probe_id.starts_with(z))
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

/// A probe index and its gene.
#[derive(Clone, Debug, Eq, Deserialize, Serialize)]
pub struct Probe {
    /// The probe ID.
    pub probe_id: String,

    /// The probe gene ID and name.
    pub gene: Gene,

    /// An optional boolean to indicate if the probe is included
    pub included: bool,

    /// The region of the probe
    pub region: Option<ProbeRegion>,
}

impl Probe {
    /// Return true if this probe ID is excluded based on its probe ID or if included
    /// is false.
    pub fn is_excluded_probe(&self) -> bool {
        is_deprecated_probe(&self.probe_id) || !self.included
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
struct ProbeSequence {
    /// LHS half-probe sequence.
    lhs: Vec<u8>,

    /// RHS half-probe sequence.
    rhs: Vec<u8>,
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
    pub fn load_from(path: &Path) -> Result<Self> {
        Ok(Self(read_csv_comment_metadata(path)?))
    }

    /// Return true if the this metadata contains probe_set_file_format.
    pub fn is_probe_set_metadata(&self) -> bool {
        self.0.contains_key("probe_set_file_format")
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

    /// The reference genome name.
    pub genome: GenomeName,

    /// The length of a probe sequence.
    /// All probe sequences must have the same length.
    probe_seq_len: usize,

    /// The map of probes to probe sequences.
    probe_to_seq: ProbeToSeqMap,

    /// The map of probes to integer indices, sorted by probe ID.
    probe_id_to_index: ProbeIDToIndex,

    /// The map of left probe sequences to probes.
    lhs_seq_to_probe: SeqToProbeMap,

    /// The map of right probe sequences to probes.
    rhs_seq_to_probe: SeqToProbeMap,
}

impl ProbeSetReference {
    /// Align a half read to the probe set reference, allowing up to one mismatch.
    /// Return no match if there is more than one possible match.
    pub fn align_half_read<'a>(
        &self,
        seq_to_probe: &'a SeqToProbeMap,
        seq: &[u8],
    ) -> Option<(&'a [Probe], i32)> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

        if seq.len() < self.probe_seq_len / 2 {
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
        if seq.len() < self.probe_seq_len / 2 {
            // Read sequence is shorter than half the probe sequence.
            return MappedProbe::new(None, None, &self.genome);
        }

        // Skip the middle base if the sequence has odd length.
        let lhs_end = self.probe_seq_len / 2;
        let rhs_start = (self.probe_seq_len + 1) / 2;
        let rhs_end = min(seq.len(), self.probe_seq_len);
        let lhs_seq = &seq[..lhs_end];
        let rhs_seq = &seq[rhs_start..rhs_end];
        let (lhs_probes, lhs_score) = self
            .align_half_read(&self.lhs_seq_to_probe, lhs_seq)
            .unwrap_or((&[], 0));
        let (rhs_probes, rhs_score) = self
            .align_half_read(&self.rhs_seq_to_probe, rhs_seq)
            .unwrap_or((&[], 0));

        match (lhs_probes, rhs_probes) {
            ([], []) => MappedProbe::new(None, None, &self.genome),
            ([..], []) => self.rescue_rhs(lhs_probes, lhs_score, rhs_seq),
            ([], [..]) => self.rescue_lhs(lhs_seq, rhs_probes, rhs_score),
            ([lhs], [rhs]) => MappedProbe::new(
                Some(MappedProbeHalf::new(lhs.clone(), lhs_seq.len() as i32)),
                Some(MappedProbeHalf::new(rhs.clone(), rhs_seq.len() as i32)),
                &self.genome,
            ),
            ([..], [..]) => {
                // Confident matches are the intersection of lhs_probes and rhs_probes.
                // Multiple matches are caused by multiple probes with identical sequence,
                // in which case use the lexicographically minimal probe ID.
                if let Some(probe) = intersect_sorted_lists(lhs_probes, rhs_probes).next() {
                    MappedProbe::new(
                        Some(MappedProbeHalf::new(probe.clone(), lhs_seq.len() as i32)),
                        Some(MappedProbeHalf::new(probe.clone(), rhs_seq.len() as i32)),
                        &self.genome,
                    )
                } else {
                    MappedProbe::new(
                        Some(MappedProbeHalf::new(
                            lhs_probes[0].clone(),
                            lhs_seq.len() as i32,
                        )),
                        Some(MappedProbeHalf::new(
                            rhs_probes[0].clone(),
                            rhs_seq.len() as i32,
                        )),
                        &self.genome,
                    )
                }
            }
        }
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
            Some(MappedProbeHalf::new(mapped_probe.clone(), score))
        } else {
            None
        }
    }

    /// Rescue an unaligned half-probe sequence.
    fn rescue_unaligned_half(
        &self,
        mapped_probes: &[Probe],
        mapped_score: i32,
        read_seq: &[u8],
        rescue_which_half: impl Fn(&ProbeSequence) -> &[u8],
    ) -> (MappedProbeHalf, Option<MappedProbeHalf>) {
        let best = mapped_probes
            .iter()
            .rev() // Return the first probe with max score.
            .filter_map(|probe| {
                let probe_seq = rescue_which_half(&self.probe_to_seq[probe]);
                self.align_unaligned_half(probe, mapped_score, probe_seq, read_seq)
            })
            .max_by_key(|x| x.score);
        if let Some(ref x) = best {
            (MappedProbeHalf::new(x.probe.clone(), mapped_score), best)
        } else {
            (
                MappedProbeHalf::new(mapped_probes[0].clone(), mapped_score),
                None,
            )
        }
    }

    /// Rescue an unaligned LHS sequence.
    fn rescue_lhs(&self, lhs_read_seq: &[u8], rhs_probes: &[Probe], rhs_score: i32) -> MappedProbe {
        let (rhs, opt_lhs) =
            self.rescue_unaligned_half(rhs_probes, rhs_score, lhs_read_seq, |x| &x.lhs);
        MappedProbe::new(opt_lhs, Some(rhs), &self.genome)
    }

    /// Rescue an unaligned RHS sequence.
    fn rescue_rhs(&self, lhs_probes: &[Probe], lhs_score: i32, rhs_read_seq: &[u8]) -> MappedProbe {
        let (lhs, opt_rhs) =
            self.rescue_unaligned_half(lhs_probes, lhs_score, rhs_read_seq, |x| &x.rhs);
        MappedProbe::new(Some(lhs), opt_rhs, &self.genome)
    }

    /// Construct a probe set refrence from a probe set CSV file.
    /// CSV header must include #probe_set_file_format.
    /// CSV header row must be gene_id,probe_seq,probe_id,included,region.
    pub fn from_path(
        target_set: &Path,
        reference_path: &Path,
        transcriptome_min_score: usize,
    ) -> Result<Self> {
        // Read the transcriptome GTF to map gene IDs to gene names.
        let gene_id_to_name: TxHashMap<_, _> = Transcriptome::from_reference_path(reference_path)?
            .genes
            .into_iter()
            .map(|x| (x.id, x.name))
            .collect();

        let metadata = ProbeSetReferenceMetadata::load_from(target_set)?;
        assert!(metadata.is_probe_set_metadata());

        // Read the probe set CSV file.
        let mut reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .from_path(target_set)
            .with_context(|| target_set.display().to_string())?;

        // Ensure that the headers are correct.
        // bait_seq and bait_id are permitted for backward compatibility.
        let header: Vec<_> = reader.headers().unwrap().iter().collect();
        assert_eq!(header[0], "gene_id");
        assert!(header[1] == "probe_seq" || header[1] == "bait_seq");
        assert!(header[2] == "probe_id" || header[2] == "bait_id");
        if let Some(&included_header) = header.get(3) {
            assert_eq!(included_header, "included");
        }

        // Parse the probe set CSV file.
        let mut probe_seq_len = None;
        let probe_to_seq: ProbeToSeqMap = reader
            .records()
            .map(|record| {
                let record = record.unwrap();
                let gene_id = record[0].to_string();
                let probe_seq = record[1].as_bytes();
                let probe_id = record[2].to_string();
                let included = record
                    .get(3)
                    .map_or(true, |x| x.to_lowercase().parse().unwrap());
                let region = record
                    .get(4)
                    .map(str::to_string)
                    .map(|r| ProbeRegion::new(&r));
                let gene_name = gene_id_to_name.get(&gene_id).unwrap_or(&gene_id).clone();
                let probe = Probe {
                    probe_id,
                    gene: Gene {
                        id: gene_id,
                        name: gene_name,
                    },
                    included,
                    region,
                };

                // Ensure that the length of the probe sequences is fixed.
                if let Some(probe_seq_len) = probe_seq_len {
                    assert_eq!(probe_seq.len(), probe_seq_len);
                } else {
                    probe_seq_len = Some(probe_seq.len());
                }

                // Skip the middle base if the sequence has odd length.
                let rhs_end = probe_seq.len() / 2;
                let lhs_start = (probe_seq.len() + 1) / 2;
                (
                    probe,
                    ProbeSequence {
                        lhs: probe_seq[..rhs_end].to_vec(),
                        rhs: probe_seq[lhs_start..].to_vec(),
                    },
                )
            })
            .collect();
        let probe_seq_len = probe_seq_len
            .unwrap_or_else(|| panic!("target_set CSV has no records: {}", target_set.display()));

        // Construct a map of probes to integer indices, sorted by probe ID.
        let probe_id_to_index: TxHashMap<_, _> = probe_to_seq
            .keys()
            .sorted()
            .enumerate()
            .map(|(i, probe)| (probe.probe_id.clone(), i))
            .collect();

        // Construct a map of probe sequences to probes.
        let mut lhs_seq_to_probe: SeqToProbeMap =
            probe_to_seq
                .iter()
                .fold(Default::default(), |mut map, (probe, seq)| {
                    map.entry(seq.lhs.clone()).or_default().push(probe.clone());
                    map
                });
        let mut rhs_seq_to_probe: SeqToProbeMap =
            probe_to_seq
                .iter()
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

        let genome = metadata.0["reference_genome"].as_str().into();
        Ok(Self {
            metadata,
            transcriptome_min_score,
            genome,
            probe_seq_len,
            probe_to_seq,
            probe_id_to_index,
            lhs_seq_to_probe,
            rhs_seq_to_probe,
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
}

impl MappedProbeHalf {
    /// Return a new MappedProbeHalf.
    pub fn new(probe: Probe, score: i32) -> Self {
        assert!(score > 0);
        MappedProbeHalf { probe, score }
    }
}

/// The result of mapping a read to a probe set reference.
#[derive(Deserialize, Serialize)]
pub struct MappedProbe {
    /// The result of mapping the left sequence.
    lhs: Option<MappedProbeHalf>,

    /// The result of mapping the right sequence.
    rhs: Option<MappedProbeHalf>,

    /// The name of the genome.
    genome: Option<GenomeName>,

    /// The read maps confidently to a probe but multimapped with STAR.
    is_rescued: bool,
}

impl MappedProbe {
    /// Create a new mapping result.
    pub fn new(
        lhs: Option<MappedProbeHalf>,
        rhs: Option<MappedProbeHalf>,
        genome: &GenomeName,
    ) -> Self {
        let genome = (lhs.is_some() || rhs.is_some()).then(|| genome.clone());
        MappedProbe {
            lhs,
            rhs,
            genome,
            is_rescued: false,
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

    /// Return the reference genome name.
    pub fn genome(&self) -> Option<&GenomeName> {
        self.genome.as_ref()
    }

    /// Return whether either left or right sequences mapped to a probe.
    pub fn is_mapped(&self) -> bool {
        self.lhs.is_some() || self.rhs.is_some()
    }

    /// Return whether both left and right sequences mapped to the same probe.
    pub fn is_conf_mapped(&self) -> bool {
        self.is_mapped() && self.lhs_probe() == self.rhs_probe()
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
                gene: Gene {
                    id: String::new(),
                    name: String::new(),
                },
                included: true,
                region: None,
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
            (Some(_), Some(_)) if self.is_conf_mapped() => HIGH_CONF_MAPQ,
            (Some(_), Some(_)) => MAPQ_SPLIT_MAPPED,
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
}

#[cfg(test)]
mod test {
    use super::*;
    use test_refdata::{refdata_available, refdata_path};

    #[test]
    fn test_intersect_sorted_lists() {
        assert!(intersect_sorted_lists(&[0, 2, 4, 6, 8], &[1, 3, 5, 7, 9])
            .next()
            .is_none());
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
        assert!(!ProbeSetReferenceMetadata::load_from(Path::new(
            "../cr_lib/test/target_panels/Immunology_targeting_hybrid.csv"
        ))?
        .is_probe_set_metadata());
        assert!(ProbeSetReferenceMetadata::load_from(Path::new(
            "../cr_lib/test/target_panels/Immunology_targeting_templated_ligation.csv"
        ))?
        .is_probe_set_metadata());
        Ok(())
    }

    #[test]
    fn test_from_path() -> Result<()> {
        if !refdata_available() {
            return Ok(());
        }
        let reference = refdata_path("GRCh38-2020-A/");

        ProbeSetReference::from_path(
            Path::new("../cr_lib/test/target_panels/Immunology_targeting_templated_ligation.csv"),
            &reference,
            0,
        )?;
        Ok(())
    }
}
