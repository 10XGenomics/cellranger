//! Identifies canonical adapter configurations in long-read sequencing data
//! and uses flanking adapters to extract 10x library constructs.
use crate::WhichRead;
use crate::adapter_trimmer::{Adapter, AdapterTrimmer};
use crate::adapters::AdapterKind;
use crate::read_pair::{
    MAX_READ_HEAD, MAX_READ_LEN, MutReadPair, PushReadError, ReadHeader, ReadPair,
    ReadPairR1OnlyWrapper, ReadPairStorage, ReadPart,
};
use anyhow::{Result, anyhow, bail};
use bio::alphabets::dna;
use bytes::BytesMut;
use fastq::{OwnedRecord, Record};
use martian_derive::{MartianStruct, MartianType};
use memchr::memmem;
use serde::{Deserialize, Serialize};
use std::fmt::Write;
use std::ops::Range;
use std::str::FromStr;

const MIN_ADAPTER_LENGTH: usize = 10;

/// Default maximum R2 length applied by [`InsertDefinition::to_read_pair`].
/// Required due to the RpRange compact representation of a read.
pub const DEFAULT_MAX_R2_LENGTH: usize = 32_000;

/// Represents the physical spatial arrangement of adapters in a 10x Genomics library.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, MartianStruct)]
pub struct AdapterLayout {
    /// Required 5' termninal adapter (adjacent to gel-bead barcode)
    pub terminal_5p: AdapterKind,
    /// Optional constant sequence affixing gel-bead oligo to probe barcode
    pub probe_prefix: Option<AdapterKind>,
    /// The sequence used to capture the target (i.e cDNA insert, probes, FB).
    /// `None` when using [`ReadSplittingStrategy::FixedOffset`].
    pub target_prefix: Option<AdapterKind>,
    /// Required 3' terminal adapter (opposite end from gel-bead barcode)
    pub terminal_3p: AdapterKind,
}

impl AdapterLayout {
    /// Returns true if the 3' terminal adapter is too short for a global search.
    pub fn terminal_3p_short_search(&self) -> bool {
        self.terminal_3p == AdapterKind::TruncatedR2
            || self.terminal_3p == AdapterKind::TruncatedR2Rc
    }
}

/// Controls how the R1/R2 split boundary is located within an insert.
/// The target_prefix adapter (if any) is not included in R1.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, MartianType)]
#[serde(tag = "type", content = "offset")]
pub enum ReadSplittingStrategy {
    /// The boundary is inferred from the observed position of the `target_prefix`
    /// adapter in the read. The `target_prefix` adapter is included in the global adapter scan
    /// to locate the split, but the split is placed at the adapter's start position.
    AdapterLocated {},
    /// The boundary is at a fixed offset (BC+UMI) from the end of the `terminal_5p` adapter.
    /// The `target_prefix` adapter (if any) is NOT included in the global scan.
    FixedOffset {
        /// The number of bases after the end of the `terminal_5p` adapter where the split occurs.
        offset: usize,
    },
}

/// A configuration that combines physical [AdapterLayout] with the
/// logic required to extract relevant sequence parts.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, MartianStruct)]
pub struct TrimmerConfig {
    /// Represents the physical spatial arrangement of adapters.
    pub layout: AdapterLayout,
    /// Determines whether to search for the intermediate adapter or use a fixed offset.
    pub splitting_strategy: ReadSplittingStrategy,
    /// If true, allow matches that ignore the final terminal adapter.
    pub allow_truncated_constructs: bool,
    /// Indicates if the library is directional.
    /// For directional libraries, we only consider the forward orientation as searchable.
    pub is_directional: bool,
    /// Aligner-driven maximum R2 length.
    pub max_r2_length: usize,
}

impl TrimmerConfig {
    /// Returns the ordered list of adapters that are searchable.
    /// For [`ReadSplittingStrategy::FixedOffset`] libraries we exclude the `target_prefix`
    /// because the split position is determined by offset, not by adapter detection.
    fn searchable_adapters(&self) -> Vec<AdapterKind> {
        let mut adapters = vec![self.layout.terminal_5p];
        if let Some(probe) = self.layout.probe_prefix {
            adapters.push(probe);
        }
        match self.splitting_strategy {
            ReadSplittingStrategy::AdapterLocated {} => {
                adapters.push(self.layout.target_prefix.expect(
                    "target_prefix must be defined with AdapterLocated splitting strategy",
                ));
            }
            ReadSplittingStrategy::FixedOffset { .. } => {
                // target_prefix is not included for FixedOffset strategy
            }
        }
        adapters.push(self.layout.terminal_3p);
        adapters
    }

    /// Generates the search targets in reverse-complement order for anti-sense reads.
    fn searchable_adapters_revcomp(&self) -> Option<Vec<AdapterKind>> {
        if self.is_directional {
            return None;
        }
        let mut adapters = Vec::new();
        for adapter_kind in self.searchable_adapters() {
            adapters.push(adapter_kind.get_revcomp_adapter());
        }
        adapters.reverse();
        Some(adapters)
    }

    fn build_stranded_adapters(&self, forward: bool) -> Vec<Adapter> {
        let mut adapters = vec![self.layout.terminal_5p.build_5p(forward)];
        if let Some(probe) = self.layout.probe_prefix {
            adapters.push(probe.build_5p(forward));
        }
        match self.splitting_strategy {
            ReadSplittingStrategy::AdapterLocated {} => {
                adapters.push(
                    self.layout
                        .target_prefix
                        .expect(
                            "target_prefix must be defined with AdapterLocated splitting strategy",
                        )
                        .build_5p(forward),
                );
            }
            ReadSplittingStrategy::FixedOffset { .. } => {
                // target_prefix is not included for FixedOffset strategy
            }
        }
        if self.layout.terminal_3p != AdapterKind::End && !self.layout.terminal_3p_short_search() {
            adapters.push(self.layout.terminal_3p.build_3p(forward));
        }
        adapters
    }

    /// Convert this configuration into a long read trimmer.
    pub fn into_trimmer(self) -> LongReadTrimmer {
        LongReadTrimmer::new(AdapterSearchStrategy {
            forward: self.searchable_adapters(),
            reverse: self.searchable_adapters_revcomp(),
            adapters: {
                let mut all = self.build_stranded_adapters(true);
                if !self.is_directional {
                    all.extend(self.build_stranded_adapters(false));
                }
                all
            },
            lib_config: self,
        })
    }
}

/// Canonical adapter configuration for a specific 10x chemistry.
#[derive(Clone, Debug)]
struct AdapterSearchStrategy {
    /// Forward orientation (gel-bead barcode adapter first)
    forward: Vec<AdapterKind>,
    /// Reverse orientation (gel-bead barcode adapter last),
    /// only considered for non-directional libraries
    reverse: Option<Vec<AdapterKind>>,
    /// All internal adapters in forward orientations plus reverse for non-directional libs.
    adapters: Vec<Adapter>,
    lib_config: TrimmerConfig,
}

impl AdapterSearchStrategy {
    fn terminal_5p_adapter(&self) -> &AdapterKind {
        self.forward.first().unwrap()
    }
    fn contains(&self, adapters: &[AdapterKind]) -> bool {
        if adapters == self.forward {
            return true;
        }
        if let Some(rev_adapters) = &self.reverse
            && adapters == rev_adapters
        {
            return true;
        }
        if self.lib_config.allow_truncated_constructs
            || self.lib_config.layout.terminal_3p_short_search()
        {
            let mut truncated_forward = self.forward.clone();
            *truncated_forward.last_mut().unwrap() = AdapterKind::End;
            if adapters == truncated_forward {
                return true;
            }
        }
        false
    }
    fn target_prefix_adapter(&self) -> Option<&AdapterKind> {
        self.lib_config.layout.target_prefix.as_ref()
    }
    fn contains_truncated(&self, adapters: &[AdapterKind]) -> bool {
        let mut truncated_forward = self.forward.clone();
        *truncated_forward.last_mut().unwrap() = AdapterKind::End;

        adapters.ends_with(&truncated_forward)
    }
}

/// Long read adapter trimmer.
#[derive(Clone)]
pub struct LongReadTrimmer {
    strategy: AdapterSearchStrategy,
    trimmers: Vec<AdapterTrimmer>,
}

impl LongReadTrimmer {
    /// Initialize a long read trimmer with the provided canonical adapter configs.
    fn new(strategy: AdapterSearchStrategy) -> Self {
        let trimmers = strategy
            .adapters
            .iter()
            .cloned()
            .map(AdapterTrimmer::new)
            .collect();
        LongReadTrimmer { strategy, trimmers }
    }

    fn trim_adapter_recursive(
        read: &[u8],
        trimmer: &AdapterTrimmer,
        adapter_pos: &mut Vec<(AdapterKind, Range<usize>)>,
    ) {
        const MAX_LIMIT: usize = 100; // Max number of adapter searches in a read
        let mut start_pos = 0;
        let mut search_count = 0;

        while search_count < MAX_LIMIT {
            let seq = &read[start_pos..];
            if let Some(trim_result) = trimmer.find(seq) {
                adapter_pos.push((
                    AdapterKind::from_str(trimmer.adapter_name()).unwrap(),
                    Range {
                        start: trim_result.adapter_range.start + start_pos,
                        end: trim_result.adapter_range.end + start_pos,
                    },
                ));
                start_pos += trim_result.adapter_range.end;
                search_count += 1;
            } else {
                break;
            }
        }
    }

    fn find_adapters(&self, read: &[u8]) -> Vec<(AdapterKind, Range<usize>)> {
        self.trimmers.iter().fold(Vec::default(), |mut acc, t| {
            Self::trim_adapter_recursive(read, t, &mut acc);
            acc
        })
    }

    /// Returns the position of the rightmost match (first found from the 3' end),
    /// or `None` if no exact match exists.
    fn reverse_scan_terminal_adapter(
        read: &[u8],
        search_region: Range<usize>,
        adapter_kind: AdapterKind,
    ) -> Option<(AdapterKind, Range<usize>)> {
        let needle = adapter_kind.sequence().as_bytes();
        let haystack = &read[search_region.clone()];
        let offset = memmem::rfind(haystack, needle)?;
        let abs_start = search_region.start + offset;
        Some((adapter_kind, abs_start..abs_start + needle.len()))
    }

    /// Detect canonical insert configurations in a read.
    pub fn boundary_detection(&self, read: &[u8]) -> (Vec<InsertDefinition>, Vec<AdapterKind>) {
        let adapters = {
            let mut a = self.find_adapters(read);
            a.sort_by_key(|(_, range)| range.start);

            // Add a symbolic End adapter to enable extracting incomplete inserts,
            // when sequencing does not reach the end of the construct.
            a.insert(0, (AdapterKind::End, 0..0));
            a.push((AdapterKind::End, read.len()..read.len()));
            a
        };

        let configs: Vec<InsertDefinition> = adapters
            .windows(self.strategy.forward.len())
            .filter_map(|window| {
                let adapter_kinds: Vec<AdapterKind> =
                    window.iter().map(|(kind, _)| *kind).collect();
                let orientation =
                    if self.strategy.terminal_5p_adapter() == &window.first().unwrap().0 {
                        Orientation::Forward
                    } else {
                        Orientation::Reverse
                    };
                let boundary = BoundaryDefinition::new(
                    window.first().unwrap().1.clone(),
                    window.last().unwrap().1.clone(),
                    orientation,
                )
                .ok()?;
                if !self.strategy.contains(&adapter_kinds) {
                    return None;
                }

                let intermediate = match &self.strategy.lib_config.splitting_strategy {
                    ReadSplittingStrategy::AdapterLocated {} => {
                        let target_prefix = self.strategy.target_prefix_adapter()?;
                        window
                            .iter()
                            .find(|(kind, _)| {
                                *kind == *target_prefix
                                    || *kind == target_prefix.get_revcomp_adapter()
                            })
                            .and_then(|(k, r)| {
                                if boundary.first.end < r.start && r.end < boundary.last.start {
                                    Some((k, r.clone()))
                                } else {
                                    None
                                }
                            })
                            .map(|(k, r)| (*k, r.clone()))
                    }
                    ReadSplittingStrategy::FixedOffset { offset } => {
                        boundary.infer_split_at_offset(*offset)
                    }
                }?;

                let mut last_adapter = window.last().unwrap().clone();

                // For short terminal adapters that were excluded from the global scan,
                // resolve the symbolic End boundary via a right-to-left exact-match search.
                if self.strategy.lib_config.layout.terminal_3p_short_search()
                    && orientation == Orientation::Forward
                {
                    assert!(window.last().unwrap().0 == AdapterKind::End);
                    let search_region = intermediate.1.end..read.len();
                    last_adapter = Self::reverse_scan_terminal_adapter(
                        read,
                        search_region,
                        self.strategy.lib_config.layout.terminal_3p,
                    )?;
                }

                Some(InsertDefinition {
                    first_adapter: window.first().unwrap().clone(),
                    last_adapter,
                    rna_prefix: intermediate,
                    orientation,
                    idx: 0, // to be filled in below
                })
            })
            .enumerate()
            .map(|(final_idx, mut insert_def)| {
                insert_def.idx = final_idx;
                insert_def
            })
            .collect();

        let raw_adapters = adapters.iter().map(|(kind, _)| *kind).collect();

        (configs, raw_adapters)
    }

    /// Trim a single-ended ReadPair into a pseudo-R1-R2 read.
    ///
    /// Return Ok(None) if the read cannot be trimmed.
    pub fn trim_read(
        &self,
        read: ReadPair,
        header_handling: ReadHeader,
        visitor: &mut impl TrimmedReadVisitor,
    ) -> Result<Option<ReadPair>> {
        let seq = read
            .get(WhichRead::R1, ReadPart::Seq)
            .ok_or_else(|| anyhow!("read has no R1 sequence"))?;
        // TODO: decide if we want to use more than the first config or not
        let (configs, raw_adapter_config) = self.boundary_detection(seq);
        let config = configs.into_iter().next();
        let is_truncated = if let Some(config) = &config {
            config.is_truncated()
        } else if self.strategy.lib_config.allow_truncated_constructs {
            // If truncated constructs are allowed then a truncated read would have generated a config,
            // so if we didn't find a config, we can infer that the read is not truncated.
            false
        } else {
            match self.strategy.lib_config.splitting_strategy {
                ReadSplittingStrategy::AdapterLocated {} => {
                    self.strategy.contains_truncated(&raw_adapter_config)
                }
                ReadSplittingStrategy::FixedOffset { .. } => todo!(
                    "Cannot infer truncated status when allow_truncated_constructs is false and the location of the R1/R2 split is inferred based on a fixed offset"
                ),
            }
        };
        visitor.visit_trimmed_read(config.as_ref(), is_truncated);
        Ok(config
            .and_then(|config| config.to_read_pair(&ReadPairR1OnlyWrapper(&read), header_handling)))
    }
}

/// A no-op visitor that can be used when metrics collection is not needed.
pub struct NoOpVisitor;
impl TrimmedReadVisitor for NoOpVisitor {
    fn visit_trimmed_read(&mut self, _configs: Option<&InsertDefinition>, _is_truncated: bool) {
        // no-op
    }
}

/// Visitor over each read trimmed by a long read trimmer.
pub trait TrimmedReadVisitor {
    /// Called once per read, with the parsed insert definition (if any) for the selected
    /// adapter configuration and a flag indicating whether that configuration is truncated.
    fn visit_trimmed_read(&mut self, configs: Option<&InsertDefinition>, is_truncated: bool);
}

/// Orientation of the insert within the original read, inferred from adapter arrangement.
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Orientation {
    /// Forward orientation
    Forward,
    /// Reverse orientation
    Reverse,
}

struct BoundaryDefinition {
    first: Range<usize>,
    last: Range<usize>,
    orientation: Orientation,
}

impl BoundaryDefinition {
    fn new(first: Range<usize>, last: Range<usize>, orientation: Orientation) -> Result<Self> {
        if last.end <= first.end {
            bail!(
                "Invalid adapter definition: last adapter ({}..{}) must end after first adapter ({}..{})",
                last.start,
                last.end,
                first.start,
                first.end
            );
        }
        if (orientation == Orientation::Forward && first.len() < MIN_ADAPTER_LENGTH)
            || (orientation == Orientation::Reverse && last.len() < MIN_ADAPTER_LENGTH)
        {
            bail!(
                "Invalid adapter definition: Adapter sequence in gel-beads must be at least {MIN_ADAPTER_LENGTH} bases"
            );
        }
        Ok(BoundaryDefinition {
            first,
            last,
            orientation,
        })
    }
    /// Returns the range of the read to retain after trimming terminal adapters.
    fn get_insert(&self) -> Range<usize> {
        match self.orientation == Orientation::Forward {
            true => Range {
                start: self.first.end,
                end: self.last.start,
            },
            false => Range {
                start: self.first.end,
                end: self.last.start,
            },
        }
    }

    /// Returns a zero-length split marker at `terminal_5p.end + offset`.
    fn infer_split_at_offset(&self, offset: usize) -> Option<(AdapterKind, Range<usize>)> {
        let split = if self.orientation == Orientation::Forward {
            self.first.end + offset
        } else {
            self.last.start.checked_sub(offset)?
        };

        // validate that the inferred split position is within the insert boundaries
        let insert = self.get_insert();
        if split < insert.start || split > insert.end {
            return None;
        }
        Some((AdapterKind::End, split..split))
    }
}

/// Represents a physical segment of within an insert.
struct SegmentDefinition {
    range: Range<usize>,
    needs_revcomp: bool,
}

impl SegmentDefinition {
    // Create a OwnedRecord from the read part definition.
    fn to_record<R: Record>(&self, record: &R, head: String) -> OwnedRecord {
        assert!(self.range.len() <= MAX_READ_LEN);
        if self.needs_revcomp {
            OwnedRecord {
                head: head.into(),
                seq: dna::revcomp(&record.seq()[self.range.clone()]),
                qual: record.qual()[self.range.clone()]
                    .iter()
                    .copied()
                    .rev()
                    .collect::<Vec<_>>(),
                sep: None,
            }
        } else {
            OwnedRecord {
                head: head.into(),
                seq: record.seq()[self.range.clone()].to_vec(),
                qual: record.qual()[self.range.clone()].to_vec(),
                sep: None,
            }
        }
    }
}

/// Represents the valid terminal adapter pairs that flank a 10x library sequence,
/// in either forward or reverse orientation.
#[derive(Debug, PartialEq, Eq)]
pub struct InsertDefinition {
    /// First terminal adapter, location coordinates are relative to original read.
    pub first_adapter: (AdapterKind, Range<usize>),
    /// Last terminal adapter, location coordinates are relative to original read.
    pub last_adapter: (AdapterKind, Range<usize>),
    /// RNA adapter, location coordinates are relative to original read.
    pub rna_prefix: (AdapterKind, Range<usize>),
    /// Orientation of the insert within the original read.
    pub orientation: Orientation,
    /// Index of this insert within the read.
    pub idx: usize,
}

impl InsertDefinition {
    /// Add a numeric suffix to the header ID.
    pub fn reformat_header_id(&self, original_id: &str) -> String {
        let mut head = String::new();
        for (i, header_entry) in original_id.split(' ').enumerate() {
            if i > 0 {
                head.push(' ');
            }
            head.push_str(header_entry);
            if i == 0 {
                write!(&mut head, "_{}", self.idx).unwrap();
            }
        }
        assert!(head.len() <= MAX_READ_HEAD);
        head
    }
    /// Identifies the two physical halves of the insert.
    /// Rna prefix is included in the second half.
    fn get_parts(&self) -> (SegmentDefinition, SegmentDefinition) {
        match self.orientation {
            Orientation::Forward => (
                SegmentDefinition {
                    range: self.first_adapter.1.end..self.rna_prefix.1.start,
                    needs_revcomp: false,
                },
                SegmentDefinition {
                    range: self.rna_prefix.1.start..self.last_adapter.1.start,
                    needs_revcomp: true,
                },
            ),
            Orientation::Reverse => (
                SegmentDefinition {
                    range: self.rna_prefix.1.end..self.last_adapter.1.start,
                    needs_revcomp: true,
                },
                SegmentDefinition {
                    range: self.first_adapter.1.end..self.rna_prefix.1.end,
                    needs_revcomp: false,
                },
            ),
        }
    }
    /// Create a ReadPair from this insert definition.
    pub fn to_read_pair<R: Record>(
        &self,
        record: &R,
        header_handling: ReadHeader,
    ) -> Option<ReadPair> {
        let head = if header_handling.skip_id() {
            String::new()
        } else {
            self.reformat_header_id(str::from_utf8(record.head()).unwrap())
        };
        let (part_a, part_b) = self.get_parts();
        let r1 = part_a.to_record(record, head.clone());
        let r2 = part_b.to_record(record, head.clone());

        let mut buffer = BytesMut::with_capacity(
            r1.seq().len() + r1.qual().len() + r2.seq().len() + r2.qual().len() + 2 * head.len(),
        );
        let mut rp = MutReadPair::empty(&mut buffer).storage(ReadPairStorage::SharedBuffer);
        rp.push_read(&r1, WhichRead::R1)
            .map_err(PushReadError::basic_id_mismatch)
            .unwrap();
        rp.push_read(&r2, WhichRead::R2)
            .map_err(PushReadError::basic_id_mismatch)
            .unwrap();
        Some(rp.freeze(header_handling))
    }
    /// Returns the range to retain the insert excluding flanking adapter sequences.
    pub fn get_retain_range(&self) -> Range<usize> {
        match self.orientation == Orientation::Forward {
            true => Range {
                start: self.first_adapter.1.end,
                end: self.last_adapter.1.start,
            },
            false => Range {
                start: self.first_adapter.1.end,
                end: self.last_adapter.1.start,
            },
        }
    }

    /// Returns true if this insert is truncated.
    pub fn is_truncated(&self) -> bool {
        self.last_adapter.0 == AdapterKind::End && self.orientation == Orientation::Forward
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_pair::{ReadPart, WhichRead};

    const RT_FRAG: TrimmerConfig = TrimmerConfig {
        layout: AdapterLayout {
            terminal_5p: AdapterKind::R1,
            probe_prefix: None,
            target_prefix: None,
            terminal_3p: AdapterKind::R2Rc,
        },
        splitting_strategy: ReadSplittingStrategy::FixedOffset { offset: 4 },
        allow_truncated_constructs: false,
        is_directional: false,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    const RT_FRAG_DIRECTIONAL: TrimmerConfig = TrimmerConfig {
        layout: AdapterLayout {
            terminal_5p: AdapterKind::R1,
            probe_prefix: None,
            target_prefix: None,
            terminal_3p: AdapterKind::R2Rc,
        },
        splitting_strategy: ReadSplittingStrategy::FixedOffset { offset: 4 },
        allow_truncated_constructs: false,
        is_directional: true,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    const RT_FRAG_NO_TERMINAL_3P: TrimmerConfig = TrimmerConfig {
        layout: AdapterLayout {
            terminal_5p: AdapterKind::R1,
            probe_prefix: None,
            target_prefix: Some(AdapterKind::PolyT),
            terminal_3p: AdapterKind::End,
        },
        splitting_strategy: ReadSplittingStrategy::FixedOffset { offset: 4 },
        allow_truncated_constructs: true,
        is_directional: false,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    const RTL_GEX: TrimmerConfig = TrimmerConfig {
        layout: AdapterLayout {
            terminal_5p: AdapterKind::R1,
            probe_prefix: None,
            target_prefix: Some(AdapterKind::Spacer),
            terminal_3p: AdapterKind::SmallRnaR2Rc,
        },
        splitting_strategy: ReadSplittingStrategy::AdapterLocated {},
        allow_truncated_constructs: false,
        is_directional: true,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    const RTL_CRISPR: TrimmerConfig = TrimmerConfig {
        layout: AdapterLayout {
            terminal_5p: AdapterKind::R1Bp12,
            probe_prefix: Some(AdapterKind::X12Capture),
            target_prefix: Some(AdapterKind::Spacer),
            terminal_3p: AdapterKind::TruncatedR2Rc,
        },
        splitting_strategy: ReadSplittingStrategy::AdapterLocated {},
        allow_truncated_constructs: false,
        is_directional: true,
        max_r2_length: DEFAULT_MAX_R2_LENGTH,
    };

    #[test]
    fn test_short_terminal_3p() {
        assert!(!RT_FRAG.layout.terminal_3p_short_search());
        assert!(!RT_FRAG_DIRECTIONAL.layout.terminal_3p_short_search());
        assert!(!RT_FRAG_NO_TERMINAL_3P.layout.terminal_3p_short_search());
        assert!(!RTL_GEX.layout.terminal_3p_short_search());
        assert!(RTL_CRISPR.layout.terminal_3p_short_search());
    }

    fn clean_read(read: &[u8]) -> Vec<u8> {
        read.iter().filter(|&&b| b != b'|').copied().collect()
    }

    fn trimmer(lib_config: TrimmerConfig) -> LongReadTrimmer {
        lib_config.into_trimmer()
    }

    fn find_configs(read: &[u8], lib_config: TrimmerConfig) -> (Vec<InsertDefinition>, String) {
        let (configs, raw_adapters) = trimmer(lib_config).boundary_detection(read);
        let raw_adapter_config = raw_adapters
            .iter()
            .map(|a| format!("{a}"))
            .collect::<Vec<_>>()
            .join(" - ");
        (configs, raw_adapter_config)
    }

    #[test]
    fn test_trim_adapter_recursive_threeprime() {
        // (4 bases)(r2_rc 34 bases)(8 bases)(r2_rc 34 bases)(6 bases)
        let read = b"ATGC|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC|ATGCATGC|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC|ATGCAT";
        let result = trimmer(RT_FRAG).find_adapters(&clean_read(read));
        assert_eq!(result[0], (AdapterKind::R2Rc, 4..38));
        assert_eq!(result[1], (AdapterKind::R2Rc, 46..80));
    }

    #[test]
    fn test_trim_adapter_recursive_fiveprime() {
        // (4 bases)(r1 22 bases)(8 bases)(r1 22 bases)(6 bases)
        let read = b"ATGC|CTACACGACGCTCTTCCGATCT|ATGCATGC|CTACACGACGCTCTTCCGATCT|ATGCAT";

        let result = trimmer(RT_FRAG).find_adapters(&clean_read(read));
        assert_eq!(result[0], (AdapterKind::R1, 4..26));
        assert_eq!(result[1], (AdapterKind::R1, 34..56));
    }

    #[test]
    fn test_no_config() {
        // (5 bases)(r1 22 bases)(8 bases)(r1 22 bases)(4 bases)(r2 34 bases)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|GAATGTCC|CTACACGACGCTCTTCCGATCT|TTCC|GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RT_FRAG);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(raw_adapter_config, "END - R1 - R1 - R2 - END");
    }

    #[test]
    fn test_no_adapter() {
        let read = b"TAGACTCTGTACTGAGTTCTATAAACGAGCCATTGCATACGAGATCGGTAGATTGATAAGGGACACAGAATATCCCCGGACGCAATAGACGGACAGCTTG";

        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RT_FRAG);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(raw_adapter_config, "END - END");
    }

    #[test]
    fn test_missing_intermediate_adapter() {
        // (5 bases)(r1 22 bases)(3 bases)(r1 22 bases)(3 bases)(r2_rc 34 bases)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|GAA|CTACACGACGCTCTTCCGATCT|TTC|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RT_FRAG);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(raw_adapter_config, "END - R1 - R1 - R2_RC - END");
    }

    #[test]
    fn test_r1_fwd_config() {
        // (5 bases)(r1 22 bases)(8 bases)(r1 22 bases)(4 bases)(poly_t 20 bases)(r2_rc 34 bases)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|GAATGTCC|CTACACGACGCTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RT_FRAG);
        let expected = InsertDefinition {
            first_adapter: (AdapterKind::R1, 35..57),
            last_adapter: (AdapterKind::R2Rc, 81..115),
            orientation: Orientation::Forward,
            rna_prefix: (AdapterKind::End, 61..61),
            idx: 0,
        };
        assert_eq!(obs_configs, [expected]);
        assert_eq!(raw_adapter_config, "END - R1 - R1 - R2_RC - END");
    }

    #[test]
    fn test_get_retain_range_r1_fwd() {
        // (5 bases)(r1 22 bases)(8 bases)(r1 22 bases)(4 bases)(poly_t 20 bases)(r2_rc 34 bases)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|GAATGTCC|CTACACGACGCTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG);
        assert_eq!(obs_configs[0].get_retain_range(), 57..81);
    }

    #[test]
    fn test_r1_rev_config() {
        // (r2 34 bases)(4 bases)(poly_a 20 bases)(4 bases)(r1_rc 22 bases)(8 bases)(r1_rc 22 bases)(5 bases)
        let read = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT|GGAA|AAAAAAAAAAAAAAAAAAAA|AGAT|CGGAAGAGCGTCGTGTAGGGACATTC|AGATCGGAAGAGCGTCGTGTAG|GACAC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RT_FRAG);
        let expected = InsertDefinition {
            first_adapter: (AdapterKind::R2, 0..34),
            last_adapter: (AdapterKind::R1Rc, 58..80),
            rna_prefix: (AdapterKind::End, 54..54),
            orientation: Orientation::Reverse,
            idx: 0,
        };
        assert_eq!(obs_configs, [expected]);
        assert_eq!(raw_adapter_config, "END - R2 - R1_RC - R1_RC - END");
    }

    #[test]
    fn test_get_retain_range_r1_rev() {
        // (r2 34 bases)(4 bases)(poly_a 20 bases)(4 bases)(r1_rc 22 bases)(8 bases)(r1_rc 22 bases)(5 bases)
        let read = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT|GGAA|AAAAAAAAAAAAAAAAAAAA|AGATCGGAAGAGCGTCGTGTAG|GGACATTC|AGATCGGAAGAGCGTCGTGTAG|GACAC";
        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG);
        assert_eq!(obs_configs[0].get_retain_range(), 34..58);
    }

    #[test]
    fn test_r1_rev_directional() {
        // (r2 34 bases)(4 bases)(poly_a 20 bases)(4 bases)(r1_rc 22 bases)(8 bases)(r1_rc 22 bases)(5 bases)
        let read = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT|GGAA|AAAAAAAAAAAAAAAAAAAA|AGATCGGAAGAGCGTCGTGTAG|GGACATTC|AGATCGGAAGAGCGTCGTGTAG|GACAC";
        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG_DIRECTIONAL);
        assert_eq!(obs_configs.len(), 0);
    }

    #[test]
    fn test_multiple_configs() {
        // (5 bases)(r1 22 bases)(8 bases)(r1 22 bases)(4 bases)(poly_t 20 bases)(r2_rc 34 bases)(r2 34 bases)(poly_a 20 bases)(4 bases)(r1_rc 22 bases)(8 bases)(r1_rc 22 bases)(5 bases)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|GAATGTCC|CTACACGACGCTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC|GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT|AAAAAAAAAAAAAAAAAAAA|GATC|AGATCGGAAGAGCGTCGTGTAG|GGACATTC|AGATCGGAAGAGCGTCGTGTAG|GACAC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RT_FRAG);
        let expected = [
            InsertDefinition {
                first_adapter: (AdapterKind::R1, 35..57),
                last_adapter: (AdapterKind::R2Rc, 81..115),
                orientation: Orientation::Forward,
                rna_prefix: (AdapterKind::End, 61..61),
                idx: 0,
            },
            InsertDefinition {
                first_adapter: (AdapterKind::R2, 115..149),
                last_adapter: (AdapterKind::R1Rc, 173..195),
                orientation: Orientation::Reverse,
                rna_prefix: (AdapterKind::End, 169..169),
                idx: 1,
            },
        ];
        assert_eq!(
            raw_adapter_config,
            "END - R1 - R1 - R2_RC - R2 - R1_RC - R1_RC - END"
        );
        assert_eq!(obs_configs.len(), 2);
        assert_eq!(obs_configs[0], expected[0]);
        assert_eq!(obs_configs[1], expected[1]);
    }

    #[test]
    fn test_to_read_pair_fwd() {
        // (5 bases)(r1 22 bases)(8 bases)(r1 22 bases)(4 bases)(poly_t 20 bases)(r2_rc 34 bases)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|GAATGTCC|CTACACGACGCTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let record = OwnedRecord {
            head: "id desc0 desc1".into(),
            seq: clean_read(read),
            qual: vec![b'I'; clean_read(read).len()],
            sep: None,
        };
        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG);
        let rp = obs_configs[0]
            .to_read_pair(&record, ReadHeader::KeepId)
            .unwrap();
        assert_eq!(str::from_utf8(rp.id().unwrap()).unwrap(), "id_0");
        // R1 = BC+UMI only (4 bytes at first_adapter.end..split)
        assert_eq!(
            str::from_utf8(rp.get(WhichRead::R1, ReadPart::Seq).unwrap()).unwrap(),
            "TTCC"
        );
        // R2 = revcomp(PolyT region) = PolyA
        assert_eq!(
            str::from_utf8(rp.get(WhichRead::R2, ReadPart::Seq).unwrap()).unwrap(),
            "AAAAAAAAAAAAAAAAAAAA"
        );
    }
    #[test]
    fn test_to_read_pair_shifted_fwd() {
        // (r1 22 bases)(8 bases)(poly_t 20 bases)(5 bases)(r2_rc 34 bases)
        let read = b"CTACACGACGCTCTTCCGATCT|TTCCTTCC|TTTTTTTTTTTTTTTTTTTT|GGGGG|AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let record = OwnedRecord {
            head: "id desc0 desc1".into(),
            seq: clean_read(read),
            qual: vec![b'I'; clean_read(read).len()],
            sep: None,
        };
        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG);
        let rp = obs_configs[0]
            .to_read_pair(&record, ReadHeader::KeepId)
            .unwrap();
        assert_eq!(rp.get(WhichRead::R1, ReadPart::Seq).unwrap(), b"TTCC");
        assert_eq!(
            rp.get(WhichRead::R2, ReadPart::Seq).unwrap(),
            b"CCCCCAAAAAAAAAAAAAAAAAAAAGGAA"
        );
    }
    #[test]
    fn test_to_read_pair_shifted_rev() {
        // (r2 34 bases)(5 bases)(poly_a 20 bases)(8 bases)(r1_rc 22 bases)
        let read = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT|CCCCC|AAAAAAAAAAAAAAAAAAAA|GGAAGGAA|AGATCGGAAGAGCGTCGTGTAG";
        let record = OwnedRecord {
            head: "id desc0 desc1".into(),
            seq: clean_read(read),
            qual: vec![b'I'; clean_read(read).len()],
            sep: None,
        };
        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG_NO_TERMINAL_3P);
        let rp = obs_configs[0]
            .to_read_pair(&record, ReadHeader::KeepId)
            .unwrap();
        assert_eq!(rp.get(WhichRead::R1, ReadPart::Seq).unwrap(), b"TTCC");

        assert_eq!(
            rp.get(WhichRead::R2, ReadPart::Seq).unwrap(),
            b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCCCCAAAAAAAAAAAAAAAAAAAAGGAA"
        );
    }

    #[test]
    fn test_truncated_read_insert() {
        // (r1 22 bases)(8 bases)(poly_t 20 bases)(25 bases)
        let read =
            b"CTACACGACGCTCTTCCGATCT|TTCCTTCC|TTTTTTTTTTTTTTTTTTTT|GATCGTACGTAGCTGATCGACTGCA";
        let (obs_configs, raw_adapter_config) =
            find_configs(&clean_read(read), RT_FRAG_NO_TERMINAL_3P);
        let expected = InsertDefinition {
            first_adapter: (AdapterKind::R1, 0..22),
            last_adapter: (AdapterKind::End, 75..75),
            rna_prefix: (AdapterKind::End, 26..26),
            orientation: Orientation::Forward,
            idx: 0,
        };
        assert_eq!(obs_configs, [expected]);
        assert_eq!(raw_adapter_config, "END - R1 - END");

        let (obs_configs, _) = find_configs(&clean_read(read), RT_FRAG);
        assert_eq!(obs_configs.len(), 0);
    }

    #[test]
    fn test_overlapping_intermediate_adapter() {
        // r1 = 16..38, spacer = 107..119, small_rna_r2_rc = 116..123
        let read = b"ACAACCGGATTGATGGCTACACGACACTCTTCCGATCTAATCGCGAGGGCACTATGTGGGCTCAGTATA|TTGCTAGGACCG|TGGGTAAGCAAACCGCGTATTCAGGTTTTCTTATATGGAATT";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RTL_GEX);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(
            raw_adapter_config,
            "END - R1 - SPACER - SMALL_RNA_R2_RC - END"
        );
    }

    #[test]
    fn test_terminal_intermediate_adapter() {
        // r1 = 28..40, x12 = 71..83, spacer = 109..121, end = 121..121
        let read = b"CAACAACCGGATTGATGGCTACACGACA|CTCTTCCGATCT|ACTCTATAGTTAACCTCCGTCGATATTCTCT|TTGCTAGGACCG|TGACAGGGAGAACCGCGTATTCAGGT|TTTCTTATATGG|";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RTL_CRISPR);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(
            raw_adapter_config,
            "END - R1_BP12 - X12_CAPTURE - SPACER - END"
        );
    }

    #[test]
    fn test_short_terminal_excluded_from_global_scan() {
        // (5 bases)(r1 22 bases)(4 bases)(poly_t 20 bases)(AGATC)(5 bases)(AGATC)
        let read = b"GTGTC|CTACACGACGCTCTTCCGATCT|TTCC|TTTTTTTTTTTTTTTTTTTT|AGATC|GGGGG|AGATC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RTL_CRISPR);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(raw_adapter_config, "END - R1_BP12 - END");
    }

    #[test]
    fn test_flex_crispr() {
        // (r1 22 bases)(GGGG)(x12 12 bases)(spacer 13 bases)(ATGCA)(truncated_r2 5 bases)
        let read = b"|CTACACGACGCTCTTCCGATCT|GGGG|TTGCTAGGACCG|TTTCTTATATGGG|ATGCA|AGATC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RTL_CRISPR);
        let expected = InsertDefinition {
            first_adapter: (AdapterKind::R1Bp12, 10..22),
            last_adapter: (AdapterKind::TruncatedR2Rc, 56..61),
            rna_prefix: (AdapterKind::Spacer, 38..51),
            orientation: Orientation::Forward,
            idx: 0,
        };
        assert_eq!(obs_configs, [expected]);
        assert_eq!(
            raw_adapter_config,
            "END - R1_BP12 - X12_CAPTURE - SPACER - END"
        );
    }

    #[test]
    fn test_flex_crispr_incomplete_layout() {
        // (r1 22 bases)(GGGG)(x12 12 bases)(spacer 13 bases)(ATGCA)(truncated_r2 5 bases)
        let read = b"|CTACACGACGCTCTTCCGATCT|GGGG|TTTCTTATATGGG|ATGCA|AGATC";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RTL_CRISPR);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(raw_adapter_config, "END - R1_BP12 - SPACER - END");
    }

    #[test]
    fn test_flex_crispr_no_terminal_rejected() {
        // Read has R1 and PolyT region but no AGATC — reverse scan fails, config rejected.
        // (r1 22 bases)(GGGG)(x12 12 bases)(spacer 13 bases)(6 bases)
        let read = b"|CTACACGACGCTCTTCCGATCT|GGGG|TTGCTAGGACCG|TTTCTTATATGGG|GGGAAA";
        let (obs_configs, raw_adapter_config) = find_configs(&clean_read(read), RTL_CRISPR);
        assert_eq!(obs_configs.len(), 0);
        assert_eq!(
            raw_adapter_config,
            "END - R1_BP12 - X12_CAPTURE - SPACER - END"
        );
    }

    #[test]
    fn test_flex_crispr_selects_rightmost_match() {
        // AGATC appears twice after the spacer region; the reverse scan must return the
        // rightmost (3'-most) occurrence as the insert boundary.
        let read = b"AAAC|CTACACGACGCTCTTCCGATCT|GGGG|TTGCTAGGACCG|TTTCTTATATGGG|AGATC|GGGGG|AGATC";
        let record = OwnedRecord {
            head: "id desc0 desc1".into(),
            seq: clean_read(read),
            qual: vec![b'I'; clean_read(read).len()],
            sep: None,
        };
        let (obs_configs, _) = find_configs(&clean_read(read), RTL_CRISPR);
        assert_eq!(obs_configs.len(), 1);
        assert_eq!(
            obs_configs[0].last_adapter,
            (AdapterKind::TruncatedR2Rc, 65..70)
        );
        let rp = obs_configs[0]
            .to_read_pair(&record, ReadHeader::KeepId)
            .unwrap();

        assert_eq!(
            rp.get(WhichRead::R1, ReadPart::Seq).unwrap(),
            b"GGGGTTGCTAGGACCG"
        );
        assert_eq!(
            rp.get(WhichRead::R2, ReadPart::Seq).unwrap(),
            b"CCCCCGATCTCCCATATAAGAAA"
        );
    }
}
