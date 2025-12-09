//! ReadPair wrapper object for RNA reads from Single Cell 3' and Single Cell 5' / VDJ ibraries.
//! Provides access to the barcode and allows for dynamic trimming.
#![expect(missing_docs)]
use crate::chemistry::{
    BarcodeExtraction, BarcodeReadComponent, ChemistryDef, UmiReadComponent, UmiTranslation,
};
use crate::sample_def::SampleDef;
use crate::serde_helpers::NumberOrStr;
use crate::types::LibraryType;
use anyhow::Result;
use arrayvec::ArrayVec;
use barcode::whitelist::find_slide_design;
use barcode::{
    Barcode, BarcodeConstruct, BarcodeSegment, BcQual, BcSegQual, BcSegSeq, BcSeq,
    SegmentedBarcode, Segments, Whitelist,
};
use fastq_set::adapter_trimmer::{AdapterTrimmer, ReadAdapterCatalog, intersect_ranges};
use fastq_set::adapters::X12_CAPTURE_SEQ;
use fastq_set::metric_utils::ILLUMINA_QUAL_OFFSET;
use fastq_set::read_pair::{ReadPair, ReadPart, RpRange, WhichRead};
use fastq_set::read_pair_iter::InputFastqs;
use fastq_set::{FastqProcessor, ProcessResult};
use itertools::{Itertools, chain, zip_eq};
use martian_derive::MartianType;
use metric::TxHashMap;
use serde::{Deserialize, Deserializer, Serialize};
use std::cmp::{max, min};
use std::ops::Range;
use umi::{SplintToUmiTranslator, Umi, UmiQual, UmiSeq};

const MAX_UMI_PARTS: usize = 4;

/// Return the position of the first match of `query` in `target` permitting up to one mismatch.
pub fn find_with_one_mismatch(target: &[u8], query: &[u8]) -> Option<usize> {
    target.windows(query.len()).position(|window| {
        zip_eq(window, query)
            .filter(|(a, b)| a != b)
            .at_most_one()
            .is_ok()
    })
}

/// `RnaChunk` represents a chunk of reads from any of our RNA products.
/// This is typically created by the `SETUP_CHUNKS` stage in our pipelines.
/// This struct knows how to interpret the raw fastq data using the `chemistry`.
/// `RnaChunk` is a `FastqProcessor`, meaning it can process the raw reads and
/// create `RnaRead`, which resolves various regions in the read such as BC, UMI etc.
/// You can also specify a subsample rate and set R1/R2 lengths.
///
/// # Example
/// See the tests
///
/// # Tests
/// * `test_rna_chunk_processor_interleaved_sc_vdj` - Test that the `RnaRead`
///   produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
///   VDJ chemistry with interleaved reads.
/// * `test_rna_chunk_processor_interleaved_sc_vdj_r2` - Test that the `RnaRead`
///   produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
///   VDJ R2-only chemistry with interleaved reads.
/// * `test_rna_chunk_processor_non_interleaved_sc_vdj` - Test that the `RnaRead`
///   produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
///   VDJ chemistry with non-interleaved reads.
/// * `test_rna_chunk_processor_non_interleaved_sc_vdj_r2` - Test that the `RnaRead`
///   produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
///   VDJ R2-only chemistry with non-interleaved reads.
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug, MartianType)]
pub struct RnaChunk {
    pub chemistry: ChemistryDef,
    pub gem_group: u16,
    pub fastqs: InputFastqs,
    pub read_group: String,
    pub subsample_rate: Option<f64>,
    pub library_type: LibraryType,
    pub read_lengths: TxHashMap<WhichRead, usize>,
    pub chunk_id: u16,
    pub library_id: u16,
    pub fastq_id: Option<String>,
    pub umi_extractor: UmiExtractor,
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianType)]
pub struct UmiExtractor {
    /// A optional translator for each of the UMI parts
    translators: ArrayVec<Option<SplintToUmiTranslator>, MAX_UMI_PARTS>,
}

impl UmiExtractor {
    fn new(umi_parts: &[UmiReadComponent]) -> Self {
        UmiExtractor {
            translators: umi_parts
                .iter()
                .map(|part| {
                    part.whitelist.as_ref().map(|wl| {
                        assert_eq!(wl.translation, UmiTranslation::SingleBase);
                        let slide_path =
                            find_slide_design(&wl.slide).expect("Failed to find slide design file");
                        SplintToUmiTranslator::new(
                            slide_design::load_oligos(&slide_path, wl.part)
                                .unwrap()
                                .iter()
                                .map(String::as_bytes),
                            1,
                        )
                    })
                })
                .collect(),
        }
    }

    fn extract_umi(
        &self,
        read: &ReadPair,
        umi_read_parts: &[UmiReadComponent],
    ) -> Result<(ArrayVec<UmiPart, MAX_UMI_PARTS>, UmiSeq)> {
        let mut umi_parts = ArrayVec::new();
        let mut umi_seq = UmiSeq::new();
        for (umi_part, translator) in zip_eq(umi_read_parts, &self.translators) {
            let umi_length = max(
                min(
                    read.len(umi_part.read_type)
                        .unwrap()
                        .saturating_sub(umi_part.offset),
                    umi_part.length,
                ),
                umi_part.min_length.unwrap_or(umi_part.length),
            );
            let range = RpRange::new(umi_part.read_type, umi_part.offset, Some(umi_length));

            read.check_range(&range, "UMI")?;

            let seq = read.get_range(range, ReadPart::Seq).unwrap();
            match translator {
                Some(t) => {
                    let base = t.translate(seq);
                    umi_seq.push_unchecked(&[base]);
                    umi_parts.push(UmiPart::SingleBaseTranslated { range, base });
                }
                None => {
                    umi_seq.push_unchecked(seq);
                    umi_parts.push(UmiPart::Untranslated { range });
                }
            };
        }
        Ok((umi_parts, umi_seq))
    }
}

impl RnaChunk {
    /// Create a new `RnaChunk`.
    pub fn new(
        chemistry: &ChemistryDef,
        sample_def: &SampleDef,
        default_library_type: LibraryType,
        fastqs: InputFastqs,
        sample_id: &str,
        library_id: u16,
        chunk_id: u16,
    ) -> Result<RnaChunk> {
        let gem_group = sample_def.gem_group.unwrap_or(1);
        let (flowcell, lane) = if fastqs.r1.is_empty() {
            // For testing.
            (String::default(), 0)
        } else {
            let ilumina_header = fastqs.get_header_info()?.unwrap_or_default();
            (ilumina_header.flowcell, ilumina_header.lane)
        };
        let read_group = format!("{sample_id}:{library_id}:{gem_group}:{flowcell}:{lane}");

        let read_lengths = [
            sample_def.r1_length.map(|x| (WhichRead::R1, x)),
            sample_def.r2_length.map(|x| (WhichRead::R2, x)),
        ]
        .into_iter()
        .flatten()
        .collect();

        Ok(RnaChunk {
            chemistry: chemistry.clone(),
            gem_group,
            fastqs,
            read_group,
            subsample_rate: sample_def.subsample_rate,
            library_type: sample_def.library_type.unwrap_or(default_library_type),
            read_lengths,
            chunk_id,
            library_id,
            fastq_id: sample_def.fastq_id.clone(),
            umi_extractor: UmiExtractor::new(&chemistry.umi),
        })
    }

    /// Set the subsample rate
    ///
    /// # Test
    /// * `test_rna_chunk_subsample()` - Make sure we get roughly as many
    ///   reads as expected after subsampling
    pub fn set_subsample_rate(&mut self, value: f64) -> &mut Self {
        self.subsample_rate = Some(value);
        self
    }
    /// Set the length to hard trim the read1 in the input fastq.
    ///
    /// # Test
    /// * `prop_test_rna_chunk_trim()` - Make sure that trimming works as
    ///   expected for arbitrary inputs
    pub fn set_illumina_r1_trim_length(&mut self, value: usize) -> &mut Self {
        self.read_lengths.insert(WhichRead::R1, value);
        self
    }
    /// Set the length to hard trim the read2 in the input fastq
    ///
    /// # Test
    /// * `prop_test_rna_chunk_trim()` - Make sure that trimming works as
    ///   expected for arbitrary inputs
    pub fn set_illumina_r2_trim_length(&mut self, value: usize) -> &mut Self {
        self.read_lengths.insert(WhichRead::R2, value);
        self
    }

    pub fn gem_group(&self) -> u16 {
        self.gem_group
    }

    pub fn library_type(&self) -> LibraryType {
        self.library_type
    }

    pub fn library_id(&self) -> u16 {
        self.library_id
    }

    pub fn library_info(&self, target_set_name: Option<String>) -> LibraryInfo {
        LibraryInfo {
            library_id: self.library_id,
            library_type: self.library_type(),
            gem_group: self.gem_group(),
            target_set_name,
        }
    }
}

#[derive(Deserialize, Serialize, Debug, PartialOrd, Ord, PartialEq, Eq)]
pub struct LibraryInfo {
    #[serde(deserialize_with = "deserialize_number_or_parse_string")]
    pub library_id: u16,
    pub library_type: LibraryType,
    pub gem_group: u16,
    #[serde(default)]
    pub target_set_name: Option<String>,
}

fn deserialize_number_or_parse_string<'de, D>(deserializer: D) -> Result<u16, D::Error>
where
    D: Deserializer<'de>,
{
    let num_or_str = NumberOrStr::deserialize(deserializer)?;
    match num_or_str {
        NumberOrStr::Number(n) => Ok(n as u16),
        NumberOrStr::Str(s) => Ok(s.parse::<u16>().unwrap()),
    }
}

/// Return a LibraryInfo for each distinct library_id.
pub fn make_library_info(
    read_chunks: &[RnaChunk],
    target_set_name: Option<&str>,
) -> Vec<LibraryInfo> {
    read_chunks
        .iter()
        .unique_by(|chunk| chunk.library_id())
        .map(|chunk| chunk.library_info(target_set_name.map(String::from)))
        .sorted()
        .collect()
}

pub struct RnaProcessor {
    chunk: RnaChunk,
    whitelist: BarcodeConstruct<Whitelist>,
    barcode_lengths: BarcodeConstruct<Range<usize>>,
}

impl RnaProcessor {
    pub fn new(chunk: RnaChunk, whitelist: BarcodeConstruct<Whitelist>) -> Self {
        RnaProcessor {
            chunk,
            barcode_lengths: whitelist.as_ref().map(Whitelist::sequence_lengths),
            whitelist,
        }
    }
}

/// Given the components of an extraction strategy, extract the barcode.
///
/// Returns the ranges from which the barcode components were extracted,
/// as well as the barcode segments.
pub fn extract_barcode(
    read: &ReadPair,
    barcode_components: BarcodeConstruct<&BarcodeReadComponent>,
    barcode_extraction: Option<&BarcodeExtraction>,
    whitelist: BarcodeConstruct<&Whitelist>,
    barcode_lengths: BarcodeConstruct<&Range<usize>>,
    gem_group: u16,
) -> Result<(BarcodeConstruct<RpRange>, SegmentedBarcode)> {
    let bc_range = match barcode_extraction {
        None => barcode_components.map(BarcodeReadComponent::rprange),
        Some(BarcodeExtraction::VariableMultiplexingBarcode {
            min_offset,
            max_offset,
        }) => {
            BarcodeConstruct::new_gel_bead_and_probe(
                // Gel bead barcode position is fixed.
                barcode_components.gel_bead().rprange(),
                get_probe_barcode_range(
                    read,
                    barcode_components,
                    whitelist.probe(),
                    *min_offset,
                    *max_offset,
                ),
            )
        }
        Some(BarcodeExtraction::JointBc1Bc2 {
            min_offset,
            max_offset,
        }) => {
            // Ensure that we are actually dealing with 2 part barcodes
            // It's more convenient to work with arrays for these two
            let component_segments = barcode_components.segments().array_vec();
            let length_segments = barcode_lengths.segments().array_vec();

            // Bunch of sanity checks
            assert!(max_offset >= min_offset);
            assert_eq!(component_segments.len(), 2);
            let which_read = component_segments[0].read_type;
            assert_eq!(which_read, component_segments[1].read_type);

            let read_len = read.len(which_read).unwrap();

            // Fallback option if we do not find any matches
            let default_range = barcode_components.map(BarcodeReadComponent::rprange);

            (*min_offset..=*max_offset)
                .flat_map(|offset| {
                    // For the given bc1 start position (offset), find the possible ranges for bc1 and bc2
                    // that could give us an exact match. To find this, consider all combinations of
                    // the lengths of both bc1 and bc2.
                    length_segments[0]
                        .clone()
                        .cartesian_product(length_segments[1].clone())
                        .filter_map(move |(len1, len2)| {
                            // Skip cases where the end of bc2 goes past the end of read
                            (offset + len1 + len2 <= read_len).then_some(
                                BarcodeConstruct::Segmented(Segments::from_iter([
                                    RpRange::new(which_read, offset, Some(len1)),
                                    RpRange::new(which_read, offset + len1, Some(len2)),
                                ])),
                            )
                        })
                })
                // Chain the default option. `max_by_key` would pick the last element if there
                // are several maximums. So we will get the default range if there are no
                // ranges with bc1 or bc2 exact matches
                .chain(std::iter::once(default_range))
                .max_by_key(|bc_range| {
                    // Find the barcode ranges that gives us most matches to the whitelist. Note
                    // that true > false, so we don't need to map the booleans to integer to pick
                    // the correct maximum
                    whitelist.zip(*bc_range).map(|(wl, range)| {
                        read.get_range(range, ReadPart::Seq)
                            .is_some_and(|seq| wl.contains(&BcSegSeq::from_bytes_unchecked(seq)))
                    })
                })
                .unwrap_or(default_range)
        }
    };

    for range in bc_range {
        read.check_range(&range, "Barcode")?;
    }

    let barcode = SegmentedBarcode::new(
        gem_group,
        bc_range.zip(whitelist).map(|(r, wl)| {
            let mut segment =
                BarcodeSegment::with_sequence_unchecked(read.get_range(r, ReadPart::Seq).unwrap());
            // The check() function changes the state of the barcode segment
            // to reflect whether it is valid before correction or not
            wl.check_and_update(&mut segment);
            segment
        }),
    );
    Ok((bc_range, barcode))
}

/// Locate the probe barcode and return its position.
/// When sequenced on R1, locate the capture sequence and return the position after it if found.
/// Otherwise return the position of the first valid probe barcode between `min_offset` and `max_offset`.
pub fn get_probe_barcode_range(
    read: &ReadPair,
    barcode_components: BarcodeConstruct<&BarcodeReadComponent>,
    whitelist: &Whitelist,
    min_offset: i64,
    max_offset: i64,
) -> RpRange {
    // Check all possible positions for a valid probe barcode.
    // Default to no offset if we don't find any.
    let default_probe_range = barcode_components.probe().rprange();

    if min_offset == 0 && max_offset == 0 {
        return default_probe_range;
    }

    // Locate the probe barcode by locating the capture sequence.
    let additional_barcode_pos = if default_probe_range.read() == WhichRead::R1 {
        let gel_bead = barcode_components.gel_bead();
        assert_eq!(gel_bead.read_type, WhichRead::R1);
        assert_eq!(gel_bead.offset(), 0);
        let r1_seq = read.get(WhichRead::R1, ReadPart::Seq).unwrap();
        // Do not return a match to the capture sequence if it overlaps with the gel-bead barcode.
        find_with_one_mismatch(r1_seq, X12_CAPTURE_SEQ)
            .filter(|&x12_pos| x12_pos >= gel_bead.length())
            .map(|x12_pos| x12_pos + X12_CAPTURE_SEQ.len())
    } else {
        None
    };

    // Return the position of the first valid probe barcode following the capture sequence
    // or between `min_offset` and `max_offset`.
    let candidate_barcode_positions = chain(
        additional_barcode_pos,
        (min_offset..=max_offset).map(|offset| {
            let new_offset = default_probe_range.offset() as i64 + offset;
            assert!(new_offset >= 0);
            new_offset as usize
        }),
    );
    match candidate_barcode_positions
        .filter_map(|position| {
            let range = RpRange::new(
                default_probe_range.read(),
                position,
                default_probe_range.len(),
            );
            read.get_range(range, ReadPart::Seq).and_then(|seq| {
                whitelist
                    .contains(&BcSegSeq::from_bytes_unchecked(seq))
                    .then_some(range)
            })
        })
        .at_most_one()
    {
        // Found no perfect match. Return the default range.
        Ok(None) => default_probe_range,
        // Found a perfect match. Return it.
        Ok(Some(range)) => range,
        // Found more than one perfect match. Return the first match.
        // TODO: Annotate the barcode as `InvalidAmbiguous`
        Err(mut iter) => iter.next().unwrap(),
    }
}

impl FastqProcessor for RnaProcessor {
    type ReadType = RnaRead;
    /// This function creates a processed `RnaRead` from a raw `ReadPair`
    /// using the chemistry stored in `self`. A rough outline of what
    /// happens inside:
    /// - Attach Barcodes and UMIs
    /// - Find the ranges for RNA read1 and RNA read2. Remember, RNA read1 could
    ///   point to illumina read2, ut it guaranteed to be empty
    /// - Optionally hard trim the illumina R1 and R2 reads
    fn process_read(&self, read: ReadPair) -> ProcessResult<RnaRead> {
        let chem = &self.chunk.chemistry;

        let (umi_parts, umi_seq) = match self.chunk.umi_extractor.extract_umi(&read, &chem.umi) {
            Ok((r, s)) => (r, s),
            Err(e) => {
                return ProcessResult::Unprocessed {
                    read,
                    reason: e.to_string(),
                };
            }
        };

        let r1_range = {
            let rna = chem.rna;
            let read_length = read.len(rna.read_type).unwrap_or(0);
            let read_rna_length = read_length.saturating_sub(rna.offset);
            let rna_length = match (rna.min_length, rna.length) {
                (None, None) => None,
                (None, Some(length)) => Some(length),
                (Some(min_length), None) => Some(read_rna_length.max(min_length)),
                (Some(min_length), Some(length)) => Some(read_rna_length.clamp(min_length, length)),
            };
            let mut range = RpRange::new(rna.read_type, rna.offset, rna_length);

            // TODO: this is no longer necessary, plumbed through FastqProcessor
            if let Some(&trim_length) = self.chunk.read_lengths.get(&rna.read_type) {
                let trim_len = min(trim_length, read_length);
                range.intersect(RpRange::new(rna.read_type, 0, Some(trim_len)));
            }

            if let Err(e) = read.check_range(&range, "First RNA read") {
                return ProcessResult::Unprocessed {
                    read,
                    reason: e.to_string(),
                };
            }
            range
        };

        let r2_range = match chem.rna2 {
            Some(rna2) => {
                let mut range: RpRange = rna2.into();
                // TODO: this is no longer necessary, plumbed through FastqProcessor
                if let Some(len) = self.chunk.read_lengths.get(&rna2.read_type) {
                    let trim_len = min(*len, read.len(rna2.read_type).unwrap_or(0));
                    range.intersect(RpRange::new(rna2.read_type, 0, Some(trim_len)));
                }

                if let Err(e) = read.check_range(&range, "Second RNA read") {
                    return ProcessResult::Unprocessed {
                        read,
                        reason: e.to_string(),
                    };
                }
                Some(range)
            }
            None => None,
        };

        let (bc_range, barcode) = match extract_barcode(
            &read,
            chem.barcode_construct(),
            chem.barcode_extraction(),
            self.whitelist.as_ref(),
            self.barcode_lengths.as_ref(),
            self.gem_group(),
        ) {
            Ok(value) => value,
            Err(e) => {
                return ProcessResult::Unprocessed {
                    read,
                    reason: e.to_string(),
                };
            }
        };

        ProcessResult::Processed(RnaRead {
            read,
            segmented_barcode: barcode,
            umi: umi_seq.into(),
            bc_range,
            umi_parts,
            r1_range,
            r2_range,
            library_type: self.chunk.library_type,
            chunk_id: self.chunk.chunk_id,
        })
    }

    /// Get the fastq files
    ///
    /// # Panics
    /// * If interleaved and r2 is not None
    /// * If not interleaved and r2 is None
    ///
    /// # Tests
    /// * `test_rna_chunk_fastq_panic_1()` - Make sure that we panic if
    ///   interleaved and r2 is not None
    /// * `test_rna_chunk_fastq_panic_2()` - Make sure that we panic if
    ///   not interleaved and r2 is None
    fn fastq_files(&self) -> InputFastqs {
        // Make sure that either
        // - r2 is None and interleaved
        // - r2 is Some and not interleaved
        if self.chunk.fastqs.r2.is_some() {
            assert!(!self.chunk.fastqs.r1_interleaved);
        }

        self.chunk.fastqs.clone()
    }

    fn bc_subsample_rate(&self) -> f64 {
        1.0
    }
    fn read_subsample_rate(&self) -> f64 {
        self.chunk.subsample_rate.unwrap_or(1.0)
    }
    fn illumina_r1_trim_length(&self) -> Option<usize> {
        self.chunk.read_lengths.get(&WhichRead::R1).copied()
    }
    fn illumina_r2_trim_length(&self) -> Option<usize> {
        self.chunk.read_lengths.get(&WhichRead::R2).copied()
    }

    fn gem_group(&self) -> u16 {
        self.chunk.gem_group
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiPart {
    SingleBaseTranslated { range: RpRange, base: u8 },
    Untranslated { range: RpRange },
}

impl UmiPart {
    pub fn range(self) -> RpRange {
        match self {
            UmiPart::SingleBaseTranslated { range, .. } => range,
            UmiPart::Untranslated { range } => range,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RnaRead {
    pub read: ReadPair,
    pub segmented_barcode: SegmentedBarcode,
    pub bc_range: BarcodeConstruct<RpRange>,
    pub umi: Umi,
    pub umi_parts: ArrayVec<UmiPart, MAX_UMI_PARTS>,
    pub r1_range: RpRange,
    pub r2_range: Option<RpRange>,
    pub library_type: LibraryType,
    pub chunk_id: u16,
}

impl RnaRead {
    pub fn barcode(&self) -> Barcode {
        self.segmented_barcode.into()
    }

    pub fn barcode_is_valid(&self) -> bool {
        self.segmented_barcode.is_valid()
    }

    pub fn raw_bc_seq(&self) -> BcSeq {
        let mut seq = BcSeq::new();
        for s in self.raw_bc_construct_seq() {
            seq.push_unchecked(s.as_bytes());
        }
        seq
    }

    pub fn raw_bc_qual(&self) -> BcQual {
        let mut qual = BcQual::new();
        for q in self.raw_bc_construct_qual() {
            qual.push_unchecked(q.as_bytes());
        }
        qual
    }
    pub fn raw_bc_construct_seq(&self) -> BarcodeConstruct<BcSegSeq> {
        self.bc_range
            .map(|r| BcSegSeq::from_bytes_unchecked(self.read.get_range(r, ReadPart::Seq).unwrap()))
    }
    pub fn raw_bc_construct_qual(&self) -> BarcodeConstruct<BcSegQual> {
        self.bc_range.map(|r| {
            BcSegQual::from_bytes_unchecked(self.read.get_range(r, ReadPart::Qual).unwrap())
        })
    }

    /// Return whether this segmented barcode has a probe barcode.
    pub fn has_probe_barcode(&self) -> bool {
        matches!(self.bc_range, BarcodeConstruct::GelBeadAndProbe(_))
    }

    pub fn umi(&self) -> Umi {
        self.umi
    }
    pub fn raw_umi(&self) -> Umi {
        self.raw_umi_seq().into()
    }
    pub fn correct_umi(&mut self, corrected: &[u8]) {
        self.umi = Umi::new(corrected);
    }
    pub fn bc_range(&self) -> BarcodeConstruct<RpRange> {
        self.bc_range
    }
    pub fn has_two_part_barcode(&self) -> bool {
        matches!(self.bc_range, BarcodeConstruct::Segmented(s) if s.is_two_part())
    }
    pub fn umi_ranges(&self) -> impl Iterator<Item = RpRange> + '_ {
        self.umi_parts.iter().map(|part| part.range())
    }
    pub fn r1_range(&self) -> RpRange {
        self.r1_range
    }
    pub fn r2_range(&self) -> Option<RpRange> {
        self.r2_range
    }
    pub fn readpair(&self) -> &ReadPair {
        &self.read
    }
    pub fn r2_exists(&self) -> bool {
        self.r2_range.is_some()
    }

    pub fn is_gene_expression(&self) -> bool {
        self.library_type == LibraryType::Gex
    }

    /// FASTQ read header
    pub fn header(&self) -> &[u8] {
        self.read.get(WhichRead::R1, ReadPart::Header).unwrap()
    }

    /// Full raw R1 sequence
    pub fn raw_illumina_read1_seq(&self) -> &[u8] {
        self.read.get(WhichRead::R1, ReadPart::Seq).unwrap()
    }

    /// Full raw R1 QVs
    pub fn raw_illumina_read1_qual(&self) -> &[u8] {
        self.read.get(WhichRead::R1, ReadPart::Qual).unwrap()
    }

    /// Full R2 sequence
    pub fn raw_illumina_read2_seq(&self) -> &[u8] {
        self.read.get(WhichRead::R2, ReadPart::Seq).unwrap()
    }

    /// Full R2 QVs
    pub fn raw_illumina_read2_qual(&self) -> &[u8] {
        self.read.get(WhichRead::R2, ReadPart::Qual).unwrap()
    }

    /// Full raw I1 sequence
    pub fn raw_illumina_i1_seq(&self) -> Option<&[u8]> {
        self.read.get(WhichRead::I1, ReadPart::Seq)
    }

    /// Full raw I1 QVs
    pub fn raw_illumina_i1_qual(&self) -> Option<&[u8]> {
        self.read.get(WhichRead::I1, ReadPart::Qual)
    }

    /// Full I2 sequence
    pub fn raw_illumina_i2_seq(&self) -> Option<&[u8]> {
        self.read.get(WhichRead::I2, ReadPart::Seq)
    }

    /// Full I2 QVs
    pub fn raw_illumina_i2_qual(&self) -> Option<&[u8]> {
        self.read.get(WhichRead::I2, ReadPart::Qual)
    }

    /// Raw, uncorrected UMI sequence
    pub fn raw_umi_seq(&self) -> UmiSeq {
        let mut umi = UmiSeq::new();
        for part in &self.umi_parts {
            match part {
                UmiPart::SingleBaseTranslated { base, .. } => {
                    umi.push_unchecked(&[*base]);
                }
                UmiPart::Untranslated { range } => {
                    umi.push_unchecked(self.read.get_range(*range, ReadPart::Seq).unwrap());
                }
            }
        }
        umi
    }

    /// Raw UMI QVs
    pub fn raw_umi_qual(&self) -> UmiQual {
        let mut qual = UmiQual::new();
        for part in &self.umi_parts {
            match part {
                UmiPart::SingleBaseTranslated { range, .. } => {
                    // We use the minimum quality among the original splint bases
                    let min_qual = self
                        .read
                        .get_range(*range, ReadPart::Qual)
                        .unwrap()
                        .iter()
                        .min()
                        .unwrap();
                    qual.push_unchecked(&[*min_qual]);
                }
                UmiPart::Untranslated { range } => {
                    qual.push_unchecked(self.read.get_range(*range, ReadPart::Qual).unwrap());
                }
            }
        }
        qual
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_seq(&self) -> &[u8] {
        self.read.get_range(self.r1_range, ReadPart::Seq).unwrap()
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_qual(&self) -> &[u8] {
        self.read.get_range(self.r1_range, ReadPart::Qual).unwrap()
    }

    /// Usable R2 bases after removal of BC and trimming
    pub fn r2_seq(&self) -> Option<&[u8]> {
        if let Some(range) = self.r2_range {
            self.read.get_range(range, ReadPart::Seq)
        } else {
            None
        }
    }

    /// Usable R2 bases after removal of BC and trimming
    pub fn r2_qual(&self) -> Option<&[u8]> {
        if let Some(range) = self.r2_range {
            self.read.get_range(range, ReadPart::Qual)
        } else {
            None
        }
    }

    /// get the RnaChunk for this read
    pub fn read_chunk<'a>(&self, read_chunks: &'a [RnaChunk]) -> &'a RnaChunk {
        &read_chunks[self.chunk_id as usize]
    }

    /// Returns the minimum barcode quality score.
    pub fn barcode_min_qual(&self) -> u8 {
        let &q = self.raw_bc_qual().as_bytes().iter().min().unwrap();
        assert!(q >= ILLUMINA_QUAL_OFFSET);
        q - ILLUMINA_QUAL_OFFSET
    }

    /// Returns the minimum UMI quality score.
    pub fn umi_min_qual(&self) -> u8 {
        let &q = self.raw_umi_qual().iter().min().unwrap();
        assert!(q >= ILLUMINA_QUAL_OFFSET);
        q - ILLUMINA_QUAL_OFFSET
    }

    /// Given an `RpRange` and a list of adapter trimmers, this function
    /// searches for all the adapters within the `RpRange` and returns a
    /// `Range<usize>`, which you need to shrink the `RpRange` to.
    /// It also returns a hashmap from the adapter name to the `RpRange`
    /// where the adapter was found. This is useful for book-keeping and
    /// computing metrics.
    fn trim_adapters_helper(
        &mut self,
        read_range: RpRange,
        trimmers: &mut [AdapterTrimmer<'_>],
    ) -> (Range<usize>, TxHashMap<String, RpRange>) {
        let seq = self.read.get_range(read_range, ReadPart::Seq).unwrap();
        let trim_results: Vec<_> = trimmers
            .iter_mut()
            .filter_map(|t| {
                t.find(seq)
                    .map(|trim_result| (t.adapter.name.clone(), trim_result))
            })
            .collect();
        let shrink_range = trim_results.iter().fold(0..seq.len(), |acc, (_, x)| {
            intersect_ranges(&acc, &x.retain_range)
        });
        let adapter_pos =
            trim_results
                .into_iter()
                .fold(TxHashMap::default(), |mut acc, (name, trim_result)| {
                    let mut this_range = read_range;
                    this_range.shrink(&trim_result.adapter_range);
                    acc.insert(name, this_range);
                    acc
                });
        (shrink_range, adapter_pos)
    }

    /// Perform adapter trimming on the `RnaRead` and return the positions
    /// of the adapters found.
    ///
    /// # Inputs
    /// * `adapter_catalog`: Packages all the adapter trimmers
    ///
    /// # Output
    /// * `TxHashMap<String, RpRange>` - where the key is the name of the adapter
    ///   and values is the `RpRange` where the adapter is found. Clearly, the adapters
    ///   which are not present in the read will not be part of the output. The `ReadAdapterCatalog`
    ///   guarantees that no two adapters share the same name, so there is no confusion.
    ///
    /// # Test
    /// * `test_rna_read_adapter_trim()`: Test that we can trim reads consistent with cutadapt.
    pub fn trim_adapters(
        &mut self,
        adapter_catalog: &mut ReadAdapterCatalog<'_>,
    ) -> TxHashMap<String, RpRange> {
        let mut result = TxHashMap::default();
        // Trim r1
        {
            let ad_trimmers = adapter_catalog.get_mut_trimmers(self.r1_range.read());
            let range = self.r1_range; // Creates a copy
            let (shrink_range, adapter_pos) = self.trim_adapters_helper(range, ad_trimmers);
            result.extend(adapter_pos);
            self.r1_range.shrink(&shrink_range);
        }

        // Trim r2
        {
            if let Some(range) = self.r2_range {
                let ad_trimmers = adapter_catalog.get_mut_trimmers(range.read());
                let (shrink_range, adapter_pos) = self.trim_adapters_helper(range, ad_trimmers);
                result.extend(adapter_pos);
                if let Some(x) = self.r2_range.as_mut() {
                    x.shrink(&shrink_range);
                }
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::{BarcodeKind, ChemistryName, UmiWhitelistSpec};
    use BarcodeConstruct::GelBeadOnly;
    use barcode::whitelist::ReqStrand;
    use barcode::{BarcodeSegmentState, Segments, WhitelistSource, WhitelistSpec};
    use itertools::assert_equal;
    use metric::{Histogram, TxHashSet, set};
    use proptest::arbitrary::any;
    use proptest::proptest;
    use serde_json;
    use std::borrow::Cow;
    use std::fs::File;
    use std::path::Path;

    #[test]
    fn test_find_seq() {
        assert_eq!(find_with_one_mismatch(b"Hello, world!", b"world"), Some(7));
        assert_eq!(find_with_one_mismatch(b"Hello, world!", b"wxrld"), Some(7));
        assert_eq!(find_with_one_mismatch(b"Hello, world!", b"wxxld"), None);
    }

    fn processor(chunk: RnaChunk) -> RnaProcessor {
        let whitelist = chunk
            .chemistry
            .barcode_whitelist_source()
            .unwrap()
            .as_ref()
            .map_result(WhitelistSource::as_whitelist)
            .unwrap();
        RnaProcessor::new(chunk, whitelist)
    }

    #[test]
    fn test_read_10bp_umis() -> Result<()> {
        // make sure we can read 10bp UMIs with the v3 chemistry setup.
        let chemistry = ChemistryDef::named(ChemistryName::ThreePrimeV3PolyA);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/small_3p_v2.fastq".into(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: true,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::Gex,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };

        let mut n = 0;
        for r in processor(chunk).iter()? {
            let umi = match r? {
                ProcessResult::Processed(rna) => rna.umi(),
                ProcessResult::Unprocessed { .. } => unreachable!(),
            };
            assert_eq!(umi.seq().len(), 10);
            n += 1;
        }

        // make sure all the records were produced
        assert_eq!(n, 1024);
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_rna_chunk_fastq_panic_1() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/r1_2k.fastq.lz4".into(),
                r2: Some("test/rna_read/r2_2k.fastq.lz4".into()),
                i1: None,
                i2: None,
                r1_interleaved: true,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let _ = processor(chunk).fastq_files();
    }

    #[test]
    #[should_panic]
    #[ignore] // this behavior is no longer panicky b/c SC5P-R1
    fn test_rna_chunk_fastq_panic_2() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/r1_1k.fasta.lz4".into(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: true,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let _ = processor(chunk).fastq_files();
    }

    fn load_expected(fname: impl AsRef<Path>) -> Vec<Vec<u8>> {
        let expected: Vec<String> = serde_json::from_reader(File::open(fname).unwrap()).unwrap();
        expected
            .into_iter()
            .map(|x| x.as_bytes().to_vec())
            .collect()
    }

    #[test]
    fn test_rna_chunk_processor_interleaved_sc_vdj() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/interleaved_2k.fastq.lz4".into(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: true,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let expected_bc = load_expected("test/rna_read/bc.json");
        let expected_bc_qual = load_expected("test/rna_read/bc_qual.json");
        let expected_umi = load_expected("test/rna_read/umi.json");
        let expected_umi_qual = load_expected("test/rna_read/umi_qual.json");
        let expected_r1 = load_expected("test/rna_read/rna_r1.json");
        let expected_r1_qual = load_expected("test/rna_read/rna_r1_qual.json");
        let expected_r2 = load_expected("test/rna_read/rna_r2.json");
        let expected_r2_qual = load_expected("test/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in processor(chunk).iter().unwrap().enumerate() {
            let ProcessResult::Processed(rna_read) = rna_read_result.unwrap() else {
                unreachable!();
            };

            // Check that the barcode data is correct
            assert_eq!(
                rna_read.bc_range,
                GelBeadOnly(RpRange::new(WhichRead::R1, 0, Some(16)))
            );
            assert_eq!(rna_read.barcode().sequence_bytes(), &expected_bc[i]);
            assert_eq!(
                rna_read.raw_bc_construct_seq(),
                GelBeadOnly(BcSegSeq::from_bytes(&expected_bc[i]))
            );
            assert_eq!(
                rna_read.raw_bc_seq(),
                BcSeq::from_bytes(expected_bc[i].as_slice())
            );
            assert_eq!(
                rna_read.raw_bc_construct_qual(),
                GelBeadOnly(BcSegQual::from_bytes(&expected_bc_qual[i]))
            );
            assert_eq!(
                rna_read.raw_bc_qual(),
                BcQual::from_bytes(expected_bc_qual[i].as_slice())
            );

            // Check that the UMI data is correct
            assert_equal(
                rna_read.umi_ranges(),
                [RpRange::new(WhichRead::R1, 16, Some(10))],
            );
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(
                rna_read.raw_umi_seq().as_bytes(),
                expected_umi[i].as_slice()
            );
            assert_eq!(
                rna_read.raw_umi_qual().as_bytes(),
                expected_umi_qual[i].as_slice()
            );

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R1, 41, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(
                rna_read.r2_range,
                Some(RpRange::new(WhichRead::R2, 0, None))
            );
            assert_eq!(rna_read.r2_seq().unwrap(), expected_r2[i].as_slice());
            assert_eq!(rna_read.r2_qual().unwrap(), expected_r2_qual[i].as_slice());
        }
    }

    #[test]
    fn test_rna_chunk_processor_interleaved_sc_vdj_r2() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjR2);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/interleaved_2k.fastq.lz4".into(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: true,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let expected_bc = load_expected("test/rna_read/bc.json");
        let expected_bc_qual = load_expected("test/rna_read/bc_qual.json");
        let expected_umi = load_expected("test/rna_read/umi.json");
        let expected_umi_qual = load_expected("test/rna_read/umi_qual.json");
        let expected_r1 = load_expected("test/rna_read/rna_r2.json");
        let expected_r1_qual = load_expected("test/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in processor(chunk).iter().unwrap().enumerate() {
            let ProcessResult::Processed(rna_read) = rna_read_result.unwrap() else {
                unreachable!();
            };

            // Check that the barcode data is correct
            assert_eq!(
                rna_read.bc_range,
                GelBeadOnly(RpRange::new(WhichRead::R1, 0, Some(16)))
            );
            assert_eq!(rna_read.barcode().sequence_bytes(), &expected_bc[i]);
            assert_eq!(
                rna_read.raw_bc_construct_seq(),
                GelBeadOnly(BcSegSeq::from_bytes(&expected_bc[i]))
            );
            assert_eq!(
                rna_read.raw_bc_seq(),
                BcSeq::from_bytes(expected_bc[i].as_slice())
            );
            assert_eq!(
                rna_read.raw_bc_construct_qual(),
                GelBeadOnly(BcSegQual::from_bytes(&expected_bc_qual[i]))
            );
            assert_eq!(
                rna_read.raw_bc_qual(),
                BcQual::from_bytes(expected_bc_qual[i].as_slice())
            );

            // Check that the UMI data is correct
            assert_equal(
                rna_read.umi_ranges(),
                [RpRange::new(WhichRead::R1, 16, Some(10))],
            );
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(
                rna_read.raw_umi_seq().as_bytes(),
                expected_umi[i].as_slice()
            );
            assert_eq!(
                rna_read.raw_umi_qual().as_bytes(),
                expected_umi_qual[i].as_slice()
            );

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R2, 0, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(rna_read.r2_range, None);
            assert_eq!(rna_read.r2_seq(), None);
            assert_eq!(rna_read.r2_qual(), None);
        }
    }

    #[test]
    fn test_rna_chunk_processor_non_interleaved_sc_vdj() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/r1_2k.fastq.lz4".into(),
                r2: Some("test/rna_read/r2_2k.fastq.lz4".into()),
                i1: None,
                i2: None,
                r1_interleaved: false,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let expected_bc = load_expected("test/rna_read/bc.json");
        let expected_bc_qual = load_expected("test/rna_read/bc_qual.json");
        let expected_umi = load_expected("test/rna_read/umi.json");
        let expected_umi_qual = load_expected("test/rna_read/umi_qual.json");
        let expected_r1 = load_expected("test/rna_read/rna_r1.json");
        let expected_r1_qual = load_expected("test/rna_read/rna_r1_qual.json");
        let expected_r2 = load_expected("test/rna_read/rna_r2.json");
        let expected_r2_qual = load_expected("test/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in processor(chunk).iter().unwrap().enumerate() {
            let ProcessResult::Processed(rna_read) = rna_read_result.unwrap() else {
                unreachable!();
            };

            // Check that the barcode data is correct
            assert_eq!(
                rna_read.bc_range,
                GelBeadOnly(RpRange::new(WhichRead::R1, 0, Some(16)))
            );
            assert_eq!(rna_read.barcode().sequence_bytes(), &expected_bc[i]);
            assert_eq!(
                rna_read.raw_bc_construct_seq(),
                GelBeadOnly(BcSegSeq::from_bytes(&expected_bc[i]))
            );
            assert_eq!(
                rna_read.raw_bc_seq(),
                BcSeq::from_bytes(expected_bc[i].as_slice())
            );
            assert_eq!(
                rna_read.raw_bc_construct_qual(),
                GelBeadOnly(BcSegQual::from_bytes(&expected_bc_qual[i]))
            );
            assert_eq!(
                rna_read.raw_bc_qual(),
                BcQual::from_bytes(expected_bc_qual[i].as_slice())
            );

            // Check that the UMI data is correct
            assert_equal(
                rna_read.umi_ranges(),
                [RpRange::new(WhichRead::R1, 16, Some(10))],
            );
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(
                rna_read.raw_umi_seq().as_bytes(),
                expected_umi[i].as_slice()
            );
            assert_eq!(
                rna_read.raw_umi_qual().as_bytes(),
                expected_umi_qual[i].as_slice()
            );

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R1, 41, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(
                rna_read.r2_range,
                Some(RpRange::new(WhichRead::R2, 0, None))
            );
            assert_eq!(rna_read.r2_seq().unwrap(), expected_r2[i].as_slice());
            assert_eq!(rna_read.r2_qual().unwrap(), expected_r2_qual[i].as_slice());
        }
    }

    #[test]
    fn test_rna_chunk_processor_non_interleaved_sc_vdj_r2() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjR2);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/r1_2k.fastq.lz4".into(),
                r2: Some("test/rna_read/r2_2k.fastq.lz4".into()),
                i1: None,
                i2: None,
                r1_interleaved: false,
            },
            read_group: "Blah".into(),
            subsample_rate: None,
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let expected_bc = load_expected("test/rna_read/bc.json");
        let expected_bc_qual = load_expected("test/rna_read/bc_qual.json");
        let expected_umi = load_expected("test/rna_read/umi.json");
        let expected_umi_qual = load_expected("test/rna_read/umi_qual.json");
        let expected_r1 = load_expected("test/rna_read/rna_r2.json");
        let expected_r1_qual = load_expected("test/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in processor(chunk).iter().unwrap().enumerate() {
            let ProcessResult::Processed(rna_read) = rna_read_result.unwrap() else {
                unreachable!();
            };

            // Check that the barcode data is correct
            assert_eq!(
                rna_read.bc_range,
                GelBeadOnly(RpRange::new(WhichRead::R1, 0, Some(16)))
            );
            assert_eq!(rna_read.barcode().sequence_bytes(), &expected_bc[i]);
            assert_eq!(
                rna_read.raw_bc_construct_seq(),
                GelBeadOnly(BcSegSeq::from_bytes(&expected_bc[i]))
            );
            assert_eq!(
                rna_read.raw_bc_seq(),
                BcSeq::from_bytes(expected_bc[i].as_slice())
            );
            assert_eq!(
                rna_read.raw_bc_construct_qual(),
                GelBeadOnly(BcSegQual::from_bytes(&expected_bc_qual[i]))
            );
            assert_eq!(
                rna_read.raw_bc_qual(),
                BcQual::from_bytes(expected_bc_qual[i].as_slice())
            );

            // Check that the UMI data is correct
            assert_equal(
                rna_read.umi_ranges(),
                [RpRange::new(WhichRead::R1, 16, Some(10))],
            );
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(
                rna_read.raw_umi_seq().as_bytes(),
                expected_umi[i].as_slice()
            );
            assert_eq!(
                rna_read.raw_umi_qual().as_bytes(),
                expected_umi_qual[i].as_slice()
            );

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R2, 0, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(rna_read.r2_range, None);
            assert_eq!(rna_read.r2_seq(), None);
            assert_eq!(rna_read.r2_qual(), None);
        }
    }

    #[test]
    fn test_rna_chunk_subsample() {
        let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
        let umi_extractor = UmiExtractor::new(&chemistry.umi);
        let chunk = RnaChunk {
            chemistry,
            umi_extractor,
            gem_group: 1,
            fastqs: InputFastqs {
                r1: "test/rna_read/interleaved_2k.fastq".into(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: true,
            },
            read_group: "Blah".into(),
            subsample_rate: Some(0.2),
            library_type: LibraryType::VdjAuto,
            read_lengths: TxHashMap::default(),
            chunk_id: 0,
            library_id: 0,
            fastq_id: None,
        };
        let processed_reads: usize = processor(chunk).iter().unwrap().map(|_| 1).sum();
        println!("{processed_reads}");
        // Expecting 400 reads since we start with 2000 reads
        assert!(processed_reads > 300); // < 1e-6 probability of this happening by chance
        assert!(processed_reads < 500); // < 1e-6 probability of this happening by chance
    }

    proptest! {
        #[test]
        fn prop_test_rna_chunk_trim(
            r1_length in 41usize..,
            trim_r1 in any::<bool>(),
            r2_length in any::<usize>(),
            trim_r2 in any::<bool>()
        ) {
            // SC-Vdj chemistry
            {
                let expected_r1_length = if trim_r1 {
                    Some(min(150-41, r1_length-41))
                } else {
                    None
                };

                let expected_r2_length = if trim_r2 {
                    Some(min(150, r2_length))
                } else {
                    None
                };

                let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
                let umi_extractor = UmiExtractor::new(&chemistry.umi);
                let mut chunk = RnaChunk {
                    chemistry,
                    umi_extractor,
                    gem_group: 1,
                    fastqs: InputFastqs{
                        r1: "test/rna_read/interleaved_2k.fastq.lz4".into(),
                        r2: None,
                        i1: None,
                        i2: None,
                        r1_interleaved: true,
                    },
                    read_group: "Blah".into(),
                    subsample_rate: None,
                    library_type: LibraryType::VdjAuto,
                    read_lengths: TxHashMap::default(),
                    chunk_id: 0,
                    library_id: 0, fastq_id: None,
                };
                if trim_r1 {
                    chunk.set_illumina_r1_trim_length(r1_length);
                }
                if trim_r2 {
                    chunk.set_illumina_r2_trim_length(r2_length);
                }

                for rna_read_result in processor(chunk).iter().unwrap().take(100) {
                    let ProcessResult::Processed(rna_read) = rna_read_result.unwrap() else {
                        unreachable!();
                    };
                    assert_eq!(rna_read.bc_range, GelBeadOnly(RpRange::new(WhichRead::R1, 0, Some(16))));
                    assert_equal(rna_read.umi_ranges(), [RpRange::new(WhichRead::R1, 16, Some(10))]);
                    assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R1, 41, expected_r1_length));
                    assert_eq!(rna_read.r2_range, Some(RpRange::new(WhichRead::R2, 0, expected_r2_length)));
                }
            }
            // SC-VDJ R2 chemistry
            {
                let expected_r1_length = if trim_r2 {
                    Some(min(150, r2_length))
                } else {
                    None
                };

                let chemistry = ChemistryDef::named(ChemistryName::VdjR2);
                let umi_extractor = UmiExtractor::new(&chemistry.umi);
                let mut chunk = RnaChunk {
                    chemistry,
                    umi_extractor,
                    gem_group: 1,
                    fastqs: InputFastqs{
                        r1: "test/rna_read/interleaved_2k.fastq.lz4".into(),
                        r2: None,
                        i1: None,
                        i2: None,
                        r1_interleaved: true,
                    },
                    read_group: "Blah".into(),
                    subsample_rate: None,
                    library_type: LibraryType::VdjAuto,
                    read_lengths: TxHashMap::default(),
                    chunk_id: 0,
                    library_id: 0, fastq_id: None,
                };
                if trim_r1 {
                    chunk.set_illumina_r1_trim_length(r1_length);
                }
                if trim_r2 {
                    chunk.set_illumina_r2_trim_length(r2_length);
                }

                for rna_read_result in processor(chunk).iter().unwrap().take(100) {
                    let ProcessResult::Processed(rna_read) = rna_read_result.unwrap() else {
                        unreachable!();
                    };
                    assert_eq!(rna_read.bc_range, GelBeadOnly(RpRange::new(WhichRead::R1, 0, Some(16))));
                    assert_equal(rna_read.umi_ranges(), [RpRange::new(WhichRead::R1, 16, Some(10))]);
                    assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R2, 0, expected_r1_length));
                    assert_eq!(rna_read.r2_range, None);
                }
            }
        }
    }

    #[test]
    fn test_rna_read_adapter_trim() {
        // The VDJ adapters are listed in `test/rna_read/vdj_adapters.json`. The reads in `test/rna_read/interleaved_2k_insert.fastq` were
        // trimmed using cutadapt and saved in `test/rna_read/interleaved_trimmed_2k.fastq` (See `test/rna_read/run_cutadapt.sh`). This
        // test makes sure that the trimming that we produce exactly matches the cutadapt outputs. There are 2000 read pairs in total and 85
        // read pairs are trimmed in total.
        use fastq_set::adapter_trimmer::{Adapter, ReadAdapterCatalog};
        use fastq_set::read_pair_iter::ReadPairIter;
        use metric::SimpleHistogram;

        let vdj_adapters: TxHashMap<WhichRead, Vec<Adapter>> =
            serde_json::from_reader(File::open("test/rna_read/vdj_adapters.json").unwrap())
                .unwrap();
        let mut ad_catalog = ReadAdapterCatalog::from(&vdj_adapters);

        let fastqs = InputFastqs {
            // This file was created by running the script `test/rna_read/run_cutadapt.sh`
            // It is be the trimmed output fastq returned by cutadapt
            r1: "test/rna_read/interleaved_trimmed_2k.fastq".into(),
            r2: None,
            i1: None,
            i2: None,
            r1_interleaved: true,
        };

        // SCVDJ Chemistry
        {
            let rp_iter = ReadPairIter::from_fastq_files(&fastqs).unwrap();

            let chemistry = ChemistryDef::named(ChemistryName::VdjPE);
            let umi_extractor = UmiExtractor::new(&chemistry.umi);
            let chunk = RnaChunk {
                chemistry,
                umi_extractor,
                gem_group: 1,
                fastqs: InputFastqs {
                    r1: "test/rna_read/interleaved_2k.fastq".into(),
                    r2: None,
                    i1: None,
                    i2: None,
                    r1_interleaved: true,
                },
                read_group: "Blah".into(),
                subsample_rate: None,
                library_type: LibraryType::VdjAuto,
                read_lengths: TxHashMap::default(),
                chunk_id: 0,
                library_id: 0,
                fastq_id: None,
            };

            let mut n_trimmed = 0;
            let mut adapter_counts: SimpleHistogram<String> = SimpleHistogram::default();
            for (rna_read_result, rp_result) in processor(chunk).iter().unwrap().zip(rp_iter) {
                let ProcessResult::Processed(mut rna_read) = rna_read_result.unwrap() else {
                    unreachable!();
                };
                adapter_counts.extend_owned(rna_read.trim_adapters(&mut ad_catalog).into_keys());
                let rp = rp_result.unwrap();
                assert_eq!(
                    rna_read.r1_seq(),
                    rp.get(WhichRead::R1, ReadPart::Seq).unwrap()
                );
                assert_eq!(
                    rna_read.r2_seq().unwrap(),
                    rp.get(WhichRead::R2, ReadPart::Seq).unwrap()
                );
                if rna_read.r1_seq().len() != 109 || rna_read.r2_seq().unwrap().len() != 150 {
                    n_trimmed += 1;
                }
            }
            println!("Trimmed {n_trimmed} sequences");
            println!("Counts {adapter_counts:#?}");

            // The counts returned by cutadapt does not completely agree with this, because they count things
            // differently when multiple adapters are found for a read.
            assert_eq!(adapter_counts.get("R2_rc"), 27);
            assert_eq!(adapter_counts.get("P7_rc"), 8);
            assert_eq!(adapter_counts.get("polyA"), 8);
            assert_eq!(adapter_counts.get("rt_primer_rc"), 2);
            assert_eq!(adapter_counts.get("spacer"), 0);
            assert_eq!(adapter_counts.get("spacer_rc"), 66);
            assert_eq!(adapter_counts.get("R1_rc"), 29);
            assert_eq!(adapter_counts.get("P5_rc"), 14);
            assert_eq!(adapter_counts.get("polyT"), 4);
            assert_eq!(adapter_counts.get("rt_primer"), 1);
        }

        // SC-VDJ R2 chemistry
        {
            let rp_iter = ReadPairIter::from_fastq_files(&fastqs).unwrap();

            let chemistry = ChemistryDef::named(ChemistryName::VdjR2);
            let umi_extractor = UmiExtractor::new(&chemistry.umi);
            let chunk = RnaChunk {
                chemistry,
                umi_extractor,
                gem_group: 1,
                fastqs: InputFastqs {
                    r1: "test/rna_read/interleaved_2k.fastq".into(),
                    r2: None,
                    i1: None,
                    i2: None,
                    r1_interleaved: true,
                },
                read_group: "Blah".into(),
                subsample_rate: None,
                library_type: LibraryType::VdjAuto,
                read_lengths: TxHashMap::default(),
                chunk_id: 0,
                library_id: 0,
                fastq_id: None,
            };

            for (rna_read_result, rp_result) in processor(chunk).iter().unwrap().zip(rp_iter) {
                let ProcessResult::Processed(mut rna_read) = rna_read_result.unwrap() else {
                    unreachable!();
                };
                rna_read.trim_adapters(&mut ad_catalog);
                let rp = rp_result.unwrap();
                assert_eq!(
                    rna_read.r1_seq(),
                    rp.get(WhichRead::R2, ReadPart::Seq).unwrap()
                );
                assert_eq!(rna_read.r2_seq(), None);
            }
        }
    }

    #[test]
    fn test_library_info() {
        assert_eq!(
            serde_json::to_string(&LibraryInfo {
                library_id: 0,
                library_type: LibraryType::Antibody,
                gem_group: 1,
                target_set_name: None
            })
            .unwrap(),
            r#"{"library_id":0,"library_type":"Antibody Capture","gem_group":1,"target_set_name":null}"#
        );
    }

    #[test]
    fn test_umi_extraction_1() {
        let read1 = fastq::OwnedRecord {
            head: b"some_name".to_vec(),
            seq: b"CTGATGGCTCAAACACAGCGACCTCGGGTGGGAACACCTTGTTCAGGT".to_vec(),
            qual: b"GGGAGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
            sep: None,
        };
        let input = [Some(read1), None, None, None];
        let umi_parts = [UmiReadComponent {
            read_type: WhichRead::R1,
            offset: 16,
            length: 12,
            min_length: Some(10),
            whitelist: None,
        }];
        let rp = ReadPair::new(input);
        let extractor = UmiExtractor::new(&umi_parts);

        let (parts, seq) = extractor.extract_umi(&rp, &umi_parts).unwrap();
        assert_eq!(seq.as_bytes(), b"AGCGACCTCGGG");
        assert_eq!(
            parts.as_slice(),
            [UmiPart::Untranslated {
                range: RpRange::new(WhichRead::R1, 16, Some(12))
            }]
        );
    }

    #[test]
    fn test_umi_extraction_2() {
        let read1 = fastq::OwnedRecord {
            head: b"some_name".to_vec(),
            seq: b"CTGATGGCTCAAACACAGCGACCTCG".to_vec(),
            qual: b"GGGAGIIIIIIIIIIIIIIIIIIIII".to_vec(),
            sep: None,
        };
        let input = [Some(read1), None, None, None];
        let umi_parts = [UmiReadComponent {
            read_type: WhichRead::R1,
            offset: 16,
            length: 12,
            min_length: Some(10),
            whitelist: None,
        }];
        let rp = ReadPair::new(input);
        let extractor = UmiExtractor::new(&umi_parts);

        let (parts, seq) = extractor.extract_umi(&rp, &umi_parts).unwrap();
        assert_eq!(seq.as_bytes(), b"AGCGACCTCG");
        assert_eq!(
            parts.as_slice(),
            [UmiPart::Untranslated {
                range: RpRange::new(WhichRead::R1, 16, Some(10))
            }]
        );
    }

    #[test]
    fn test_umi_extraction_3() {
        let translator =
            |seqs: [&str; 4]| SplintToUmiTranslator::new(seqs.iter().map(|x| x.as_bytes()), 1);
        let extractor = UmiExtractor {
            translators: [
                Some(translator(["AATA", "CGCC", "GCGG", "TTAT"])),
                Some(translator(["AGCA", "CATC", "GTAG", "TCGT"])),
                Some(translator(["ACTC", "CGAT", "GTGA", "TACG"])),
                None,
            ]
            .into(),
        };
        let umi_parts = [
            UmiReadComponent {
                read_type: WhichRead::R1,
                offset: 7,
                length: 4,
                min_length: None,
                whitelist: Some(UmiWhitelistSpec {
                    slide: "test".to_string(),
                    part: slide_design::OligoPart::Bc1,
                    translation: UmiTranslation::SingleBase,
                }),
            },
            UmiReadComponent {
                read_type: WhichRead::R1,
                offset: 18,
                length: 4,
                min_length: None,
                whitelist: Some(UmiWhitelistSpec {
                    slide: "test".to_string(),
                    part: slide_design::OligoPart::Bc1,
                    translation: UmiTranslation::SingleBase,
                }),
            },
            UmiReadComponent {
                read_type: WhichRead::R1,
                offset: 29,
                length: 4,
                min_length: None,
                whitelist: Some(UmiWhitelistSpec {
                    slide: "test".to_string(),
                    part: slide_design::OligoPart::Bc1,
                    translation: UmiTranslation::SingleBase,
                }),
            },
            UmiReadComponent {
                read_type: WhichRead::R1,
                offset: 40,
                length: 8,
                min_length: None,
                whitelist: None,
            },
        ];

        {
            let rp = ReadPair::new([
                Some(fastq::OwnedRecord {
                    head: b"some_name".to_vec(),
                    seq: concat!(
                        "CTGATGG", "AATA", // UMI base 'A'
                        "AACACAG", "TCGA", // UMI base 'T'
                        "CTCGGGT", "CGAT", // UMI base 'C'
                        "ACACCTT", "GTTCAGGT", // Last 8 UMI bases
                    )
                    .into(),
                    qual: b"GGGAGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
                    sep: None,
                }),
                None,
                None,
                None,
            ]);

            let (parts, seq) = extractor.extract_umi(&rp, &umi_parts).unwrap();

            assert_eq!(seq.as_bytes(), b"ATCGTTCAGGT");
            assert_eq!(
                parts.as_slice(),
                [
                    UmiPart::SingleBaseTranslated {
                        range: RpRange::new(WhichRead::R1, 7, Some(4)),
                        base: b'A'
                    },
                    UmiPart::SingleBaseTranslated {
                        range: RpRange::new(WhichRead::R1, 18, Some(4)),
                        base: b'T'
                    },
                    UmiPart::SingleBaseTranslated {
                        range: RpRange::new(WhichRead::R1, 29, Some(4)),
                        base: b'C'
                    },
                    UmiPart::Untranslated {
                        range: RpRange::new(WhichRead::R1, 40, Some(8)),
                    },
                ]
            );
        }

        {
            let rp = ReadPair::new([
                Some(fastq::OwnedRecord {
                    head: b"some_name".to_vec(),
                    seq: [
                        "CTGATGG", "ACGA", // UMI base 'N'
                        "AACACAG", "TCGA", // UMI base 'T'
                        "CTCGGGT", "CGAT", // UMI base 'C'
                        "ACACCTT", "GTTCAGGT", // Last 8 UMI bases
                    ]
                    .iter()
                    .flat_map(|x| x.as_bytes())
                    .copied()
                    .collect(),
                    qual: b"GGGAGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
                    sep: None,
                }),
                None,
                None,
                None,
            ]);

            let (parts, seq) = extractor.extract_umi(&rp, &umi_parts).unwrap();

            assert_eq!(seq.as_bytes(), b"NTCGTTCAGGT");
            assert_eq!(
                parts.as_slice(),
                [
                    UmiPart::SingleBaseTranslated {
                        range: RpRange::new(WhichRead::R1, 7, Some(4)),
                        base: b'N'
                    },
                    UmiPart::SingleBaseTranslated {
                        range: RpRange::new(WhichRead::R1, 18, Some(4)),
                        base: b'T'
                    },
                    UmiPart::SingleBaseTranslated {
                        range: RpRange::new(WhichRead::R1, 29, Some(4)),
                        base: b'C'
                    },
                    UmiPart::Untranslated {
                        range: RpRange::new(WhichRead::R1, 40, Some(8))
                    },
                ]
            );
        }
    }

    #[test]
    fn test_library_info_deserialize() {
        let library_info: LibraryInfo = serde_json::from_str(
            r#"{"gem_group": 1, "library_id": "0", "library_type": "Gene Expression"}"#,
        )
        .unwrap();
        dbg!(library_info);
    }

    #[test]
    fn test_barcode_extraction_variable_multiplexing_barcode() {
        let gb_bc = "CGTAGCTGCTAGA";
        let probe_bc = "ACTGCTGA";

        let whitelist =
            make_seqs_into_whitelists(BarcodeConstruct::new_gel_bead_and_probe(gb_bc, probe_bc));

        let components = BarcodeConstruct::new_gel_bead_and_probe(
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::GelBead,
                offset: 0,
                length: gb_bc.len(),
                whitelist: dummy_whitelist_spec(),
            },
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::GelBead,
                offset: gb_bc.len(),
                length: probe_bc.len(),
                whitelist: dummy_whitelist_spec(),
            },
        );

        let extract = |read: &ReadPair, extraction: Option<&BarcodeExtraction>| {
            let (ranges, bc) = extract_barcode(
                read,
                components.as_ref(),
                extraction,
                whitelist.as_ref(),
                whitelist.as_ref().map(Whitelist::sequence_lengths).as_ref(),
                0,
            )
            .unwrap();
            (ranges, bc.segments())
        };

        // Create a single read with an extra base between the two barcodes.
        let read = make_read(&format!("{gb_bc}T{probe_bc}"));

        // Default extraction should fail.
        let (ranges, segments) = extract(&read, None);
        assert!(segments.gel_bead().is_valid());
        assert_eq!(components.as_ref().gel_bead().rprange(), ranges.gel_bead());
        // should have an invalid probe barcode with default range
        assert!(!segments.probe().is_valid());
        assert_eq!(components.as_ref().probe().rprange(), ranges.probe());

        // Extraction allowing for a 1 base offset should succeed.
        let (ranges, segments) = extract(
            &read,
            Some(&BarcodeExtraction::VariableMultiplexingBarcode {
                min_offset: 0,
                max_offset: 1,
            }),
        );
        assert!(segments.gel_bead().is_valid());
        assert_eq!(components.as_ref().gel_bead().rprange(), ranges.gel_bead());
        assert!(segments.probe().is_valid());
        assert_eq!(
            components.as_ref().probe().rprange().offset() + 1,
            ranges.probe().offset()
        );
    }

    fn make_read(seq: &str) -> ReadPair {
        let read1 = fastq::OwnedRecord {
            head: b"some_name".to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: vec![b'I'; seq.len()],
            sep: None,
        };
        let input = [Some(read1), None, None, None];
        ReadPair::new(input)
    }

    fn dummy_whitelist_spec() -> WhitelistSpec {
        WhitelistSpec::TxtFile {
            name: "custom".into(),
            translation: false,
            strand: ReqStrand::Forward,
        }
    }

    /// Create whitelists with a single item for each provided seq.
    fn make_seqs_into_whitelists(seqs: BarcodeConstruct<&str>) -> BarcodeConstruct<Whitelist> {
        seqs.map(|bc| Whitelist::Plain(set![BcSegSeq::from_bytes(bc.as_bytes())]))
    }

    fn test_barcode_extraction_joint(
        prefix: String,
        bc1: String,
        bc2: String,
        use_bc1: bool,
        use_bc2: bool,
    ) {
        let bc1_fill = if use_bc1 {
            Cow::Borrowed(bc1.as_str())
        } else {
            (0..bc1.len()).map(|_| 'N').collect()
        };
        let bc2_fill = if use_bc2 {
            Cow::Borrowed(bc2.as_str())
        } else {
            (0..bc2.len()).map(|_| 'N').collect()
        };
        let read = make_read(&format!("{prefix}{bc1_fill}{bc2_fill}TTCAGGT"));

        let whitelist = make_seqs_into_whitelists(BarcodeConstruct::Segmented(
            [bc1.as_str(), bc2.as_str()].into_iter().collect(),
        ));

        let components = BarcodeConstruct::Segmented(Segments {
            segment1: BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::SpotSegment,
                offset: 9,
                length: 19,
                whitelist: dummy_whitelist_spec(),
            },
            segment2: BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::SpotSegment,
                offset: 28,
                length: 18,
                whitelist: dummy_whitelist_spec(),
            },
            segment3: None,
            segment4: None,
        });
        let (bc_range, barcode) = extract_barcode(
            &read,
            components.as_ref(),
            Some(&BarcodeExtraction::JointBc1Bc2 {
                min_offset: 8,
                max_offset: 11,
            }),
            whitelist.as_ref(),
            whitelist.as_ref().map(Whitelist::sequence_lengths).as_ref(),
            1,
        )
        .unwrap();

        let default_range = components.as_ref().map(BarcodeReadComponent::rprange);

        match (use_bc1, use_bc2) {
            (true, true) => {
                assert_eq!(
                    barcode.segments(),
                    BarcodeConstruct::Segmented(Segments::from_iter([
                        BarcodeSegment::with_sequence(
                            bc1.as_bytes(),
                            BarcodeSegmentState::ValidBeforeCorrection
                        ),
                        BarcodeSegment::with_sequence(
                            bc2.as_bytes(),
                            BarcodeSegmentState::ValidBeforeCorrection
                        )
                    ]))
                );

                assert_eq!(
                    bc_range,
                    BarcodeConstruct::Segmented(Segments {
                        segment1: RpRange::new(WhichRead::R1, prefix.len(), Some(bc1.len())),
                        segment2: RpRange::new(
                            WhichRead::R1,
                            prefix.len() + bc1.len(),
                            Some(bc2.len())
                        ),
                        segment3: None,
                        segment4: None
                    })
                );
            }
            (false, false) => {
                assert_eq!(bc_range, default_range);
                assert_eq!(
                    barcode.segments(),
                    default_range.map(|r| BarcodeSegment::with_sequence(
                        read.get_range(r, ReadPart::Seq).unwrap(),
                        BarcodeSegmentState::Invalid
                    ))
                );
            }
            (true, false) => {
                assert_eq!(
                    bc_range.segments().segment1,
                    RpRange::new(WhichRead::R1, prefix.len(), Some(bc1.len()))
                );
                assert_eq!(
                    barcode.segments().segments().segment1,
                    BarcodeSegment::with_sequence(
                        bc1.as_bytes(),
                        BarcodeSegmentState::ValidBeforeCorrection
                    ),
                );
                assert_eq!(
                    barcode.segments().segments().segment2.state,
                    BarcodeSegmentState::Invalid,
                );
            }
            (false, true) => {
                assert_eq!(
                    bc_range.segments().segment2,
                    RpRange::new(WhichRead::R1, prefix.len() + bc1.len(), Some(bc2.len()))
                );
                assert_eq!(
                    barcode.segments().segments().segment2,
                    BarcodeSegment::with_sequence(
                        bc2.as_bytes(),
                        BarcodeSegmentState::ValidBeforeCorrection
                    ),
                );
                assert_eq!(
                    barcode.segments().segments().segment1.state,
                    BarcodeSegmentState::Invalid,
                );
            }
        }
    }

    proptest! {
        #[test]
        fn prop_test_join_barcode_extraction(
            umi in "[ACGT]{8, 11}",
            bc1 in "[ACGT]{15, 16}",
            bc2 in "[ACGT]{15, 16}",
            use_bc1 in any::<bool>(),
            use_bc2 in any::<bool>(),
        ) {
            test_barcode_extraction_joint(
                umi,
                format!("{bc1}GGG"),
                format!("{bc2}GGG"),
                use_bc1,
                use_bc2,
            );
        }
    }
}
