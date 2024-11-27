use anyhow::{Context, Result};
use cr_types::reference::genome_of_chrom::genome_of_chrom;
use cr_types::{utils, GenomeName, ReqStrand};
use fastq_set::WhichEnd;
use martian_filetypes::tabular_file::TsvFileNoHeader;
use martian_filetypes::FileTypeRead;
use metric::AsMetricPrefix;
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::{cmp, str};
use strum_macros::{Display, IntoStaticStr};
use transcriptome::{Gene, TranscriptIndex, Transcriptome};

#[derive(Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct PairAnnotationData {
    pub genes: Vec<Gene>,
}

impl PairAnnotationData {
    /// Annotate a pair of alignments
    /// Take the intersection of the non-empty gene sets of the mates
    pub fn from_pair(anno1: &AnnotationData, anno2: &AnnotationData) -> PairAnnotationData {
        let genes = match (!anno1.genes.is_empty(), !anno2.genes.is_empty()) {
            (true, false) => anno1.genes.clone(),
            (false, true) => anno2.genes.clone(),
            _ if anno1.genome == anno2.genome => utils::intersect_vecs(&anno1.genes, &anno2.genes),
            _ => vec![],
        };
        PairAnnotationData { genes }
    }

    pub fn make_gx_gn_tags(&self) -> Option<(String, String)> {
        format_gx_gn_tags(&self.genes)
    }
}

#[derive(Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct AnnotationData {
    pub transcripts: Vec<TranscriptAlignment>,
    pub antisense: Vec<TranscriptAlignment>,
    pub genes: Vec<Gene>,
    pub region: AnnotationRegion,
    pub rescued: bool,
    /// Name of the genome this annotation corresponds to
    pub genome: GenomeName,
}

impl AnnotationData {
    fn new(genome: GenomeName) -> AnnotationData {
        AnnotationData {
            transcripts: Vec::new(),
            antisense: Vec::new(),
            genes: Vec::new(),
            region: AnnotationRegion::Intergenic,
            rescued: false,
            genome,
        }
    }

    pub fn is_antisense(&self) -> bool {
        self.transcripts.is_empty() && !self.antisense.is_empty()
    }

    pub fn is_sense(&self) -> bool {
        !self.transcripts.is_empty() && self.antisense.is_empty()
    }

    pub fn make_tx_tag(&self) -> Option<String> {
        if !self.transcripts.is_empty() {
            Some(
                self.transcripts
                    .iter()
                    .map(TranscriptAlignment::to_tag)
                    .collect::<Vec<String>>()
                    .join(";"),
            )
        } else {
            None
        }
    }

    pub fn make_an_tag(&self) -> Option<String> {
        if !self.antisense.is_empty() {
            Some(
                self.antisense
                    .iter()
                    .map(TranscriptAlignment::to_tag)
                    .collect::<Vec<String>>()
                    .join(";"),
            )
        } else {
            None
        }
    }

    pub fn make_re_tag(&self) -> Option<&'static str> {
        self.region.to_tag()
    }

    pub fn make_gx_gn_tags(&self) -> Option<(String, String)> {
        format_gx_gn_tags(&self.genes)
    }

    pub fn make_mm_tag(&self) -> Option<i32> {
        if self.rescued {
            Some(1)
        } else {
            None
        }
    }
}

fn format_gx_gn_tags(genes: &[Gene]) -> Option<(String, String)> {
    if !genes.is_empty() {
        let mut gx = Vec::new();
        let mut gn = Vec::new();
        for gene in genes {
            gx.push(gene.id.clone());
            gn.push(gene.name.clone());
        }
        Some((gx.join(";"), gn.join(";")))
    } else {
        None
    }
}

// These quantities are well defined for a valid transcriptomic alignment
#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Serialize, Deserialize)]
pub struct TxAlignProperties {
    pub id: String,
    pub pos: i64,
    pub cigar: String,
    pub alen: i64,
    pub se_insert_size: u32,
}

/*
Describes the transcript alignment of a non-intergenic genomic alignment
gene:     gene that is associated with the transcript
strand:   whether the genomic alignment is in the same orientation as the
          transcript or in reverse
tx_align: defined when the alignment is compatible with the transcript splicing
          pattern. An alignment is not compatible if it is:
          (1) intronic
          (2) exonic but does not agree with annotated isoforms
*/

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Serialize, Deserialize)]
pub struct TranscriptAlignment {
    pub gene: Gene,
    pub strand: ReqStrand,
    pub tx_align: Option<TxAlignProperties>,
}

impl TranscriptAlignment {
    // Note: when intron counting mode is enabled the TX (and AN) tags for an intronic alignment
    // consist of the gene ID and the strand. While for transcriptomic alignments they are in the
    // usual format of transcript ID, strand, position, cigar
    fn to_tag(&self) -> String {
        let strand = if self.strand == ReqStrand::Forward {
            '+'
        } else {
            '-'
        };
        match &self.tx_align {
            Some(tal) => format!("{},{strand}{},{}", tal.id, tal.pos, tal.cigar),
            None => format!("{},{strand}", self.gene.id),
        }
    }
    // in intron mode, identifies alignments that are still transcript compatible, from those
    // that are intronic/exonic
    pub fn is_txomic(&self) -> bool {
        self.tx_align.is_some()
    }
}

#[derive(
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Debug,
    Clone,
    Copy,
    Hash,
    Display,
    IntoStaticStr,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "snake_case")]
pub enum AnnotationRegion {
    Exonic,
    Intronic,
    Intergenic,
}

impl AsMetricPrefix for AnnotationRegion {
    fn as_metric_prefix(&self) -> Option<&'static str> {
        Some(self.into())
    }
}

impl AnnotationRegion {
    fn to_tag(self) -> Option<&'static str> {
        match self {
            AnnotationRegion::Exonic => Some("E"),
            AnnotationRegion::Intronic => Some("N"),
            AnnotationRegion::Intergenic => Some("I"),
        }
    }
}

/// A STAR transcript.
/// The order of these fields must match the order of the columns in the file
/// transcriptInfo.tab produced by STAR, which has no header.
#[derive(Deserialize, Serialize)]
pub struct StarTranscript {
    pub(crate) id: String,
    pub(crate) start: i64,
    pub(crate) end: i64, // inclusive
    pub(crate) max_end: i64,
    pub(crate) strand: u8,
    pub(crate) num_exons: usize,
    pub(crate) break_idx: usize, // TODO rename
}

/// A STAR exon.
/// The order of these fields must match the order of the columns in the file
/// exonInfo.tab produced by STAR, which has no header.
#[derive(Deserialize, Serialize)]
struct StarExon {
    start: i64,
    end: i64,     // inclusive
    cum_len: i64, // cumulative length of previous exons, not including this one
}

impl StarExon {
    fn len(&self) -> i64 {
        self.end - self.start + 1
    }
}

#[derive(Clone)]
pub struct AnnotationParams {
    pub chemistry_strandedness: ReqStrand,
    pub chemistry_endedness: WhichEnd,
    pub intergenic_trim_bases: i64,
    pub intronic_trim_bases: i64,
    pub junction_trim_bases: i64,
    pub region_min_overlap: f64,
    pub include_exons: bool,
    pub include_introns: bool,
}

#[derive(Debug)]
struct SpliceSegment {
    start: i64,
    end: i64,
    cigar: Vec<Cigar>,
}

pub struct TranscriptAnnotator {
    params: AnnotationParams,
    transcript_info: Vec<StarTranscript>,
    exon_info: Vec<StarExon>,
    chrom_starts: Vec<i64>,
    transcript_index: TranscriptIndex,
    /// The genome corresponding to each chromosome. Given a read alignment, we can readily
    /// find the genome by using the `read.tid()` index into this vector
    genome_of_tid: Vec<GenomeName>,
}

/// Read a STAR tab TSV file with no header. Skip the first line, which is the number of records.
fn read_star_tab<T: DeserializeOwned + Serialize>(path: &Path) -> Result<Vec<T>> {
    let mut reader = BufReader::new(File::open(path).with_context(|| path.display().to_string())?);
    let mut line = String::new();
    reader.read_line(&mut line)?;
    TsvFileNoHeader::read_from(reader)
}

impl TranscriptAnnotator {
    pub fn new(reference_path: &Path, params: AnnotationParams) -> Result<TranscriptAnnotator> {
        let transcript_info: Vec<StarTranscript> =
            read_star_tab(&reference_path.join("star/transcriptInfo.tab"))?;
        let exon_info: Vec<StarExon> = read_star_tab(&reference_path.join("star/exonInfo.tab"))?;
        let chrom_starts = cr_types::utils::load_txt(&reference_path.join("star/chrStart.txt"))?;
        let genome_of_tid = genome_of_chrom(reference_path)?;
        let txome = Transcriptome::from_reference_path(reference_path)?;
        let transcript_index = TranscriptIndex::from_transcriptome(&txome);

        Ok(TranscriptAnnotator {
            params,
            transcript_info,
            exon_info,
            chrom_starts,
            transcript_index,
            genome_of_tid,
        })
    }

    pub fn get_params(&self) -> &AnnotationParams {
        &self.params
    }

    pub fn annotate_alignment(&self, read: &Record) -> AnnotationData {
        assert!(
            !read.is_unmapped(),
            "Unmapped alignments cannot be annotated"
        );
        let genome = self.genome_of_tid[read.tid() as usize].clone();
        let mut annotation_data = AnnotationData::new(genome);

        let ref_offset = self.chrom_starts[read.tid() as usize];
        let clipped_read_start = ref_offset + read.pos();
        let clipped_read_end = clipped_read_start + cr_bam::bam::alen(read) - 1;

        // find maximum index of transcripts that could possibly overlap the read
        let mut tx_idx = cr_types::utils::bisect(
            &self.transcript_info,
            clipped_read_end,
            &|tx| tx.start,
            utils::BisectDirection::Right,
        );
        //let mut tx_idx = transcript_info.binary_search_by_key(&clipped_read_end, |tx| tx.start).unwrap();

        if tx_idx == 0 {
            // read ends before the first transcript start
            return annotation_data;
        }
        // first transcript on the left
        tx_idx -= 1;

        // cycle back through overlapping transcripts
        let mut curr_transcript = &self.transcript_info[tx_idx];
        let mut any_exonic = false;
        let mut any_intronic = false;
        let mut alignments = Vec::new();
        while clipped_read_start <= cmp::max(curr_transcript.end, curr_transcript.max_end) {
            if clipped_read_start <= curr_transcript.end {
                // transcript overlaps the read
                // if we found a transcript alignment, then skip processing of intronic alignments
                let (curr_aln, curr_region) =
                    self.align_to_transcript(read, ref_offset, curr_transcript);
                match curr_region {
                    AnnotationRegion::Exonic => {
                        any_exonic = true;
                        if self.params.include_exons {
                            alignments.push(curr_aln.unwrap());
                        }
                    }
                    AnnotationRegion::Intronic => {
                        any_intronic = true;
                        if self.params.include_introns {
                            alignments.push(curr_aln.unwrap());
                        };
                    }
                    AnnotationRegion::Intergenic => (),
                };
            }

            if tx_idx == 0 {
                break;
            }
            tx_idx -= 1;
            curr_transcript = &self.transcript_info[tx_idx];
        }
        let mut seen_genes = HashSet::new();
        let mut transcripts = BTreeMap::new();
        let mut antisense = BTreeMap::new();
        if any_exonic {
            annotation_data.region = AnnotationRegion::Exonic;
            // Check if there are transcriptome compatible alignments

            let txome_exists = alignments.iter().any(|x| x.tx_align.is_some());
            let count_intronic = !txome_exists && self.params.include_introns;
            for aln in alignments.into_iter().rev() {
                match (aln.strand, &aln.tx_align, count_intronic) {
                    (ReqStrand::Forward, Some(tx_align), _) => {
                        // Transcript sense alignment
                        seen_genes.insert(aln.gene.clone());
                        transcripts.insert(tx_align.id.clone(), aln);
                    }
                    (ReqStrand::Reverse, Some(tx_align), _) => {
                        // Transcript anti-sense alignment
                        antisense.insert(tx_align.id.clone(), aln);
                    }
                    // Exonic but not transcriptomic, only count if no transcriptomic alignment
                    // exists and we are in include_introns mode
                    (ReqStrand::Forward, None, true) => {
                        seen_genes.insert(aln.gene.clone());
                        transcripts.insert(aln.gene.id.clone(), aln);
                    }
                    (ReqStrand::Reverse, None, true) => {
                        antisense.insert(aln.gene.id.clone(), aln);
                    }
                    _ => (),
                }
            }
        } else if any_intronic {
            annotation_data.region = AnnotationRegion::Intronic;
            if self.params.include_introns {
                for aln in alignments.into_iter().rev() {
                    assert!(
                        aln.tx_align.is_none(),
                        "any_exonic was False, but encountered tx alignment"
                    );
                    match aln.strand {
                        ReqStrand::Forward => {
                            seen_genes.insert(aln.gene.clone());
                            transcripts.insert(aln.gene.id.clone(), aln);
                        }
                        ReqStrand::Reverse => {
                            antisense.insert(aln.gene.id.clone(), aln);
                        }
                    }
                }
            }
        } else {
            annotation_data.region = AnnotationRegion::Intergenic;
        };
        annotation_data.transcripts = transcripts.into_values().collect::<Vec<_>>();
        annotation_data.antisense = antisense.into_values().collect::<Vec<_>>();
        annotation_data.genes = seen_genes.into_iter().collect::<Vec<Gene>>();
        // Sorting this makes life easier later.
        annotation_data.genes.sort_unstable();

        annotation_data
    }

    fn align_to_transcript(
        &self,
        read: &Record,
        ref_offset: i64,
        transcript_info: &StarTranscript,
    ) -> (Option<TranscriptAlignment>, AnnotationRegion) {
        // figure out coordinates
        let tx_start = transcript_info.start;
        let tx_end = transcript_info.end;
        let tx_first_exon = transcript_info.break_idx;
        let tx_num_exons = transcript_info.num_exons;
        let tx_last_exon = tx_first_exon + tx_num_exons - 1; // inclusive
        let star_genomic_start = ref_offset + read.pos();
        let star_genomic_end = star_genomic_start + cr_bam::bam::alen(read) - 1; // inclusive
        let clipped_read_start = star_genomic_start - tx_start; // convert to transcript coordinates
        let clipped_read_end = star_genomic_end - tx_start;

        // compute region
        let (left_clip, right_clip, splice_segments) = get_cigar_segments(read, clipped_read_start);
        let is_exonic = is_read_exonic(
            &splice_segments,
            &self.exon_info,
            tx_first_exon,
            tx_last_exon,
            self.params.region_min_overlap,
        );
        let is_intronic = !is_exonic
            && get_overlap(star_genomic_start, star_genomic_end, tx_start, tx_end) >= 1.0;
        let region = if is_exonic {
            AnnotationRegion::Exonic
        } else if is_intronic {
            AnnotationRegion::Intronic
        } else {
            // If neither exonic nor intronic => intergenic.
            return (None, AnnotationRegion::Intergenic);
        };
        // compute strand
        let tx_reverse_strand = transcript_info.strand == 2;
        let mut read_reverse_strand = read.is_reverse();
        if read.is_paired() && read.is_last_in_template() {
            read_reverse_strand = !read_reverse_strand;
        };
        let is_antisense = match self.params.chemistry_strandedness {
            ReqStrand::Forward => tx_reverse_strand != read_reverse_strand,
            ReqStrand::Reverse => tx_reverse_strand == read_reverse_strand,
        };

        let tx_strand = if is_antisense {
            ReqStrand::Reverse
        } else {
            ReqStrand::Forward
        };

        // Prior to aligning to the transcript we have a gene alignment
        let gene = self
            .transcript_index
            .get_gene_from_transcript(&transcript_info.id)
            .clone();
        let gene_aln = TranscriptAlignment {
            gene: gene.clone(),
            strand: tx_strand,
            tx_align: None,
        };
        // if intronic return region and strand
        if is_intronic {
            return (Some(gene_aln), region);
        };

        // must be exonic to continue
        assert!(is_exonic);

        // find the beginning / ending exons
        let Some((ex_start, ex_end)) = find_exons(
            &self.exon_info,
            clipped_read_start,
            clipped_read_end,
            tx_first_exon,
            tx_last_exon,
            self.params.intergenic_trim_bases,
            self.params.intronic_trim_bases,
        ) else {
            return (Some(gene_aln), region);
        };

        // compute offsets
        let ex_offset = cmp::max(clipped_read_start - self.exon_info[ex_start].start, 0);
        let mut tx_offset = self.exon_info[ex_start].cum_len + ex_offset;
        let tx_length = self.exon_info[tx_last_exon].cum_len + self.exon_info[tx_last_exon].len();

        // align the read to the exons
        let tx_alignment = align_junctions(
            &left_clip,
            &right_clip,
            &splice_segments,
            &self.exon_info[ex_start..ex_end + 1],
            self.params.junction_trim_bases,
        );

        // discard misaligned reads
        let Some((mut tx_cigar, tx_aligned_bases)) = tx_alignment else {
            return (Some(gene_aln), region);
        };

        // flip reverse strand
        if tx_reverse_strand {
            tx_offset = tx_length - (tx_offset + tx_aligned_bases);
            tx_cigar.reverse();
        };

        // Apparent insert size, assuming a SE read coming from this transcript
        let se_insert_size = match self.params.chemistry_endedness {
            WhichEnd::FivePrime => tx_offset + tx_aligned_bases,
            WhichEnd::ThreePrime => {
                let tx_len = self
                    .transcript_index
                    .get_transcript_length(&transcript_info.id);
                tx_len - tx_offset
            }
        } as u32;

        let alignment = TranscriptAlignment {
            gene,
            strand: tx_strand,
            tx_align: Some(TxAlignProperties {
                id: transcript_info.id.clone(),
                pos: tx_offset,
                cigar: CigarString(tx_cigar).to_string(),
                alen: tx_aligned_bases,
                se_insert_size,
            }),
        };

        (Some(alignment), region)
    }
}

fn is_read_exonic(
    splice_segments: &[SpliceSegment],
    exon_info: &[StarExon],
    first_exon: usize,
    last_exon: usize,
    min_overlap_frac: f64,
) -> bool {
    for segment in splice_segments {
        // find first exon that ends to the right of the segment start
        let ex_idx = utils::bisect(
            &exon_info[first_exon..last_exon + 1],
            segment.start,
            &|ex| ex.end,
            utils::BisectDirection::Right,
        ) + first_exon;
        //let ex_idx = exon_info.binary_search_by_key(&segment.start, |ex| ex.end - 1).unwrap() as u64;
        if ex_idx > last_exon {
            // read is out of bounds
            return false;
        }
        let exon = &exon_info[ex_idx];
        if get_overlap(segment.start, segment.end, exon.start, exon.end) < min_overlap_frac {
            return false;
        }
    }
    true
}

fn find_exons(
    exon_info: &[StarExon],
    read_start: i64,
    read_end: i64,
    first_exon: usize,
    last_exon: usize,
    intergenic_trim_bases: i64,
    intronic_trim_bases: i64,
) -> Option<(usize, usize)> {
    let ex_start = cr_types::utils::bisect(
        &exon_info[first_exon..last_exon + 1],
        read_start,
        &|ex| ex.end,
        utils::BisectDirection::Right,
    ) + first_exon; // find first exon that ends to the right of the read start
    let ex_end = cr_types::utils::bisect(
        &exon_info[first_exon..last_exon + 1],
        read_end,
        &|ex| ex.start,
        cr_types::utils::BisectDirection::Left,
    ) - 1
        + first_exon; // find first exon that starts to the left of the read end
                      //let ex_start = exon_info.binary_search_by_key(&read_start, |ex| ex.end - 1).unwrap() as u64; // find first exon that ends to the right of the read start
                      //let ex_end = exon_info.binary_search_by_key(&read_end, |ex| ex.start - 1).unwrap() as u64 - 1; // find first exon that starts to the left of the read end

    if (ex_start > last_exon) | (ex_end < first_exon) {
        // read is out of bounds
        return None;
    }

    let starting_exon = &exon_info[ex_start];
    let ending_exon = &exon_info[ex_end];

    if read_start < starting_exon.start {
        // read overhangs exon on the left
        let overhang = starting_exon.start - read_start;
        let trim_bases = if ex_start < first_exon {
            intergenic_trim_bases
        } else {
            intronic_trim_bases
        };
        if overhang > trim_bases {
            // too much overhang
            return None;
        };
    }

    if read_end > ending_exon.end {
        // read overhangs exon on the right
        let overhang = read_end - ending_exon.end;
        let trim_bases = if ex_end > last_exon {
            intergenic_trim_bases
        } else {
            intronic_trim_bases
        };
        if overhang > trim_bases {
            // too much overhang
            return None;
        };
    }

    Some((ex_start, ex_end))
}

/// Fraction of read interval covered by ref interval
fn get_overlap(read_start: i64, read_end: i64, ref_start: i64, ref_end: i64) -> f64 {
    let mut overlap_bases = cmp::min(ref_end, read_end) - cmp::max(ref_start, read_start);
    if overlap_bases < 0 {
        overlap_bases = 0;
    }
    (overlap_bases as f64) / ((read_end - read_start) as f64)
}

/// Segment the cigar string into the three segments: the initial clipping, the spliced alignment, and the final clipping.
/// The spliced alignment is returned as `Vec<SpliceSegment>`. Each `SpliceSegment` represents a block of cigar operations
/// not containing any RefSkip sections, and is annotated with the reference start/end position of the block.
fn get_cigar_segments(
    read: &Record,
    alignment_start: i64,
) -> (Vec<Cigar>, Vec<Cigar>, Vec<SpliceSegment>) {
    let mut left_clip: Vec<Cigar> = Vec::new();
    let mut right_clip: Vec<Cigar> = Vec::new();
    let mut splice_segments: Vec<SpliceSegment> = Vec::new();
    let mut seen_nonclips = false; // whether we've seen non-clip bases yet
    let mut curr_segment = SpliceSegment {
        start: alignment_start,
        end: alignment_start,
        cigar: Vec::new(),
    };

    for c in &read.cigar() {
        match c {
            Cigar::HardClip(_) | Cigar::SoftClip(_) => {
                if seen_nonclips {
                    right_clip.push(*c);
                } else {
                    left_clip.push(*c);
                }
            }
            Cigar::RefSkip(_) => {
                seen_nonclips = true;
                let next_start = curr_segment.end + c.len() as i64;
                splice_segments.push(curr_segment);
                curr_segment = SpliceSegment {
                    start: next_start,
                    end: next_start,
                    cigar: Vec::new(),
                };
            }
            Cigar::Ins(_) => {
                seen_nonclips = true;
                curr_segment.cigar.push(*c);
            }
            Cigar::Match(_) | Cigar::Del(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                seen_nonclips = true;
                curr_segment.end += c.len() as i64;
                curr_segment.cigar.push(*c);
            }
            Cigar::Pad(_) => unreachable!(),
        }
    }

    splice_segments.push(curr_segment);
    (left_clip, right_clip, splice_segments)
}

fn align_junctions(
    left_clip: &[Cigar],
    right_clip: &[Cigar],
    splice_segments: &[SpliceSegment],
    exon_info: &[StarExon],
    tolerance: i64,
) -> Option<(Vec<Cigar>, i64)> {
    if splice_segments.len() != exon_info.len() {
        return None;
    }
    let mut full_cigar = Vec::new();
    let mut aligned_bases = 0;
    full_cigar.extend_from_slice(left_clip);
    for i in 0..splice_segments.len() {
        let curr_segment = &splice_segments[i];
        let curr_exon = &exon_info[i];
        aligned_bases += curr_exon.len();
        let mut tmp_cigar = Vec::new();
        tmp_cigar.extend_from_slice(&curr_segment.cigar);

        // align the start
        let start_diff = curr_exon.start - curr_segment.start;
        if i == 0 {
            // first segment
            if start_diff > 0 {
                // overhang -> softclip
                tmp_cigar =
                    mask_read_bases(&mut tmp_cigar, Cigar::SoftClip(start_diff as u32), false);
            } else if start_diff < 0 {
                // underhang -> decrement aligned bases
                aligned_bases -= start_diff.abs();
            }
        } else if start_diff.abs() > tolerance {
            return None; // can't align properly
        } else if start_diff > 0 {
            // overhang -> mark as insertion
            tmp_cigar = mask_read_bases(&mut tmp_cigar, Cigar::Ins(start_diff as u32), false);
        } else if start_diff < 0 {
            // underhang -> mark as deletion
            tmp_cigar = mark_deleted_ref_bases(
                &mut tmp_cigar,
                start_diff.unsigned_abs().try_into().unwrap(),
                false,
            );
        }

        // align the end
        let end_diff = curr_segment.end - curr_exon.end - 1;
        if i == splice_segments.len() - 1 {
            // last segment
            if end_diff > 0 {
                // overhang -> softclip
                tmp_cigar = mask_read_bases(&mut tmp_cigar, Cigar::SoftClip(end_diff as u32), true);
            } else if end_diff < 0 {
                // underhang -> decrement aligned bases
                aligned_bases -= end_diff.abs();
            }
        } else if end_diff.abs() > tolerance {
            return None; // can't align properly
        } else if end_diff > 0 {
            // overhang -> mark as insertion
            tmp_cigar = mask_read_bases(&mut tmp_cigar, Cigar::Ins(end_diff as u32), true);
        } else if end_diff < 0 {
            // underhang -> mark as deletion
            tmp_cigar = mark_deleted_ref_bases(
                &mut tmp_cigar,
                end_diff.unsigned_abs().try_into().unwrap(),
                true,
            );
        }

        // extend
        full_cigar = extend_cigar(&full_cigar, &tmp_cigar);
    }
    full_cigar = extend_cigar(&full_cigar, right_clip);

    Some((full_cigar, aligned_bases))
}

fn extend_cigar(old_cigar: &[Cigar], new_cigar: &[Cigar]) -> Vec<Cigar> {
    // extend list of cigar ops, checking the ends to see if they should be merged
    let old_len = old_cigar.len();
    let new_len = new_cigar.len();
    let mut merged_cigar = Vec::new();
    if (old_len > 0) && (new_len > 0) {
        let old_tail = old_cigar[old_len - 1];
        let new_head = new_cigar[0];
        if old_tail.char() == new_head.char() {
            let merged_ends = set_cigar_len(&new_head, old_tail.len() + new_head.len());
            merged_cigar.extend_from_slice(&old_cigar[..old_len - 1]);
            merged_cigar.push(merged_ends);
            merged_cigar.extend_from_slice(&new_cigar[1..]);
        } else {
            merged_cigar.extend_from_slice(old_cigar);
            merged_cigar.extend_from_slice(new_cigar);
        }
    } else {
        merged_cigar.extend_from_slice(old_cigar);
        merged_cigar.extend_from_slice(new_cigar);
    }
    merged_cigar
}

fn mark_deleted_ref_bases(cigar: &mut [Cigar], del_len: u32, reverse: bool) -> Vec<Cigar> {
    let mut new_cigar = Vec::new();
    let del = Cigar::Del(del_len);
    if reverse {
        new_cigar.push(del);
        new_cigar.extend_from_slice(cigar);
    } else {
        new_cigar.extend_from_slice(cigar);
        new_cigar.push(del);
    };
    new_cigar
}

fn set_cigar_len(cig: &Cigar, len: u32) -> Cigar {
    use rust_htslib::bam::record::Cigar::{
        Del, Diff, Equal, HardClip, Ins, Match, Pad, RefSkip, SoftClip,
    };

    match cig {
        Match(_) => Match(len),
        Ins(_) => Ins(len),
        Del(_) => Del(len),
        RefSkip(_) => RefSkip(len),
        SoftClip(_) => SoftClip(len),
        HardClip(_) => HardClip(len),
        Pad(_) => Pad(len),
        Equal(_) => Equal(len),
        Diff(_) => Diff(len),
    }
}

fn mask_read_bases(cigar: &mut [Cigar], mask: Cigar, reverse: bool) -> Vec<Cigar> {
    // NOTE: this assumes that refskips have been removed
    let mut new_cigar = Vec::new();
    let mask_len = mask.len();
    let mut consumed_bases = 0;
    new_cigar.push(mask);
    if reverse {
        cigar.reverse();
    }
    for c in cigar {
        if consumed_bases < mask_len {
            // this op should be masked
            let read_bases = match c {
                Cigar::Del(_) => 0, // deletions don't consume read bases
                _ => c.len(),
            };
            if consumed_bases + read_bases >= mask_len {
                let truncated = set_cigar_len(c, read_bases + consumed_bases - mask_len);
                new_cigar.push(truncated);
            };
            consumed_bases += read_bases;
        } else {
            // just copy the op
            new_cigar.push(*c);
        };
    }
    if reverse {
        new_cigar.reverse();
    }
    new_cigar
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use std::collections::HashMap;

    #[allow(dead_code)]
    struct TranscriptomeTest {
        chrom_starts: Vec<i64>,
        transcript_info: Vec<StarTranscript>,
        exon_info: Vec<StarExon>,
        transcript_index: TranscriptIndex,
    }

    fn make_test_annotator(params: AnnotationParams) -> TranscriptAnnotator {
        /*
        Build a transcriptome. Things to test:
        - transcripts w/ single or multiple exons
        - transcripts that start / end at same pos
        - transcript that where max_end > end
        - pos / neg strand

        Transcript layout:
            00...05...10...15...20...25...30...35...40...45...50...55
        tx0 >>>>>>>>>>
        tx1 >>>>>>>>>>------------------------->>>>>>>>>>>>>>>
        tx2                          <<<<<-----<<<<<<<<<<

        NOTE: that STAR uses closed intervals (so transcript / exon ends are inclusive)
        */
        let transcript_info = vec![
            StarTranscript {
                id: "tx0".into(),
                start: 0,
                end: 9,
                max_end: 9,
                strand: 1,
                num_exons: 1,
                break_idx: 0,
            },
            StarTranscript {
                id: "tx1".into(),
                start: 0,
                end: 49,
                max_end: 49,
                strand: 1,
                num_exons: 2,
                break_idx: 1,
            },
            StarTranscript {
                id: "tx2".into(),
                start: 25,
                end: 44,
                max_end: 49,
                strand: 2,
                num_exons: 2,
                break_idx: 3,
            },
        ];

        let exon_info = vec![
            // tx0
            StarExon {
                start: 0,
                end: 9,
                cum_len: 0,
            },
            // tx1
            StarExon {
                start: 0,
                end: 9,
                cum_len: 0,
            },
            StarExon {
                start: 35,
                end: 49,
                cum_len: 10,
            },
            // tx2
            StarExon {
                start: 0,
                end: 4,
                cum_len: 0,
            },
            StarExon {
                start: 10,
                end: 19,
                cum_len: 5,
            },
        ];

        let chrom_starts = vec![0];
        let genome_of_tid = vec!["GRCh38".into()];

        let mut transcript_genes = HashMap::new();
        transcript_genes.insert(
            "tx0".into(),
            Gene {
                id: "gx0".into(),
                name: "gene0".into(),
            },
        );
        transcript_genes.insert(
            "tx1".into(),
            Gene {
                id: "gx1".into(),
                name: "gene1".into(),
            },
        );
        transcript_genes.insert(
            "tx2".into(),
            Gene {
                id: "gx2".into(),
                name: "gene2".into(),
            },
        );

        let mut transcript_lengths = HashMap::new();
        transcript_lengths.insert("tx0".into(), 1000);
        transcript_lengths.insert("tx1".into(), 2000);
        transcript_lengths.insert("tx2".into(), 2500);

        let transcript_index = TranscriptIndex {
            transcript_genes,
            transcript_lengths,
        };

        TranscriptAnnotator {
            chrom_starts,
            transcript_info,
            exon_info,
            transcript_index,
            params,
            genome_of_tid,
        }
    }

    fn check_annotation(
        rec: &Record,
        annotator: &TranscriptAnnotator,
        expected_transcripts: usize,
        expected_antisense: usize,
        expected_genes: Option<Vec<&str>>,
        expected_region: AnnotationRegion,
    ) {
        let result = annotator.annotate_alignment(rec);
        println!("{result:?}");
        assert_eq!(result.transcripts.len(), expected_transcripts);
        assert_eq!(result.antisense.len(), expected_antisense);
        if let Some(exp_genes) = expected_genes {
            assert_eq!(result.genes.len(), exp_genes.len());
            for (i, gene) in result.genes.iter().enumerate() {
                assert_eq!(gene.id, exp_genes[i]);
            }
        }
        assert_eq!(result.region, expected_region);
    }

    fn default_params() -> AnnotationParams {
        AnnotationParams {
            chemistry_strandedness: ReqStrand::Forward,
            chemistry_endedness: WhichEnd::ThreePrime,
            intergenic_trim_bases: 0,
            intronic_trim_bases: 0,
            junction_trim_bases: 0,
            region_min_overlap: 0.5,
            include_exons: true,
            include_introns: false,
        }
    }

    fn include_introns_params() -> AnnotationParams {
        AnnotationParams {
            chemistry_strandedness: ReqStrand::Forward,
            chemistry_endedness: WhichEnd::ThreePrime,
            intergenic_trim_bases: 0,
            intronic_trim_bases: 0,
            junction_trim_bases: 0,
            region_min_overlap: 0.5,
            include_exons: true,
            include_introns: true,
        }
    }

    #[test]
    fn test_include_introns() {
        // annotator without introns
        let tx_def = make_test_annotator(default_params());
        // annotator with introns
        let tx_int = make_test_annotator(include_introns_params());
        /*
        Transcript layout:
                    00...05...10...15...20...25...30...35...40...45...50...55
                tx0 >>>>>>>>>>
                tx1 >>>>>>>>>>------------------------->>>>>>>>>>>>>>>
                tx2                          <<<<<-----<<<<<<<<<<
        */
        let mut read = Record::new();
        let qname = "ClippedRead".as_bytes();
        let cigar = CigarString(vec![
            Cigar::SoftClip(2),
            Cigar::Match(5),
            Cigar::SoftClip(2),
        ]);
        let seq = "AAAAAAAAA".as_bytes();
        let qual = "IIIIIIIII".as_bytes();
        read.unset_flags();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(qname, Some(&cigar), seq, qual);

        let gx1 = Some(vec!["gx1"]);
        let gx2 = Some(vec!["gx2"]);
        let no_gene = Some(Vec::new());
        // Intron sense alignment
        read.set_pos(15);
        read.unset_reverse();
        check_annotation(&read, &tx_def, 0, 0, None, AnnotationRegion::Intronic);
        check_annotation(
            &read,
            &tx_int,
            1,
            0,
            gx1.clone(),
            AnnotationRegion::Intronic,
        );

        // Intron antisense
        read.set_pos(15);
        read.set_reverse();
        check_annotation(&read, &tx_def, 0, 0, None, AnnotationRegion::Intronic);
        check_annotation(
            &read,
            &tx_int,
            0,
            1,
            no_gene.clone(),
            AnnotationRegion::Intronic,
        );

        // Exon sense of tx2 vs intron antisense of tx1: choose exon
        read.set_pos(25);
        read.set_reverse();
        check_annotation(&read, &tx_def, 1, 0, None, AnnotationRegion::Exonic);
        check_annotation(&read, &tx_int, 1, 0, gx2.clone(), AnnotationRegion::Exonic);

        // Exon anti-sense of tx2 vs intron sense of tx1: choose exon
        read.unset_reverse();
        check_annotation(&read, &tx_def, 0, 1, None, AnnotationRegion::Exonic);
        check_annotation(
            &read,
            &tx_int,
            0,
            1,
            no_gene.clone(),
            AnnotationRegion::Exonic,
        );

        // Intron sense (tx1) vs intron anti-sense (tx2)
        read.set_pos(30);
        read.unset_reverse();
        check_annotation(&read, &tx_def, 0, 0, None, AnnotationRegion::Intronic);
        check_annotation(
            &read,
            &tx_int,
            1,
            1,
            gx1.clone(),
            AnnotationRegion::Intronic,
        );

        // 2 exons sense
        read.set_pos(0);
        read.unset_reverse();
        check_annotation(&read, &tx_def, 2, 0, None, AnnotationRegion::Exonic);
        check_annotation(
            &read,
            &tx_int,
            2,
            0,
            Some(vec!["gx0", "gx1"]),
            AnnotationRegion::Exonic,
        );

        read.set_pos(35);
        read.unset_reverse();
        check_annotation(&read, &tx_def, 1, 1, None, AnnotationRegion::Exonic);
        check_annotation(&read, &tx_int, 1, 1, gx1.clone(), AnnotationRegion::Exonic);
    }

    #[test]
    fn test_transcriptome_basic() {
        let mut txome = make_test_annotator(default_params());

        let mut read = Record::new();
        let qname = "ClippedRead".as_bytes();
        let cigar = CigarString(vec![
            Cigar::SoftClip(2),
            Cigar::Match(8),
            Cigar::SoftClip(2),
        ]);
        let seq = "AAAAAAAAAAAA".as_bytes();
        let qual = "IIIIIIIIIIII".as_bytes();
        read.unset_flags();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(qname, Some(&cigar), seq, qual);

        // test different positions
        read.set_pos(0); // non-clipped portion entirely exonic (tx0/tx1)
        check_annotation(&read, &txome, 2, 0, None, AnnotationRegion::Exonic);
        read.set_pos(5); // non-clipped portion mostly exonic, but is not consistent with exon boundaries
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Exonic);
        read.set_pos(8); // mostly intronic (overlap exon on left)
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intronic);
        read.set_pos(10); // entirely intronic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intronic);
        read.set_pos(20); // mostly intronic (overlap exon on right)
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intronic);
        read.set_pos(23); // overlaps intron / intergenic on either side, but mostly exonic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Exonic);
        read.set_pos(25); // mostly exonic, antisense to tx2
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Exonic);
        read.set_pos(28); // overlaps exons on either side, but mostly intronic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intronic);
        read.set_pos(35); // entirely exonic (tx1/tx2, tx2 is antisense)
        check_annotation(&read, &txome, 1, 1, None, AnnotationRegion::Exonic);
        read.set_pos(45); // mostly exonic (hanging off end)
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Exonic);
        read.set_pos(48); // mostly intergenic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intergenic);
        read.set_pos(50); // totally intergenic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intergenic);

        // test different overlap threshold
        read.set_pos(5);
        txome.params.region_min_overlap = 0.8;
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intronic);

        // test unmapped
    }

    #[test]
    fn test_transcriptome_strand() {
        let mut txome = make_test_annotator(default_params());

        let mut read = Record::new();
        let qname = "ClippedRead".as_bytes();
        let cigar = CigarString(vec![
            Cigar::SoftClip(2),
            Cigar::Match(8),
            Cigar::SoftClip(2),
        ]);
        let seq = "AAAAAAAAAAAA".as_bytes();
        let qual = "IIIIIIIIIIII".as_bytes();
        read.unset_flags();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(qname, Some(&cigar), seq, qual);
        read.set_pos(0);
        read.set_paired();

        // positive chemistry, R2 forward strand
        read.set_last_in_template();
        check_annotation(&read, &txome, 0, 2, None, AnnotationRegion::Exonic);

        // positive chemistry, R2 negative strand
        read.set_reverse();
        check_annotation(&read, &txome, 2, 0, None, AnnotationRegion::Exonic);

        // negative chemistry, R2 negative strand
        txome.params.chemistry_strandedness = ReqStrand::Reverse;
        check_annotation(&read, &txome, 0, 2, None, AnnotationRegion::Exonic);
    }

    #[test]
    fn test_transcriptome_splice() {
        let txome = make_test_annotator(default_params());

        let mut read = Record::new();
        let qname = "SplicedRead".as_bytes();
        let cigar = CigarString(vec![Cigar::Match(5), Cigar::RefSkip(25), Cigar::Match(7)]);
        let seq = "AAAAAAAAAAAA".as_bytes();
        let qual = "IIIIIIIIIIII".as_bytes();
        read.unset_flags();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(qname, Some(&cigar), seq, qual);

        read.set_pos(0); // first segment is exonic, but second segment is intronic -> intronic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intronic);

        read.set_pos(5); // exonic w/ correct splice junctions for tx1
        check_annotation(&read, &txome, 1, 0, None, AnnotationRegion::Exonic);

        read.set_pos(6); // misaligned splice junction, but each segment is exonic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Exonic);

        read.set_pos(25); // one segment intronic, other intergenic -> intergenic
        check_annotation(&read, &txome, 0, 0, None, AnnotationRegion::Intergenic);
    }
}
