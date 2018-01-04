//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use rust_htslib::bam::record::{Cigar, Record, Aux};
use std::path::Path;
use std::str;
use std::cmp;
use std::collections::{HashSet, BTreeMap};
use utils;

use reference::{TranscriptIndex, Gene};

// These only exist in BAM land; just use bytes.
const TRANSCRIPT_TAG: &'static [u8]  = b"TX";
const GENE_ID_TAG: &'static [u8]     = b"GX";
const GENE_NAME_TAG: &'static [u8]   = b"GN";
const REGION_TAG: &'static [u8]      = b"RE";
const MULTIMAPPER_TAG: &'static [u8] = b"MM";
const ANTISENSE_TAG: &'static [u8]   = b"AN"; // equivalent to TX, but for antisense alignments
const EXTRA_FLAGS_TAG: &'static [u8] = b"xf";
// These are set to the original single-read annotations
// if a read-pair had gene disagreement.
const UNPAIRED_GENE_ID_TAG: &'static [u8] = b"gX";
const UNPAIRED_GENE_NAME_TAG: &'static [u8] = b"gN";

bitflags! {
    #[derive(Default)]
    pub struct ExtraFlags: u32 {
        // Confidently mapped to transcriptome
        const CONF_MAPPED = 1u32;
        const LOW_SUPPORT_UMI = 2u32;
        // Mates mapped to incompatible sets of genes
        const GENE_DISCORDANT = 4u32;
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct AnnotationData {
    pub transcripts:    BTreeMap<String, TranscriptAlignment>,
    pub antisense:      BTreeMap<String, TranscriptAlignment>,
    pub genes:          Vec<Gene>,
    pub region:         AnnotationRegion,
    pub rescued:        bool,
}

#[derive(Eq, PartialEq, Debug)]
pub struct PairAnnotationData {
    pub genes: Vec<Gene>,
}

fn format_gx_gn_tags(genes: &Vec<Gene>) -> Option<(String, String)> {
    if genes.len() > 0 {
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

impl PairAnnotationData {
    /// Annotate a pair of alignments
    /// Take the intersection of the non-empty gene sets of the mates
    pub fn from_pair(anno1: &AnnotationData, anno2: &AnnotationData)
                     -> PairAnnotationData {
        let genes = match (anno1.genes.len() > 0, anno2.genes.len() > 0) {
            (true, false) => anno1.genes.clone(),
            (false, true) => anno2.genes.clone(),
            _ => utils::intersect_vecs(&anno1.genes, &anno2.genes),
        };
        PairAnnotationData { genes: genes }
    }

    pub fn make_gx_gn_tags(&self) -> Option<(String, String)> {
        format_gx_gn_tags(&self.genes)
    }
}

impl AnnotationData {
    fn new() -> AnnotationData {
        AnnotationData {
            transcripts:    BTreeMap::new(),
            antisense:      BTreeMap::new(),
            genes:          Vec::new(),
            region:         AnnotationRegion::Unmapped,
            rescued:        false,
        }
    }

    pub fn is_antisense(&self) -> bool {
        return (self.transcripts.len() == 0) && (self.antisense.len() > 0)
    }

    fn make_tx_tag(&self) -> Option<String> {
        if self.transcripts.len() > 0 {
            Some(self.transcripts.values().map(|aln| aln.to_tag()).collect::<Vec<String>>().join(";"))
        } else {
            None
        }
    }

    fn make_an_tag(&self) -> Option<String> {
        if self.antisense.len() > 0 {
            Some(self.antisense.values().map(|aln| aln.to_tag()).collect::<Vec<String>>().join(";"))
        } else {
            None
        }
    }

    fn make_re_tag(&self) -> Option<String> {
        self.region.to_tag()
    }

    fn make_gx_gn_tags(&self) -> Option<(String, String)> {
        format_gx_gn_tags(&self.genes)
    }

    fn make_mm_tag(&self) -> Option<i32> {
        if self.rescued {
            Some(1)
        } else {
            None
        }
    }

    /// Add tags to a BAM record.
    /// Set is_conf_mapped to true if the qname is confidently mapped to
    /// the transcriptome.
    pub fn attach_tags(&mut self, record: &mut Record, is_conf_mapped: bool,
                       is_gene_discordant: bool, pair_anno: Option<&PairAnnotationData>) {
        if let Some(tag) = self.make_tx_tag() {
            record.push_aux(TRANSCRIPT_TAG, &Aux::String(tag.as_bytes()));
        }

        if let Some(pair) = pair_anno {
            // Write the pair's annotation
            if let Some((tag_gx, tag_gn)) = pair.make_gx_gn_tags() {
                record.push_aux(GENE_ID_TAG, &Aux::String(tag_gx.as_bytes()));
                record.push_aux(GENE_NAME_TAG, &Aux::String(tag_gn.as_bytes()));
            }

            if self.genes != pair.genes {
                // This record disagrees with the pair, so store its single-end
                // gene annotations separately.
                if let Some((tag_gx, tag_gn)) = self.make_gx_gn_tags() {
                    record.push_aux(UNPAIRED_GENE_ID_TAG, &Aux::String(tag_gx.as_bytes()));
                    record.push_aux(UNPAIRED_GENE_NAME_TAG, &Aux::String(tag_gn.as_bytes()));
                }
            }

        } else {
            // Unpaired case
            if let Some((tag_gx, tag_gn)) = self.make_gx_gn_tags() {
                record.push_aux(GENE_ID_TAG, &Aux::String(tag_gx.as_bytes()));
                record.push_aux(GENE_NAME_TAG, &Aux::String(tag_gn.as_bytes()));
            }
        }

        if let Some(tag) = self.make_re_tag() {
            record.push_aux(REGION_TAG, &Aux::Char(tag.as_bytes()[0]));
        }
        if let Some(tag) = self.make_mm_tag() {
            record.push_aux(MULTIMAPPER_TAG, &Aux::Integer(tag));
        }
        if let Some(tag) = self.make_an_tag() {
            record.push_aux(ANTISENSE_TAG, &Aux::String(tag.as_bytes()));
        }

        // Note: only attach these flags to primary alignment
        if !record.is_secondary() {
            let mut flags: ExtraFlags = Default::default();
            if is_conf_mapped {
                flags |= CONF_MAPPED;
            }
            if is_gene_discordant {
                flags |= GENE_DISCORDANT;
            }
            record.push_aux(EXTRA_FLAGS_TAG, &Aux::Integer(flags.bits() as i32));
        }
    }
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct TranscriptAlignment {
    pub id: String,
    pub strand: Strand,
    pub pos:    i64,
    pub cigar:  String,
    pub alen:   i64,
}

impl TranscriptAlignment {
    fn to_tag(&self) -> String {
        format!("{},{}{},{}", self.id, self.strand.to_string(), self.pos, self.cigar)
    }
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Copy, Hash)]
pub enum AnnotationRegion {
    Exonic,
    Intronic,
    Intergenic,
    Unmapped,
}

impl AnnotationRegion {
    fn to_tag(&self) -> Option<String> {
        match *self {
            AnnotationRegion::Exonic          => Some("E".to_string()),
            AnnotationRegion::Intronic        => Some("N".to_string()),
            AnnotationRegion::Intergenic      => Some("I".to_string()),
            AnnotationRegion::Unmapped        => None,
        }
    }
}

#[derive(Deserialize, Debug)]
struct StarTranscript {
    id:         String,
    start:      i64,
    end:        i64, // inclusive
    max_end:    i64,
    strand:     u8,
    num_exons:  usize,
    break_idx:  usize, // TODO rename
}

#[derive(Deserialize, Debug)]
struct StarExon {
    start:      i64,
    end:        i64, // inclusive
    cum_len:    i64, // cumulative length of previous exons, not including this one
}

impl StarExon {
    fn len(&self) -> i64 {
        self.end - self.start + 1
    }
}

pub struct AnnotationParams {
    pub chemistry_strandedness:     Strandedness,
    pub chemistry_fiveprime:        bool,
    pub intergenic_trim_bases:      i64,
    pub intronic_trim_bases:        i64,
    pub junction_trim_bases:        i64,
    pub region_min_overlap:         f64,
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    fn to_string(&self) -> String {
        match self {
            &Strand::Forward => "+".to_string(),
            &Strand::Reverse => "-".to_string(),
        }
    }
}

pub enum Strandedness {
    Forward,
    Reverse,
    Mixed,
}

impl Strandedness {
    pub fn from_string(strandedness: String) -> Strandedness {
        match strandedness.as_ref() {
            "+"     => Strandedness::Forward,
            "-"     => Strandedness::Reverse,
            "mixed" => Strandedness::Mixed,
            _       => panic!("Invalid strand arg: {}", strandedness),
        }
    }
}

#[derive(Debug)]
struct SpliceSegment {
    start:      i64,
    end:        i64,
    cigar:      Vec<FlatCigar>, // TODO could be a slice of the full cigar to avoid copying
}

#[derive(Eq, PartialEq, Debug, Clone, Copy)]
enum CigarType {
    Match,
    Ins,
    Del,
    RefSkip,
    SoftClip,
    HardClip,
    Pad,
    Equal,
    Diff,
    Back,
}

// TODO: rust_htslib now has an improved CIGAR API, so we should try using that if it's simpler
#[derive(Eq, PartialEq, Debug, Clone, Copy)]
struct FlatCigar {
    op:     CigarType,
    len:    u32,
}

impl FlatCigar {
    fn to_string(&self) -> String {
        let symbol = match self.op {
            CigarType::Match    => 'M',
            CigarType::Ins      => 'I',
            CigarType::Del      => 'D',
            CigarType::RefSkip  => 'N',
            CigarType::SoftClip => 'S',
            CigarType::HardClip => 'H',
            CigarType::Pad      => 'P',
            CigarType::Equal    => '=',
            CigarType::Diff     => 'X',
            CigarType::Back     => 'B',
        };
        return format!("{}{}", self.len, symbol)
    }
}

fn flatten_cigar(cig: &Cigar) -> FlatCigar {
    return match cig {
        &Cigar::Match(len)    => FlatCigar { op: CigarType::Match, len: len },
        &Cigar::Ins(len)      => FlatCigar { op: CigarType::Ins, len: len },
        &Cigar::Del(len)      => FlatCigar { op: CigarType::Del, len: len },
        &Cigar::RefSkip(len)  => FlatCigar { op: CigarType::RefSkip, len: len },
        &Cigar::SoftClip(len) => FlatCigar { op: CigarType::SoftClip, len: len },
        &Cigar::HardClip(len) => FlatCigar { op: CigarType::HardClip, len: len },
        &Cigar::Pad(len)      => FlatCigar { op: CigarType::Pad, len: len },
        &Cigar::Equal(len)    => FlatCigar { op: CigarType::Equal, len: len },
        &Cigar::Diff(len)     => FlatCigar { op: CigarType::Diff, len: len },
        &Cigar::Back(len)     => FlatCigar { op: CigarType::Back, len: len },
    }
}

fn align_to_transcriptome(read: &Record, chrom_starts: &[i64], transcript_info: &[StarTranscript],
    exon_info: &[StarExon], transcript_index: &TranscriptIndex, params: &AnnotationParams) -> AnnotationData {
    let mut annotation_data = AnnotationData::new();

    if read.is_unmapped() {
        return annotation_data
    } else {
        annotation_data.region = AnnotationRegion::Intergenic;
    }

    let ref_offset = chrom_starts[read.tid() as usize];
    let clipped_read_start = ref_offset + read.pos() as i64;
    let clipped_read_end = clipped_read_start + utils::alen(read) - 1;

    // find maximum index of transcripts that could possibly overlap the read
    let mut tx_idx = utils::bisect(transcript_info, clipped_read_end, &|tx| tx.start, utils::BisectDirection::Right);
    //let mut tx_idx = transcript_info.binary_search_by_key(&clipped_read_end, |tx| tx.start).unwrap();

    if tx_idx == 0 {
        // read ends before the first transcript start
        return annotation_data;
    } else {
        tx_idx -= 1; // first transcript on the left
    }

    // cycle back through overlapping transcripts
    let mut seen_genes = HashSet::new();
    let mut curr_transcript = &transcript_info[tx_idx];
    let mut any_exonic = false;
    let mut any_intronic = false;
    while clipped_read_start <= cmp::max(curr_transcript.end, curr_transcript.max_end) {
        if clipped_read_start <= curr_transcript.end { // transcript overlaps the read
            let (curr_alignment, curr_region, curr_strand) = align_to_transcript(read, ref_offset, curr_transcript, exon_info, params);
            match curr_region {
                AnnotationRegion::Exonic => any_exonic = true,
                AnnotationRegion::Intronic => any_intronic = true,
                AnnotationRegion::Intergenic => (),
                AnnotationRegion::Unmapped => (), // TODO should never happen - should we panic?
            };
            match curr_alignment {
                Some(aln) => match curr_strand {
                    Strand::Forward => {
                        seen_genes.insert(transcript_index.get_gene_from_transcript(&aln.id).to_owned());
                        annotation_data.transcripts.insert(aln.id.clone(), aln);
                    },
                    Strand::Reverse => {
                        annotation_data.antisense.insert(aln.id.clone(), aln);
                    }
                },
                None => (),
            }
        }

        if tx_idx == 0 {
            break
        } else {
            tx_idx -= 1;
            curr_transcript = &transcript_info[tx_idx];
        }
    }

    annotation_data.region = if any_exonic { AnnotationRegion::Exonic } else if any_intronic { AnnotationRegion::Intronic } else { AnnotationRegion::Intergenic };

    annotation_data.genes = seen_genes.into_iter().collect::<Vec<Gene>>();
    // Sorting this makes life easier later.
    annotation_data.genes.sort_unstable();

    return annotation_data
}

fn align_to_transcript(read: &Record, ref_offset: i64, transcript_info: &StarTranscript, exon_info: &[StarExon], params: &AnnotationParams) -> (Option<TranscriptAlignment>, AnnotationRegion, Strand) {
    // figure out coordinates
    let tx_start = transcript_info.start;
    let tx_end = transcript_info.end;
    let tx_first_exon = transcript_info.break_idx;
    let tx_num_exons = transcript_info.num_exons;
    let tx_last_exon = tx_first_exon + tx_num_exons - 1; // inclusive
    let star_genomic_start = ref_offset + read.pos() as i64;
    let star_genomic_end = star_genomic_start + utils::alen(read) - 1; // inclusive
    let clipped_read_start = star_genomic_start - tx_start; // convert to transcript coordinates
    let clipped_read_end = star_genomic_end - tx_start;

    // compute region
    let (left_clip, right_clip, splice_segments) = get_cigar_segments(read, clipped_read_start);
    let is_exonic = is_read_exonic(&splice_segments, exon_info, tx_first_exon, tx_last_exon, params.region_min_overlap);
    let is_intronic = !is_exonic && get_overlap(star_genomic_start, star_genomic_end, tx_start, tx_end) >= 1.0;
    let region = if is_exonic { AnnotationRegion::Exonic } else if is_intronic { AnnotationRegion::Intronic } else { AnnotationRegion::Intergenic };

    // discard non-exonic reads
    if !is_exonic {
        return (None, region, Strand::Forward)
    };

    // compute strand
    let tx_reverse_strand = transcript_info.strand == 2;
    let mut read_reverse_strand = read.is_reverse();
    if read.is_paired() && read.is_last_in_template() { read_reverse_strand = !read_reverse_strand };
    let is_antisense = match params.chemistry_strandedness {
        Strandedness::Forward => tx_reverse_strand != read_reverse_strand,
        Strandedness::Reverse => tx_reverse_strand == read_reverse_strand,
        Strandedness::Mixed   => false,
    };

    let tx_strand = if is_antisense { Strand::Reverse } else { Strand::Forward };

    // find the beginning / ending exons
    let (ex_start, ex_end) = match find_exons(exon_info, clipped_read_start, clipped_read_end, tx_first_exon, tx_last_exon, params.intergenic_trim_bases, params.intronic_trim_bases) {
        Some((s,e)) => (s,e),
        None => return (None, region, tx_strand),
    };

    // compute offsets
    let ex_offset = cmp::max(clipped_read_start - exon_info[ex_start].start, 0);
    let mut tx_offset = exon_info[ex_start].cum_len + ex_offset;
    let tx_length = exon_info[tx_last_exon].cum_len + exon_info[tx_last_exon].len();

    // align the read to the exons
    let tx_alignment = align_junctions(&left_clip, &right_clip, &splice_segments, &exon_info[ex_start..ex_end+1], params.junction_trim_bases);

    // discard misaligned reads
    let (mut tx_cigar, tx_aligned_bases) = match tx_alignment {
        Some((cig, bases)) => (cig, bases),
        None => return (None, region, tx_strand),
    };

    // flip reverse strand
    if tx_reverse_strand {
        tx_offset = tx_length - (tx_offset + tx_aligned_bases);
        tx_cigar.reverse();
    };

    let tx_cigar_string = tx_cigar.iter().map(|c| c.to_string()).collect::<Vec<String>>().join("");

    let alignment = TranscriptAlignment {
        id: transcript_info.id.clone(),
        strand: tx_strand,
        pos: tx_offset,
        cigar: tx_cigar_string,
        alen: tx_aligned_bases,
    };

    return (Some(alignment), region, tx_strand)
}

fn is_read_exonic(splice_segments: &[SpliceSegment], exon_info: &[StarExon], first_exon: usize, last_exon: usize, min_overlap_frac: f64) -> bool {
    for segment in splice_segments {
        // find first exon that ends to the right of the segment start
        let ex_idx = utils::bisect(&exon_info[first_exon..last_exon+1], segment.start, &|ex| ex.end, utils::BisectDirection::Right) + first_exon;
        //let ex_idx = exon_info.binary_search_by_key(&segment.start, |ex| ex.end - 1).unwrap() as u64;
        if ex_idx > last_exon {
            // read is out of bounds
            return false
        } else {
            let exon = &exon_info[ex_idx];
            if get_overlap(segment.start, segment.end, exon.start, exon.end) < min_overlap_frac {
                return false
            }
        }
    };
    return true
}

fn find_exons(exon_info: &[StarExon], read_start: i64, read_end: i64, first_exon: usize, last_exon: usize, intergenic_trim_bases: i64, intronic_trim_bases: i64) -> Option<(usize, usize)> {
    let ex_start = utils::bisect(&exon_info[first_exon..last_exon+1], read_start, &|ex| ex.end, utils::BisectDirection::Right) + first_exon; // find first exon that ends to the right of the read start
    let ex_end = utils::bisect(&exon_info[first_exon..last_exon+1], read_end, &|ex| ex.start, utils::BisectDirection::Left) - 1 + first_exon; // find first exon that starts to the left of the read end
    //let ex_start = exon_info.binary_search_by_key(&read_start, |ex| ex.end - 1).unwrap() as u64; // find first exon that ends to the right of the read start
    //let ex_end = exon_info.binary_search_by_key(&read_end, |ex| ex.start - 1).unwrap() as u64 - 1; // find first exon that starts to the left of the read end

    if (ex_start > last_exon) | (ex_end < first_exon) {
        // read is out of bounds
        return None
    }

    let starting_exon = &exon_info[ex_start];
    let ending_exon = &exon_info[ex_end];

    if read_start < starting_exon.start {
        // read overhangs exon on the left
        let overhang = starting_exon.start - read_start;
        let trim_bases = if ex_start < first_exon { intergenic_trim_bases } else { intronic_trim_bases };
        if overhang > trim_bases {
            // too much overhang
            return None
        };
    }

    if read_end > ending_exon.end {
        // read overhangs exon on the right
        let overhang = read_end - ending_exon.end;
        let trim_bases = if ex_end > last_exon { intergenic_trim_bases } else { intronic_trim_bases };
        if overhang > trim_bases {
            // too much overhang
            return None
        };
    }

    return Some((ex_start, ex_end))
}

fn get_overlap(read_start: i64, read_end: i64, ref_start: i64, ref_end: i64) -> f64 {
    let mut overlap_bases = cmp::min(ref_end, read_end) - cmp::max(ref_start, read_start);
    if overlap_bases < 0 {
        overlap_bases = 0;
    }
    return (overlap_bases as f64) / ((read_end - read_start) as f64);
}

fn get_cigar_segments(read: &Record, alignment_start: i64) -> (Vec<FlatCigar>, Vec<FlatCigar>, Vec<SpliceSegment>) {
    let mut left_clip: Vec<FlatCigar> = Vec::new();
    let mut right_clip: Vec<FlatCigar> = Vec::new();
    let mut splice_segments: Vec<SpliceSegment> = Vec::new();
    let mut seen_nonclips = false; // whether we've seen non-clip bases yet
    let mut curr_segment = SpliceSegment {
        start: alignment_start,
        end: alignment_start,
        cigar: Vec::new(),
    };

    for c in read.cigar().iter().map(|c| flatten_cigar(c)) {
        match c.op {
            CigarType::HardClip | CigarType::SoftClip => {
                if seen_nonclips { right_clip.push(c) }
                else { left_clip.push(c) }
            },
            CigarType::RefSkip => {
                seen_nonclips = true;
                let next_start = curr_segment.end + c.len as i64;
                splice_segments.push(curr_segment);
                curr_segment = SpliceSegment {
                    start: next_start,
                    end: next_start,
                    cigar: Vec::new(),
                };
            },
            CigarType::Ins => {
                seen_nonclips = true;
                curr_segment.cigar.push(c);
            },
            CigarType::Match | CigarType::Del | CigarType::Equal | CigarType::Diff => {
                seen_nonclips = true;
                curr_segment.end += c.len as i64;
                curr_segment.cigar.push(c);
            },
            _ => ()
        }
    }
    splice_segments.push(curr_segment);
    return (left_clip, right_clip, splice_segments)
}

fn align_junctions(left_clip: &[FlatCigar], right_clip: &[FlatCigar], splice_segments: &[SpliceSegment], exon_info: &[StarExon], tolerance: i64) -> Option<(Vec<FlatCigar>, i64)> {
    if splice_segments.len() != exon_info.len() {
        return None
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
        if i == 0 { // first segment
            if start_diff > 0 { // overhang -> softclip
                tmp_cigar = mask_read_bases(&mut tmp_cigar, CigarType::SoftClip, start_diff as u32, false);
            } else if start_diff < 0 { // underhang -> decrement aligned bases
                aligned_bases -= start_diff.abs();
            }
        } else {
            if start_diff.abs() > tolerance {
                return None // can't align properly
            } else if start_diff > 0 { // overhang -> mark as insertion
                tmp_cigar = mask_read_bases(&mut tmp_cigar, CigarType::Ins, start_diff as u32, false);
            } else if start_diff < 0 { // underhang -> mark as deletion
                tmp_cigar = mark_deleted_ref_bases(&mut tmp_cigar, start_diff.abs() as u32, false);
            }
        }

        // align the end
        let end_diff = curr_segment.end - curr_exon.end - 1;
        if i == splice_segments.len() - 1 { // last segment
            if end_diff > 0 { // overhang -> softclip
                tmp_cigar = mask_read_bases(&mut tmp_cigar, CigarType::SoftClip, end_diff as u32, true);
            } else if end_diff < 0 { // underhang -> decrement aligned bases
                aligned_bases -= end_diff.abs();
            }
        } else {
            if end_diff.abs() > tolerance {
                return None // can't align properly
            } else if end_diff > 0 { // overhang -> mark as insertion
                tmp_cigar = mask_read_bases(&mut tmp_cigar, CigarType::Ins, end_diff as u32, true);
            } else if end_diff < 0 { // underhang -> mark as deletion
                tmp_cigar = mark_deleted_ref_bases(&mut tmp_cigar, end_diff.abs() as u32, true);
            }
        }

        // extend
        full_cigar = extend_cigar(&full_cigar, &tmp_cigar);
    }
    full_cigar = extend_cigar(&full_cigar, right_clip);

    return Some((full_cigar, aligned_bases))
}

fn extend_cigar(old_cigar: &[FlatCigar], new_cigar: &[FlatCigar]) -> Vec<FlatCigar> {
    // extend list of cigar ops, checking the ends to see if they should be merged
    let old_len = old_cigar.len();
    let new_len = new_cigar.len();
    let mut merged_cigar = Vec::new();
    if (old_len > 0) && (new_len > 0) {
        let old_tail = old_cigar[old_len - 1];
        let new_head = new_cigar[0];
        if old_tail.op == new_head.op {
            let merged_ends = FlatCigar { op: new_head.op, len: old_tail.len + new_head.len };
            merged_cigar.extend_from_slice(&old_cigar[..old_len-1]);
            merged_cigar.push(merged_ends);
            merged_cigar.extend_from_slice(&new_cigar[1..]);
        } else {
            merged_cigar.extend_from_slice(&old_cigar);
            merged_cigar.extend_from_slice(&new_cigar);
        }
    } else {
        merged_cigar.extend_from_slice(&old_cigar);
        merged_cigar.extend_from_slice(&new_cigar);
    }
    return merged_cigar
}

fn mark_deleted_ref_bases(cigar: &mut [FlatCigar], del_len: u32, reverse: bool) -> Vec<FlatCigar> {
    let mut new_cigar = Vec::new();
    let del = FlatCigar { op: CigarType::Del, len: del_len};
    if reverse {
        new_cigar.push(del);
        new_cigar.extend_from_slice(&cigar);
    } else {
        new_cigar.extend_from_slice(&cigar);
        new_cigar.push(del);
    };
    return new_cigar
}

fn mask_read_bases(cigar: &mut [FlatCigar], mask_op: CigarType, mask_len: u32, reverse: bool) -> Vec<FlatCigar> {
    // NOTE: this assumes that refskips have been removed
    let mut new_cigar = Vec::new();
    let mask = FlatCigar { op: mask_op, len: mask_len };
    let mut consumed_bases = 0;
    new_cigar.push(mask);
    if reverse { cigar.reverse(); }
    for c in cigar {
        if consumed_bases < mask_len {
            // this op should be masked
            let read_bases = match c.op {
                CigarType::Del => 0, // deletions don't consume read bases
                _ => c.len,
            };
            if consumed_bases + read_bases >= mask_len {
                let truncated = FlatCigar { op: c.op, len: read_bases + consumed_bases - mask_len };
                new_cigar.push(truncated);
            };
            consumed_bases += read_bases;
        } else {
            // just copy the op
            new_cigar.push(*c);
        };
    };
    if reverse { new_cigar.reverse(); }
    return new_cigar
}

pub struct TranscriptAnnotator {
    params:             AnnotationParams,
    transcript_info:    Vec<StarTranscript>,
    exon_info:          Vec<StarExon>,
    chrom_starts:       Vec<i64>,
    transcript_index:   TranscriptIndex,
}

impl TranscriptAnnotator {
    pub fn new(reference_path: &str, transcript_index_tabfile: &str, params: AnnotationParams) -> TranscriptAnnotator {
        // TODO less hacky paths
        let transcript_info: Vec<StarTranscript> = utils::load_tabular(Path::new(reference_path).join("star/transcriptInfo.tab").as_path().to_str().unwrap(), true);
        let exon_info: Vec<StarExon> = utils::load_tabular(Path::new(reference_path).join("star/exonInfo.tab").as_path().to_str().unwrap(), true);
        let chrom_starts: Vec<i64> = utils::load_txt(Path::new(reference_path).join("star/chrStart.txt").as_path().to_str().unwrap());
        let transcript_index = TranscriptIndex::new(transcript_index_tabfile);
        return TranscriptAnnotator {
            params:             params,
            transcript_info:    transcript_info,
            exon_info:          exon_info,
            chrom_starts:       chrom_starts,
            transcript_index:   transcript_index,
        }
    }

    pub fn annotate_genomic_alignment(&self, read: &Record) -> AnnotationData {
        align_to_transcriptome(read, &self.chrom_starts, &self.transcript_info, &self.exon_info, &self.transcript_index, &self.params)
    }

    pub fn get_transcript_index(&self) -> &TranscriptIndex {
        return &self.transcript_index
    }

    pub fn get_params(&self) -> &AnnotationParams {
        return &self.params
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use rust_htslib::bam::record::CigarString;

    struct TranscriptomeTest {
        chrom_starts:       Vec<i64>,
        transcript_info:    Vec<StarTranscript>,
        exon_info:          Vec<StarExon>,
        transcript_index:   TranscriptIndex,
    }

    pub fn set_forward(record: &mut Record) {
        // NOTE: rust_htslib doesn't have a method for this either
        if record.is_reverse() {
            record.inner_mut().core.flag -= 16u16;
        }
    }

    impl TranscriptomeTest {
        fn new() -> TranscriptomeTest {
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
                StarTranscript { id: "tx0".into(), start: 0, end: 9, max_end: 9, strand: 1, num_exons: 1, break_idx: 0 },
                StarTranscript { id: "tx1".into(), start: 0, end: 49, max_end: 49, strand: 1, num_exons: 2, break_idx: 1 },
                StarTranscript { id: "tx2".into(), start: 25, end: 44, max_end: 49, strand: 2, num_exons: 2, break_idx: 3 },
            ];

            let exon_info = vec![
                // tx0
                StarExon { start: 0, end: 9, cum_len: 0 },
                // tx1
                StarExon { start: 0, end: 9, cum_len: 0 },
                StarExon { start: 35, end: 49, cum_len: 10 },
                // tx2
                StarExon { start: 0, end: 4, cum_len: 0 },
                StarExon { start: 10, end: 19, cum_len: 5 },
            ];

            let chrom_starts = vec![0];

            let mut transcript_genes = HashMap::new();
            transcript_genes.insert("tx0".into(), Gene { id: "gx0".into(), name: "gene0".into() });
            transcript_genes.insert("tx1".into(), Gene { id: "gx1".into(), name: "gene1".into() });
            transcript_genes.insert("tx2".into(), Gene { id: "gx2".into(), name: "gene2".into() });

            let transcript_lengths = HashMap::new(); // doesn't need to be populated

            let transcript_index = TranscriptIndex {
                transcript_genes: transcript_genes,
                transcript_lengths: transcript_lengths,
            };

            return TranscriptomeTest {
                chrom_starts:       chrom_starts,
                transcript_info:    transcript_info,
                exon_info:          exon_info,
                transcript_index:   transcript_index,
            }
        }

        fn check_annotation(&self, rec: &Record, params: &AnnotationParams, expected_transcripts: usize,
                            expected_antisense: usize, expected_region: AnnotationRegion) {
            let result = align_to_transcriptome(rec, &self.chrom_starts, &self.transcript_info,
                                                &self.exon_info, &self.transcript_index, params);
            println!("{:?}", result);
            assert_eq!(result.transcripts.len(), expected_transcripts);
            assert_eq!(result.antisense.len(), expected_antisense);
            assert_eq!(result.region, expected_region);
        }
    }

    fn default_params() -> AnnotationParams {
        AnnotationParams {
            chemistry_strandedness:     Strandedness::Forward,
            chemistry_fiveprime:        false,
            intergenic_trim_bases:      0,
            intronic_trim_bases:        0,
            junction_trim_bases:        0,
            region_min_overlap:         0.5,
        }
    }

    #[test]
    fn test_transcriptome_basic() {
        let txome = TranscriptomeTest::new();
        let mut params = default_params();

        let mut read = Record::new();
        let qname = "ClippedRead".as_bytes();
        let cigar = CigarString(vec![Cigar::SoftClip(2), Cigar::Match(8), Cigar::SoftClip(2)]);
        let seq = "AAAAAAAAAAAA".as_bytes();
        let qual = "IIIIIIIIIIII".as_bytes();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(&qname, &cigar, &seq, &qual);

        // test different positions
        read.set_pos(0); // non-clipped portion entirely exonic (tx0/tx1)
        txome.check_annotation(&read, &params, 2, 0, AnnotationRegion::Exonic);
        read.set_pos(5); // non-clipped portion mostly exonic, but is not consistent with exon boundaries
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Exonic);
        read.set_pos(8); // mostly intronic (overlap exon on left)
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intronic);
        read.set_pos(10); // entirely intronic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intronic);
        read.set_pos(20); // mostly intronic (overlap exon on right)
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intronic);
        read.set_pos(23); // overlaps intron / intergenic on either side, but mostly exonic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Exonic);
        read.set_pos(25); // mostly exonic, antisense to tx2
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Exonic);
        read.set_pos(28); // overlaps exons on either side, but mostly intronic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intronic);
        read.set_pos(35); // entirely exonic (tx1/tx2, tx2 is antisense)
        txome.check_annotation(&read, &params, 1, 1, AnnotationRegion::Exonic);
        read.set_pos(45); // mostly exonic (hanging off end)
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Exonic);
        read.set_pos(48); // mostly intergenic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intergenic);
        read.set_pos(50); // totally intergenic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intergenic);

        // test different overlap threshold
        read.set_pos(5);
        params.region_min_overlap = 0.8;
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intronic);

        // test unmapped
        read.set_unmapped();
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Unmapped);
    }

    #[test]
    fn test_transcriptome_strand() {
        let txome = TranscriptomeTest::new();
        let mut params = default_params();

        let mut read = Record::new();
        let qname = "ClippedRead".as_bytes();
        let cigar = CigarString(vec![Cigar::SoftClip(2), Cigar::Match(8), Cigar::SoftClip(2)]);
        let seq = "AAAAAAAAAAAA".as_bytes();
        let qual = "IIIIIIIIIIII".as_bytes();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(&qname, &cigar, &seq, &qual);
        read.set_pos(0);
        read.set_paired();

        // positive chemistry, R2 forward strand
        read.set_last_in_template();
        txome.check_annotation(&read, &params, 0, 2, AnnotationRegion::Exonic);

        // positive chemistry, R2 negative strand
        read.set_reverse();
        txome.check_annotation(&read, &params, 2, 0, AnnotationRegion::Exonic);

        // negative chemistry, R2 negative strand
        params.chemistry_strandedness = Strandedness::Reverse;
        txome.check_annotation(&read, &params, 0, 2, AnnotationRegion::Exonic);

        // mixed-strand chemistry
        params.chemistry_strandedness = Strandedness::Mixed;
        read.set_reverse();
        txome.check_annotation(&read, &params, 2, 0, AnnotationRegion::Exonic);
        set_forward(&mut read);
        txome.check_annotation(&read, &params, 2, 0, AnnotationRegion::Exonic);
    }

    #[test]
    fn test_transcriptome_splice() {
        let txome = TranscriptomeTest::new();
        let params = default_params();

        let mut read = Record::new();
        let qname = "SplicedRead".as_bytes();
        let cigar = CigarString(vec![Cigar::Match(5), Cigar::RefSkip(25), Cigar::Match(7)]);
        let seq = "AAAAAAAAAAAA".as_bytes();
        let qual = "IIIIIIIIIIII".as_bytes();
        read.set_mapq(255);
        read.set_tid(0);
        read.set(&qname, &cigar, &seq, &qual);

        read.set_pos(0); // first segment is exonic, but second segment is intronic -> intronic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intronic);

        read.set_pos(5); // exonic w/ correct splice junctions for tx1
        txome.check_annotation(&read, &params, 1, 0, AnnotationRegion::Exonic);

        read.set_pos(6); // misaligned splice junction, but each segment is exonic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Exonic);

        read.set_pos(25); // one segment intronic, other intergenic -> intergenic
        txome.check_annotation(&read, &params, 0, 0, AnnotationRegion::Intergenic);
    }
}
