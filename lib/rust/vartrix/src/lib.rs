//! Crate for VarTrix.
//!
//! Computes allele counts given alignment shards and a barcode list.
#![deny(missing_docs)]
use anyhow::{Context, Result, bail};
use barcode::Barcode;
use bio::alignment::Alignment;
use bio::alignment::pairwise::banded::Aligner;
use bio::io::fasta;
use counter::Counter;
use cr_types::{BamContigInfo, BarcodeIndex};
use itertools::Itertools;
use log::{debug, warn};
use metric::{TxHashMap, TxHashSet};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{IndexedReader as BamReader, Read as BamRead, Record as BamRecord};
use rust_htslib::bcf::{Read as BcfRead, Reader as BcfReader, Record as BcfRecord};
use serde::{Deserialize, Serialize};
use std::cmp::{max, min};
use std::fmt::{Debug, Display};
use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use std::path::Path;
use std::str;
use std::str::FromStr;

/// Alignment scores below this are ignored
const MIN_SCORE: i32 = 25;
/// Minimum fraction of a UMI's reads that must agree on the `AlleleType` for it to be
/// considered not `Unknown`
const INTRA_UMI_CONS_THRESHOLD: f64 = 0.75;
/// kmer match length
const K: usize = 6;
/// Window size for creating the band
const W: usize = 20;
/// Match score
const MATCH: i32 = 1;
/// Mismatch score
const MISMATCH: i32 = -5;
/// Gap open score
const GAP_OPEN: i32 = -5;
/// Gap extend score
const GAP_EXTEND: i32 = -1;
/// UMI BAM tag
const UMI_BAM_TAG: &str = "UB";
/// Cell barcode BAM tag
const CB_BAM_TAG: &str = "CB";

/// Chunk of VCF records.
pub type VcfChunk = Vec<Snp>;
/// DNA/RNA sequence represented as a vector of bytes.
pub type Seq = Vec<u8>;
/// Index of a cell barcode in CellBarcodeIndex.
pub type CellBarcodeIdx = usize;
/// Index of a VCF record in the VCF file.
pub type VcfRecordIdx = usize;

#[derive(thiserror::Error, Debug, PartialEq)]
enum VartrixError {
    #[error("record {:?} does not have a rid", rec)]
    VcfRecordRidNotFound { rec: String },
    #[error("record {:?} does not have a name for its rid {}", rec, rid)]
    VcfRecordRidNameNotFound { rec: String, rid: u32 },
    #[error("record {:?} contig name is not valid UTF-8", rec)]
    VcfRecordContigNameNotValidUtf8 { rec: String },
    #[error("FASTA fetch error: {}", err_msg)]
    FastaFetchError { err_msg: String },
    #[error("FASTA reading error: {}", err_msg)]
    FastaReadError { err_msg: String },
    #[error("BAM record {} has invalid barcode {}", rec, bc)]
    InvalidBamRecordBarcode { rec: String, bc: String },
}

/// SNP record.
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct Snp {
    /// Index of the record in the VCF file.
    pub idx: VcfRecordIdx,
    /// Contig target ID
    pub tid: u32,
    /// Locus start
    pub start: u64,
    /// Locus end
    pub end: u64,
    /// Alternative allele
    pub alt: Vec<u8>,
}

/// Represents an allele.
#[derive(Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub enum AlleleType {
    /// Indicates a reference allele.
    Ref,
    /// Indicates an alternative allele.
    Alt,
    /// Indicates an unknown allele.
    Unknown,
}

impl AlleleType {
    /// Convert the allele to an u8 integer.
    pub fn to_u8(&self) -> u8 {
        match self {
            AlleleType::Ref => 0,
            AlleleType::Alt => 1,
            AlleleType::Unknown => 2,
        }
    }
}

/// Vartrix settings.
pub struct VartrixSettings {
    /// Specifies whether the operation is primary.
    primary_alignments_only: bool,
    /// The mapping quality threshold.
    mapq: u8,
    /// Specifies whether duplicates should be considered.
    skip_duplicate_alignments: bool,
    /// Specifies whether to use unique molecular identifiers (UMIs).
    use_umi: bool,
    /// The BAM tag to use for cell barcodes.
    cb_bam_tag: String,
    /// The BAM tag to use for UMIs.
    umi_bam_tag: String,
    /// The set of valid characters.
    valid_chars: TxHashSet<u8>,
    /// The padding value.
    padding: u64,
}

/// Default settings for Vartrix. Based on souporcell_pipeline.py
impl Default for VartrixSettings {
    fn default() -> Self {
        VartrixSettings {
            primary_alignments_only: false,
            mapq: 30,
            skip_duplicate_alignments: false,
            use_umi: true,
            cb_bam_tag: CB_BAM_TAG.to_string(),
            umi_bam_tag: UMI_BAM_TAG.to_string(),
            valid_chars: b"ATGCatgc".iter().copied().collect(),
            padding: 100,
        }
    }
}

/// Represents the metrics for a set of reads.
#[derive(Default, Serialize, Deserialize, Debug)]
pub struct VcfMetrics {
    /// Alignment use status statistics.
    pub aln_use_counter: Counter<AlignmentUseStatus>,
    /// The record in invalid data.
    pub has_invalid_chars: bool,
    /// The record has multiple alleles.
    pub is_multiallelic: bool,
}

/// Alignment use status
#[derive(Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub enum AlignmentUseStatus {
    /// Alignment used to compute allele counts.
    Useful,
    /// Alignment skipped due to low mapping quality.
    LowMappingQuality,
    /// Alignment skipped due to not being a primary alignment.
    NonPrimary,
    /// Alignment skipped due to being a duplicate.
    Duplicate,
    /// Alignment skipped due to not overlapping with the locus.
    NoOverlap,
    /// Alignment skipped due to not having a cell barcode.
    NoCellBarcode,
    /// Alignment skipped due to not having a UMI.
    NoUmi,
    /// Alignment skipped due to its CIGAR having invalid characters.
    InvalidChars,
    /// Alignment skipped due to being multi-allelic.
    MultiAllelic,
}

/// Results of evaluating alignments.
#[derive(Default, Serialize, Deserialize, Debug)]
pub struct SnpResults {
    /// VCF record index.
    pub rec_idx: VcfRecordIdx,
    /// Metrics for evaluating alignments.
    pub metrics: VcfMetrics,
    /// Cell barcode calls
    pub cell_barcode_calls: Vec<CellBarcodeCalls>,
}

#[derive(Debug, Clone)]
/// Represents a locus in a genome.
struct Locus<'a> {
    /// The chromosome name.
    contig: &'a str,
    /// The start position.
    start: u64,
    /// The end position.
    end: u64,
}

impl Display for Locus<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
    }
}

/// Represents a collection of variant haplotypes.
#[derive(Debug)]
struct VariantHaps<'a> {
    locus: Locus<'a>,
    ref_hap: Seq,
    alt_hap: Seq,
}

/// Represents the scores for a cell barcode and UMI.
#[derive(Default, Debug, Serialize, Deserialize)]
struct Scores {
    /// Index of the cell barcode.
    cb_idx: CellBarcodeIdx,
    /// Unique Molecular Identifier (UMI) vector.
    umi: Seq,
    /// Reference score.
    ref_score: i32,
    /// Alternative score.
    alt_score: i32,
}

/// Cell barcode index and its computed allele counts.
#[derive(Serialize, Deserialize, Debug)]
pub struct CellBarcodeCalls {
    /// The index of the cell barcode
    pub cb_idx: CellBarcodeIdx,
    /// The allele counts for the cell barcode
    pub allele_counts: Counter<AlleleType>,
}

/// Parse VCF records and return an iterator of SNP records.
pub fn snp_records<'a>(
    vcf_reader: &'a mut BcfReader,
    contig_info: &'a TxHashMap<String, BamContigInfo>,
) -> Result<impl Iterator<Item = Result<Snp>> + 'a> {
    let records = vcf_reader
        .records()
        .enumerate()
        .map(|(idx, r)| {
            let rec = r?;
            let Some(rid) = rec.rid() else {
                bail!(VartrixError::VcfRecordRidNotFound {
                    rec: format!("{rec:?}")
                })
            };
            let Ok(name) = rec.header().rid2name(rid) else {
                bail!(VartrixError::VcfRecordRidNameNotFound {
                    rec: format!("{rec:?}"),
                    rid
                })
            };
            let name = str::from_utf8(name)?;
            let Some(contig_info) = contig_info.get(name) else {
                bail!(VartrixError::VcfRecordRidNameNotFound {
                    rec: format!("{rec:?}"),
                    rid
                })
            };

            let alt = match rec.alleles().as_slice() {
                [_ref] => Vec::new(),
                [_ref, alt] => alt.to_vec(),
                _ => return Ok(None), // skip multi-allelic variants
            };

            Ok(Some(Snp {
                idx,
                tid: contig_info.tid,
                start: rec.pos() as u64,
                end: rec.end() as u64,
                alt,
            }))
        })
        .filter_map(Result::transpose);
    Ok(records)
}

/// Write the variants to a file.
///
/// # Arguments
///
/// * `out_variants` - The path to the output file.
/// * `vcf_file` - The path to the VCF file.
///
/// # Returns
///
/// Returns a `Result` indicating success or failure.
pub fn write_variants(out_variants: &Path, vcf_file: &Path) -> Result<()> {
    // write the variants to a TSV file for easy loading into Seraut
    let mut rdr = BcfReader::from_path(vcf_file).context("error opening vcf file")?;
    let mut of = BufWriter::new(File::create(out_variants)?);
    for rec in rdr.records() {
        let rec = rec?;
        writeln!(of, "{chr}_{pos}", chr = get_contig(&rec)?, pos = rec.pos())?;
    }

    Ok(())
}

fn get_contig(rec: &BcfRecord) -> Result<&str> {
    let Some(rid) = rec.rid() else {
        bail!(VartrixError::VcfRecordRidNotFound {
            rec: format!("{rec:?}")
        })
    };
    let Ok(name) = rec.header().rid2name(rid) else {
        bail!(VartrixError::VcfRecordRidNameNotFound {
            rec: format!("{rec:?}"),
            rid
        })
    };
    let Ok(contig) = str::from_utf8(name) else {
        bail!(VartrixError::VcfRecordContigNameNotValidUtf8 {
            rec: format!("{rec:?}")
        })
    };
    Ok(contig)
}

/// Evaluate a SNP record returning its cell barcode calls and metrics.
pub fn evaluate_snp(
    snp: &Snp,
    bam_reader: &mut BamReader,
    fasta_reader: &mut fasta::IndexedReader<File>,
    contig_info_by_tid: &TxHashMap<u32, BamContigInfo>,
    cb_index: &BarcodeIndex,
    settings: &VartrixSettings,
) -> Result<SnpResults> {
    let contig_info = &contig_info_by_tid[&snp.tid];

    let locus = Locus {
        contig: contig_info.name.as_str(),
        start: snp.start,
        end: snp.end,
    };

    bam_reader.fetch((locus.contig.as_bytes(), locus.start, locus.end))?;
    let records = bam_reader.records();

    let haps = construct_haplotypes(
        fasta_reader,
        contig_info.length,
        locus,
        &snp.alt,
        settings.padding,
    )?;

    // make sure our alt is sane; if it is not, bail before alignment and warn the user
    for c in &haps.alt_hap {
        if !settings.valid_chars.contains(c) {
            warn!(
                "Variant at {} has invalid alternative characters. This record will be ignored.",
                haps.locus
            );
            return Ok(SnpResults {
                rec_idx: snp.idx,
                metrics: VcfMetrics {
                    has_invalid_chars: true,
                    ..Default::default()
                },
                ..Default::default()
            });
        }
    }

    let (aln_use_counter, alignment_scores) =
        records.process_results(|iter| evaluate_alns(iter, &haps, cb_index, settings))?;

    Ok(SnpResults {
        rec_idx: snp.idx,
        cell_barcode_calls: parse_scores(&alignment_scores, settings.use_umi),
        metrics: VcfMetrics {
            aln_use_counter,
            ..Default::default()
        },
    })
}

/// Return an `Option` containing the cell barcode index if found, or `None` if not found.
fn get_cell_barcode(
    rec: &BamRecord,
    cb_index: &BarcodeIndex,
    bam_tag: &str,
) -> Result<Option<CellBarcodeIdx>> {
    match rec.aux(bam_tag.as_bytes()) {
        Ok(Aux::String(hp)) => {
            let cb = Barcode::from_str(hp).map_err(|_| VartrixError::InvalidBamRecordBarcode {
                rec: format!("{rec:?}"),
                bc: hp.to_string(),
            })?;
            Ok(cb_index.get(&cb))
        }
        _ => Ok(None),
    }
}

/// Return an `Option` containing the UMI if found, or `None` if not found.
fn get_umi<'a>(rec: &'a BamRecord, bam_tag: &'a str) -> Option<&'a [u8]> {
    match rec.aux(bam_tag.as_bytes()) {
        Ok(Aux::String(umi)) => Some(umi.as_bytes()),
        _ => None,
    }
}

/// Check if an alignment is useful for evaluation.
/// Filters alignments to ensure that they truly overlap the region of interest.
/// Overlap is defined as having an aligned base anywhere in the locus
fn is_overlapping_alignment(start: u64, end: u64, rec: &BamRecord) -> bool {
    let cigar = rec.cigar();
    for i in start..=end {
        let pos = cigar
            .read_pos(i as u32, false, true)
            .with_context(|| format!("BAM record {rec:?} has invalid CIGAR"))
            .unwrap();
        if pos.is_some() {
            return true;
        }
    }
    false
}

fn get_qname(rec: &BamRecord) -> String {
    String::from_utf8(rec.qname().to_vec()).unwrap()
}

fn compute_alignment_scores(
    rec: BamRecord,
    cb_index: &BarcodeIndex,
    settings: &VartrixSettings,
    haps: &VariantHaps<'_>,
) -> (Option<Scores>, AlignmentUseStatus) {
    let qname = get_qname(&rec);
    let skip_read_log = |msg: &str| {
        debug!("{} skipping read {} due to {}", haps.locus, qname, msg);
    };

    if rec.mapq() < settings.mapq {
        skip_read_log("low mapping quality");
        return (None, AlignmentUseStatus::LowMappingQuality);
    }
    if settings.primary_alignments_only && (rec.is_secondary() | rec.is_supplementary()) {
        skip_read_log("not being the primary alignment");
        return (None, AlignmentUseStatus::NonPrimary);
    }
    if settings.skip_duplicate_alignments && rec.is_duplicate() {
        skip_read_log("being a duplicate alignment");
        return (None, AlignmentUseStatus::Duplicate);
    }
    if !is_overlapping_alignment(haps.locus.start, haps.locus.end, &rec) {
        skip_read_log("not overlapping the locus");
        return (None, AlignmentUseStatus::NoOverlap);
    }

    let Some(cb_idx) = get_cell_barcode(&rec, cb_index, &settings.cb_bam_tag).unwrap() else {
        skip_read_log("not having a cell barcode");
        return (None, AlignmentUseStatus::NoCellBarcode);
    };

    let umi = get_umi(&rec, &settings.umi_bam_tag).map(<[u8]>::to_vec);
    if (settings.use_umi) && umi.is_none() {
        skip_read_log("not having a UMI");
        return (None, AlignmentUseStatus::NoUmi);
    }
    // if no UMIs in this dataset, just plug in dummy UMI
    let umi = if !settings.use_umi {
        vec![1_u8]
    } else {
        umi.unwrap()
    };

    let seq = &rec.seq().as_bytes();

    let (ref_alignment, alt_alignment) = align(seq, &haps.ref_hap, &haps.alt_hap);

    debug!(
        "{} {} ref_aln:\n{}",
        haps.locus,
        qname,
        ref_alignment.pretty(seq, &haps.ref_hap, 100)
    );
    debug!(
        "{} {} alt_aln:\n{}",
        haps.locus,
        qname,
        alt_alignment.pretty(seq, &haps.alt_hap, 100)
    );
    debug!(
        "{} {} ref_score: {} alt_score: {}",
        haps.locus, qname, ref_alignment.score, alt_alignment.score
    );
    (
        Some(Scores {
            cb_idx,
            umi,
            ref_score: ref_alignment.score,
            alt_score: alt_alignment.score,
        }),
        AlignmentUseStatus::Useful,
    )
}

fn align(seq: &[u8], ref_hap: &[u8], alt_hap: &[u8]) -> (Alignment, Alignment) {
    let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
    let mut aligner = Aligner::new(GAP_OPEN, GAP_EXTEND, score, K, W);
    let ref_alignment = aligner.local(seq, ref_hap);
    let alt_alignment = aligner.local(seq, alt_hap);

    (ref_alignment, alt_alignment)
}

/// Construct haplotypes for a given locus.
/// This function reads the reference and alternative haplotypes from a FASTA reader.
fn evaluate_alns(
    bam_records: impl Iterator<Item = BamRecord>,
    haps: &VariantHaps<'_>,
    cb_index: &BarcodeIndex,
    settings: &VartrixSettings,
) -> (Counter<AlignmentUseStatus>, Vec<Scores>) {
    // loop over all alignments in the region of interest
    // if the alignments are useful (aligned over this region)
    // perform Smith-Waterman against both haplotypes
    // and report the scores

    debug!("Evaluating record {}", haps.locus);
    let (scores, use_counter): (Vec<_>, Counter<_>) = bam_records
        .map(|rec| compute_alignment_scores(rec, cb_index, settings, haps))
        .collect();
    // filter out the None scores and turn them into a Vec<Scores>
    let scores: Vec<Scores> = scores
        .into_iter()
        .flatten()
        .sorted_by_key(|s| s.cb_idx)
        .collect();
    (use_counter, scores)
}

/// Read the sequence for a given locus from a FASTA reader.
///
/// # Arguments
///
/// * `fasta_reader` - The FASTA reader.
/// * `contig_len` - The length of the contig.
/// * `locus` - The locus to read.
/// * `pad_left` - The number of bases to pad on the left side of the locus.
/// * `pad_right` - The number of bases to pad on the right side of the locus.
///
/// # Returns
///
/// Returns a tuple containing the sequence
fn read_locus(
    fasta_reader: &mut fasta::IndexedReader<File>,
    contig_len: u64,
    locus: &Locus<'_>,
    pad_left: u64,
    pad_right: u64,
) -> Result<Seq> {
    let mut seq = Vec::new();

    let new_start = max(0, locus.start as i64 - pad_left as i64) as u64;
    let new_end = min(locus.end + pad_right, contig_len);
    let contig = locus.contig;

    fasta_reader
        .fetch(contig, new_start, new_end)
        .map_err(|e| VartrixError::FastaFetchError {
            err_msg: format!("tried to fetch {contig}:{new_start}-{new_end} but got error: {e}"),
        })?;

    fasta_reader
        .read(&mut seq)
        .map_err(|e| VartrixError::FastaReadError {
            err_msg: format!("tried to read {contig}:{new_start}-{new_end} but got error: {e}"),
        })?;
    assert!(seq.len() as u64 == new_end - new_start);
    Ok(seq.as_mut_slice().to_ascii_uppercase())
}

/// Construct haplotypes for a given locus.
///
/// # Arguments
///
/// * `fasta_reader` - The FASTA reader.
/// * `contig_lens` - The map of contig lengths.
/// * `locus` - The locus to construct haplotypes for.
/// * `alt` - The alternative allele sequence.
/// * `padding` - Number of bases to pad on either side of the locus.
///
/// # Returns
///
/// Returns a tuple containing the reference haplotype and the alternative haplotype.
fn construct_haplotypes<'a>(
    fasta_reader: &mut fasta::IndexedReader<File>,
    contig_len: u64,
    locus: Locus<'a>,
    alt: &[u8],
    padding: u64,
) -> Result<VariantHaps<'a>> {
    let alt_hap: Result<Seq> = {
        let prefix_loc = Locus {
            contig: locus.contig,
            start: locus.start.saturating_sub(padding),
            end: locus.start,
        };
        let prefix = read_locus(fasta_reader, contig_len, &prefix_loc, 0, 0)?;
        let suffix_loc = Locus {
            contig: prefix_loc.contig,
            start: locus.end,
            end: min(locus.end + padding, contig_len),
        };
        let suffix = read_locus(fasta_reader, contig_len, &suffix_loc, 0, 0)?;
        let alt_hap = prefix
            .into_iter()
            .chain(alt.iter().copied())
            .chain(suffix)
            .collect();
        Ok(alt_hap)
    };

    let alt_hap = alt_hap?;

    let ref_hap = read_locus(fasta_reader, contig_len, &locus, padding, padding)?;
    debug!(
        "{locus} -- ref: {} alt: {}",
        String::from_utf8(ref_hap.clone()).unwrap(),
        String::from_utf8(alt_hap.clone()).unwrap()
    );
    Ok(VariantHaps {
        locus,
        ref_hap,
        alt_hap,
    })
}

fn evaluate_scores(ref_score: i32, alt_score: i32) -> Option<AlleleType> {
    if (ref_score < MIN_SCORE) && (alt_score < MIN_SCORE) {
        return None;
    }
    if ref_score > alt_score {
        return Some(AlleleType::Ref);
    }
    if alt_score > ref_score {
        return Some(AlleleType::Alt);
    }
    assert!(ref_score == alt_score);
    Some(AlleleType::Unknown)
}

fn parse_scores(scores: &[Scores], use_umi: bool) -> Vec<CellBarcodeCalls> {
    // parse the score vector into collapsed calls
    scores
        .iter()
        .chunk_by(|s| s.cb_idx)
        .into_iter()
        .map(|(cell_barcode_index, cell_scores)| {
            if use_umi {
                let allele_counts = count_alleles_using_umi(cell_scores);
                CellBarcodeCalls {
                    allele_counts,
                    cb_idx: cell_barcode_index,
                }
            } else {
                let allele_counts = count_alleles_not_using_umi(cell_scores);
                CellBarcodeCalls {
                    allele_counts,
                    cb_idx: cell_barcode_index,
                }
            }
        })
        .collect()
}

fn count_alleles_not_using_umi<'a>(
    cell_scores: impl Iterator<Item = &'a Scores>,
) -> Counter<AlleleType> {
    // Map of CB to Score objects. This is basically trivial in the non-UMI case
    cell_scores
        .filter_map(|score| evaluate_scores(score.ref_score, score.alt_score))
        .collect()
}

fn count_alleles_using_umi<'a>(
    cell_scores: impl Iterator<Item = &'a Scores>,
) -> Counter<AlleleType> {
    // Map of UMI to Score objects; keep track of all Scores for a given CB/UMI pair
    let parsed_scores: TxHashMap<_, Vec<_>> = cell_scores
        .filter_map(|score| {
            let eval = evaluate_scores(score.ref_score, score.alt_score)?;
            Some((&score.umi, eval))
        })
        .fold(TxHashMap::default(), |mut acc, (umi, eval)| {
            acc.entry(umi).or_default().push(eval);
            acc
        });
    // collapse each UMI into a consensus value
    parsed_scores
        .into_values()
        .map(|alleles| {
            let counts: Counter<_> = alleles.iter().collect();
            let denom = counts.values().sum::<usize>();
            assert!(denom > 0);
            let denom = denom as f64;
            let ref_frac = counts[&AlleleType::Ref] as f64 / denom;
            let alt_frac = counts[&AlleleType::Alt] as f64 / denom;
            if (ref_frac < INTRA_UMI_CONS_THRESHOLD) && (alt_frac < INTRA_UMI_CONS_THRESHOLD) {
                AlleleType::Unknown
            } else if alt_frac >= INTRA_UMI_CONS_THRESHOLD {
                AlleleType::Alt
            } else {
                assert!(ref_frac >= INTRA_UMI_CONS_THRESHOLD);
                AlleleType::Ref
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam;
    use std::string::ToString;

    const FASTA_CONTENT: &str = "\
>chr1
ACGTACGTACGT
";

    #[ctor::ctor]
    fn init() {
        // this ensures insta knows where to find its snap tests
        let cwd = std::env::current_dir().unwrap();
        let workspace_root = cwd.parent().unwrap();
        unsafe { std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root) }
    }

    #[test]
    fn test_evaluate_scores() {
        let high_score = MIN_SCORE + 1;
        let low_score = MIN_SCORE;
        let below_min_score = MIN_SCORE - 1;
        assert_eq!(
            evaluate_scores(high_score, low_score),
            Some(AlleleType::Ref)
        );
        assert_eq!(
            evaluate_scores(low_score, high_score),
            Some(AlleleType::Alt)
        );
        assert_eq!(evaluate_scores(below_min_score, below_min_score), None);
    }

    #[test]
    fn test_collapse_calls_not_using_umi() {
        let scores = [
            Scores {
                // Alt: Alt > Ref >= MIN_SCORE
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE + 1,
                ..Default::default()
            },
            Scores {
                // Ref: Ref > Alt >= MIN_SCORE
                ref_score: MIN_SCORE + 1,
                alt_score: MIN_SCORE,
                ..Default::default()
            },
            Scores {
                // Unknown: Ref = Alt >= MIN_SCORE
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE,
                ..Default::default()
            },
            Scores {
                // Skipped because both scores are below MIN_SCORE
                ref_score: MIN_SCORE - 1,
                alt_score: MIN_SCORE - 1,
                ..Default::default()
            },
            Scores {
                // Ref: Ref >= MIN_SCORE and Ref < Alt
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE - 1,
                ..Default::default()
            },
        ];

        let act_allele_counts = count_alleles_not_using_umi(scores.iter());
        let exp_allele_counts: Counter<_> = vec![
            (AlleleType::Ref, 2),
            (AlleleType::Alt, 1),
            (AlleleType::Unknown, 1),
        ]
        .into_iter()
        .collect();
        assert_eq!(act_allele_counts, exp_allele_counts);
    }

    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_valid_consensus_threshold() {
        assert!(
            INTRA_UMI_CONS_THRESHOLD > 0.5,
            "Consensus threshold should be above 0.5"
        );
        assert!(
            INTRA_UMI_CONS_THRESHOLD <= 1.0,
            "Consensus threshold should be below 1.0"
        );
    }

    #[test]
    fn test_collapse_calls_using_umi() {
        let mut scores = Vec::new();
        // First UMI is >= INTRA_UMI_CONS_THRESHOLD ref --> Ref++
        let total = 5;
        let dominant_allele_count = 1 + (total as f64 * INTRA_UMI_CONS_THRESHOLD) as usize;
        for _ in 0..dominant_allele_count {
            scores.push(Scores {
                ref_score: MIN_SCORE + 1,
                alt_score: MIN_SCORE,
                umi: vec![1],
                ..Default::default()
            });
        }
        for _ in 0..(total - dominant_allele_count) {
            scores.push(Scores {
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE + 1,
                umi: vec![1],
                ..Default::default()
            });
        }
        // Second UMI is >= INTRA_UMI_CONS_THRESHOLD alt --> Alt++
        for _ in 0..dominant_allele_count {
            scores.push(Scores {
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE + 1,
                umi: vec![2],
                ..Default::default()
            });
        }
        for _ in 0..(total - dominant_allele_count) {
            scores.push(Scores {
                ref_score: MIN_SCORE + 1,
                alt_score: MIN_SCORE,
                umi: vec![2],
                ..Default::default()
            });
        }
        // Third UMI has no dominant allele (1:1) --> Unknown++
        for _ in 0..(total / 2) {
            scores.push(Scores {
                ref_score: MIN_SCORE + 1,
                alt_score: MIN_SCORE,
                umi: vec![3],
                ..Default::default()
            });
        }
        for _ in 0..(total - (total / 2)) {
            scores.push(Scores {
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE + 1,
                umi: vec![3],
                ..Default::default()
            });
        }
        // Fourth UMI has no dominant allele (1:1:1) --> Unknown++
        for _ in 0..(total / 3) {
            scores.push(Scores {
                ref_score: MIN_SCORE + 1,
                alt_score: MIN_SCORE,
                umi: vec![4],
                ..Default::default()
            });
        }
        for _ in 0..(total - (total / 3)) {
            scores.push(Scores {
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE + 1,
                umi: vec![4],
                ..Default::default()
            });
        }
        for _ in 0..(total / 3) {
            scores.push(Scores {
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE,
                umi: vec![4],
                ..Default::default()
            });
        }
        // Fifth UMI is >= INTRA_UMI_CONS_THRESHOLD Unknown --> Unknown++
        for _ in 0..dominant_allele_count {
            scores.push(Scores {
                ref_score: MIN_SCORE,
                alt_score: MIN_SCORE,
                umi: vec![5],
                ..Default::default()
            });
        }
        for _ in 0..(total - dominant_allele_count) {
            scores.push(Scores {
                ref_score: MIN_SCORE + 1,
                alt_score: MIN_SCORE,
                umi: vec![5],
                ..Default::default()
            });
        }
        // Sixth UMI has only low scores --> Skipped
        for _ in 0..total {
            scores.push(Scores {
                ref_score: MIN_SCORE - 1,
                alt_score: MIN_SCORE - 1,
                umi: vec![6],
                ..Default::default()
            });
        }

        let act_allele_counts = count_alleles_using_umi(scores.iter());
        let exp_allele_counts: Counter<_> = vec![
            (AlleleType::Ref, 1),
            (AlleleType::Alt, 1),
            (AlleleType::Unknown, 3),
        ]
        .into_iter()
        .collect();
        assert_eq!(act_allele_counts, exp_allele_counts);
    }

    #[test]
    fn test_read_locus() {
        let mut fasta = tempfile::NamedTempFile::with_suffix(".fasta").unwrap();
        fasta.write_all(FASTA_CONTENT.as_bytes()).unwrap();
        fasta.flush().unwrap();
        rust_htslib::faidx::build(fasta.path()).expect("Failed to index FASTA file");
        let mut fasta_reader = fasta::IndexedReader::from_file(&fasta).unwrap();

        let contig_len = fasta_reader
            .index
            .sequences()
            .iter()
            .find(|s| s.name == "chr1")
            .unwrap()
            .len;
        let locus = Locus {
            contig: "chr1",
            start: 4,
            end: 6,
        };
        let seq = read_locus(&mut fasta_reader, contig_len, &locus, 2, 2).unwrap();
        assert_eq!(str::from_utf8(&seq).unwrap(), "GTACGT");

        let locus = Locus {
            contig: "chr2", // Invalid contig
            start: 4,
            end: 8,
        };
        let act_err = read_locus(&mut fasta_reader, contig_len, &locus, 0, 0)
            .unwrap_err()
            .downcast::<VartrixError>()
            .unwrap();
        let exp_err = VartrixError::FastaFetchError {
            err_msg: "tried to fetch chr2:4-8 but got error: Unknown sequence name: chr2."
                .to_string(),
        };
        assert_eq!(act_err, exp_err);

        let locus = Locus {
            contig: "chr1",
            start: 2,
            end: 1000,
        };
        let seq = read_locus(&mut fasta_reader, contig_len, &locus, 10, 10).unwrap();
        assert_eq!(str::from_utf8(&seq).unwrap(), "ACGTACGTACGT");

        let locus = Locus {
            contig: "chr1",
            start: 5,
            end: 2,
        };
        let act_err = read_locus(&mut fasta_reader, contig_len, &locus, 0, 0)
            .unwrap_err()
            .downcast::<VartrixError>()
            .unwrap();
        let exp_err = VartrixError::FastaReadError {
            err_msg: "tried to read chr1:5-2 but got error: Invalid query interval".to_string(),
        };
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn test_get_cell_barcode() {
        let mut record = BamRecord::new();
        record.push_aux(b"CB", Aux::String("ACGT-1")).unwrap();
        let record_no_cb = BamRecord::new();
        let mut record_invalid = BamRecord::new();
        record_invalid.push_aux(b"CB", Aux::String("ACGT")).unwrap();

        let cb_index = BarcodeIndex::from_sorted(vec![
            Barcode::from_str("ACGT-1").unwrap(),
            Barcode::from_str("TGCA-1").unwrap(),
        ]);

        // Test valid cell barcode
        assert_eq!(get_cell_barcode(&record, &cb_index, "CB").unwrap(), Some(0));
        // Test incorrect bam_tag
        assert_eq!(get_cell_barcode(&record, &cb_index, "XX").unwrap(), None);
        // Test record with missing cell barcode
        assert_eq!(
            get_cell_barcode(&record_no_cb, &cb_index, "CB").unwrap(),
            None
        );
        // Test record with invalid cell barcode
        let act_err = get_cell_barcode(&record_invalid, &cb_index, "CB")
            .unwrap_err()
            .downcast::<VartrixError>()
            .unwrap();
        let exp_err = VartrixError::InvalidBamRecordBarcode {
            rec: format!("{record_invalid:?}"),
            bc: "ACGT".to_string(),
        };
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn test_get_umi() {
        let mut record = BamRecord::new();
        record
            .push_aux(UMI_BAM_TAG.as_bytes(), Aux::String("ACGT"))
            .unwrap();

        let umi = get_umi(&record, UMI_BAM_TAG);
        assert_eq!(umi, Some(b"ACGT".as_ref()));

        let record_no_umi = BamRecord::new();
        let umi_none = get_umi(&record_no_umi, UMI_BAM_TAG);
        assert_eq!(umi_none, None);
    }

    #[test]
    fn test_is_overlapping_alignment() {
        let mut rec = BamRecord::new();
        let qname = b"read1";
        let seq = b"ACCGT";
        let qual = b"KKKKK";
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(5)]);
        rec.set(qname, Some(&cigar), seq, qual);
        let start = 10;
        let end = 10;

        // Before the locus
        rec.set_pos(5);
        assert!(!is_overlapping_alignment(start, end, &rec));
        // Overlaps the locus
        for i in 6..=10 {
            rec.set_pos(i);
            assert!(is_overlapping_alignment(start, end, &rec));
        }
        // After the locus
        rec.set_pos(11);
        assert!(!is_overlapping_alignment(start, end, &rec));
    }

    #[test]
    fn test_construct_haplotypes_success() {
        let mut fasta = tempfile::NamedTempFile::with_suffix(".fasta").unwrap();
        fasta.write_all(FASTA_CONTENT.as_bytes()).unwrap();
        fasta.flush().unwrap();
        rust_htslib::faidx::build(fasta.path()).expect("Failed to index FASTA file");
        let mut fasta_reader = fasta::IndexedReader::from_file(&fasta).unwrap();
        let locus = Locus {
            contig: "chr1",
            start: 4,
            end: 8,
        };
        let contig_lens: TxHashMap<_, _> = fasta_reader
            .index
            .sequences()
            .into_iter()
            .map(|s| (s.name, s.len))
            .collect();

        let alt = b"TT";
        let haps = construct_haplotypes(
            &mut fasta_reader,
            contig_lens["chr1"],
            locus.clone(),
            alt,
            2,
        )
        .unwrap();
        assert_eq!(str::from_utf8(&haps.ref_hap).unwrap(), "GTACGTAC");
        assert_eq!(str::from_utf8(&haps.alt_hap).unwrap(), "GTTTAC");

        let alt = b"";
        let haps =
            construct_haplotypes(&mut fasta_reader, contig_lens["chr1"], locus, alt, 2).unwrap();
        assert_eq!(str::from_utf8(&haps.ref_hap).unwrap(), "GTACGTAC");
        assert_eq!(str::from_utf8(&haps.alt_hap).unwrap(), "GTAC");
    }
}
