use crate::bam_tags::{
    ExtraFlags, EXTRA_FLAGS_TAG, FEATURE_IDS_TAG, PROC_BC_SEQ_TAG, PROC_UMI_SEQ_TAG,
    RAW_BARCODE_SEQ_TAG, RAW_UMI_SEQ_TAG,
};
use crate::constants::{ALN_BC_DISK_CHUNK_SZ, ALN_BC_GIB};
use barcode::{Barcode, BarcodeContent};
use rust_htslib::bam::record::{Aux, Cigar, Record};
use shardio::{ShardReader, SortKey, SHARD_ITER_SZ as SHARD_SZ};
use std::borrow::Cow;
use std::hash::{Hash, Hasher};
use umi::UmiSeq;

// a gibibyte is 1024**3 bytes

/* BAM */

pub fn set_primary(record: &mut Record) {
    // NOTE: rust_htslib doesn't have a method for marking an alignment as primary :/
    if record.is_secondary() {
        record.inner_mut().core.flag -= 256u16;
    }
}

pub fn alen(read: &Record) -> i64 {
    // NOTE: end_pos in rust_htslib was recently fixed to account for deletions, we could update to that instead
    let mut alen = 0;

    for c in &read.cigar() {
        match *c {
            Cigar::Match(l)
            | Cigar::Del(l)
            | Cigar::RefSkip(l)
            | Cigar::Equal(l)
            | Cigar::Diff(l) => alen += l as i64,
            _ => (),
        }
    }
    alen
}

pub trait AuxExt {
    fn integer(&self) -> i32;
    fn str(&self) -> &str;
    fn str_bytes(&self) -> &[u8] {
        self.str().as_bytes()
    }
}

impl<'a> AuxExt for Aux<'a> {
    fn integer(&self) -> i32 {
        // From the spec: <https://samtools.github.io/hts-specs/SAMv1.pdf>
        // While all single (i.e., non-array) integer types are stored in SAM
        // as ‘i’, in BAM any of ‘cCsSiI’ may be used together with the
        // correspondingly-sized binary integer value, chosen according to
        // the field value’s magnitude
        //
        // So the type we set need not be equal to the type we read :/
        match self {
            Aux::I8(i) => *i as i32,
            Aux::U8(i) => *i as i32,
            Aux::I16(i) => *i as i32,
            Aux::U16(i) => *i as i32,
            Aux::I32(i) => *i,
            Aux::U32(i) => *i as i32,
            a => panic!("Expected integer, got {a:?}"),
        }
    }
    fn str(&self) -> &str {
        match self {
            Aux::String(s) => s,
            a => panic!("Expected Aux::String, got {a:?}"),
        }
    }
}

pub trait CrRecord {
    fn raw_barcode(&self) -> Option<BarcodeContent>;
    fn processed_barcode(&self) -> Option<Barcode>;

    fn raw_umi(&self) -> Option<UmiSeq>;
    fn processed_umi(&self) -> Option<UmiSeq>;
    fn is_umi_count(&self) -> bool;

    fn cr_extra_flags(&self) -> ExtraFlags;
    fn is_dup_candidate(&self) -> bool;
}

impl CrRecord for Record {
    fn raw_barcode(&self) -> Option<BarcodeContent> {
        self.aux(RAW_BARCODE_SEQ_TAG)
            .map(|x| {
                x.str()
                    .parse::<BarcodeContent>()
                    .unwrap_or_else(|e| panic!("Invalid raw barcode {}: {e}", x.str()))
            })
            .ok()
    }

    fn processed_barcode(&self) -> Option<Barcode> {
        self.aux(PROC_BC_SEQ_TAG)
            .map(|x| Barcode::from_bytes(x.str_bytes()).unwrap())
            .ok()
    }

    fn raw_umi(&self) -> Option<UmiSeq> {
        self.aux(RAW_UMI_SEQ_TAG)
            .map(|x| UmiSeq::from_bytes(x.str_bytes()))
            .ok()
    }

    fn processed_umi(&self) -> Option<UmiSeq> {
        self.aux(PROC_UMI_SEQ_TAG)
            .map(|x| UmiSeq::from_bytes(x.str_bytes()))
            .ok()
    }

    fn is_umi_count(&self) -> bool {
        self.cr_extra_flags().intersects(ExtraFlags::UMI_COUNT)
    }

    /// Get an alignment record's extra flags
    fn cr_extra_flags(&self) -> ExtraFlags {
        self.aux(EXTRA_FLAGS_TAG).map_or(Default::default(), |x| {
            ExtraFlags::from_bits_truncate(x.integer() as u32)
        })
    }

    /// Return true if a read should be considered for PCR-duplicate marking
    fn is_dup_candidate(&self) -> bool {
        let is_primary_alignment = !self.is_secondary();
        let extra_flags = self.cr_extra_flags();

        let is_conf_mapped = extra_flags.intersects(ExtraFlags::CONF_FEATURE)
            || extra_flags.intersects(ExtraFlags::CONF_MAPPED);
        is_primary_alignment && is_conf_mapped
    }
}

/// Get an alignment record's cell barcode
pub fn get_read_barcode(record: &Record) -> Option<Vec<u8>> {
    record
        .aux(PROC_BC_SEQ_TAG)
        .map(|x| x.str_bytes().to_owned())
        .ok()
}

/// Get an alignment record's processed UMI sequence
pub fn get_read_umi(record: &Record) -> Option<Vec<u8>> {
    record
        .aux(PROC_UMI_SEQ_TAG)
        .map(|x| x.str_bytes().to_owned())
        .ok()
}

/// Produce a processed BC sequence.
pub fn get_processed_bc(corrected_seq: &str, gem_group: u8) -> String {
    format!("{corrected_seq}-{gem_group}")
}

pub fn is_read_umi_count(read: &Record) -> bool {
    read.cr_extra_flags().intersects(ExtraFlags::UMI_COUNT)
}

pub fn is_read_low_support_umi(read: &Record) -> bool {
    read.cr_extra_flags()
        .intersects(ExtraFlags::LOW_SUPPORT_UMI)
}

pub fn is_read_conf_mapped_to_transcriptome(read: &Record, high_conf_mapq: u8) -> bool {
    if (read.is_unmapped()) || (read.mapq() < high_conf_mapq) {
        false
    } else {
        get_read_gene_ids(read).is_some_and(|l| l.len() == 1)
    }
}

pub fn is_read_conf_mapped_to_feature(read: &Record) -> bool {
    read.cr_extra_flags().intersects(ExtraFlags::CONF_FEATURE)
}

pub fn get_read_gene_ids(record: &Record) -> Option<Vec<String>> {
    record
        .aux(RAW_UMI_SEQ_TAG)
        .map(|x| x.str().split(';').map(String::from).collect())
        .ok()
}

/// Get an alignment record's processed UMI sequence
pub fn get_read_raw_umi(record: &Record) -> Option<Vec<u8>> {
    record
        .aux(RAW_UMI_SEQ_TAG)
        .map(|x| x.str_bytes().to_owned())
        .ok()
}

/// Return true if a read should be considered for PCR-duplicate marking
pub fn is_read_dup_candidate(record: &Record) -> bool {
    let is_primary_alignment = !record.is_secondary();
    let extra_flags = record.cr_extra_flags();
    let is_conf_mapped = extra_flags.intersects(ExtraFlags::CONF_FEATURE)
        || extra_flags.intersects(ExtraFlags::CONF_MAPPED);
    is_primary_alignment && is_conf_mapped
}

/// Get the metric prefix for a given library type.
/// Some are hardcoded for historical reasons.
pub fn get_library_type_metric_prefix(lib_type: &str) -> String {
    match lib_type {
        "Gene Expression" => String::new(),
        "CRISPR Guide Capture" => "CRISPR_".to_owned(),
        "Antibody Capture" => "ANTIBODY_".to_owned(),
        _ => format!("{lib_type}_"),
    }
}

/// Marker trait to sort BAM records by their position, used for ShardIO
pub struct BamPosSort;

// fxhash has issues with the least bytes of byte arrays, use siphash
fn hash32<T: Hash + ?Sized>(v: &T) -> u32 {
    let h = hash64(v);
    let mask = u32::MAX as u64;
    (h & mask) as u32 ^ ((h >> 32) & mask) as u32
}

pub fn hash64<T: Hash + ?Sized>(v: &T) -> u64 {
    let mut s = wyhash::WyHash::with_seed(0);
    v.hash(&mut s);
    s.finish()
}

impl SortKey<Record> for BamPosSort {
    type Key = (u32, i64, bool, u64, u32, u16);

    fn sort_key(rec: &Record) -> Cow<'_, Self::Key> {
        Cow::Owned((
            rec.tid() as u32,
            rec.pos(),
            // We parse unmapped reads in REPORT_UNMAPPED_READS_PD stage downstream
            // This allows us to avoid decompresssing blocks of reads that
            // are mapped to the genome or feature barcodes
            rec.is_unmapped() && rec.aux(FEATURE_IDS_TAG).is_err(),
            hash64(rec.qname()),
            hash32(rec.raw_cigar()),
            rec.flags(),
        ))
    }
}

/// Return the estimated memory usage in GiB of reading this shard range.
pub fn estimate_memory_for_range_gib(
    range: &shardio::Range<<BamPosSort as SortKey<Record>>::Key>,
    reader: &ShardReader<Record, BamPosSort>,
) -> f64 {
    // BufReader + lz4 buffer + 1KB record
    let n_blocks = (reader.est_len_range(range) as f64 / ALN_BC_DISK_CHUNK_SZ as f64).ceil();
    (SHARD_SZ + 1_024) as f64 * n_blocks / ALN_BC_GIB
}
