//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::collections::{HashMap, HashSet};
use std::io::{BufReader, BufRead};
use std::iter::FromIterator;
use std::string::String;
use std::cmp;
use std::path::{Path, PathBuf};
use std::str;
use std::f64;
use failure::Error;

use rust_htslib::bam::record::{Record, Aux};

use utils;

const RAW_BC_SEQ_TAG: &'static str   = "CR";
const RAW_BC_QUAL_TAG: &'static str  = "CY";
const PROC_BC_SEQ_TAG: &'static str  = "CB";
const RAW_UMI_SEQ_TAG: &'static str  = "UR";
const RAW_UMI_QUAL_TAG: &'static str = "UY";
const PROC_UMI_SEQ_TAG: &'static str = "UB";
const SI_SEQ_TAG: &'static str       = "BC";
const SI_QUAL_TAG: &'static str      = "QT";

const ILLUMINA_QUAL_OFFSET: u8 = 33;

// TODO constants below should be configurable

const BC_CONF_THRESHOLD: f64 = 0.975;
const BC_MAX_QV: u8 = 33; // cap observed QVs at this value
const BC_MIN_QV: u8 = 10;

const UMI_POLYT_SUFFIX_LENGTH: usize = 5;
const UMI_MIN_QV: u8 = 10;

#[derive(Debug)]
pub struct SampleIndexData {
    pub seq: String,
    pub qual: Vec<u8>,
}

#[derive(Debug)]
pub struct BarcodeData {
    pub raw_seq:        String,
    pub qual:           Vec<u8>,
    pub processed_seq:  Option<String>,
    pub on_whitelist:   bool,
    pub low_min_qual:   bool,
    pub has_n:          bool,
    pub is_homopolymer: bool,
}

#[derive(Debug)]
pub struct UmiData {
    pub raw_seq:        String,
    pub qual:           Vec<u8>,
    pub is_valid:       bool,
    pub low_min_qual:   bool,
    pub has_n:          bool,
    pub is_homopolymer: bool,
    pub has_primer:     bool,
    pub has_polyt:      bool,
}

// TODO move these sequence functions into utils?

fn has_n(seq: &[u8]) -> bool {
    return seq.contains(&('N' as u8))
}

fn has_polyt(seq: &[u8]) -> bool {
    let end_idx = if UMI_POLYT_SUFFIX_LENGTH <= seq.len() { seq.len() - UMI_POLYT_SUFFIX_LENGTH } else { seq.len() };
    for base in seq[..end_idx].iter() {
        if *base != 'T' as u8 {
            return false
        }
    }
    return true
}

fn has_primer(seq: &[u8]) -> bool {
    // TODO implement. for now it doesn't matter because we don't list any primer sequences when checking UMIs in the Python code
    let _ = seq;
    return false
}

fn is_homopolymer(seq: &[u8]) -> bool {
    // TODO this would call NNNNN as a homopolymer
    for i in 1..seq.len() {
        if seq[i-1] != seq[i] {
            return false
        }
    }
    return true
}

fn low_min_qual(qual: &[u8], min_qv: u8) -> bool {
    for qv in qual {
        if qv - ILLUMINA_QUAL_OFFSET < min_qv {
            return true
        }
    }
    return false
}

fn load_barcode_dist(bc_counts_json: &str,
                     bc_whitelist_txt: &str,
                     gem_group: &u8,
                     proportions: bool,
                     library_type: &str) -> HashMap<String, f64> {
    // TODO untangle this mess and handle errors
    let bc_counts: Vec<f64> = utils::read_json(bc_counts_json)
        .get(library_type)
        .expect(&format!("Failed to find library type {} in barcode counts JSON", library_type))
        .as_array()
        .expect("Failed to load barcode counts JSON data as array")
        .to_vec().iter().map(|elt| elt.as_f64().unwrap()).collect();
    let bc_whitelist: Vec<String> = utils::load_txt_maybe_gz(bc_whitelist_txt);
    let start_idx = (gem_group - 1) as usize * bc_whitelist.len();
    let end_idx = *gem_group as usize * bc_whitelist.len();
    let mut bc_dist = HashMap::new();
    let mut total_counts = 0_f64;
    for (bc, value) in bc_whitelist.iter().zip(bc_counts[start_idx..end_idx].iter()) {
        total_counts += *value;
        bc_dist.insert(bc.clone(), *value);
    }
    if proportions {
        for val in bc_dist.values_mut() {
            *val /= total_counts;
        }
    }
    return bc_dist;
}

fn load_barcode_whitelist(bc_whitelist_txt: &str) -> HashSet<String> {
    let bc_whitelist: Vec<String> = utils::load_txt_maybe_gz(bc_whitelist_txt);
    return HashSet::from_iter(bc_whitelist.into_iter())
}

/// Given a barcode whitelist path, return the expected
/// path of its translation file (e.g., xyz/abc.txt => xyz/translation/abc.txt)
fn get_translation_path<P: AsRef<Path>>(whitelist_path: P) -> PathBuf {
    whitelist_path.as_ref().parent().unwrap_or(Path::new(&"."))
        .join("translation")
        .join(whitelist_path.as_ref().file_name().unwrap())
}

fn load_barcode_translation<P: AsRef<Path>>(whitelist_path: P,
                                            whitelist: &HashSet<String>) ->
    Option<HashMap<String, String>> {
        let path = get_translation_path(whitelist_path);
        match path.exists() {
            false => None,
            true => {
                let reader = BufReader::new(utils::open_maybe_compressed(path));
                let mut translation = HashMap::new();

                for maybe_line in reader.lines() {
                    let line = maybe_line.ok().expect("Failed to read line");
                    let mut row = line.split('\t');
                    let src = row.next().expect("Expected first column containing source barcode");
                    let dest = row.next().expect("Expected second column containing destination barcode");
                    if !whitelist.contains(src) {
                        panic!("Barcode translation contains sequence not on whitelist: {}", src);
                    }
                    if !whitelist.contains(dest) {
                        panic!("Barcode translation contains sequence not on whitelist: {}", dest)
                    }
                    translation.insert(src.to_owned(), dest.to_owned());
                }
                Some(translation)
            },
        }
}

pub struct BarcodeUmiData {
    pub sample_index_data: Option<SampleIndexData>,
    pub barcode_data:      Option<BarcodeData>,
    pub umi_data:          Option<UmiData>,
}

pub struct BarcodeUmiChecker {
    // NOTE: barcode_whitelist is redundant if barcode_dist is defined, but we should try not to couple them
    barcode_whitelist: HashSet<String>,
    barcode_translation: Option<HashMap<String, String>>,
    barcode_dist: HashMap<String, f64>,
    barcodes_checkable: bool,
}

impl BarcodeUmiChecker {
    pub fn new(bc_counts_json: &str,
               bc_whitelist_txt: &str,
               gem_group: &u8,
               translate: bool,
               library_type: &str) -> BarcodeUmiChecker {
        let mut barcodes_checkable = false;
        let mut barcode_whitelist = HashSet::new();
        let mut barcode_translation = None;
        let mut barcode_dist = HashMap::new();
        if bc_whitelist_txt != "null" {
            barcodes_checkable = true;
            barcode_whitelist = load_barcode_whitelist(bc_whitelist_txt);
            if translate {
                barcode_translation = load_barcode_translation(bc_whitelist_txt,
                                                               &barcode_whitelist);
            }
            barcode_dist = load_barcode_dist(bc_counts_json, bc_whitelist_txt,
                                             &gem_group, true, library_type);
        }
        return BarcodeUmiChecker {
            barcode_whitelist: barcode_whitelist,
            barcode_translation: barcode_translation,
            barcode_dist: barcode_dist,
            barcodes_checkable: barcodes_checkable,
        }
    }

    fn correct_bc_error(&self, seq: &[u8], qual: &[u8]) -> Option<String> {
        if !self.barcodes_checkable {
            return None
        }
        let nucs = "ACGT".as_bytes();
        let mut test_seq = Vec::new();
        test_seq.extend_from_slice(seq);
        let mut whitelist_likelihoods = HashMap::new();
        let mut likelihood_sum = 0.0;

        // compute likelihoods for Hamming distance 1 sequences
        for i in 0..test_seq.len() {
            let orig_base = test_seq[i];
            for base in nucs {
                if *base != orig_base {
                    test_seq[i] = *base;
                    let test_str = str::from_utf8(&test_seq).unwrap().to_string();
                    match self.barcode_dist.get(&test_str) {
                        Some(p_whitelist) => {
                            let qv = cmp::min(qual[i] - ILLUMINA_QUAL_OFFSET, BC_MAX_QV) as f64;
                            let p_edit = 10.0_f64.powf(-qv / 10.0);
                            let likelihood = p_whitelist * p_edit;
                            whitelist_likelihoods.insert(test_str.clone(), likelihood);
                            likelihood_sum += likelihood;
                        }
                        None => (),
                    }
                }
            }
            test_seq[i] = orig_base;
        }

        // find maximum likelihood
        let mut max_likelihood = -1.0;
        let mut best_whitelist_bc = String::new();
        for (whitelist_bc, likelihood) in whitelist_likelihoods {
            if likelihood > max_likelihood {
                max_likelihood = likelihood;
                best_whitelist_bc = whitelist_bc;
            }
        }

        if max_likelihood / likelihood_sum >= BC_CONF_THRESHOLD {
            return Some(best_whitelist_bc)
        } else {
            return None
        }
    }

    fn is_whitelisted(&self, seq: &[u8]) -> bool {
        if !self.barcodes_checkable {
            return true
        }
        let seq_string = str::from_utf8(seq).unwrap().to_string();
        return self.barcode_whitelist.contains(&seq_string)
    }

    fn check_barcode(&self, seq: &[u8], qual: &[u8]) -> BarcodeData {
        let on_whitelist = self.is_whitelisted(seq);
        let low_min_qual = low_min_qual(qual, BC_MIN_QV);
        let is_homopolymer = is_homopolymer(seq);
        let has_n = has_n(seq);
        let raw_seq_str = str::from_utf8(seq).unwrap().to_string();
        let corrected_seq = if on_whitelist {
            Some(raw_seq_str.clone())
        } else {
            self.correct_bc_error(seq, qual)
        };

        // Translate barcode if a translation is specified
        let processed_seq = match self.barcode_translation {
            // No translation map; don't translate
            None => corrected_seq,
            // There is a translation map; attempt translation if bc on whitelist
            Some(ref tmap) =>
                corrected_seq.map(|corr_seq| tmap.get(&corr_seq)
                                  .expect(&format!("Whitelisted barcode not in translation map: {}.", &corr_seq))).cloned(),
        };

        return BarcodeData {
            raw_seq:        raw_seq_str,
            qual:           qual.to_owned(),
            processed_seq:  processed_seq,
            on_whitelist:   on_whitelist,
            low_min_qual:   low_min_qual,
            has_n:          has_n,
            is_homopolymer: is_homopolymer,
        }
    }

    fn check_umi(&self, seq: &[u8], qual: &[u8]) -> UmiData {
        let has_n = has_n(seq);
        let has_polyt = has_polyt(seq);
        let has_primer = has_primer(seq);
        let is_homopolymer = is_homopolymer(seq);
        let low_min_qual = low_min_qual(qual, UMI_MIN_QV);
        let is_valid = !(has_n || is_homopolymer || low_min_qual);
        let raw_seq_str = str::from_utf8(seq).unwrap().to_string();
        UmiData {
            raw_seq:        raw_seq_str,
            qual:           qual.to_owned(),
            is_valid:       is_valid,
            low_min_qual:   low_min_qual,
            has_n:          has_n,
            is_homopolymer: is_homopolymer,
            has_primer:     has_primer,
            has_polyt:      has_polyt
        }
    }

    fn get_sample_index_data(&self, bc_tags: &HashMap<String, String>) -> Option<SampleIndexData> {
        let seq = bc_tags.get(SI_SEQ_TAG.into());
        let qual = bc_tags.get(SI_QUAL_TAG.into());
        match (seq, qual) {
            (Some(seq), Some(qual)) => Some(SampleIndexData{
                seq: seq.to_owned(),
                qual: qual.as_bytes().to_owned(),
            }),
            _ => None,
        }
    }

    fn get_barcode_data(&self, bc_tags: &HashMap<String, String>) -> Option<BarcodeData> {
        let bc_seq = bc_tags.get(RAW_BC_SEQ_TAG.into());
        let bc_qual = bc_tags.get(RAW_BC_QUAL_TAG.into());
        match (bc_seq, bc_qual) {
            (Some(seq), Some(qual)) => Some(self.check_barcode(seq.as_bytes(), qual.as_bytes())),
            _ => None,
        }
    }

    fn get_umi_data(&self, bc_tags: &HashMap<String, String>) -> Option<UmiData> {
        let umi_seq = bc_tags.get(RAW_UMI_SEQ_TAG.into());
        let umi_qual = bc_tags.get(RAW_UMI_QUAL_TAG.into());
        match (umi_seq, umi_qual) {
            (Some(seq), Some(qual)) => Some(self.check_umi(seq.as_bytes(), qual.as_bytes())),
            _ => None,
        }
    }

    pub fn process_barcodes_and_umis(&self, bc_tags: &HashMap<String, String>) -> BarcodeUmiData {
        let si_data = self.get_sample_index_data(&bc_tags);
        let bc_data = self.get_barcode_data(&bc_tags);
        let umi_data = self.get_umi_data(&bc_tags);

        BarcodeUmiData {
            sample_index_data: si_data,
            barcode_data:      bc_data,
            umi_data:          umi_data,
        }
    }
}

impl BarcodeUmiData {
    /// Add tags to a BAM record
    pub fn attach_tags(&self, record: &mut Record, gem_group: &u8) -> Result<(), Error> {
        // Attach sample index
        if let Some(ref si_data) = self.sample_index_data {
            record.push_aux(&SI_SEQ_TAG.as_bytes(),
                            &Aux::String(si_data.seq.as_bytes()))?;
            record.push_aux(&SI_QUAL_TAG.as_bytes(),
                            &Aux::String(&si_data.qual))?;
        };

        if let Some(ref bc_data) = self.barcode_data {
            record.push_aux(&RAW_BC_SEQ_TAG.as_bytes(),
                            &Aux::String(bc_data.raw_seq.as_bytes()))?;
            record.push_aux(&RAW_BC_QUAL_TAG.as_bytes(),
                            &Aux::String(&bc_data.qual))?;
            if let Some(ref corrected_seq) = bc_data.processed_seq {
                let processed_bc = utils::get_processed_bc(&corrected_seq, gem_group);
                record.push_aux(&PROC_BC_SEQ_TAG.as_bytes(),
                                &Aux::String(processed_bc.as_bytes()))?;
            }
        }

        // Attach UMI
        if let Some(ref umi_data) = self.umi_data {
            record.push_aux(&RAW_UMI_SEQ_TAG.as_bytes(),
                            &Aux::String(umi_data.raw_seq.as_bytes()))?;
            record.push_aux(&RAW_UMI_QUAL_TAG.as_bytes(),
                            &Aux::String(&umi_data.qual))?;
            if umi_data.is_valid {
                record.push_aux(&PROC_UMI_SEQ_TAG.as_bytes(),
                                &Aux::String(umi_data.raw_seq.as_bytes()))?;
            }
        };

        Ok(())
    }
}
