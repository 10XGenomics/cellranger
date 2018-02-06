//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::io::{BufRead, BufReader, Read, Write};
use std::str::FromStr;
use std::fs::File;
use std::collections::{HashMap};
use std::path::Path;

use rust_htslib::bam::record::{Record, Cigar};
use serde::de::DeserializeOwned;
use serde_json;
use serde_json::Value;
use csv;
use flate2::read::{GzDecoder};

pub const NUCLEOTIDES: &'static [u8] = b"ACGT";
pub const ILLUMINA_QUAL_OFFSET: u8 = 33;

/* File IO */

pub fn load_txt<T: FromStr>(file: &str) -> Vec<T> {
    let rdr = BufReader::new(File::open(file).unwrap());
    let rows = rdr.lines().map(|l| l.ok().and_then(|s| s.parse::<T>().ok()).unwrap()).collect();
    return rows
}

pub fn load_txt_gz<T: FromStr>(file: &str) -> Vec<T> {
    let gz = GzDecoder::new(File::open(file).unwrap()).unwrap();
    let rdr = BufReader::new(gz);
    let rows = rdr.lines().map(|l| l.ok().and_then(|s| s.parse::<T>().ok()).unwrap()).collect();
    return rows
}

pub fn load_txt_maybe_gz<T: FromStr>(file: &str) -> Vec<T> {
    match file.ends_with(".gz") {
        true => load_txt_gz(file),
        false => load_txt(file),
    }
}

pub fn save_txt(vec: Vec<&String>, file: &str) {
     let mut writer = File::create(file).unwrap();
     for entry in vec {
         writer.write_all(format!("{}\n", entry).as_bytes()).unwrap();
     }
}

pub fn load_tabular<T: DeserializeOwned>(file: &str, skip_row: bool) -> Vec<T> {
    // TODO: serde-based CSV parsing appears to be slower than rustc_serialize; figure out how to optimize this
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).flexible(true).from_path(Path::new(file)).unwrap();
    let mut iter = rdr.deserialize::<T>();
    if skip_row { let _ = iter.next(); }
    let rows = iter.collect::<csv::Result<Vec<T>>>().unwrap();
    return rows
}

pub fn read_json(path: &str) -> Value {
    let mut file = File::open(path).unwrap();
    let mut data = String::new();
    file.read_to_string(&mut data).unwrap();
    let json = serde_json::from_str(&data).unwrap();
    return json
}

pub fn write_json(data: Value, path: &str) {
    let writer = File::create(path).unwrap();
    serde_json::to_writer_pretty(writer, &data).expect("Failed to serialize JSON");
}

pub fn get_reference_genomes(reference_path: &str) -> Vec<String> {
    let metadata_path = Path::new(reference_path).join("reference.json").as_path().to_str().unwrap().to_owned();
    let metadata = read_json(&metadata_path);
    let genomes = metadata.as_object().unwrap().get("genomes").unwrap().as_array().unwrap().to_owned();
    return genomes.iter().map(|x| x.as_str().unwrap().to_string()).collect()
}

/* Miscellaneous */

#[derive(Eq, PartialEq, Debug, Clone, Copy)]
pub enum BisectDirection {
    Left,
    Right,
}

pub fn bisect<T,U,F>(foo: &[T], bar: U, key_func: &F, direction: BisectDirection) -> usize where U: Ord, F: Fn(&T) -> U {
    // NOTE: assumes that input is sorted by key_func
    match foo.binary_search_by_key(&bar, key_func) {
        Ok(init_idx) => { // found an exact match -> seek to left / right if necessary
            let mut idx = init_idx as i32;
            let inc = match direction {
                BisectDirection::Right => 1,
                BisectDirection::Left => -1,
            };
            while (idx >= 0) && (idx < foo.len() as i32) {
                if key_func(&foo[idx as usize]) != bar { break; }
                idx += inc;
            };
            if direction == BisectDirection::Left { idx += 1; }
            idx as usize
        },
        Err(idx) => idx, // no exact match -> this is the only correct insertion point
    }
}

pub fn median(vec: &[i64]) -> i64 {
    // NOTE: assumes input is sorted
    let n = vec.len();
    if n == 0 {
        return 0
    } else if n % 2 == 0 {
        let i = n / 2;
        return (vec[i-1] + vec[i]) / 2
    } else {
        return vec[n/2]
    }
}

/* FASTQ */

// TODO: this was copied from vdj_asm, we should move it into a shared library
pub struct CellrangerFastqHeader {
    pub qname: String,
    pub header: String,
    pub tags: HashMap<String, String>,
}

impl CellrangerFastqHeader {
    pub fn new(header: String) -> CellrangerFastqHeader {
        let parts = header.clone().split("|||").map(|x| x.to_string()).collect::<Vec<String>>();

        let mut tags = HashMap::new();

        if parts.len() > 1 {
            for idx in 1..parts.len() {
                if idx % 2 == 0 { continue; }
                if parts[idx + 1].len() == 0 { continue; } // skip empty tags
                tags.insert(parts[idx].to_string(), parts[idx + 1].to_string());
            }
        }

        CellrangerFastqHeader {
            qname: parts[0].to_string(),
            header: header,
            tags: tags,
        }
    }
}

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
    for c in read.cigar().into_iter() {
        match *c {
            Cigar::Match(l) | Cigar::Del(l) | Cigar::RefSkip(l) | Cigar::Equal(l) | Cigar::Diff(l) => alen += l as i64,
            Cigar::Back(l) => alen -= l as i64,
            _ => ()
        }
    }
    return alen
}

/// Produce a processed BC sequence.
pub fn get_processed_bc(corrected_seq: &str, gem_group: &u8) -> String {
    format!("{}-{}", corrected_seq, gem_group)
}

/// Find intersection of two sorted, deduped vectors
pub fn intersect_vecs<T: PartialOrd + Eq + Clone>(x: &Vec<T>, y: &Vec<T>)
                                                  -> Vec<T> {
    let mut i = 0;
    let mut j = 0;
    let mut result = Vec::new();

    while i < x.len() && j < y.len() {
        if x[i] < y[j] {
            i += 1;
        } else if y[j] < x[i] {
            j += 1
        } else {
            result.push(x[i].clone());
            i += 1;
            j += 1;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisect() {
        let v = vec![-1, -1, 0, 2, 4];
        let id = |x: &i32| *x;
        assert_eq!(bisect(&v, -2, &id, BisectDirection::Left), 0);
        assert_eq!(bisect(&v, -2, &id, BisectDirection::Right), 0);
        assert_eq!(bisect(&v, -1, &id, BisectDirection::Left), 0);
        assert_eq!(bisect(&v, -1, &id, BisectDirection::Right), 2);
        assert_eq!(bisect(&v, 0, &id, BisectDirection::Left), 2);
        assert_eq!(bisect(&v, 0, &id, BisectDirection::Right), 3);
        assert_eq!(bisect(&v, 3, &id, BisectDirection::Left), 4);
        assert_eq!(bisect(&v, 3, &id, BisectDirection::Right), 4);
        assert_eq!(bisect(&v, 5, &id, BisectDirection::Left), 5);
        assert_eq!(bisect(&v, 5, &id, BisectDirection::Right), 5);
        // NOTE: input is assumed to be sorted by key_func, so only apply functions that maintain sort
        let minus = |x: &i32| (*x)-1;
        assert_eq!(bisect(&v, 0, &minus, BisectDirection::Left), 3);
        assert_eq!(bisect(&v, 0, &minus, BisectDirection::Right), 3);

    }

    #[test]
    fn test_median() {
        assert_eq!(median(&vec![]), 0);
        assert_eq!(median(&vec![1]), 1);
        assert_eq!(median(&vec![1,2]), 1);
        assert_eq!(median(&vec![2,2]), 2);
        assert_eq!(median(&vec![1,2,3]), 2);
    }

    #[test]
    fn test_intersect() {
        assert_eq!(intersect_vecs(&vec![0u64,1], &vec![2,3]), Vec::<u64>::new());
        assert_eq!(intersect_vecs(&vec![0,2,4], &vec![1,2,5]), vec![2]);
    }
}
