extern crate flate2;

use std::path::Path;
use std::fs::{File};
use std::io::{Read, BufRead, BufReader, Lines};
use std::boxed::Box;
use flate2::read::MultiGzDecoder;

/// Head, Seq, Qual from single FASTQ
pub type FqRec = (Vec<u8>, Vec<u8>, Vec<u8>);

/// R1, R2, optional SI
pub type RawReadSet = (FqRec, FqRec, Option<FqRec>);

/// Read a single FqRec from a line 
pub fn get_fastq_item<R: BufRead>(lines: &mut Lines<R>) -> Option<FqRec>{
    match lines.next() {
        Some(head) => {
            let r1 = lines.next().unwrap().unwrap().into_bytes();
            let _  = lines.next().unwrap().unwrap().into_bytes();
            let q1 = lines.next().unwrap().unwrap().into_bytes();

            // trim after first space of head
            //let head_str = head.unwrap();
            //let mut split = head_str.split_whitespace();
            //let trim_head = split.next().unwrap().to_string().into_bytes();

            Some((head.unwrap().into_bytes(), r1, q1))
        },
        None => None
    }
}

pub struct FastqPairIter {
    r1: Lines<BufReader<Box<Read>>>,
    r2: Lines<BufReader<Box<Read>>>,
    si: Option<Lines<BufReader<Box<Read>>>>,
}

pub fn open_w_gz<P: AsRef<Path>>(p: P) -> Box<Read> {
    let r = File::open(p.as_ref()).unwrap();

    if p.as_ref().extension().unwrap() == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(32*1024, gz);
        Box::new(buf_reader)
    } else {
        let buf_reader = BufReader::with_capacity(32*1024, r);
        Box::new(buf_reader)
    }
}


pub fn open_fastq_pair_iter<P: AsRef<Path>>(r1: P, r2: P, si: Option<P>) -> Box<Iterator<Item=RawReadSet>> {
    Box::new(
    FastqPairIter::init(
        open_w_gz(r1),
        open_w_gz(r2),
        si.map(|si| open_w_gz(si)),
    ))
}

impl FastqPairIter {
    pub fn init(r1: Box<Read>, r2: Box<Read>, si: Option<Box<Read>>) -> FastqPairIter {

        FastqPairIter {
            r1: BufReader::new(r1).lines(),
            r2: BufReader::new(r2).lines(),
            si: si.map(|x| BufReader::new(x).lines()),
        }
    }
}

impl Iterator for FastqPairIter {
    type Item = RawReadSet;

    fn next(&mut self) -> Option<RawReadSet> {

        match get_fastq_item(&mut self.r1) {
            Some(f_r1) => {
                let f_r2 = get_fastq_item(&mut self.r2).unwrap();
                
                let f_si = self.si.as_mut().map(|_si| get_fastq_item(_si).unwrap());
                Some((f_r1, f_r2, f_si))
            }
            None => None,
        }
    }
}

pub struct InterleavedFastqPairIter {
    ra: Lines<BufReader<Box<Read>>>,
    si: Option<Lines<BufReader<Box<Read>>>>,
}

pub fn open_interleaved_fastq_pair_iter<P: AsRef<Path>>(ra: P, si: Option<P>) -> Box<Iterator<Item=RawReadSet>> {
    Box::new(
    InterleavedFastqPairIter::init(
        open_w_gz(ra),
        si.map(|si| open_w_gz(si)),
    ))
}

impl InterleavedFastqPairIter {

    pub fn init(ra: Box<Read>, si: Option<Box<Read>>) -> InterleavedFastqPairIter {
        InterleavedFastqPairIter {
            ra: BufReader::new(ra).lines(),
            si: si.map(|x| BufReader::new(x).lines()),
        }
    }
}

impl Iterator for InterleavedFastqPairIter {
    type Item = RawReadSet;

    fn next(&mut self) -> Option<RawReadSet> {

        match get_fastq_item(&mut self.ra) {
            Some(f_r1) => {
                let f_r2 = get_fastq_item(&mut self.ra).unwrap();
                let f_si = self.si.as_mut().map(|_si| get_fastq_item(_si).unwrap());
                Some((f_r1, f_r2, f_si))
            }
            None => None,
        }
    }
}