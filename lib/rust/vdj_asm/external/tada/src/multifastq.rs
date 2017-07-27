//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Reader for the FASTQ-ish format produced by SORT_FASTQ_BY_BARCODE pipeline

use std::fs::File;
use std::path::Path;
use std::io::{BufRead, BufReader, Lines};
use std::collections::HashMap;
use flate2::read::GzDecoder;

/// Single record from a FASTQ file
#[derive(Debug)]
pub struct Fastq {
    pub read: String,
    pub qual: String,
}

pub struct InputRead {
    pub header: String,
    pub r1: Fastq,
    pub r2: Fastq,
    pub bc_id: u32,
    pub read_id: u16,
}

pub struct MultiFastqIter<'a> {
    lines: Lines<BufReader<GzDecoder<File>>>,
    current_bc: Option<u32>,
    bc_map: &'a HashMap<String, u32>,
    default_bc: u32,
    bc_count: usize,
}

impl<'a> MultiFastqIter<'a> {
    /// Iterate reads from an (optionally compressed) FASTQ file
    pub fn new(name: &Path,
               bc_map: &'a HashMap<String, u32>,
               default_bc: u32)
               -> MultiFastqIter<'a> {

        println!("{:?}", name);
        let f = File::open(name).unwrap();

        match name.extension() {
            None => panic!("Not a gz file"),
            Some(ext) => {
                if ext == "gz" {
                    let gz = GzDecoder::new(f).unwrap();
                    let br = BufReader::new(gz);
                    let lines = br.lines();
                    MultiFastqIter {
                        lines: lines,
                        current_bc: None,
                        bc_count: 0,
                        bc_map: bc_map,
                        default_bc: default_bc,
                    }
                } else {
                    panic!("Not a gz file")
                }
            }
        }
    }
}

impl<'a> Iterator for MultiFastqIter<'a> {
    type Item = InputRead;

    fn next(&mut self) -> Option<InputRead> {
        let ref mut lines = self.lines;

        match lines.next() {
            Some(head) => {

                let r1 = lines.next().unwrap().unwrap();
                let q1 = lines.next().unwrap().unwrap();
                let read1 = Fastq {
                    read: r1,
                    qual: q1,
                };

                let r2 = lines.next().unwrap().unwrap();
                let q2 = lines.next().unwrap().unwrap();
                let read2 = Fastq {
                    read: r2,
                    qual: q2,
                };

                let bc = lines.next().unwrap().unwrap();
                let _ = lines.next().unwrap().unwrap();

                let _ = lines.next().unwrap().unwrap();
                let _ = lines.next().unwrap().unwrap();


                // Process barcode
                let bc_id =
                    match bc.find(",") {
                        Some(_) => {
                            let seq = bc.split(",").nth(0).unwrap();
                            *self.bc_map.get(seq).unwrap_or(&self.default_bc)
                        },
                        None => 0,
                    };

                if Some(bc_id) == self.current_bc && bc_id > 0 {
                    self.bc_count = self.bc_count + 1;
                } else {
                    self.current_bc = Some(bc_id);
                    self.bc_count = 0;
                }

                Some(InputRead {
                    header: head.unwrap(),
                    r1: read1,
                    r2: read2,
                    bc_id: bc_id,
                    read_id: self.bc_count as u16,
                })
            }
            None => None,
        }
    }
}
