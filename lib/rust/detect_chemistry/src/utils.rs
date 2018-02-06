//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::path::Path;
use std::io::{Read, BufReader};
use std::fs::File;
use std::ffi::OsStr;

use flate2::read::{MultiGzDecoder};

pub fn open_maybe_gzip<P: AsRef<Path>> (path: &P) -> Box<Read> {
    let file = File::open(path).expect("Failed to open file");
    let ext = path.as_ref().extension();
    if ext.unwrap_or(OsStr::new("")) == "gz" {
        let gz = MultiGzDecoder::new(file).expect("Failed to open gz stream");
        Box::new(BufReader::with_capacity(32*1024, gz))
    } else {
        Box::new(BufReader::with_capacity(32*1024, file))
    }
}
