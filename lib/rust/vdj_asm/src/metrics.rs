//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

extern crate rustc_serialize;
use rustc_serialize::{json, Encodable};
use std::fs::File;
use std::io::{Write};
use std::io;
use std::collections::HashMap;

#[derive(RustcEncodable, RustcDecodable, Default)]
pub struct AssemblyMetrics {
    pub assemblable_read_pairs_by_bc: HashMap<String, u64>,
}

pub fn write_summary<T: Encodable>(file: &mut File, metrics: &T) -> io::Result<()> {
    writeln!(file, "{}", json::as_pretty_json(metrics))
}
