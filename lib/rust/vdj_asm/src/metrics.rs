//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use serde::{Serialize};
use serde_json;
use std::fs::File;
use std::io::{Write};
use std::io;
use std::collections::HashMap;

#[derive(Serialize, Deserialize, Default)]
pub struct AssemblyMetrics {
    pub assemblable_read_pairs_by_bc: HashMap<String, u64>,
}

pub fn write_summary<T: Serialize>(file: &mut File, metrics: &T) -> io::Result<()> {
    writeln!(file, "{}", serde_json::to_string(metrics)?)
}
