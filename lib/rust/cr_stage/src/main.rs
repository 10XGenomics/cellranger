//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

extern crate log;
extern crate env_logger;
extern crate serde;
#[macro_use] extern crate serde_derive;
#[macro_use] extern crate serde_json;
extern crate bincode;
extern crate rust_htslib;
extern crate martian;
extern crate itertools;
extern crate byteorder;
extern crate tenkit;
extern crate docopt;
extern crate chrono;
#[macro_use] extern crate bitflags;
extern crate failure;

use std::collections::{HashMap};
use std::io::Write;
use std::env;
use std::str;
use docopt::Docopt;
use itertools::Itertools;
use rust_htslib::bam;
use log::LevelFilter;
use env_logger::Builder;
use chrono::Local;
use martian::*;

mod cmd_mark_dups;
mod utils;

const USAGE: &'static str = "
Cell Ranger Rust stages
Usage:
  cr_stage martian <adapter>...
  cr_stage (-h | --help)
Options:
  -h --help            Show this screen.
";

#[derive(Debug, Deserialize)]
struct Args {
    // Martian interface
    cmd_martian: bool,
    arg_adapter: Vec<String>,
}

fn main() {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);

    // Martian does its own logging -- don't configure if we're being called by martian
    if !args.cmd_martian {
        // Setup logging
        Builder::new()
            .format(|buf, record| {
                write!(buf, "{} [{}] - {}",
                        Local::now().format("%Y-%m-%dT%H:%M:%S"),
                        record.level(),
                        record.args())
            })
            .filter(None, LevelFilter::Info)
            .init();

    } else if args.cmd_martian {

        let mut stage_registry : HashMap<String, Box<MartianStage>> = HashMap::new();
        stage_registry.insert("mark_duplicates".to_string(), Box::new(cmd_mark_dups::MarkDuplicatesStage));

        println!("{}", env::args().join(" "));

        // Run the built-in martian adapter
        martian::martian_main(args.arg_adapter, stage_registry);
    }
}
