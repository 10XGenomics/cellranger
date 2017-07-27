//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! # tada, assembler components in Rust
//! tada implements a front-end DeBruijn graph construction using the MSP strategy.

#![crate_name = "tada"]
#![allow(dead_code)]

extern crate bincode;
extern crate rustc_serialize;
extern crate bio;
extern crate flate2;
extern crate docopt;
extern crate rand;
extern crate chrono;
extern crate time;
extern crate itertools;
extern crate linked_hash_map;
extern crate csv;
extern crate martian;
extern crate rayon;
extern crate crossbeam;

#[macro_use]
extern crate log;
extern crate env_logger;

pub mod cmd_msp;
pub mod cmd_shard_asm;
pub mod cmd_main_asm;
pub mod cmd_graph_stats;

pub mod kmer;
pub mod msp;
pub mod debruijn;

pub mod bitenc;
mod bwt;

pub mod shard_buf;
pub mod load_graph;
pub mod external;
pub mod multifastq;
pub mod utils;
pub mod sim_tests;


use docopt::Docopt;
use std::path::{Path, PathBuf};
use std::env;
use log::{LogRecord, LogLevelFilter};
use env_logger::LogBuilder;
use std::collections::HashMap;
use chrono::*;
use martian::MartianStage;

const USAGE: &'static str = "
Top-down assembler

Usage:
  tada exp <graph> [<edge-db>]
  tada msp <whitelist> <permutation> <output> <files>...
  tada shard-asm <output> <files>...
  tada martian <adapter>...
  tada (-h | --help)
  tada --version

Options:
  -h --help     Show this screen.
  --version     Show version.
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_whitelist: String,
    arg_permutation: String,
    arg_output: String,
    arg_files: Vec<String>,
    arg_edge_db: Option<String>,

    arg_graph: String,
    //arg_sequence: String,

    cmd_exp: bool,
    cmd_msp: bool,
    cmd_shard_asm: bool,

    // Martian interface mode
    cmd_martian: bool,
    arg_adapter: Vec<String>,
}

fn main() {


    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.decode())
                         .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);
    info!("{:?}", args);

    // Martian does it's own logging -- don't configure if we're being called by martian
    if !args.cmd_martian {

        // setup logging
        let format = |record: &LogRecord| {
            let time = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
            format!("[{}] {} - {}", time, record.level(), record.args())
        };

        let mut builder = LogBuilder::new();
        builder.format(format).filter(None, LogLevelFilter::Info);

        if env::var("RUST_LOG").is_ok() {
           builder.parse(&env::var("RUST_LOG").unwrap());
        }

        builder.init().unwrap();
    }

    if args.cmd_msp {
        let mut paths: Vec<PathBuf> = Vec::new();
        for p in args.arg_files {
            paths.push(PathBuf::from(p));
        }

        cmd_msp::main_msp(
            Path::new(&args.arg_whitelist), paths,
            Path::new(&args.arg_permutation),
            Path::new(&args.arg_output));

    } else if args.cmd_shard_asm {

        let mut paths: Vec<PathBuf> = Vec::new();
        for p in args.arg_files {
            paths.push(PathBuf::from(p));
        }

        cmd_shard_asm::main_shard_asm(paths, Path::new(&args.arg_output), Path::new("FIXME"));
    } else if args.cmd_martian {

        let mut stage_registry : HashMap<String, Box<MartianStage>> = HashMap::new();
        stage_registry.insert("msp".to_string(), Box::new(cmd_msp::MspMartian));
        stage_registry.insert("shard-asm".to_string(), Box::new(cmd_shard_asm::ShardAsmMartian));
        stage_registry.insert("main-asm".to_string(), Box::new(cmd_main_asm::MainAsmMartian));

        // Run the built-in martian adapter
        martian::martian_main(args.arg_adapter, stage_registry);
    } else if args.cmd_exp {
        debruijn::do_explore(args.arg_graph, args.arg_edge_db);
    }
}
