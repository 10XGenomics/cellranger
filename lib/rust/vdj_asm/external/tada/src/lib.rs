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
extern crate crossbeam;
extern crate rayon;

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
pub mod phased_graph;
pub mod ec;
pub mod bc_sets;


pub mod bitenc;
mod bwt;

pub mod shard_buf;
pub mod load_graph;
pub mod external;
pub mod multifastq;
pub mod utils;
pub mod sim_tests;
