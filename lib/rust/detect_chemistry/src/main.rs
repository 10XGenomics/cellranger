//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::alloc::System;

#[global_allocator]
static A: System = System;

extern crate docopt;
#[macro_use] extern crate serde;
extern crate regex;
#[macro_use] extern crate lazy_static;
extern crate bio;
extern crate debruijn;
extern crate fnv;
extern crate itertools;
extern crate bincode;
extern crate serde_json;
extern crate boomphf;
extern crate pdqsort;
#[macro_use] extern crate log;
extern crate env_logger;
extern crate chrono;
extern crate rand;
extern crate flate2;

use docopt::Docopt;
use std::path::Path;
use serde::{Serialize};
use std::fs;
use std::fs::File;
use std::io::{BufWriter, BufReader, Read};
use std::hash::Hash;

use bio::io::{fasta,fastq};

use debruijn::{Kmer};
use debruijn::kmer::{IntKmer};

mod mylog;
mod transcript;
mod index;
mod map;
mod tests;
mod utils;

use index::{KmerIndexType,load_index,BBHashKmerIndex,HashSetKmerIndex,KmerPresenceQuery};

const USAGE: &'static str = "
Usage:
  detect_chemistry get-transcripts <genome-fa> <in-gtf> <out-fa> [--all-transcripts]
  detect_chemistry index-transcripts <in-fa> <out-idx> [--skip=N] [--hash=XXX] [--gamma=F]
  detect_chemistry map-reads <in-idx> <in-fq> <out-json> [--skip=N] [--min-kmers=N] [--initial-reads=N] [--fastqs=F]
  detect_chemistry (-h | --help)

Options:
  -h --help            Show this screen.
  --all-transcripts    Use all transcripts (not just 1st encountered per gene)
  --skip=N             Skip bases between kmers [default: 0]
  --hash=XXX           Hash method (bbhash|hashset) [default: bbhash]
  --gamma=F            BBHash gamma parameter for time/space tradeoff [default: 2.0]
  --min-kmers=N        Min matching kmers to consider a read mapped [default: 1]
  --initial-reads=N    Process initial N reads
  --fastqs=F           Path to JSON file w/ info for each input fastq (ignores <in-fq>)
";

#[derive(Debug, Deserialize, Clone)]
struct Args {
    cmd_get_transcripts:    bool,
    cmd_index_transcripts:  bool,
    cmd_map_reads:  bool,

    // get_transcripts
    arg_genome_fa:       Option<String>,
    arg_in_gtf:          Option<String>,
    arg_out_fa:          Option<String>,
    flag_all_transcripts:bool,

    // index_transcripts
    arg_in_fa:           Option<String>,
    arg_out_idx:         Option<String>,
    flag_skip:           usize,
    flag_hash:           Option<String>,
    flag_gamma:          f64,

    // map
    arg_in_idx:          Option<String>,
    arg_in_fq:           Option<String>,
    arg_out_json:        Option<String>,
    flag_min_kmers:      usize,
    flag_initial_reads:  Option<usize>,
    flag_fastqs:         Option<String>,
}

type MyKmer = IntKmer<u64>;

fn main() {
    mylog::init_log();
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    if args.cmd_get_transcripts {
        get_transcripts(args);
    } else if args.cmd_index_transcripts {
        match args.clone().flag_hash.unwrap().as_str() {
            "bbhash" => index_transcripts_mphf::<MyKmer>(args),
            "hashset" => index_transcripts_hashset::<MyKmer>(args),
            _ => panic!("Unexpected value for --hash"),
        }
    } else if args.cmd_map_reads {
        map_reads(args);
    }
    info!("Done");
}


fn get_transcripts(args: Args) {
    info!("Parsing GTF");
    let gtf_reader = BufReader::new(File::open(&args.arg_in_gtf.expect("No GTF file specified"))
                                .expect("Failed to open GTF file"));
    let transcripts = transcript::parse_gtf(gtf_reader, !args.flag_all_transcripts);

    info!("Fetching sequences");
    let mut fasta = fasta::IndexedReader::from_file(
        &Path::new(&args.arg_genome_fa.expect("No genome FASTA file specified").to_owned()))
        .expect("Failed to open genome FASTA file");

    let mut out_fasta = fasta::Writer::to_file(
        &Path::new(&args.arg_out_fa.expect("No output FASTA file specified").to_owned()))
        .expect("Failed to open output FASTA file");

    for (tx_id, transcript) in transcripts {
        let seq = transcript::fetch_transcript(&mut fasta, transcript);

        out_fasta.write(&tx_id, None, &seq)
            .expect("Failed to write sequence to output FASTA file");
    }
}

fn index_transcripts_hashset<K: Kmer + Serialize + Hash + Eq>(args: Args) {
    let fa_reader = fasta::Reader::from_file(&Path::new(&args.arg_in_fa.unwrap()))
        .expect("Failed to load transcript FASTA file");

    let idx: HashSetKmerIndex<K> = index::index_transcripts_hashset(fa_reader,
                                                                           args.flag_skip);

    info!("Writing hashset");
    let mut writer = BufWriter::new(File::create(&Path::new(&args.arg_out_idx.unwrap()))
                                    .expect("Failed to create index file"));

    bincode::serialize_into(&mut writer, &idx)
        .expect("Failed to serialize index");
}

fn index_transcripts_mphf<K: Kmer + Serialize>(args: Args) {
    let in_fa_str = &args.arg_in_fa.unwrap();
    let fa_path = Path::new(in_fa_str);

    // Get length so we can avoid allocating extra mem when storing kmers
    // Approximate number of kmers as the filesize in bytes
    let metadata = fs::metadata(&fa_path)
        .expect("Unable to get metadata of FASTA file");
    let fa_size = metadata.len();

    let fa_reader = fasta::Reader::from_file(&fa_path)
        .expect("Failed to load FASTA file");

    let idx: BBHashKmerIndex<K> = index::index_transcripts_mphf(fa_reader, fa_size as usize,
                                                                       args.flag_skip,
                                                                       args.flag_gamma);

    let mut writer = BufWriter::new(File::create(&Path::new(&args.arg_out_idx.unwrap()))
                                    .expect("Failed to create index file"));

    info!("Writing index");
    bincode::serialize_into(&mut writer, &idx)
        .expect("Failed to serialize index");
}

#[derive(Serialize, Deserialize)]
struct FastqDef {
    read_type: String,
    input: Vec<String>,
    interleaved: bool,
}
type FastqDefList = Vec<FastqDef>;

fn map_reads(args: Args) {
    // Get the kmer index type
    let idx_type: KmerIndexType;
    {
        let mut reader = BufReader::new(File::open(&Path::new(&args.clone().arg_in_idx.unwrap()))
                                        .expect("Failed to open index file for detecting type"));
        idx_type = index::get_index_type(&mut reader);
    }
    let mut idx_reader = BufReader::new(File::open(&Path::new(&args.clone().arg_in_idx.unwrap()))
                                    .expect("Failed to open index file for deserialization"));

    println!("Loading index");
    // Use Box here to make the match arms consistent
    let idx = match idx_type {
        KmerIndexType::BBHashIndex => {
            Box::new(load_index::<BBHashKmerIndex<MyKmer>, BufReader<File>>(&mut idx_reader)) as Box<KmerPresenceQuery<MyKmer>>
        },
        KmerIndexType::HashSetIndex => {
            Box::new(load_index::<HashSetKmerIndex<MyKmer>, BufReader<File>>(&mut idx_reader)) as Box<KmerPresenceQuery<MyKmer>>
        },
    };

    let min_kmers = args.flag_min_kmers;
    let skip_bases = args.flag_skip;
    let initial_reads = args.flag_initial_reads.unwrap_or(std::usize::MAX);
    let flag_fastqs = args.flag_fastqs;
    let arg_in_fq = args.arg_in_fq.clone();
    let arg_out_json = args.arg_out_json;

    // Prep fastq definitions
    let mut fq_defs = Vec::new();
    if flag_fastqs.is_some() {
        let mut f = File::open(flag_fastqs.unwrap()).expect("Failed to open fastq-defining JSON file for reading");
        let mut buf = String::new();
        f.read_to_string(&mut buf).expect("Failed to read fastq-defining JSON file");
        let fq_def_list: FastqDefList = serde_json::from_str(&buf).expect("Failed to parse fastq-defining JSON file");
        for fq_def in fq_def_list {
            fq_defs.push(fq_def)
        }
    } else {
        fq_defs.push(FastqDef{read_type: String::from("R1"), input: vec![arg_in_fq.unwrap()], interleaved: false});
    }

    let mut metrics: map::Metrics = Default::default();
    for fq_def in &fq_defs {
        for fq_filename in &fq_def.input {
            let reader = utils::open_maybe_gzip(&fq_filename);
            let fq_reader = fastq::Reader::new(reader);
            let fq_metrics = map::map_reads(&*idx, fq_reader, min_kmers, skip_bases, initial_reads, fq_def.interleaved, Some(fq_def.read_type.clone()));
            map::add_metrics(&mut metrics, &fq_metrics);
        }
    }
    serde_json::to_writer_pretty(File::create(&Path::new(&arg_out_json.unwrap()))
                                 .expect("Couldn't create output JSON file"), &metrics)
        .expect("Failed to write JSON output");
}
