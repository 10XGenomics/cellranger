//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Command to generate barcoded MSP substrings from FASTQ data
//!
//! For each input read, compute the MSP substrings. Each substring is at most 2K-P long, and
//! is stored in a fixed-size struct Bsp, along with the barcode id, read id, and read position,
//! and an Exts bitfield indicating the neighboring bases (if any).

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::collections::{BTreeMap, HashMap};
use std::env;

use utils::load_barcode_whitelist;
use shard_buf::{Shard, ShardProducer};
use multifastq::{MultiFastqIter, InputRead};
use msp::simple_scan;
use martian::{JsonDict, MartianStage, obj_decode, obj_encode, set_file_handle_limit};
use rustc_serialize::json::{Json, ToJson};
use crossbeam;

use kmer::{self, Pmer, Bsp, base_to_bits};
use utils;


const K: usize = kmer::K;
const P: usize = 8;

#[inline(never)]
pub fn main_msp(bc_wl: &Path, fastq: Vec<PathBuf>, permutation: &Path, out_path: &Path) -> Vec<PathBuf> {
    let shards = 1024;
    set_file_handle_limit(shards + 100);

    let bc_wl = &load_barcode_whitelist(bc_wl);
    let permutation : &Vec<usize> = &utils::read_obj(permutation).expect("couldn't read permutation");

    let s = Shard::new(out_path, shards, 2);

    crossbeam::scope(|scope| {
        for f in fastq.chunks(2)
        {
            let chnk = Vec::from(f);
            let producer = ShardProducer::new(&s);

            scope.spawn(move ||
            {
                let fq_src = chnk.iter().flat_map(|p| MultiFastqIter::new(&p, &bc_wl, 0));
                process_fastq_set(K, P, &permutation, fq_src, producer);
            });
        }
    });

    s.finish()
}

pub fn msp_get_pmer_permutation(fastq: Vec<PathBuf>, bc_wl: &Path, permutation_path: &Path) {
    let bc_wl = load_barcode_whitelist(bc_wl);
    let fq_src_small = fastq.iter().flat_map(|p| MultiFastqIter::new(&p, &bc_wl, 0).take(25000));
    let pmer_counts = count_pmers_fastq(fq_src_small, P);
    let permutation = pmer_inverse_freq_permutation(P, &pmer_counts);

    info!("Writing pmer permutation: {:?}", permutation_path);
    utils::write_obj(&permutation, permutation_path).expect("write failed");
}



#[inline(never)]
fn msp_read(k: usize,
            p: usize,
            partition: u32,
            permutation: &Vec<usize>,
            read: u16,
            seq: &[u8],
            shard_writer: &mut ShardProducer<Bsp>) {
    let msp_parts = simple_scan(k, p, seq, permutation, true);

    for (msp_key, _, start_pos, slice_length) in msp_parts {
        assert!(slice_length >= k);
        let b = Bsp::new(seq, start_pos, slice_length, partition, read);


        if (b.len() as usize) < K || (b.len() as usize) != b.sequence.len()
        {
            println!("bad bsp: {:?}", b);
            panic!("bad bsp: {:?}", b);
        }

        shard_writer.send(msp_key as usize, b);
    }
}


fn process_fastq_set<T>(k: usize,
                        p: usize,
                        permutation: &Vec<usize>,
                        input: T,
                        mut shard_writer: ShardProducer<Bsp>)
    where T: IntoIterator<Item = InputRead>
{
    let mut iter = input.into_iter();
    loop {
        match iter.next() {
            Some(inp) => {
                // R1 is numbered with an even number, R2 is the next value
                let r1b: Vec<u8> = inp.r1.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r1_slice = &r1b[..];
                msp_read(k,
                         p,
                         inp.bc_id,
                         permutation,
                         inp.read_id,
                         r1_slice,
                         &mut shard_writer);

                let r2b: Vec<u8> = inp.r2.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r2_slice = &r2b[..];
                msp_read(k,
                         p,
                         inp.bc_id,
                         permutation,
                         inp.read_id + 1,
                         r2_slice,
                         &mut shard_writer);
            }
            None => break,
        }
    }
}

fn pmer_inverse_freq_permutation(p: usize, counts: &HashMap<Pmer, usize>) -> Vec<usize> {
    let mut freq_vec: Vec<(usize, Pmer)> = Vec::new();

    for pmer in Pmer::enumerate(p) {
        match counts.get(&pmer) {
            Some(count) => freq_vec.push((*count, pmer)),
            None => freq_vec.push((0, pmer)),
        }
    }

    freq_vec.sort();

    let mut perm: Vec<(usize, usize)> = freq_vec.into_iter().enumerate().map(|(idx, (_, pmer))| (pmer.value(), idx)).collect();
    perm.sort();
    perm.into_iter().map(|(_, idx)| idx).collect()
    /*
    for (idx, (count, pmer)) in Iterator::enumerate(freq_vec.iter().cloned())
    {
        println!("pmer: {:?}, count: {:?}, idx: {:?}", pmer, count, idx);
    }
    freq_vec.iter().map(|&(_, pmer)| pmer.value()).collect()
    */
}


fn count_pmers_fastq<T: Iterator<Item = InputRead>>(mut input: T,
                                                    p: usize)
                                                    -> HashMap<Pmer, usize> {
    let mut counter = HashMap::new();
    loop {
        match input.next() {
            Some(inp) => {

                // R1 is numbered with an even number, R2 is the next value
                let r1b: Vec<u8> = inp.r1.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r1_slice = &r1b[..];
                let pmers = Pmer::pmers_from_string(r1_slice, p);

                for pmer in pmers {
                    let (min_pmer, _) = pmer.min_rc();
                    let e = counter.entry(min_pmer).or_insert_with(|| 0);
                    *e = *e + 1;
                }

                let r2b: Vec<u8> = inp.r2.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r2_slice = &r2b[..];
                let pmers = Pmer::pmers_from_string(r2_slice, p);

                for pmer in pmers {
                    let (min_pmer, _) = pmer.min_rc();
                    let e = counter.entry(min_pmer).or_insert_with(|| 0);
                    *e = *e + 1;
                }

            }
            None => break,
        }
    }

    counter
}

pub struct MspMartian;


impl MartianStage for MspMartian {
    fn split(&self, args: JsonDict) -> JsonDict {

        let fastqs : Vec<String> = args["fastqs"].as_array().unwrap().into_iter().map(|x| x.as_string().unwrap().to_string()).collect();
        //let fastqs_fofn = args["fastqs"].as_string().unwrap();
        //let fastqs = utils::read_fofn(Path::new(fastqs_fofn));

        let _whitelist = args["barcode_whitelist"].as_string().unwrap();
        let whitelist  = Path::new(_whitelist);

        let mut permutation_file = env::current_dir().unwrap();
        permutation_file.push("permutation.perm");

        // Compute the permutation
        let paths = fastqs.clone().iter().map(PathBuf::from).collect();
        msp_get_pmer_permutation(paths, &whitelist, &permutation_file);

        let mut chunks : Vec<Json> = Vec::new();
        for in_files in fastqs.chunks(4)
        {
            let fastqs = in_files.iter().map(|c| Json::String(c.clone())).collect();

            let mut bt = BTreeMap::new();
            bt.insert("chunk".to_string(), Json::Array(fastqs));
            let p2s = permutation_file.clone().into_os_string().into_string().unwrap();
            bt.insert("permutation".to_string(), Json::String(p2s));
            bt.insert("__mem_gb".to_string(), Json::F64(4.0));
            bt.insert("__threads".to_string(), Json::I64(8));
            chunks.push(Json::Object(bt));
        }


        let chunk_array = Json::Array(chunks);
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), chunk_array);
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {
        let chunk = args["chunk"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();
        //let chunk  = vec![PathBuf::from(_chunk)];

        let _whitelist = args["barcode_whitelist"].as_string().unwrap();
        let whitelist  = Path::new(_whitelist);

        let _permutation = args["permutation"].as_string().unwrap();
        let permutation = Path::new(_permutation);

        let _out_fn = {outs.get("chunks").unwrap().as_string().unwrap().clone()};
        let mut out_fn = PathBuf::from(_out_fn);
        out_fn.pop();
        out_fn.push("shards");
        let out_dir = out_fn.as_path();

        let files = main_msp(&whitelist, chunk, &permutation, out_dir);
        let mut final_outs = outs.clone();
        final_outs.insert("msp".to_string(), Json::Array(files.iter().map(|x| Json::String(x.to_str().unwrap().to_string())).collect()));
        final_outs
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

            let mut msps_vec: Vec<MspChunk> = Vec::new();

            for c in chunk_outs {
                let msps =  obj_decode(c);
                msps_vec.push(msps);
            }

            let mut obj = BTreeMap::new();
            obj.insert("chunks".to_string(), obj_encode(&msps_vec));
            obj
    }
}


#[derive(RustcDecodable, RustcEncodable)]
pub struct MspChunk  {
    pub msp: Vec<String>,
}

// Specify encoding method manually
impl ToJson for MspChunk {
    fn to_json(&self) -> Json {
        let mut d = BTreeMap::new();
        d.insert("msp".to_string(), self.msp.to_json());
        Json::Object(d)
    }
}
